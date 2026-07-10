"""M6B Event-A -> Event-B transfer: a frozen Event-A teacher as one Event-B feature.

The teacher is trained ONLY on Event-A labels, using the same length-agnostic M6A core
features, then scores every Event-B candidate (long SLPs included). Its per-candidate
score is added as ``event_a_teacher_score`` -- the single variable that separates the
M6B candidate model from the frozen M6A baseline. Labels are never merged.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier

from epicurus_neo.m6.dataset import completeness_report
from epicurus_neo.m6.evaluate import macro_paired_delta
from epicurus_neo.m6.features import build_feature_matrix, feature_columns
from epicurus_neo.m6.loso import loso_folds
from epicurus_neo.m6.models import _logistic, fit_predict
from epicurus_neo.m6.ranking import classification_metrics, patient_rank_vectors

TEACHER_SCORE_COLUMN = "event_a_teacher_score"
HARM_MARGIN = -0.05
TRANSFER_QUESTION = (
    "Does Event-A information improve Event-B transfer beyond the Event-B-only model?"
)


@dataclass(frozen=True)
class FrozenTeacher:
    """A teacher fit on Event-A only; ``feature_cols`` pins the training column order."""

    model: object
    feature_cols: tuple[str, ...]
    tier: str


def _teacher_model(model_name: str, seed: int):
    if model_name == "logistic":
        return _logistic(seed)
    if model_name == "boosting":
        return HistGradientBoostingClassifier(
            max_depth=3, max_iter=200, learning_rate=0.05, random_state=seed
        )
    raise ValueError(f"unknown teacher model: {model_name!r}")


def train_frozen_teacher(
    event_a_frame: pd.DataFrame, *, tier: str = "core", model_name: str = "logistic", seed: int = 17
) -> FrozenTeacher:
    """Fit the frozen Event-A teacher on the shared M6A feature tier (Event-A labels only)."""
    matrix = build_feature_matrix(event_a_frame, tier)
    cols = feature_columns(matrix)
    model = _teacher_model(model_name, seed)
    model.fit(matrix[cols], matrix.label.to_numpy())
    return FrozenTeacher(model=model, feature_cols=tuple(cols), tier=tier)


def add_teacher_score(event_b_frame: pd.DataFrame, teacher: FrozenTeacher) -> pd.DataFrame:
    """Append ``event_a_teacher_score`` for every Event-B row using the frozen teacher.

    The Event-B feature matrix is reindexed to the teacher's training columns; any column
    the teacher never saw is left NaN and filled by the teacher's own fitted imputer
    (logistic) or treated as missing (boosting), so scoring is defined for every candidate.
    """
    matrix = build_feature_matrix(event_b_frame, teacher.tier)
    aligned = matrix.reindex(columns=list(teacher.feature_cols))
    scores = teacher.model.predict_proba(aligned)[:, 1].astype(float)
    out = event_b_frame.copy()
    out[TEACHER_SCORE_COLUMN] = scores
    return out


def assert_no_event_a_leakage(event_a_frame: pd.DataFrame, event_b_frame: pd.DataFrame) -> None:
    """Raise if any (mutant_peptide, hla_allele) pair is shared across Event-A and Event-B."""
    event_a_pairs = set(zip(event_a_frame.mutant_peptide, event_a_frame.hla_allele, strict=False))
    event_b_pairs = set(zip(event_b_frame.mutant_peptide, event_b_frame.hla_allele, strict=False))
    shared = event_a_pairs & event_b_pairs
    if shared:
        raise AssertionError(
            f"Event-A/Event-B peptide-HLA leakage: {len(shared)} shared pair(s), "
            f"e.g. {sorted(str(pair) for pair in shared)[:3]}"
        )


def _transfer_verdict(
    macro_auroc_delta: float, folds_improved: int, n_folds: int, macro_hits: dict, macro_p1: dict
) -> str:
    """The declared M6B gate (frozen in the pre-registration): AUROC-primary, >=3/4 folds."""
    if n_folds == 0 or not np.isfinite(macro_auroc_delta):
        return "CONSISTENT_WITH_NO_EFFECT_TRANSFER"

    def _lo(entry: dict) -> float:
        ci = entry.get("delta_ci", [None, None])
        return ci[0] if ci and ci[0] is not None else float("nan")

    harm = (_lo(macro_hits) < HARM_MARGIN) or (_lo(macro_p1) < HARM_MARGIN)
    improved = macro_auroc_delta > 0 and folds_improved >= 3
    worsened = macro_auroc_delta < 0 and (n_folds - folds_improved) >= 3
    if improved and not harm:
        return "ACCEPT_TRANSFER"
    if worsened or harm:
        return "REJECT_TRANSFER"
    return "CONSISTENT_WITH_NO_EFFECT_TRANSFER"


def evaluate_transfer_track(
    event_b_frame: pd.DataFrame,
    event_a_frame: pd.DataFrame,
    *,
    tier: str = "core",
    model_name: str = "logistic",
    k_cap: int = 20,
    seed: int = 17,
    bootstrap_n: int = 20_000,
) -> dict:
    """M6B: does a frozen Event-A teacher feature improve the Event-B-only model under LOSO?

    Candidate = ``logistic(core + event_a_teacher_score)``; baseline = ``logistic(core)`` (frozen
    M6A). The teacher feature is the sole difference. Verdict follows the pre-registered gate.
    """
    assert_no_event_a_leakage(event_a_frame, event_b_frame)
    teacher = train_frozen_teacher(event_a_frame, tier=tier, model_name=model_name, seed=seed)
    scored_frame = add_teacher_score(event_b_frame, teacher)

    report = completeness_report(scored_frame, k_cap=k_cap)
    rankable_patients = set(
        report.loc[report.denominator_type == "HAS_TESTED_NEGATIVE", "patient_id"]
    )

    per_fold: dict[str, dict] = {}
    hits_per_study: dict[str, tuple] = {}
    p1_per_study: dict[str, tuple] = {}
    informative_total = 0
    for fold in loso_folds(scored_frame):
        base_train = build_feature_matrix(fold.train, tier).reset_index(drop=True)
        base_eval = build_feature_matrix(fold.evaluation, tier).reset_index(drop=True)
        core_cols = [c for c in feature_columns(base_train) if c != TEACHER_SCORE_COLUMN]
        base_train[TEACHER_SCORE_COLUMN] = fold.train[TEACHER_SCORE_COLUMN].to_numpy()
        base_eval[TEACHER_SCORE_COLUMN] = fold.evaluation[TEACHER_SCORE_COLUMN].to_numpy()
        base_train["label"] = fold.train.label.to_numpy()
        cand_cols = [*core_cols, TEACHER_SCORE_COLUMN]
        y_eval = fold.evaluation.label.to_numpy()
        base_scores = fit_predict(model_name, base_train, base_eval, core_cols, seed=seed)
        cand_scores = fit_predict(model_name, base_train, base_eval, cand_cols, seed=seed)
        base_auroc = classification_metrics(y_eval, base_scores)["auroc"]
        cand_auroc = classification_metrics(y_eval, cand_scores)["auroc"]
        delta = (
            float(cand_auroc - base_auroc)
            if np.isfinite(base_auroc) and np.isfinite(cand_auroc)
            else float("nan")
        )
        per_fold[fold.held_out_study] = {
            "baseline_auroc": float(base_auroc),
            "candidate_auroc": float(cand_auroc),
            "auroc_delta": delta,
            "n_eval": int(len(fold.evaluation)),
        }
        eval_scored = fold.evaluation.assign(_cand=cand_scores, _base=base_scores)
        rankable = eval_scored[eval_scored.patient_id.isin(rankable_patients)]
        if not rankable.empty:
            cand_v = patient_rank_vectors(rankable, "_cand", k_cap=k_cap)
            base_v = patient_rank_vectors(rankable, "_base", k_cap=k_cap)
            hits_per_study[fold.held_out_study] = (cand_v["hits_at_k"], base_v["hits_at_k"])
            p1_per_study[fold.held_out_study] = (cand_v["p_at_least_1"], base_v["p_at_least_1"])
            informative_total += int((rankable.groupby("patient_id").size() > k_cap).sum())

    macro_hits = macro_paired_delta(hits_per_study, seed=seed, n=bootstrap_n)
    macro_p1 = macro_paired_delta(p1_per_study, seed=seed, n=bootstrap_n)
    finite = [e["auroc_delta"] for e in per_fold.values() if np.isfinite(e["auroc_delta"])]
    macro_auroc_delta = float(np.mean(finite)) if finite else float("nan")
    folds_improved = int(sum(1 for d in finite if d > 0))
    n_folds = len(finite)
    verdict = _transfer_verdict(macro_auroc_delta, folds_improved, n_folds, macro_hits, macro_p1)
    return {
        "track": "event_a_transfer",
        "question": TRANSFER_QUESTION,
        "verdict": verdict,
        "macro_auroc_delta": macro_auroc_delta,
        "folds_improved": folds_improved,
        "n_folds_scored": n_folds,
        "per_fold": per_fold,
        "macro_delta_hits_at_k": macro_hits,
        "macro_delta_p_at_least_1": macro_p1,
        "ranking_informative_patients": informative_total,
        "teacher": {
            "model": model_name,
            "tier": tier,
            "n_event_a": int(len(event_a_frame)),
            "n_event_a_positive": int((event_a_frame.label == 1).sum()),
        },
    }
