from __future__ import annotations

import numpy as np
import pandas as pd

from benchmark.scorecard import pre_registered_verdict
from benchmark.stats import paired_bootstrap
from epicurus_neo.m6.dataset import completeness_report
from epicurus_neo.m6.features import build_feature_matrix, feature_columns
from epicurus_neo.m6.loso import loso_folds
from epicurus_neo.m6.models import fit_predict
from epicurus_neo.m6.ranking import classification_metrics, patient_rank_vectors

_GUARD_METRICS = ("hits_at_k", "capture_fraction", "p_at_least_1")


def _entry(candidate: np.ndarray, baseline: np.ndarray, *, n: int = 20_000) -> dict:
    comparison = paired_bootstrap(candidate, baseline, n=n, seed=17)
    return {
        "delta_vs_baseline": comparison.delta,
        "delta_ci": [comparison.lo, comparison.hi],
        "p_better": comparison.p_better,
        "n": comparison.n,
    }


def _finite_diffs(candidate: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    diff = np.asarray(candidate, dtype=float) - np.asarray(baseline, dtype=float)
    return diff[np.isfinite(diff)]


def macro_paired_delta(
    per_study: dict[str, tuple[np.ndarray, np.ndarray]], *, seed: int = 17, n: int = 20_000
) -> dict:
    """Equal-weight-per-study delta with a study-stratified patient bootstrap CI."""
    studies = sorted(per_study)
    diffs = {study: _finite_diffs(*per_study[study]) for study in studies}
    populated = [study for study in studies if len(diffs[study])]
    point = (
        float(np.mean([diffs[study].mean() for study in populated])) if populated else float("nan")
    )
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(n):
        study_means = [
            diffs[study][rng.integers(0, len(diffs[study]), size=len(diffs[study]))].mean()
            for study in populated
        ]
        if study_means:
            draws.append(float(np.mean(study_means)))
    if not draws:
        return {"delta": point, "delta_ci": [float("nan"), float("nan")], "per_study": {}}
    lo, hi = np.quantile(draws, [0.025, 0.975])
    return {
        "delta": point,
        "delta_ci": [float(lo), float(hi)],
        "per_study": {
            study: float(diffs[study].mean()) if len(diffs[study]) else float("nan")
            for study in studies
        },
    }


def _macro_entry(macro: dict) -> dict:
    return {"delta_vs_baseline": macro["delta"], "delta_ci": macro["delta_ci"]}


def _macro_classification(per_fold: dict, pooled_true: list, pooled_score: list) -> dict:
    """Macro-average per-held-out-study classification; pooled OOF is demoted/caveated.

    Pooling out-of-fold scores across studies with 10%->83% prevalence conflates
    cross-study calibration shift with discrimination and reads far worse than any
    single fold, so the macro-average is the honest headline.
    """
    folds = [pf["classification"] for pf in per_fold.values() if "classification" in pf]

    def _macro(metric: str) -> float:
        values = [f[metric] for f in folds if np.isfinite(f[metric])]
        return float(np.mean(values)) if values else float("nan")

    pooled = classification_metrics(np.concatenate(pooled_true), np.concatenate(pooled_score))
    return {
        "auroc": _macro("auroc"),
        "brier": _macro("brier"),
        "average_precision": _macro("average_precision"),
        "per_fold_auroc": {
            study: pf["classification"]["auroc"]
            for study, pf in per_fold.items()
            if "classification" in pf
        },
        "pooled_out_of_fold_auroc": pooled["auroc"],
        "pooled_note": (
            "pooled OOF AUROC conflates cross-study calibration shift with "
            "discrimination; macro/per-fold is the honest read"
        ),
        "calibration": pooled["calibration"],
    }


def evaluate_track(
    frame: pd.DataFrame,
    *,
    model_name: str,
    baseline_name: str,
    track: str,
    tier: str = "core",
    k_cap: int = 20,
    seed: int = 17,
    bootstrap_n: int = 20_000,
) -> dict:
    """Run LOSO for one (model vs baseline) comparison and assemble the registered report."""
    report = completeness_report(frame, k_cap=k_cap)
    rankable_patients = set(
        report.loc[report.denominator_type == "HAS_TESTED_NEGATIVE", "patient_id"]
    )

    per_study: dict[str, dict[str, tuple[np.ndarray, np.ndarray]]] = {m: {} for m in _GUARD_METRICS}
    micro: dict[str, dict[str, list]] = {m: {"cand": [], "base": []} for m in _GUARD_METRICS}
    pooled_true, pooled_score = [], []
    per_fold = {}
    informative_total = 0
    for fold in loso_folds(frame):
        train = build_feature_matrix(fold.train, tier)
        evaluation = build_feature_matrix(fold.evaluation, tier)
        cols = feature_columns(train)
        model_scores = fit_predict(model_name, train, evaluation, cols, seed=seed)
        base_scores = fit_predict(baseline_name, train, evaluation, cols, seed=seed)
        scored = evaluation.assign(_model=model_scores, _base=base_scores)
        # Pooled classification uses every candidate; primary top-k excludes
        # patients with no tested negative (nothing to rank against).
        pooled_true.append(scored.label.to_numpy())
        pooled_score.append(model_scores)
        fold_class = classification_metrics(scored.label.to_numpy(), model_scores)
        fold_report = {
            "classification": {k: fold_class[k] for k in ("auroc", "average_precision", "brier")}
        }
        rankable = scored[scored.patient_id.isin(rankable_patients)]
        if rankable.empty:
            per_fold[fold.held_out_study] = {
                "n_patients": 0,
                "n_ranking_informative": 0,
                **fold_report,
            }
            continue
        cand = patient_rank_vectors(rankable, "_model", k_cap=k_cap)
        base = patient_rank_vectors(rankable, "_base", k_cap=k_cap)
        for metric in _GUARD_METRICS:
            per_study[metric][fold.held_out_study] = (cand[metric], base[metric])
            micro[metric]["cand"].append(cand[metric])
            micro[metric]["base"].append(base[metric])
        informative = int((rankable.groupby("patient_id").size() > k_cap).sum())
        informative_total += informative
        per_fold[fold.held_out_study] = {
            "hits_at_k": _entry(cand["hits_at_k"], base["hits_at_k"], n=bootstrap_n),
            "n_patients": int(rankable.patient_id.nunique()),
            "n_ranking_informative": informative,
            **fold_report,
        }
    macro = {
        metric: macro_paired_delta(per_study[metric], seed=seed, n=bootstrap_n)
        for metric in _GUARD_METRICS
    }
    micro_hits = _entry(
        np.concatenate(micro["hits_at_k"]["cand"]),
        np.concatenate(micro["hits_at_k"]["base"]),
        n=bootstrap_n,
    )
    classification = _macro_classification(per_fold, pooled_true, pooled_score)
    verdict = pre_registered_verdict(
        _macro_entry(macro["hits_at_k"]),
        _macro_entry(macro["capture_fraction"]),
        _macro_entry(macro["p_at_least_1"]),
    )
    return {
        "track": track,
        "model": model_name,
        "baseline": baseline_name,
        "macro_hits_at_k": macro["hits_at_k"],
        "macro_capture": macro["capture_fraction"],
        "macro_p_at_least_1": macro["p_at_least_1"],
        "micro_hits_at_k": micro_hits,
        "per_fold": per_fold,
        "classification": classification,
        "verdict": verdict,
        "ranking_informative_patients": informative_total,
    }
