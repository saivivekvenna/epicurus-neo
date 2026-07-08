from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from epicurus_neo.features import add_baseline_scores
from epicurus_neo.leakage import LeakageReport, detect_exact_leakage
from epicurus_neo.metrics import group_metrics, summarize_group_metrics
from epicurus_neo.model import FittedEpicurusRanker, fit_ranker
from epicurus_neo.pairwise_ranker import fit_pairwise_ranker


@dataclass(frozen=True)
class BenchmarkResult:
    score_col: str
    summary: dict[str, float]


@dataclass(frozen=True)
class TrainEvaluateResult:
    leakage: LeakageReport
    feature_columns: tuple[str, ...]
    benchmark_results: tuple[BenchmarkResult, ...]
    scored_test: pd.DataFrame


def add_groupwise_ensemble_scores(
    frame: pd.DataFrame,
    *,
    group_col: str,
    component_cols: list[str],
    output_col: str = "epicurus_blend_score",
) -> pd.DataFrame:
    out = frame.copy()
    available = [col for col in component_cols if col in out.columns]
    if not available:
        return out

    rank_components: list[pd.Series] = []
    for col in available:
        ranks = out.groupby(group_col)[col].rank(method="average", pct=True)
        rank_components.append(ranks)
    out[output_col] = sum(rank_components) / len(rank_components)
    return out


def add_weighted_groupwise_score(
    frame: pd.DataFrame,
    *,
    group_col: str,
    weights: dict[str, float],
    output_col: str,
) -> pd.DataFrame:
    out = frame.copy()
    usable = {col: weight for col, weight in weights.items() if col in out.columns and weight > 0}
    if not usable:
        return out

    total = sum(usable.values())
    score = 0.0
    for col, weight in usable.items():
        percentile_rank = out.groupby(group_col)[col].rank(method="average", pct=True)
        score = score + (weight / total) * percentile_rank
    out[output_col] = score
    return out


def add_transferable_presentation_score(
    frame: pd.DataFrame,
    *,
    group_col: str,
    output_col: str = "epicurus_transfer_score",
) -> pd.DataFrame:
    """Add a peptide/HLA transfer score that does not rely on source-specific predictors."""
    required = {
        "mhcflurry_presentation_score",
        "seq_hydrophobicity_mean",
        "seq_cysteine_fraction",
        "seq_aromatic_fraction",
    }
    if not required.issubset(frame.columns):
        return frame.copy()

    out = frame.copy()
    presentation = out.groupby(group_col)["mhcflurry_presentation_score"].rank(
        method="average", pct=True
    )
    hydrophilic = 1.0 - out.groupby(group_col)["seq_hydrophobicity_mean"].rank(
        method="average", pct=True
    )
    low_cysteine = 1.0 - out.groupby(group_col)["seq_cysteine_fraction"].rank(
        method="average", pct=True
    )
    low_aromatic = 1.0 - out.groupby(group_col)["seq_aromatic_fraction"].rank(
        method="average", pct=True
    )
    out[output_col] = (
        0.70 * presentation
        + 0.15 * hydrophilic
        + 0.10 * low_cysteine
        + 0.05 * low_aromatic
    )
    return out


def add_retrieval_score(
    frame: pd.DataFrame,
    *,
    output_col: str = "epicurus_retrieval_score",
) -> pd.DataFrame:
    if "retrieval_max_positive_similarity" not in frame.columns:
        return frame.copy()
    out = frame.copy()
    out[output_col] = pd.to_numeric(out["retrieval_max_positive_similarity"], errors="coerce")
    return out


def evaluate_score_columns(
    frame: pd.DataFrame,
    *,
    group_col: str,
    score_columns: list[str],
    k: int = 20,
) -> list[BenchmarkResult]:
    results: list[BenchmarkResult] = []
    for score_col in score_columns:
        if score_col not in frame.columns:
            continue
        per_group = group_metrics(frame, group_col=group_col, score_col=score_col, k=k)
        results.append(
            BenchmarkResult(
                score_col=score_col,
                summary=summarize_group_metrics(per_group),
            )
        )
    return results


def train_and_evaluate(
    train: pd.DataFrame,
    test: pd.DataFrame,
    *,
    group_col: str = "patient_id",
    feature_columns: list[str] | None = None,
    k: int = 20,
    allow_exact_leakage: bool = False,
    include_shared_studies_as_leakage: bool = True,
    uncertainty_ensemble_size: int = 1,
    uncertainty_penalty: float = 1.0,
) -> TrainEvaluateResult:
    leakage = detect_exact_leakage(train, test)
    if (
        leakage.has_blocking_leakage(include_shared_studies=include_shared_studies_as_leakage)
        and not allow_exact_leakage
    ):
        raise ValueError(f"Exact train/test leakage detected: {leakage}")

    train_features = add_baseline_scores(train)
    test_features = add_baseline_scores(test)
    ranker: FittedEpicurusRanker = fit_ranker(
        train_features,
        feature_columns=feature_columns,
        uncertainty_ensemble_size=uncertainty_ensemble_size,
        uncertainty_penalty=uncertainty_penalty,
    )
    scored_test = ranker.predict_scores(test_features)
    scored_test = add_groupwise_ensemble_scores(
        scored_test,
        group_col=group_col,
        component_cols=[
            "epicurus_immunogenicity_prob",
            "baseline_gartner_nmer_score",
            "baseline_netmhcpan_el_score",
            "baseline_mhcflurry_score",
        ],
    )
    scored_test = add_weighted_groupwise_score(
        scored_test,
        group_col=group_col,
        weights={
            "baseline_gartner_nmer_score": 0.9,
            "baseline_netmhcpan_el_score": 0.1,
        },
        output_col="epicurus_hits20_score",
    )
    scored_test = add_transferable_presentation_score(scored_test, group_col=group_col)
    scored_test = add_retrieval_score(scored_test)
    try:
        pairwise_ranker = fit_pairwise_ranker(train_features, group_col=group_col)
        scored_test = pairwise_ranker.predict_scores(scored_test)
    except ValueError:
        pass

    score_columns = [
        "epicurus_retrieval_score",
        "epicurus_transfer_score",
        "epicurus_hits20_score",
        "epicurus_blend_score",
        "epicurus_pairwise_score",
        "epicurus_lower_confidence_score",
        "epicurus_score",
        "epicurus_immunogenicity_prob",
        "baseline_gartner_nmer_score",
        "baseline_netmhcpan_el_score",
        "baseline_mhcflurry_score",
        "baseline_pvac_style_score",
        "baseline_presentation_score",
        "baseline_binding_score",
    ]
    results = evaluate_score_columns(
        scored_test, group_col=group_col, score_columns=score_columns, k=k
    )

    return TrainEvaluateResult(
        leakage=leakage,
        feature_columns=ranker.feature_columns,
        benchmark_results=tuple(results),
        scored_test=scored_test,
    )
