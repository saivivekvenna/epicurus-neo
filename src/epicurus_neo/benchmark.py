from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from epicurus_neo.features import add_baseline_scores
from epicurus_neo.leakage import LeakageReport, detect_exact_leakage
from epicurus_neo.metrics import group_metrics, summarize_group_metrics
from epicurus_neo.model import FittedEpicurusRanker, fit_ranker


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
) -> TrainEvaluateResult:
    leakage = detect_exact_leakage(train, test)
    if (
        leakage.has_blocking_leakage(include_shared_studies=include_shared_studies_as_leakage)
        and not allow_exact_leakage
    ):
        raise ValueError(f"Exact train/test leakage detected: {leakage}")

    train_features = add_baseline_scores(train)
    test_features = add_baseline_scores(test)
    ranker: FittedEpicurusRanker = fit_ranker(train_features, feature_columns=feature_columns)
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

    score_columns = [
        "epicurus_blend_score",
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
