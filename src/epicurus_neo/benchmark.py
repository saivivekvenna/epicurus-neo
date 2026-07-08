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
) -> TrainEvaluateResult:
    leakage = detect_exact_leakage(train, test)
    if leakage.has_leakage and not allow_exact_leakage:
        raise ValueError(f"Exact train/test leakage detected: {leakage}")

    train_features = add_baseline_scores(train)
    test_features = add_baseline_scores(test)
    ranker: FittedEpicurusRanker = fit_ranker(train_features, feature_columns=feature_columns)
    scored_test = ranker.predict_scores(test_features)

    score_columns = [
        "epicurus_score",
        "epicurus_immunogenicity_prob",
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

