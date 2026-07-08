from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from epicurus_neo.benchmark import BenchmarkResult, train_and_evaluate
from epicurus_neo.splits import leave_group_out_splits


@dataclass(frozen=True)
class FoldResult:
    name: str
    status: str
    test_groups: tuple[str, ...]
    feature_columns: tuple[str, ...]
    benchmark_results: tuple[BenchmarkResult, ...]
    reason: str = ""


def grouped_cross_validate(
    frame: pd.DataFrame,
    *,
    group_col: str,
    metric_group_col: str = "patient_id",
    feature_columns: list[str] | None = None,
    k: int = 20,
    max_splits: int | None = None,
) -> list[FoldResult]:
    """Run leave-group-out evaluation while blocking exact leakage."""
    folds: list[FoldResult] = []
    for split in leave_group_out_splits(frame, group_col=group_col, max_splits=max_splits):
        test_groups = tuple(sorted(split.test[group_col].dropna().astype(str).unique()))
        if split.train.empty or split.test.empty:
            folds.append(
                FoldResult(
                    name=split.name,
                    status="skipped",
                    test_groups=test_groups,
                    feature_columns=(),
                    benchmark_results=(),
                    reason="empty train or test split",
                )
            )
            continue
        if split.leakage.has_leakage:
            folds.append(
                FoldResult(
                    name=split.name,
                    status="leakage_blocked",
                    test_groups=test_groups,
                    feature_columns=(),
                    benchmark_results=(),
                    reason=str(split.leakage),
                )
            )
            continue

        result = train_and_evaluate(
            split.train,
            split.test,
            group_col=metric_group_col,
            feature_columns=feature_columns,
            k=k,
        )
        folds.append(
            FoldResult(
                name=split.name,
                status="ok",
                test_groups=test_groups,
                feature_columns=result.feature_columns,
                benchmark_results=result.benchmark_results,
            )
        )
    return folds


def summarize_cross_validation(folds: list[FoldResult]) -> dict[str, object]:
    ok_folds = [fold for fold in folds if fold.status == "ok"]
    blocked = [fold for fold in folds if fold.status == "leakage_blocked"]
    summaries: dict[str, list[dict[str, float]]] = {}

    for fold in ok_folds:
        for benchmark in fold.benchmark_results:
            summaries.setdefault(benchmark.score_col, []).append(benchmark.summary)

    aggregate: dict[str, dict[str, float]] = {}
    for score_col, items in summaries.items():
        keys = sorted({key for item in items for key in item})
        aggregate[score_col] = {
            key: float(sum(item.get(key, 0.0) for item in items) / len(items))
            for key in keys
        }

    return {
        "folds": len(folds),
        "ok_folds": len(ok_folds),
        "leakage_blocked_folds": len(blocked),
        "aggregate": aggregate,
    }

