from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from epicurus_neo.benchmark import BenchmarkResult, train_and_evaluate
from epicurus_neo.leakage import detect_exact_leakage, purge_train_overlaps
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
    purge_exact_overlaps: bool = True,
) -> list[FoldResult]:
    """Run leave-group-out evaluation while blocking exact leakage."""
    folds: list[FoldResult] = []
    include_shared_studies = group_col == "study_id"
    for split in leave_group_out_splits(frame, group_col=group_col, max_splits=max_splits):
        test_groups = tuple(sorted(split.test[group_col].dropna().astype(str).unique()))
        train = split.train
        if purge_exact_overlaps:
            train = purge_train_overlaps(train, split.test)
            split_leakage = detect_exact_leakage(train, split.test)
        else:
            split_leakage = split.leakage

        if train.empty or split.test.empty:
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
        if split_leakage.has_blocking_leakage(include_shared_studies=include_shared_studies):
            folds.append(
                FoldResult(
                    name=split.name,
                    status="leakage_blocked",
                    test_groups=test_groups,
                    feature_columns=(),
                    benchmark_results=(),
                    reason=str(split_leakage),
                )
            )
            continue

        result = train_and_evaluate(
            train,
            split.test,
            group_col=metric_group_col,
            feature_columns=feature_columns,
            k=k,
            include_shared_studies_as_leakage=include_shared_studies,
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
