"""Candidate-list random baselines, leakage canaries, and reranking headroom."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from benchmark.metrics import _positive_mask, hits_at_k, identity_tiebreak


def random_expectation(
    df: pd.DataFrame,
    *,
    group_col: str = "patient_id",
    label_col: str = "label",
    k: int = 20,
) -> float:
    """Analytic expected mean hits@k under uniform random rankings."""
    values: list[float] = []
    for _, group in df.groupby(group_col, sort=True, dropna=False):
        positives = int(_positive_mask(group[label_col]).sum())
        values.append(min(k, len(group)) * positives / len(group) if len(group) else 0.0)
    return float(np.mean(values)) if values else float("nan")


def source_order_hits(
    df: pd.DataFrame,
    *,
    group_col: str = "patient_id",
    label_col: str = "label",
    k: int = 20,
) -> float:
    """Unsafe diagnostic only: expose the exact source-order leakage canary."""
    values = [
        float(_positive_mask(group[label_col]).astype(int)[:k].sum())
        for _, group in df.groupby(group_col, sort=True, dropna=False)
    ]
    return float(np.mean(values))


def tie_break_canary(
    df: pd.DataFrame,
    *,
    group_col: str = "patient_id",
    label_col: str = "label",
    k: int = 20,
    leaked_value: float = 2.4714285714285715,
) -> dict[str, float]:
    """Fail if a constant score reproduces source-order top-k membership."""
    work = df.copy()
    work["_constant_score"] = 1.0
    deterministic = float(
        hits_at_k(
            work,
            group_col=group_col,
            score_col="_constant_score",
            label_col=label_col,
            k=k,
        ).mean()
    )
    source = source_order_hits(work, group_col=group_col, label_col=label_col, k=k)
    if np.isclose(deterministic, leaked_value) or np.isclose(deterministic, source):
        raise AssertionError(
            f"Tie-break leakage detected: constant-score hits@{k}={deterministic:.4f}"
        )
    return {
        "constant_md5": deterministic,
        "source_order": source,
        "random_expectation": random_expectation(
            work, group_col=group_col, label_col=label_col, k=k
        ),
    }


def headroom_table(
    df: pd.DataFrame,
    *,
    base_score_col: str,
    group_col: str = "patient_id",
    label_col: str = "label",
    slates: tuple[int | None, ...] = (20, 50, 100, 200, None),
    k: int = 20,
    ascending: bool = False,
) -> list[dict[str, Any]]:
    """Report base hits, oracle within each base-ranked slate, and unreachable groups."""
    work = df.copy()
    work["_tiebreak"] = identity_tiebreak(work)
    rows: list[dict[str, Any]] = []
    for slate in slates:
        base_values: list[float] = []
        oracle_values: list[float] = []
        zero_positive = 0
        for _, group in work.groupby(group_col, sort=True, dropna=False):
            ranked = group.sort_values(
                [base_score_col, "_tiebreak"],
                ascending=[ascending, True],
                kind="mergesort",
            )
            selected = ranked if slate is None else ranked.head(slate)
            labels = _positive_mask(selected[label_col])
            positives = int(labels.sum())
            base_values.append(float(labels[:k].sum()))
            oracle_values.append(float(min(positives, k)))
            zero_positive += int(positives == 0)
        rows.append(
            {
                "slate": "all" if slate is None else slate,
                "base": float(np.mean(base_values)),
                "oracle": float(np.mean(oracle_values)),
                "headroom": float(np.mean(oracle_values) - np.mean(base_values)),
                "zero_positive_patients": zero_positive,
                "total_patients": len(base_values),
            }
        )
    return rows
