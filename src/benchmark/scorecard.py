"""The sole Milestone-1 reporting path for ranked candidate experiments."""

from __future__ import annotations

from collections.abc import Callable
from hashlib import sha256
from typing import Any

import numpy as np
import pandas as pd

from benchmark.metrics import (
    capture_fraction,
    group_positive_counts,
    hits_at_k,
    mrr,
    p_at_least_one,
    precision_at_k,
)
from benchmark.headroom import random_expectation
from benchmark.metrics import identity_tiebreak
from benchmark.stats import bootstrap_ci, mde, paired_bootstrap


Metric = Callable[..., np.ndarray]


def _entry(candidate: np.ndarray, baseline: np.ndarray) -> dict[str, Any]:
    estimate = bootstrap_ci(candidate)
    comparison = paired_bootstrap(candidate, baseline)
    return {
        "value": estimate.mean,
        "ci": [estimate.lo, estimate.hi],
        "delta_vs_baseline": comparison.delta,
        "delta_ci": [comparison.lo, comparison.hi],
        "p_better": comparison.p_better,
        "n": comparison.n,
    }


def scorecard(
    df: pd.DataFrame,
    score_col: str,
    baseline_col: str,
    group_col: str = "patient_id",
    k: int = 20,
    *,
    label_col: str = "label",
    ascending: bool = False,
    baseline_ascending: bool = False,
) -> dict[str, Any]:
    """Compute the complete pre-registered scorecard, including bad news."""
    metrics: list[tuple[str, Metric]] = [
        (f"hits@{k}", hits_at_k),
        ("capture_fraction", capture_fraction),
        ("p_at_least_one", p_at_least_one),
        (f"precision@{k}", precision_at_k),
        ("mrr", mrr),
    ]
    report: dict[str, Any] = {"score_col": score_col, "baseline_col": baseline_col, "k": k}
    vectors: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for name, function in metrics:
        candidate = function(
            df,
            group_col=group_col,
            score_col=score_col,
            label_col=label_col,
            k=k,
            ascending=ascending,
        )
        baseline = function(
            df,
            group_col=group_col,
            score_col=baseline_col,
            label_col=label_col,
            k=k,
            ascending=baseline_ascending,
        )
        vectors[name] = (candidate, baseline)
        report[name] = _entry(candidate, baseline)

    positive_counts = group_positive_counts(df, group_col=group_col, label_col=label_col)
    unreachable = int((positive_counts == 0).sum())
    report["unreachable_patients"] = {
        "count": unreachable,
        "total": int(len(positive_counts)),
        "display": f"{unreachable} / {len(positive_counts)}",
    }
    identities = np.sort(identity_tiebreak(df).to_numpy(dtype="uint64"))
    # A compact, deterministic candidate-universe identifier. It prevents a
    # selection-conditioned result from silently changing its denominator.
    fingerprint = sha256(identities.tobytes()).hexdigest()[:16]
    report["candidate_list"] = {
        "rows": int(len(df)),
        "groups": int(len(positive_counts)),
        "identity_sha256": fingerprint,
    }
    report["random_baseline"] = random_expectation(
        df, group_col=group_col, label_col=label_col, k=k
    )
    report["candidate_recall"] = None  # Milestone 3.
    primary_name = f"hits@{k}"
    try:
        report["mde_at_current_n"] = float(mde(*vectors[primary_name]))
    except ValueError:
        report["mde_at_current_n"] = None

    primary = report[primary_name]
    co_primary = report["capture_fraction"]
    clinical = report["p_at_least_one"]
    primary_significant = primary["delta_ci"][0] > 0.0
    no_co_primary_regression = co_primary["delta_ci"][1] >= 0.0
    no_clinical_regression = clinical["delta_ci"][1] >= 0.0
    if primary_significant and no_co_primary_regression and no_clinical_regression:
        verdict = "ACCEPT"
    elif (
        primary["delta_vs_baseline"] > 0.0
        and primary["delta_ci"][0] <= 0.0 <= primary["delta_ci"][1]
    ):
        verdict = "CONSISTENT_WITH_NO_EFFECT"
    else:
        verdict = "REJECT"
    report["verdict"] = verdict
    return report
