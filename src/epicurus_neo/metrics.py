from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class GroupMetrics:
    group_id: str
    positives: int
    evaluated: int
    hits_at_k: int
    precision_at_k: float
    recall_at_k: float
    ndcg_at_k: float
    mrr: float


def _dcg(labels: np.ndarray) -> float:
    return float(sum((2**rel - 1) / math.log2(idx + 2) for idx, rel in enumerate(labels)))


def _binary_labels(labels: pd.Series) -> np.ndarray:
    return labels.map({"positive": 1, "negative": 0}).astype(int).to_numpy()


def group_metrics(
    frame: pd.DataFrame,
    *,
    group_col: str,
    score_col: str,
    label_col: str = "label",
    k: int = 20,
) -> list[GroupMetrics]:
    """Compute top-k metrics per patient/case group."""
    metrics: list[GroupMetrics] = []
    labeled = frame[frame[label_col].isin(["positive", "negative"])].copy()

    for group_id, group in labeled.groupby(group_col, sort=True):
        ranked = group.sort_values(score_col, ascending=False, kind="mergesort")
        labels = _binary_labels(ranked[label_col])
        positives = int(labels.sum())
        top = labels[:k]
        hits = int(top.sum())
        evaluated = int(min(k, len(ranked)))
        precision = hits / evaluated if evaluated else 0.0
        recall = hits / positives if positives else 0.0

        ideal = np.sort(labels)[::-1][:k]
        ideal_dcg = _dcg(ideal)
        ndcg = _dcg(top) / ideal_dcg if ideal_dcg > 0 else 0.0

        positive_positions = np.flatnonzero(labels == 1)
        mrr = 1.0 / float(positive_positions[0] + 1) if len(positive_positions) else 0.0

        metrics.append(
            GroupMetrics(
                group_id=str(group_id),
                positives=positives,
                evaluated=evaluated,
                hits_at_k=hits,
                precision_at_k=precision,
                recall_at_k=recall,
                ndcg_at_k=ndcg,
                mrr=mrr,
            )
        )
    return metrics


def summarize_group_metrics(metrics: list[GroupMetrics]) -> dict[str, float]:
    if not metrics:
        return {
            "groups": 0.0,
            "mean_hits_at_k": 0.0,
            "mean_precision_at_k": 0.0,
            "mean_recall_at_k": 0.0,
            "mean_ndcg_at_k": 0.0,
            "mean_mrr": 0.0,
        }

    return {
        "groups": float(len(metrics)),
        "mean_hits_at_k": float(np.mean([m.hits_at_k for m in metrics])),
        "mean_precision_at_k": float(np.mean([m.precision_at_k for m in metrics])),
        "mean_recall_at_k": float(np.mean([m.recall_at_k for m in metrics])),
        "mean_ndcg_at_k": float(np.mean([m.ndcg_at_k for m in metrics])),
        "mean_mrr": float(np.mean([m.mrr for m in metrics])),
    }


def expected_calibration_error(
    y_true: np.ndarray,
    y_prob: np.ndarray,
    *,
    bins: int = 10,
) -> float:
    """Compute binary expected calibration error with equal-width probability bins."""
    if len(y_true) == 0:
        return 0.0
    y_true = np.asarray(y_true, dtype=float)
    y_prob = np.asarray(y_prob, dtype=float)
    boundaries = np.linspace(0.0, 1.0, bins + 1)
    ece = 0.0

    for lower, upper in zip(boundaries[:-1], boundaries[1:]):
        if upper == 1.0:
            mask = (y_prob >= lower) & (y_prob <= upper)
        else:
            mask = (y_prob >= lower) & (y_prob < upper)
        if not np.any(mask):
            continue
        confidence = float(np.mean(y_prob[mask]))
        accuracy = float(np.mean(y_true[mask]))
        ece += float(np.mean(mask)) * abs(accuracy - confidence)
    return ece

