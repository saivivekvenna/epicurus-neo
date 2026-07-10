from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

from benchmark.metrics import identity_tiebreak


def _ranked_labels(df: pd.DataFrame, score_col: str, label_col: str, ascending: bool):
    work = df.copy()
    work["_score"] = pd.to_numeric(work[score_col], errors="coerce")
    work["_tiebreak"] = identity_tiebreak(work)
    for _patient, group in work.groupby("patient_id", sort=True):
        ranked = group.sort_values(
            ["_score", "_tiebreak"],
            ascending=[ascending, True],
            kind="mergesort",
            na_position="last",
        )
        yield ranked[label_col].to_numpy(dtype=float)


def patient_rank_vectors(
    df: pd.DataFrame,
    score_col: str,
    *,
    k_cap: int = 20,
    label_col: str = "label",
    ascending: bool = False,
) -> dict[str, np.ndarray]:
    """Per-patient ranking metrics with k = min(k_cap, n_candidates)."""
    keys = (
        "hits_at_k",
        "precision_at_k",
        "capture_fraction",
        "p_at_least_1",
        "p_at_least_2",
        "p_at_least_4",
    )
    out: dict[str, list[float]] = {key: [] for key in keys}
    for labels in _ranked_labels(df, score_col, label_col, ascending):
        k = min(k_cap, len(labels))
        top = labels[:k]
        hits = float(top.sum())
        positives = int(labels.sum())
        out["hits_at_k"].append(hits)
        out["precision_at_k"].append(hits / k if k else float("nan"))
        out["capture_fraction"].append(hits / min(positives, k) if positives else float("nan"))
        for threshold in (1, 2, 4):
            out[f"p_at_least_{threshold}"].append(
                float(hits >= threshold) if positives else float("nan")
            )
    return {key: np.asarray(values, dtype=float) for key, values in out.items()}


def classification_metrics(y_true: np.ndarray, y_score: np.ndarray, *, bins: int = 5) -> dict:
    """Threshold-free classification metrics over pooled out-of-fold predictions."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    finite = np.isfinite(y_score)
    y_true, y_score = y_true[finite], y_score[finite]
    edges = np.linspace(0.0, 1.0, bins + 1)
    which = np.clip(np.digitize(y_score, edges[1:-1]), 0, bins - 1)
    calibration = []
    for b in range(bins):
        mask = which == b
        if mask.any():
            calibration.append(
                {
                    "bin": b,
                    "n": int(mask.sum()),
                    "mean_score": float(y_score[mask].mean()),
                    "mean_label": float(y_true[mask].mean()),
                }
            )
    return {
        "auroc": float(roc_auc_score(y_true, y_score))
        if len(set(y_true.tolist())) > 1
        else float("nan"),
        "average_precision": float(average_precision_score(y_true, y_score))
        if y_true.any()
        else float("nan"),
        "brier": float(brier_score_loss(y_true, np.clip(y_score, 0.0, 1.0))),
        "calibration": calibration,
    }
