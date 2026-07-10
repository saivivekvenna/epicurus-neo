"""Model-free label-frequency accounting for Positive-Unlabeled datasets."""

from __future__ import annotations

import numpy as np


def estimate_label_frequency(known_positive: np.ndarray, labeled_positive: np.ndarray) -> float:
    """Estimate c=P(labeled|positive) from an ascertainment/audit subset.

    This is intentionally an estimator over known audit labels; it does not fit
    a recognition model in the instrumentation milestone.
    """
    positive = np.asarray(known_positive, dtype=bool)
    labeled = np.asarray(labeled_positive, dtype=bool)
    if positive.shape != labeled.shape:
        raise ValueError("known_positive and labeled_positive must have the same shape")
    if np.any(labeled & ~positive):
        raise ValueError("a labeled positive cannot be outside the known-positive set")
    denominator = int(positive.sum())
    if denominator == 0:
        raise ValueError("cannot estimate c without known positives")
    return float((positive & labeled).sum() / denominator)
