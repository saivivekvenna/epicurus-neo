from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.model_selection import StratifiedGroupKFold

from epicurus_neo.m6.features import build_feature_matrix, feature_columns


def prevalence_by_study(frame: pd.DataFrame) -> pd.DataFrame:
    table = frame.groupby("study_id").agg(
        n=("label", "size"), positives=("label", "sum")
    ).reset_index()
    table["positive_rate"] = table.positives / table.n
    return table[["study_id", "n", "positive_rate"]]


def study_only_classifier(frame: pd.DataFrame, *, seed: int = 17) -> dict:
    """Can pre-vaccine core features predict which study a candidate came from?"""
    matrix = build_feature_matrix(frame, "core")
    cols = feature_columns(matrix)
    x = matrix[cols].to_numpy()
    y = matrix.study_id.to_numpy()
    groups = matrix.patient_id.to_numpy()
    splitter = StratifiedGroupKFold(n_splits=4, shuffle=True, random_state=seed)
    predictions = np.empty(len(y), dtype=object)
    for train_idx, test_idx in splitter.split(x, y, groups):
        model = HistGradientBoostingClassifier(max_depth=3, random_state=seed)
        model.fit(x[train_idx], y[train_idx])
        predictions[test_idx] = model.predict(x[test_idx])
    correct = predictions == y
    per_study = {study: float(correct[y == study].mean()) for study in sorted(set(y.tolist()))}
    _, counts = np.unique(y, return_counts=True)
    return {
        "accuracy": float(correct.mean()),
        "majority_rate": float(counts.max() / counts.sum()),
        "per_study": per_study,
    }
