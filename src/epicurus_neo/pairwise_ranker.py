from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from epicurus_neo.features import infer_numeric_feature_columns
from epicurus_neo.schema import supervised_rows


@dataclass(frozen=True)
class FittedPairwiseRanker:
    feature_columns: tuple[str, ...]
    model: Any

    def predict_scores(self, frame: pd.DataFrame, *, output_col: str = "epicurus_pairwise_score") -> pd.DataFrame:
        out = frame.copy()
        x = out.reindex(columns=list(self.feature_columns))
        out[output_col] = self.model.decision_function(x)
        return out


def _pairwise_training_rows(
    labeled: pd.DataFrame,
    *,
    feature_columns: list[str],
    group_col: str,
    max_negatives_per_positive: int,
) -> tuple[pd.DataFrame, np.ndarray]:
    diffs: list[pd.Series] = []
    labels: list[int] = []
    for _, group in labeled.groupby(group_col, dropna=False):
        positives = group[group["label"] == "positive"]
        negatives = group[group["label"] == "negative"]
        if positives.empty or negatives.empty:
            continue
        negatives = negatives.sort_values(feature_columns, ascending=False, kind="mergesort")
        for _, positive in positives.iterrows():
            sampled_negatives = negatives.head(max_negatives_per_positive)
            positive_features = positive.loc[feature_columns]
            for _, negative in sampled_negatives.iterrows():
                negative_features = negative.loc[feature_columns]
                diffs.append(positive_features - negative_features)
                labels.append(1)
                diffs.append(negative_features - positive_features)
                labels.append(0)
    if not diffs:
        return pd.DataFrame(columns=feature_columns), np.array([], dtype=int)
    return pd.DataFrame(diffs).reset_index(drop=True), np.array(labels, dtype=int)


def fit_pairwise_ranker(
    frame: pd.DataFrame,
    *,
    group_col: str,
    feature_columns: list[str] | None = None,
    max_negatives_per_positive: int = 20,
    random_state: int = 17,
) -> FittedPairwiseRanker:
    if group_col not in frame.columns:
        raise ValueError(f"Missing group column: {group_col}")
    if max_negatives_per_positive < 1:
        raise ValueError("max_negatives_per_positive must be at least 1")

    labeled = supervised_rows(frame)
    if feature_columns is None:
        feature_columns = infer_numeric_feature_columns(labeled)
    if not feature_columns:
        raise ValueError("No numeric feature columns available for pairwise ranking.")

    x_pairs, y_pairs = _pairwise_training_rows(
        labeled,
        feature_columns=feature_columns,
        group_col=group_col,
        max_negatives_per_positive=max_negatives_per_positive,
    )
    if len(y_pairs) == 0:
        raise ValueError("No positive/negative pairs available for pairwise ranking.")

    model = Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            ("scaler", StandardScaler()),
            (
                "classifier",
                LogisticRegression(
                    class_weight="balanced",
                    max_iter=1000,
                    random_state=random_state,
                ),
            ),
        ]
    )
    model.fit(x_pairs, y_pairs)
    return FittedPairwiseRanker(feature_columns=tuple(feature_columns), model=model)
