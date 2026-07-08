from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline

from epicurus_neo.features import infer_numeric_feature_columns
from epicurus_neo.schema import supervised_rows


class ConstantProbabilityModel:
    def __init__(self, probability: float):
        self.probability = float(probability)

    def predict_proba(self, x: pd.DataFrame | np.ndarray) -> np.ndarray:
        count = len(x)
        positive = np.full(count, self.probability)
        return np.column_stack([1.0 - positive, positive])


@dataclass(frozen=True)
class FittedEpicurusRanker:
    feature_columns: tuple[str, ...]
    immunogenicity_model: Any
    dud_model: Any

    def predict_scores(self, frame: pd.DataFrame) -> pd.DataFrame:
        x = frame.reindex(columns=list(self.feature_columns))
        immunogenicity_prob = self.immunogenicity_model.predict_proba(x)[:, 1]
        dud_prob = self.dud_model.predict_proba(x)[:, 1]

        out = frame.copy()
        out["epicurus_immunogenicity_prob"] = immunogenicity_prob
        out["epicurus_dud_prob"] = dud_prob
        out["epicurus_score"] = immunogenicity_prob * (1.0 - dud_prob)
        return out


def _make_classifier(random_state: int) -> Pipeline:
    return Pipeline(
        steps=[
            ("imputer", SimpleImputer(strategy="median")),
            (
                "classifier",
                HistGradientBoostingClassifier(
                    learning_rate=0.06,
                    max_iter=200,
                    l2_regularization=0.05,
                    class_weight="balanced",
                    random_state=random_state,
                ),
            ),
        ]
    )


def _fit_binary_model(
    x: pd.DataFrame,
    y: np.ndarray,
    sample_weight: np.ndarray,
    *,
    random_state: int,
) -> Any:
    classes = set(y.tolist())
    if classes == {0}:
        return ConstantProbabilityModel(0.0)
    if classes == {1}:
        return ConstantProbabilityModel(1.0)

    model = _make_classifier(random_state)
    model.fit(x, y, classifier__sample_weight=sample_weight)
    return model


def fit_ranker(
    frame: pd.DataFrame,
    *,
    feature_columns: list[str] | None = None,
    random_state: int = 17,
) -> FittedEpicurusRanker:
    labeled = supervised_rows(frame)
    if labeled.empty:
        raise ValueError("No supervised rows available; positive/negative labels are required.")

    if feature_columns is None:
        feature_columns = infer_numeric_feature_columns(labeled)
    if not feature_columns:
        raise ValueError("No numeric feature columns available for training.")

    x = labeled.loc[:, feature_columns]
    y_immunogenic = (labeled["label"] == "positive").astype(int).to_numpy()
    y_dud = (labeled["label"] == "negative").astype(int).to_numpy()
    weights = pd.to_numeric(labeled.get("label_weight", 1.0), errors="coerce").fillna(1.0).to_numpy()

    immunogenicity_model = _fit_binary_model(
        x, y_immunogenic, weights, random_state=random_state
    )
    dud_model = _fit_binary_model(
        x, y_dud, weights, random_state=random_state + 1
    )

    return FittedEpicurusRanker(
        feature_columns=tuple(feature_columns),
        immunogenicity_model=immunogenicity_model,
        dud_model=dud_model,
    )
