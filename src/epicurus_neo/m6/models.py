from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

REQUIRED_MODELS = ("prevalence", "logistic", "boosting")


def _logistic(seed: int) -> Pipeline:
    return Pipeline(
        [
            ("impute", SimpleImputer(strategy="median")),
            ("scale", StandardScaler()),
            ("clf", LogisticRegression(max_iter=1000, random_state=seed)),
        ]
    )


def fit_predict(
    model_name: str,
    train: pd.DataFrame,
    evaluation: pd.DataFrame,
    feature_cols: list[str],
    *,
    seed: int = 17,
) -> np.ndarray:
    """Fit ``model_name`` on train and return one POSITIVE-likelihood score per eval row."""
    if model_name == "prevalence":
        return np.full(len(evaluation), float(train.label.mean()), dtype=float)
    if model_name == "presentation":
        if "presentation_score" not in evaluation.columns:
            raise ValueError("presentation model requires a presentation_score column")
        return pd.to_numeric(evaluation.presentation_score, errors="coerce").to_numpy(dtype=float)

    x_train, y_train = train[feature_cols], train.label.to_numpy()
    if model_name == "logistic":
        model = _logistic(seed)
    elif model_name == "boosting":
        model = HistGradientBoostingClassifier(
            max_depth=3, max_iter=200, learning_rate=0.05, random_state=seed
        )
    else:
        raise ValueError(f"unknown model: {model_name!r}")
    model.fit(x_train, y_train)
    return model.predict_proba(evaluation[feature_cols])[:, 1].astype(float)
