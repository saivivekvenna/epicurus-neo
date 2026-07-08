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
    ensemble_immunogenicity_models: tuple[Any, ...] = ()
    ensemble_dud_models: tuple[Any, ...] = ()
    uncertainty_penalty: float = 1.0

    def predict_scores(self, frame: pd.DataFrame) -> pd.DataFrame:
        x = frame.reindex(columns=list(self.feature_columns))
        immunogenicity_models = self.ensemble_immunogenicity_models or (
            self.immunogenicity_model,
        )
        dud_models = self.ensemble_dud_models or (self.dud_model,)
        immunogenicity_matrix = np.column_stack(
            [model.predict_proba(x)[:, 1] for model in immunogenicity_models]
        )
        dud_matrix = np.column_stack([model.predict_proba(x)[:, 1] for model in dud_models])
        score_matrix = immunogenicity_matrix * (1.0 - dud_matrix)

        immunogenicity_prob = immunogenicity_matrix.mean(axis=1)
        dud_prob = dud_matrix.mean(axis=1)
        score = score_matrix.mean(axis=1)
        score_std = score_matrix.std(axis=1)

        out = frame.copy()
        out["epicurus_immunogenicity_prob"] = immunogenicity_prob
        out["epicurus_immunogenicity_prob_std"] = immunogenicity_matrix.std(axis=1)
        out["epicurus_dud_prob"] = dud_prob
        out["epicurus_dud_prob_std"] = dud_matrix.std(axis=1)
        out["epicurus_score"] = score
        out["epicurus_score_std"] = score_std
        out["epicurus_lower_confidence_score"] = np.clip(
            score - self.uncertainty_penalty * score_std,
            0.0,
            1.0,
        )
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


def _bootstrap_indices(
    labeled: pd.DataFrame,
    *,
    random_state: int,
) -> np.ndarray:
    rng = np.random.default_rng(random_state)
    positive_idx = labeled.index[labeled["label"] == "positive"].to_numpy()
    negative_idx = labeled.index[labeled["label"] == "negative"].to_numpy()
    if len(positive_idx) == 0 or len(negative_idx) == 0:
        return rng.choice(labeled.index.to_numpy(), size=len(labeled), replace=True)
    positive_sample = rng.choice(positive_idx, size=len(positive_idx), replace=True)
    negative_sample = rng.choice(negative_idx, size=len(negative_idx), replace=True)
    sampled = np.concatenate([positive_sample, negative_sample])
    rng.shuffle(sampled)
    return sampled


def fit_ranker(
    frame: pd.DataFrame,
    *,
    feature_columns: list[str] | None = None,
    random_state: int = 17,
    uncertainty_ensemble_size: int = 1,
    uncertainty_penalty: float = 1.0,
) -> FittedEpicurusRanker:
    labeled = supervised_rows(frame)
    if labeled.empty:
        raise ValueError("No supervised rows available; positive/negative labels are required.")
    if uncertainty_ensemble_size < 1:
        raise ValueError("uncertainty_ensemble_size must be at least 1")
    if uncertainty_penalty < 0:
        raise ValueError("uncertainty_penalty must be non-negative")

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
    ensemble_immunogenicity_models: list[Any] = []
    ensemble_dud_models: list[Any] = []
    if uncertainty_ensemble_size > 1:
        for offset in range(uncertainty_ensemble_size):
            sampled_idx = _bootstrap_indices(labeled, random_state=random_state + 1000 + offset)
            boot = labeled.loc[sampled_idx]
            boot_x = boot.loc[:, feature_columns]
            boot_y_immunogenic = (boot["label"] == "positive").astype(int).to_numpy()
            boot_y_dud = (boot["label"] == "negative").astype(int).to_numpy()
            boot_weights = (
                pd.to_numeric(boot.get("label_weight", 1.0), errors="coerce")
                .fillna(1.0)
                .to_numpy()
            )
            ensemble_immunogenicity_models.append(
                _fit_binary_model(
                    boot_x,
                    boot_y_immunogenic,
                    boot_weights,
                    random_state=random_state + 2000 + offset,
                )
            )
            ensemble_dud_models.append(
                _fit_binary_model(
                    boot_x,
                    boot_y_dud,
                    boot_weights,
                    random_state=random_state + 3000 + offset,
                )
            )

    return FittedEpicurusRanker(
        feature_columns=tuple(feature_columns),
        immunogenicity_model=immunogenicity_model,
        dud_model=dud_model,
        ensemble_immunogenicity_models=tuple(ensemble_immunogenicity_models),
        ensemble_dud_models=tuple(ensemble_dud_models),
        uncertainty_penalty=uncertainty_penalty,
    )
