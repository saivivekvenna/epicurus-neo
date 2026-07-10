"""M6B Event-A -> Event-B transfer: a frozen Event-A teacher as one Event-B feature.

The teacher is trained ONLY on Event-A labels, using the same length-agnostic M6A core
features, then scores every Event-B candidate (long SLPs included). Its per-candidate
score is added as ``event_a_teacher_score`` -- the single variable that separates the
M6B candidate model from the frozen M6A baseline. Labels are never merged.
"""

from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier

from epicurus_neo.m6.features import build_feature_matrix, feature_columns
from epicurus_neo.m6.models import _logistic

TEACHER_SCORE_COLUMN = "event_a_teacher_score"


@dataclass(frozen=True)
class FrozenTeacher:
    """A teacher fit on Event-A only; ``feature_cols`` pins the training column order."""

    model: object
    feature_cols: tuple[str, ...]
    tier: str


def _teacher_model(model_name: str, seed: int):
    if model_name == "logistic":
        return _logistic(seed)
    if model_name == "boosting":
        return HistGradientBoostingClassifier(
            max_depth=3, max_iter=200, learning_rate=0.05, random_state=seed
        )
    raise ValueError(f"unknown teacher model: {model_name!r}")


def train_frozen_teacher(
    event_a_frame: pd.DataFrame, *, tier: str = "core", model_name: str = "logistic", seed: int = 17
) -> FrozenTeacher:
    """Fit the frozen Event-A teacher on the shared M6A feature tier (Event-A labels only)."""
    matrix = build_feature_matrix(event_a_frame, tier)
    cols = feature_columns(matrix)
    model = _teacher_model(model_name, seed)
    model.fit(matrix[cols], matrix.label.to_numpy())
    return FrozenTeacher(model=model, feature_cols=tuple(cols), tier=tier)


def add_teacher_score(event_b_frame: pd.DataFrame, teacher: FrozenTeacher) -> pd.DataFrame:
    """Append ``event_a_teacher_score`` for every Event-B row using the frozen teacher.

    The Event-B feature matrix is reindexed to the teacher's training columns; any column
    the teacher never saw is left NaN and filled by the teacher's own fitted imputer
    (logistic) or treated as missing (boosting), so scoring is defined for every candidate.
    """
    matrix = build_feature_matrix(event_b_frame, teacher.tier)
    aligned = matrix.reindex(columns=list(teacher.feature_cols))
    scores = teacher.model.predict_proba(aligned)[:, 1].astype(float)
    out = event_b_frame.copy()
    out[TEACHER_SCORE_COLUMN] = scores
    return out


def assert_no_event_a_leakage(event_a_frame: pd.DataFrame, event_b_frame: pd.DataFrame) -> None:
    """Raise if any (mutant_peptide, hla_allele) pair is shared across Event-A and Event-B."""
    event_a_pairs = set(zip(event_a_frame.mutant_peptide, event_a_frame.hla_allele, strict=False))
    event_b_pairs = set(zip(event_b_frame.mutant_peptide, event_b_frame.hla_allele, strict=False))
    shared = event_a_pairs & event_b_pairs
    if shared:
        raise AssertionError(
            f"Event-A/Event-B peptide-HLA leakage: {len(shared)} shared pair(s), "
            f"e.g. {sorted(str(pair) for pair in shared)[:3]}"
        )
