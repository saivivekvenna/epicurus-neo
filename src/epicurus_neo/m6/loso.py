from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class LosoFold:
    held_out_study: str
    train: pd.DataFrame
    evaluation: pd.DataFrame


def loso_folds(frame: pd.DataFrame) -> list[LosoFold]:
    """Yield one leave-one-study-out fold per study, in sorted order."""
    studies = sorted(frame.study_id.unique())
    if len(studies) < 2:
        raise ValueError("LOSO requires at least two studies")
    folds: list[LosoFold] = []
    for study in studies:
        is_eval = frame.study_id == study
        folds.append(
            LosoFold(
                held_out_study=study,
                train=frame[~is_eval].reset_index(drop=True),
                evaluation=frame[is_eval].reset_index(drop=True),
            )
        )
    return folds
