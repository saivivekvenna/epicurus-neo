from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import pandas as pd

from epicurus_neo.leakage import LeakageReport, detect_exact_leakage


@dataclass(frozen=True)
class SplitFrames:
    name: str
    train: pd.DataFrame
    test: pd.DataFrame
    leakage: LeakageReport


def assign_holdout_split(
    frame: pd.DataFrame,
    *,
    group_col: str,
    holdout_values: Iterable[str],
    split_col: str = "split",
) -> pd.DataFrame:
    holdouts = {str(value) for value in holdout_values}
    out = frame.copy()
    out[split_col] = out[group_col].astype(str).map(lambda value: "test" if value in holdouts else "train")
    return out


def split_from_column(
    frame: pd.DataFrame,
    *,
    split_col: str = "split",
    train_value: str = "train",
    test_value: str = "test",
    name: str = "split",
) -> SplitFrames:
    train = frame[frame[split_col] == train_value].copy()
    test = frame[frame[split_col] == test_value].copy()
    leakage = detect_exact_leakage(train, test)
    return SplitFrames(name=name, train=train, test=test, leakage=leakage)


def leave_group_out_splits(
    frame: pd.DataFrame,
    *,
    group_col: str,
    max_splits: int | None = None,
) -> list[SplitFrames]:
    splits: list[SplitFrames] = []
    values = sorted(frame[group_col].dropna().astype(str).unique())
    if max_splits is not None:
        values = values[:max_splits]

    for value in values:
        assigned = assign_holdout_split(frame, group_col=group_col, holdout_values=[value])
        split = split_from_column(assigned, name=f"leave-{group_col}-{value}-out")
        splits.append(split)
    return splits

