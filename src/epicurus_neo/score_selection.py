from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from epicurus_neo.metrics import group_metrics, summarize_group_metrics


@dataclass(frozen=True)
class ScoreSelection:
    default_score_col: str
    group_score_cols: dict[str, str]
    validation_summary: dict[str, dict[str, float]]


def _positive_count(frame: pd.DataFrame) -> int:
    return int((frame["label"] == "positive").sum())


def _score_key(summary: dict[str, float]) -> tuple[float, float, float, float]:
    return (
        summary["mean_hits_at_k"],
        summary["mean_recall_at_k"],
        summary["mean_ndcg_at_k"],
        summary["mean_mrr"],
    )


def best_score_column(
    frame: pd.DataFrame,
    *,
    group_col: str,
    score_columns: list[str],
    k: int = 20,
) -> tuple[str, dict[str, float]]:
    best_col = ""
    best_summary: dict[str, float] = {}
    best_key: tuple[float, float, float, float] | None = None
    for score_col in score_columns:
        if score_col not in frame.columns:
            continue
        per_group = group_metrics(frame, group_col=group_col, score_col=score_col, k=k)
        summary = summarize_group_metrics(per_group)
        key = _score_key(summary)
        if best_key is None or key > best_key:
            best_key = key
            best_col = score_col
            best_summary = summary
    if not best_col:
        raise ValueError("No usable score columns found")
    return best_col, best_summary


def select_score_columns_by_group(
    validation: pd.DataFrame,
    *,
    group_col: str,
    score_columns: list[str],
    k: int = 20,
    min_positive: int = 1,
) -> ScoreSelection:
    default_col, default_summary = best_score_column(
        validation,
        group_col=group_col,
        score_columns=score_columns,
        k=k,
    )
    group_score_cols: dict[str, str] = {}
    validation_summary: dict[str, dict[str, float]] = {"__default__": default_summary}

    for group_value, group in validation.groupby(group_col):
        if _positive_count(group) < min_positive:
            group_score_cols[str(group_value)] = default_col
            continue
        best_col, summary = best_score_column(
            group,
            group_col=group_col,
            score_columns=score_columns,
            k=k,
        )
        group_score_cols[str(group_value)] = best_col
        validation_summary[str(group_value)] = summary

    return ScoreSelection(
        default_score_col=default_col,
        group_score_cols=group_score_cols,
        validation_summary=validation_summary,
    )


def apply_score_selection(
    frame: pd.DataFrame,
    selection: ScoreSelection,
    *,
    group_col: str,
    output_col: str = "epicurus_selected_score",
) -> pd.DataFrame:
    out = frame.copy()
    out[output_col] = pd.NA
    out["epicurus_selected_score_source"] = ""
    for group_value, index in out.groupby(group_col).groups.items():
        score_col = selection.group_score_cols.get(str(group_value), selection.default_score_col)
        out.loc[index, output_col] = pd.to_numeric(out.loc[index, score_col], errors="coerce")
        out.loc[index, "epicurus_selected_score_source"] = score_col
    out[output_col] = pd.to_numeric(out[output_col], errors="coerce")
    return out


def apply_score_selection_files(
    validation_path: str | Path,
    target_path: str | Path,
    output_path: str | Path,
    *,
    group_col: str,
    score_columns: list[str],
    k: int = 20,
    min_positive: int = 1,
) -> tuple[Path, ScoreSelection]:
    validation = pd.read_csv(validation_path)
    target = pd.read_csv(target_path)
    selection = select_score_columns_by_group(
        validation,
        group_col=group_col,
        score_columns=score_columns,
        k=k,
        min_positive=min_positive,
    )
    out = apply_score_selection(target, selection, group_col=group_col)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output, selection
