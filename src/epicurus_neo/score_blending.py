from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import pandas as pd

from epicurus_neo.metrics import group_metrics, summarize_group_metrics
from epicurus_neo.score_selection import _score_key


@dataclass(frozen=True)
class BlendSelection:
    default_weights: dict[str, float]
    group_weights: dict[str, dict[str, float]]
    validation_summary: dict[str, dict[str, float]]


def _positive_count(frame: pd.DataFrame) -> int:
    return int((frame["label"] == "positive").sum())


def _available_score_columns(frame: pd.DataFrame, score_columns: list[str]) -> list[str]:
    return [col for col in score_columns if col in frame.columns]


def _candidate_weight_sets(score_columns: list[str]) -> list[dict[str, float]]:
    weights: list[dict[str, float]] = [{col: 1.0} for col in score_columns]
    for left, right in combinations(score_columns, 2):
        for left_weight in (0.25, 0.5, 0.75):
            weights.append({left: left_weight, right: 1.0 - left_weight})
    return weights


def blend_name(weights: dict[str, float]) -> str:
    return "+".join(f"{weight:.2f}*{col}" for col, weight in sorted(weights.items()))


def add_rank_blend_score(
    frame: pd.DataFrame,
    *,
    group_col: str,
    weights: dict[str, float],
    output_col: str,
) -> pd.DataFrame:
    out = frame.copy()
    blended = pd.Series(0.0, index=out.index, dtype=float)
    for score_col, weight in weights.items():
        if score_col not in out.columns:
            raise ValueError(f"Missing score column for blend: {score_col}")
        numeric = pd.to_numeric(out[score_col], errors="coerce")
        ranks = numeric.groupby(out[group_col], sort=False).rank(method="average", pct=True)
        blended += float(weight) * ranks.fillna(0.0)
    out[output_col] = blended
    return out


def best_blend(
    frame: pd.DataFrame,
    *,
    group_col: str,
    score_columns: list[str],
    k: int = 20,
    objective: str = "hits",
) -> tuple[dict[str, float], dict[str, float]]:
    available = _available_score_columns(frame, score_columns)
    if not available:
        raise ValueError("No usable score columns found")

    best_weights: dict[str, float] = {}
    best_summary: dict[str, float] = {}
    best_key: tuple[float, float, float, float] | None = None
    for weights in _candidate_weight_sets(available):
        scored = add_rank_blend_score(
            frame,
            group_col=group_col,
            weights=weights,
            output_col="__epicurus_blend_candidate",
        )
        per_group = group_metrics(
            scored,
            group_col=group_col,
            score_col="__epicurus_blend_candidate",
            k=k,
        )
        summary = summarize_group_metrics(per_group)
        key = _score_key(summary, objective=objective)
        if best_key is None or key > best_key:
            best_key = key
            best_weights = weights
            best_summary = summary
    return best_weights, best_summary


def select_blends_by_group(
    validation: pd.DataFrame,
    *,
    group_col: str,
    score_columns: list[str],
    k: int = 20,
    min_positive: int = 1,
    objective: str = "hits",
) -> BlendSelection:
    default_weights, default_summary = best_blend(
        validation,
        group_col=group_col,
        score_columns=score_columns,
        k=k,
        objective=objective,
    )
    group_weights: dict[str, dict[str, float]] = {}
    validation_summary: dict[str, dict[str, float]] = {"__default__": default_summary}

    for group_value, group in validation.groupby(group_col):
        if _positive_count(group) < min_positive:
            group_weights[str(group_value)] = default_weights
            continue
        weights, summary = best_blend(
            group,
            group_col=group_col,
            score_columns=score_columns,
            k=k,
            objective=objective,
        )
        group_weights[str(group_value)] = weights
        validation_summary[str(group_value)] = summary

    return BlendSelection(
        default_weights=default_weights,
        group_weights=group_weights,
        validation_summary=validation_summary,
    )


def apply_blend_selection(
    frame: pd.DataFrame,
    selection: BlendSelection,
    *,
    group_col: str,
    output_col: str = "epicurus_blend_score",
) -> pd.DataFrame:
    out = frame.copy()
    out[output_col] = pd.NA
    out["epicurus_blend_score_source"] = ""
    for group_value, index in out.groupby(group_col).groups.items():
        weights = selection.group_weights.get(str(group_value), selection.default_weights)
        scored_group = add_rank_blend_score(
            out.loc[index],
            group_col=group_col,
            weights=weights,
            output_col=output_col,
        )
        out.loc[index, output_col] = scored_group[output_col]
        out.loc[index, "epicurus_blend_score_source"] = blend_name(weights)
    out[output_col] = pd.to_numeric(out[output_col], errors="coerce")
    return out


def apply_blend_selection_files(
    validation_path: str | Path,
    target_path: str | Path,
    output_path: str | Path,
    *,
    group_col: str,
    score_columns: list[str],
    k: int = 20,
    min_positive: int = 1,
    objective: str = "hits",
) -> tuple[Path, BlendSelection]:
    validation = pd.read_csv(validation_path)
    target = pd.read_csv(target_path)
    selection = select_blends_by_group(
        validation,
        group_col=group_col,
        score_columns=score_columns,
        k=k,
        min_positive=min_positive,
        objective=objective,
    )
    out = apply_blend_selection(target, selection, group_col=group_col)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output, selection
