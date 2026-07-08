from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


@dataclass(frozen=True)
class PrecisionThreshold:
    score_col: str
    threshold: float
    target_precision: float
    validation_selected: int
    validation_hits: int
    validation_precision: float
    validation_recall: float
    achieved_target: bool


@dataclass(frozen=True)
class GroupedPrecisionThreshold:
    score_col: str
    group_col: str
    default_threshold: PrecisionThreshold
    group_thresholds: dict[str, PrecisionThreshold]
    target_precision: float
    min_group_positives: int
    min_group_selected: int


def _positive_count(frame: pd.DataFrame) -> int:
    return int((frame["label"] == "positive").sum())


def _threshold_summary(
    frame: pd.DataFrame,
    *,
    score_col: str,
    threshold: float,
) -> dict[str, float | int]:
    labeled = frame[frame["label"].isin(["positive", "negative"])].copy()
    scores = pd.to_numeric(labeled[score_col], errors="coerce")
    selected = labeled[scores >= threshold]
    hits = _positive_count(selected)
    selected_count = len(selected)
    positives = _positive_count(labeled)
    return {
        "selected": selected_count,
        "hits": hits,
        "precision": hits / selected_count if selected_count else 0.0,
        "recall": hits / positives if positives else 0.0,
    }


def calibrate_precision_threshold(
    validation: pd.DataFrame,
    *,
    score_col: str,
    target_precision: float = 0.5,
    min_selected: int = 1,
) -> PrecisionThreshold:
    if score_col not in validation.columns:
        raise ValueError(f"Missing score column: {score_col}")
    if target_precision <= 0.0 or target_precision > 1.0:
        raise ValueError("target_precision must be in (0, 1]")
    if min_selected < 1:
        raise ValueError("min_selected must be at least 1")

    labeled = validation[validation["label"].isin(["positive", "negative"])].copy()
    if labeled.empty:
        raise ValueError("Validation table has no positive/negative labels")
    scores = pd.to_numeric(labeled[score_col], errors="coerce").dropna()
    if scores.empty:
        raise ValueError(f"Validation score column has no numeric values: {score_col}")

    candidates: list[tuple[float, dict[str, float | int]]] = []
    for threshold in sorted(scores.unique(), reverse=True):
        summary = _threshold_summary(labeled, score_col=score_col, threshold=float(threshold))
        if int(summary["selected"]) >= min_selected:
            candidates.append((float(threshold), summary))
    if not candidates:
        raise ValueError("No threshold selected enough validation rows")

    passing = [
        (threshold, summary)
        for threshold, summary in candidates
        if float(summary["precision"]) >= target_precision
    ]
    achieved_target = bool(passing)
    if passing:
        threshold, summary = max(
            passing,
            key=lambda item: (int(item[1]["selected"]), float(item[1]["precision"])),
        )
    else:
        threshold, summary = max(
            candidates,
            key=lambda item: (float(item[1]["precision"]), int(item[1]["selected"])),
        )

    return PrecisionThreshold(
        score_col=score_col,
        threshold=threshold,
        target_precision=target_precision,
        validation_selected=int(summary["selected"]),
        validation_hits=int(summary["hits"]),
        validation_precision=float(summary["precision"]),
        validation_recall=float(summary["recall"]),
        achieved_target=achieved_target,
    )


def calibrate_grouped_precision_threshold(
    validation: pd.DataFrame,
    *,
    group_col: str,
    score_col: str,
    target_precision: float = 0.5,
    min_selected: int = 1,
    min_group_positives: int = 2,
    min_group_selected: int = 1,
) -> GroupedPrecisionThreshold:
    if group_col not in validation.columns:
        raise ValueError(f"Missing group column: {group_col}")
    if min_group_positives < 1:
        raise ValueError("min_group_positives must be at least 1")
    if min_group_selected < 1:
        raise ValueError("min_group_selected must be at least 1")

    default = calibrate_precision_threshold(
        validation,
        score_col=score_col,
        target_precision=target_precision,
        min_selected=min_selected,
    )

    group_thresholds: dict[str, PrecisionThreshold] = {}
    labeled = validation[validation["label"].isin(["positive", "negative"])].copy()
    for group_value, group in labeled.groupby(group_col, dropna=True):
        positives = _positive_count(group)
        if positives < min_group_positives:
            continue
        if len(group) < min_group_selected:
            continue
        try:
            threshold = calibrate_precision_threshold(
                group,
                score_col=score_col,
                target_precision=target_precision,
                min_selected=min_group_selected,
            )
        except ValueError:
            continue
        if threshold.achieved_target:
            group_thresholds[str(group_value)] = threshold

    return GroupedPrecisionThreshold(
        score_col=score_col,
        group_col=group_col,
        default_threshold=default,
        group_thresholds=group_thresholds,
        target_precision=target_precision,
        min_group_positives=min_group_positives,
        min_group_selected=min_group_selected,
    )


def apply_precision_threshold(
    frame: pd.DataFrame,
    threshold: PrecisionThreshold,
    *,
    output_col: str = "epicurus_precision_selected",
) -> pd.DataFrame:
    if threshold.score_col not in frame.columns:
        raise ValueError(f"Missing score column: {threshold.score_col}")
    out = frame.copy()
    scores = pd.to_numeric(out[threshold.score_col], errors="coerce")
    out[output_col] = scores >= threshold.threshold
    return out


def apply_grouped_precision_threshold(
    frame: pd.DataFrame,
    threshold: GroupedPrecisionThreshold,
    *,
    output_col: str = "epicurus_precision_selected",
) -> pd.DataFrame:
    if threshold.score_col not in frame.columns:
        raise ValueError(f"Missing score column: {threshold.score_col}")
    if threshold.group_col not in frame.columns:
        raise ValueError(f"Missing group column: {threshold.group_col}")

    out = frame.copy()
    scores = pd.to_numeric(out[threshold.score_col], errors="coerce")
    out[output_col] = scores >= threshold.default_threshold.threshold
    out["epicurus_precision_threshold"] = threshold.default_threshold.threshold
    out["epicurus_precision_threshold_source"] = "default"

    groups = out[threshold.group_col].astype("string")
    for group_value, group_threshold in threshold.group_thresholds.items():
        mask = groups == group_value
        out.loc[mask, output_col] = scores.loc[mask] >= group_threshold.threshold
        out.loc[mask, "epicurus_precision_threshold"] = group_threshold.threshold
        out.loc[mask, "epicurus_precision_threshold_source"] = "group"
    return out


def precision_selection_summary(
    frame: pd.DataFrame,
    *,
    selected_col: str = "epicurus_precision_selected",
) -> dict[str, float | int]:
    if selected_col not in frame.columns:
        raise ValueError(f"Missing selected column: {selected_col}")
    selected = frame[frame[selected_col].astype(bool)]
    selected_count = len(selected)
    payload: dict[str, float | int] = {"selected": selected_count}
    if "label" in frame.columns and frame["label"].isin(["positive", "negative"]).any():
        hits = _positive_count(selected)
        positives = _positive_count(frame)
        payload.update(
            {
                "hits": hits,
                "precision": hits / selected_count if selected_count else 0.0,
                "recall": hits / positives if positives else 0.0,
            }
        )
    return payload


def apply_precision_threshold_files(
    validation_path: str | Path,
    target_path: str | Path,
    output_path: str | Path,
    *,
    score_col: str,
    target_precision: float = 0.5,
    min_selected: int = 1,
    group_col: str | None = None,
    min_group_positives: int = 2,
    min_group_selected: int = 1,
) -> tuple[Path, PrecisionThreshold | GroupedPrecisionThreshold, dict[str, float | int]]:
    validation = pd.read_csv(validation_path)
    target = pd.read_csv(target_path)
    if group_col:
        threshold = calibrate_grouped_precision_threshold(
            validation,
            group_col=group_col,
            score_col=score_col,
            target_precision=target_precision,
            min_selected=min_selected,
            min_group_positives=min_group_positives,
            min_group_selected=min_group_selected,
        )
        out = apply_grouped_precision_threshold(target, threshold)
    else:
        threshold = calibrate_precision_threshold(
            validation,
            score_col=score_col,
            target_precision=target_precision,
            min_selected=min_selected,
        )
        out = apply_precision_threshold(target, threshold)
    summary = precision_selection_summary(out)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output, threshold, summary
