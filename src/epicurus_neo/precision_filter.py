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
) -> tuple[Path, PrecisionThreshold, dict[str, float | int]]:
    validation = pd.read_csv(validation_path)
    target = pd.read_csv(target_path)
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
