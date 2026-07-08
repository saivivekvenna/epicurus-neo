from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


METRIC_COLUMNS = [
    "mean_hits_at_k",
    "mean_precision_at_k",
    "mean_recall_at_k",
    "mean_ndcg_at_k",
    "mean_mrr",
]


def compare_metric_reports(
    paths: list[str | Path],
    *,
    sort_by: str = "mean_hits_at_k",
) -> pd.DataFrame:
    if sort_by not in METRIC_COLUMNS:
        raise ValueError(f"Unsupported sort metric: {sort_by}")

    rows: list[dict[str, object]] = []
    for path_like in paths:
        path = Path(path_like)
        payload = json.loads(path.read_text())
        for benchmark in payload.get("benchmarks", []):
            summary = benchmark.get("summary", {})
            row: dict[str, object] = {
                "report": str(path),
                "table": payload.get("table", ""),
                "score_col": benchmark.get("score_col", ""),
            }
            for metric in METRIC_COLUMNS:
                row[metric] = summary.get(metric, 0.0)
            rows.append(row)

    frame = pd.DataFrame(rows)
    if frame.empty:
        return pd.DataFrame(columns=["report", "table", "score_col", *METRIC_COLUMNS])
    return frame.sort_values(sort_by, ascending=False, kind="mergesort").reset_index(drop=True)


def compare_metric_reports_file(
    paths: list[str | Path],
    output_path: str | Path,
    *,
    sort_by: str = "mean_hits_at_k",
) -> Path:
    frame = compare_metric_reports(paths, sort_by=sort_by)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    return output
