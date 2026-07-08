import json

import pandas as pd
import pytest

from epicurus_neo.cli import build_parser, cmd_compare_metrics
from epicurus_neo.metric_compare import compare_metric_reports, compare_metric_reports_file


def _write_report(path, score_col: str, hits: float, ndcg: float) -> None:
    payload = {
        "table": "scored.csv",
        "group_col": "hla_allele",
        "k": 20,
        "benchmarks": [
            {
                "score_col": score_col,
                "summary": {
                    "mean_hits_at_k": hits,
                    "mean_precision_at_k": hits / 20.0,
                    "mean_recall_at_k": 0.5,
                    "mean_ndcg_at_k": ndcg,
                    "mean_mrr": 0.25,
                },
            }
        ],
    }
    path.write_text(json.dumps(payload) + "\n")


def test_compare_metric_reports_sorts_by_metric(tmp_path):
    low = tmp_path / "low.json"
    high = tmp_path / "high.json"
    _write_report(low, "low_score", hits=1.0, ndcg=0.9)
    _write_report(high, "high_score", hits=2.0, ndcg=0.1)

    frame = compare_metric_reports([low, high], sort_by="mean_hits_at_k")

    assert frame["score_col"].tolist() == ["high_score", "low_score"]


def test_compare_metric_reports_rejects_unknown_sort_metric(tmp_path):
    report = tmp_path / "report.json"
    _write_report(report, "score", hits=1.0, ndcg=0.5)

    with pytest.raises(ValueError, match="Unsupported sort metric"):
        compare_metric_reports([report], sort_by="not_a_metric")


def test_compare_metric_reports_file_writes_csv(tmp_path):
    report = tmp_path / "report.json"
    output = tmp_path / "compare.csv"
    _write_report(report, "score", hits=1.0, ndcg=0.5)

    assert compare_metric_reports_file([report], output).exists()
    assert pd.read_csv(output).loc[0, "score_col"] == "score"


def test_compare_metrics_cli_writes_output(tmp_path):
    report = tmp_path / "report.json"
    output = tmp_path / "compare.csv"
    _write_report(report, "score", hits=1.0, ndcg=0.5)

    parser = build_parser()
    args = parser.parse_args(
        [
            "compare-metrics",
            str(report),
            "--output",
            str(output),
        ]
    )

    assert cmd_compare_metrics(args) == 0
    assert "mean_hits_at_k" in output.read_text()
