import pandas as pd

from epicurus_neo.cli import build_parser, cmd_precision_filter
from epicurus_neo.precision_filter import (
    apply_grouped_precision_threshold,
    apply_precision_threshold,
    calibrate_grouped_precision_threshold,
    calibrate_precision_threshold,
    precision_selection_summary,
)


def test_calibrate_precision_threshold_selects_largest_passing_set():
    validation = pd.DataFrame(
        {
            "label": ["positive", "negative", "positive", "negative"],
            "score": [0.9, 0.8, 0.7, 0.1],
        }
    )

    threshold = calibrate_precision_threshold(
        validation,
        score_col="score",
        target_precision=0.5,
    )

    assert threshold.threshold == 0.1
    assert threshold.validation_selected == 4
    assert threshold.validation_hits == 2
    assert threshold.achieved_target


def test_calibrate_precision_threshold_falls_back_to_best_available_precision():
    validation = pd.DataFrame(
        {
            "label": ["negative", "positive", "negative"],
            "score": [0.9, 0.8, 0.7],
        }
    )

    threshold = calibrate_precision_threshold(
        validation,
        score_col="score",
        target_precision=0.8,
    )

    assert threshold.threshold == 0.8
    assert threshold.validation_precision == 0.5
    assert not threshold.achieved_target


def test_apply_precision_threshold_and_summary():
    target = pd.DataFrame(
        {
            "label": ["positive", "negative", "negative"],
            "score": [0.9, 0.8, 0.1],
        }
    )
    threshold = calibrate_precision_threshold(target, score_col="score", target_precision=0.5)

    out = apply_precision_threshold(target, threshold)
    summary = precision_selection_summary(out)

    assert out["epicurus_precision_selected"].tolist() == [True, True, False]
    assert summary["selected"] == 2
    assert summary["hits"] == 1
    assert summary["precision"] == 0.5


def test_grouped_precision_threshold_uses_group_threshold_when_supported():
    validation = pd.DataFrame(
        {
            "label": ["positive", "negative", "negative", "positive", "negative", "positive"],
            "score": [0.9, 0.8, 0.2, 0.7, 0.6, 0.1],
            "hla": ["A", "A", "A", "B", "B", "B"],
        }
    )

    threshold = calibrate_grouped_precision_threshold(
        validation,
        group_col="hla",
        score_col="score",
        target_precision=0.5,
        min_selected=1,
        min_group_positives=1,
        min_group_selected=1,
    )
    target = pd.DataFrame(
        {
            "label": ["negative", "positive", "positive"],
            "score": [0.8, 0.7, 0.05],
            "hla": ["A", "B", "C"],
        }
    )

    out = apply_grouped_precision_threshold(target, threshold)

    assert set(threshold.group_thresholds) == {"A", "B"}
    assert out["epicurus_precision_selected"].tolist() == [True, True, False]
    assert out["epicurus_precision_threshold_source"].tolist() == ["group", "group", "default"]


def test_grouped_precision_threshold_skips_low_evidence_groups():
    validation = pd.DataFrame(
        {
            "label": ["positive", "negative", "positive", "negative"],
            "score": [0.9, 0.8, 0.7, 0.6],
            "hla": ["A", "A", "B", "B"],
        }
    )

    threshold = calibrate_grouped_precision_threshold(
        validation,
        group_col="hla",
        score_col="score",
        target_precision=0.5,
        min_group_positives=2,
    )

    assert threshold.group_thresholds == {}
    assert threshold.default_threshold.achieved_target


def test_precision_filter_cli_writes_outputs(tmp_path):
    validation = tmp_path / "validation.csv"
    target = tmp_path / "target.csv"
    output = tmp_path / "selected.csv"
    report = tmp_path / "report.json"
    pd.DataFrame(
        {
            "label": ["positive", "negative", "positive", "negative"],
            "score": [0.9, 0.8, 0.7, 0.1],
        }
    ).to_csv(validation, index=False)
    pd.DataFrame(
        {
            "label": ["positive", "negative"],
            "score": [0.9, 0.1],
        }
    ).to_csv(target, index=False)

    parser = build_parser()
    args = parser.parse_args(
        [
            "precision-filter",
            "--validation",
            str(validation),
            "--target",
            str(target),
            "--output",
            str(output),
            "--report-output",
            str(report),
            "--score-col",
            "score",
            "--target-precision",
            "0.5",
        ]
    )

    assert cmd_precision_filter(args) == 0
    assert "epicurus_precision_selected" in pd.read_csv(output).columns
    assert "target_summary" in report.read_text()


def test_precision_filter_cli_writes_grouped_outputs(tmp_path):
    validation = tmp_path / "validation.csv"
    target = tmp_path / "target.csv"
    output = tmp_path / "selected.csv"
    report = tmp_path / "report.json"
    pd.DataFrame(
        {
            "label": ["positive", "negative", "positive", "negative"],
            "score": [0.9, 0.8, 0.7, 0.1],
            "hla": ["A", "A", "B", "B"],
        }
    ).to_csv(validation, index=False)
    pd.DataFrame(
        {
            "label": ["positive", "negative"],
            "score": [0.9, 0.1],
            "hla": ["A", "C"],
        }
    ).to_csv(target, index=False)

    parser = build_parser()
    args = parser.parse_args(
        [
            "precision-filter",
            "--validation",
            str(validation),
            "--target",
            str(target),
            "--output",
            str(output),
            "--report-output",
            str(report),
            "--score-col",
            "score",
            "--group-col",
            "hla",
            "--target-precision",
            "0.5",
        ]
    )

    assert cmd_precision_filter(args) == 0
    scored = pd.read_csv(output)
    assert "epicurus_precision_threshold_source" in scored.columns
    assert "group_thresholds" in report.read_text()
