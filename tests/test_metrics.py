import numpy as np
import pandas as pd

from epicurus_neo.metrics import expected_calibration_error, group_metrics, summarize_group_metrics


def test_group_metrics_focus_on_top_k():
    frame = pd.DataFrame(
        {
            "patient_id": ["p1", "p1", "p1", "p2", "p2"],
            "label": ["negative", "positive", "positive", "positive", "negative"],
            "score": [0.9, 0.8, 0.1, 0.7, 0.6],
        }
    )

    metrics = group_metrics(frame, group_col="patient_id", score_col="score", k=2)
    summary = summarize_group_metrics(metrics)

    assert metrics[0].hits_at_k == 1
    assert metrics[0].precision_at_k == 0.5
    assert metrics[0].recall_at_k == 0.5
    assert metrics[1].hits_at_k == 1
    assert summary["groups"] == 2.0
    assert summary["mean_hits_at_k"] == 1.0
    assert summary["oracle_mean_hits_at_k"] == 1.5
    assert summary["oracle_capture_rate"] == 2 / 3


def test_expected_calibration_error():
    y_true = np.array([0, 0, 1, 1])
    y_prob = np.array([0.1, 0.2, 0.8, 0.9])
    assert expected_calibration_error(y_true, y_prob, bins=2) < 0.2


def test_group_metrics_ties_do_not_depend_on_source_row_order():
    frame = pd.DataFrame(
        {
            "patient_id": ["p1"] * 4,
            "mutant_peptide": ["AAAAAAAAA", "CCCCCCCCC", "DDDDDDDDD", "EEEEEEEEE"],
            "hla_allele": ["HLA-A*02:01"] * 4,
            "label": ["positive", "negative", "positive", "negative"],
            "constant_score": [1.0] * 4,
        }
    )

    original = group_metrics(
        frame,
        group_col="patient_id",
        score_col="constant_score",
        k=2,
    )
    shuffled = group_metrics(
        frame.sample(frac=1.0, random_state=19),
        group_col="patient_id",
        score_col="constant_score",
        k=2,
    )

    assert original == shuffled
