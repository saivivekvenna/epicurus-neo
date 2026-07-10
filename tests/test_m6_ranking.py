import numpy as np
import pandas as pd

from epicurus_neo.m6.ranking import classification_metrics, patient_rank_vectors


def _frame():
    # Patient P1: 3 candidates, k=2. Patient P2: perfect ranking.
    return pd.DataFrame(
        {
            "patient_id": ["P1", "P1", "P1", "P2", "P2", "P2"],
            "mutant_peptide": list("ABCDEF"),
            "hla_allele": ["HLA-A*02:01"] * 6,
            "label": [1, 0, 0, 1, 1, 0],
            "score": [0.9, 0.1, 0.2, 0.9, 0.8, 0.1],
        }
    )


def test_k_patient_uses_min_of_cap_and_length():
    vectors = patient_rank_vectors(_frame(), "score", k_cap=2)
    # P1 has 3 candidates -> k=2; top-2 by score = A(1), C(0) -> 1 hit, precision 1/2.
    # P2 has 3 candidates -> k=2; top-2 = D(1), E(1) -> 2 hits, precision 2/2.
    assert vectors["hits_at_k"].tolist() == [1.0, 2.0]
    assert vectors["precision_at_k"].tolist() == [0.5, 1.0]
    assert vectors["p_at_least_2"].tolist() == [0.0, 1.0]


def test_classification_metrics_on_separable_scores():
    metrics = classification_metrics(np.array([0, 0, 1, 1]), np.array([0.1, 0.2, 0.8, 0.9]))
    assert metrics["auroc"] == 1.0
    assert 0.0 <= metrics["brier"] <= 0.25
    assert len(metrics["calibration"]) >= 1
