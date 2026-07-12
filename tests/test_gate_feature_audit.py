"""Tests for the isolated gate-feature-audit helpers.

These cover the pure computational core: identifying the high-presentation
stratum (the candidates a presentation gate would keep / that outrank
positives), and quantifying whether an ORTHOGONAL feature can separate the
top-ranked TESTED_NEGATIVE decoys from POSITIVES on that stratum. This is the
exact discrimination the dynamic gate could not achieve with presentation
features alone.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from benchmark.gate_feature_audit import (
    POSITIVE,
    TESTED_NEGATIVE,
    UNTESTED,
    conditional_auroc,
    feature_coverage,
    grouped_oof_auroc,
    high_presentation_mask,
    within_patient_variation,
)


def test_high_presentation_mask_top_fraction_within_patient():
    frame = pd.DataFrame(
        {
            "patient_id": ["A", "A", "A", "A", "B", "B"],
            "score": [1.0, 2.0, 3.0, 4.0, 10.0, 20.0],
        }
    )
    mask = high_presentation_mask(frame, "score", higher_better=True, quantile=0.5)
    # Patient A top-2 by score -> rows with 3,4 ; patient B top-1 -> row with 20
    assert list(mask) == [False, False, True, True, False, True]


def test_high_presentation_mask_top_k_per_patient():
    frame = pd.DataFrame(
        {
            "patient_id": ["A"] * 5 + ["B"] * 3,
            "score": [5.0, 4.0, 3.0, 2.0, 1.0, 9.0, 8.0, 7.0],
        }
    )
    mask = high_presentation_mask(frame, "score", higher_better=True, top_k=2)
    # top-2 per patient by score
    assert list(mask) == [True, True, False, False, False, True, True, False]


def test_high_presentation_mask_lower_is_better_uses_ascending():
    frame = pd.DataFrame(
        {
            "patient_id": ["A", "A", "A", "A"],
            "rank": [0.1, 0.2, 5.0, 9.0],  # percentile rank: lower = better presenter
        }
    )
    mask = high_presentation_mask(frame, "rank", higher_better=False, quantile=0.5)
    # Best presenters are the two smallest ranks
    assert list(mask) == [True, True, False, False]


def test_conditional_auroc_orientation_and_untested_drop():
    frame = pd.DataFrame(
        {
            "label": [POSITIVE, POSITIVE, TESTED_NEGATIVE, TESTED_NEGATIVE, UNTESTED],
            "feat": [0.9, 0.8, 0.2, 0.1, 0.5],
        }
    )
    hi = conditional_auroc(frame, "feat", higher_better=True, min_per_class=2)
    assert hi["auroc"] == 1.0
    assert hi["n_pos"] == 2 and hi["n_neg"] == 2  # UNTESTED never a negative
    lo = conditional_auroc(frame, "feat", higher_better=False, min_per_class=2)
    assert lo["auroc"] == 0.0


def test_conditional_auroc_nan_excluded_from_scoring_but_counted_in_coverage():
    frame = pd.DataFrame(
        {
            "label": [POSITIVE, POSITIVE, TESTED_NEGATIVE, TESTED_NEGATIVE],
            "feat": [0.9, np.nan, 0.1, 0.2],
        }
    )
    out = conditional_auroc(frame, "feat", higher_better=True, min_per_class=1)
    assert out["coverage"] == 0.75  # 3/4 non-null among POS/TN rows
    assert out["n_pos"] == 1 and out["n_neg"] == 2  # NaN positive dropped from AUROC


def test_conditional_auroc_degenerate_returns_none():
    frame = pd.DataFrame({"label": [POSITIVE, POSITIVE], "feat": [0.5, 0.6]})
    out = conditional_auroc(frame, "feat", higher_better=True)
    assert out["auroc"] is None


def test_grouped_oof_auroc_recovers_real_signal_and_rejects_noise():
    rng = np.random.default_rng(0)
    n = 240
    groups = np.repeat(np.arange(12), n // 12)
    y = rng.integers(0, 2, size=n)
    signal = y + rng.normal(0, 0.3, size=n)  # feature carries the label
    noise = rng.normal(0, 1, size=n)
    frame = pd.DataFrame(
        {
            "patient_id": groups,
            "label": np.where(y == 1, POSITIVE, TESTED_NEGATIVE),
            "signal": signal,
            "noise": noise,
        }
    )
    real = grouped_oof_auroc(frame, ["signal"], n_splits=4, seed=0)
    rand = grouped_oof_auroc(frame, ["noise"], n_splits=4, seed=0)
    assert real["oof_auroc"] > 0.85
    assert 0.35 < rand["oof_auroc"] < 0.65


def test_within_patient_variation_separates_candidate_varying_from_constant():
    frame = pd.DataFrame(
        {
            "patient_id": ["A", "A", "A", "B", "B"],
            "candidate_varying": [1.0, 2.0, 3.0, 9.0, 8.0],  # varies in both patients
            "patient_constant": [5.0, 5.0, 5.0, 7.0, 7.0],  # constant within each patient
        }
    )
    # candidate_varying varies within both A and B -> fraction 1.0
    assert within_patient_variation(frame, "candidate_varying") == 1.0
    # patient_constant is constant within each patient -> fraction 0.0
    assert within_patient_variation(frame, "patient_constant") == 0.0
    # absent column -> 0.0, no crash
    assert within_patient_variation(frame, "missing") == 0.0


def test_feature_coverage_reports_nonnull_fraction():
    frame = pd.DataFrame(
        {"a": [1.0, 2.0, np.nan, 4.0], "b": [np.nan, np.nan, np.nan, np.nan]}
    )
    cov = feature_coverage(frame, ["a", "b", "missing_col"])
    assert cov["a"] == 0.75
    assert cov["b"] == 0.0
    assert cov["missing_col"] == 0.0  # absent column -> zero coverage, no crash
