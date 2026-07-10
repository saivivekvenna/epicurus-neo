import numpy as np

from epicurus_neo.m6.evaluate import macro_paired_delta


def test_macro_delta_equal_weights_studies():
    # Study A: delta +1 per patient (2 patients). Study B: delta 0 (4 patients).
    per_study = {
        "A": (np.array([2.0, 2.0]), np.array([1.0, 1.0])),
        "B": (np.array([1.0, 1.0, 1.0, 1.0]), np.array([1.0, 1.0, 1.0, 1.0])),
    }
    result = macro_paired_delta(per_study, seed=17)
    # Macro = mean of per-study mean deltas = mean(+1, 0) = 0.5 (not patient-weighted 1/3).
    assert abs(result["delta"] - 0.5) < 1e-9
    assert result["delta_ci"][0] <= result["delta"] <= result["delta_ci"][1]
