import pandas as pd
import pytest

from epicurus_neo.plm_finetune import (
    fit_numeric_stats,
    make_hard_pairs,
    numeric_matrix,
)


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "hla_allele": ["A0201"] * 5 + ["B0702"] * 3,
            "label": [
                "positive",
                "negative",
                "negative",
                "negative",
                "negative",
                "positive",
                "negative",
                "negative",
            ],
            "mhcflurry_affinity": [10, 20, 30, 40, 50, 15, 25, 35],
            "mhcflurry_processing_score": [0.5] * 8,
            "mhcflurry_presentation_score": [0.1, 0.9, 0.8, 0.7, 0.2, 0.1, 0.9, 0.2],
            "mhcflurry_presentation_percentile": [1, 2, 3, 4, 5, 1, 2, 3],
        }
    )


def test_make_hard_pairs_samples_high_presentation_negatives_by_group():
    frame = _frame()

    positives, negatives = make_hard_pairs(
        frame,
        group_col="hla_allele",
        negatives_per_positive=2,
        hard_negative_pool=2,
        random_state=17,
    )

    assert len(positives) == 4
    assert set(positives) == {0, 5}
    assert set(negatives).issubset({1, 2, 6, 7})
    assert all(
        frame.loc[p, "hla_allele"] == frame.loc[n, "hla_allele"]
        for p, n in zip(positives, negatives, strict=True)
    )


def test_numeric_features_are_imputed_and_standardized():
    frame = _frame()
    frame.loc[0, "mhcflurry_affinity"] = None
    means, stds = fit_numeric_stats(frame)
    matrix = numeric_matrix(frame, means, stds)

    assert matrix.shape == (8, 4)
    assert not pd.isna(matrix).any()


def test_make_hard_pairs_requires_both_classes():
    frame = _frame()
    frame["label"] = "negative"

    with pytest.raises(ValueError, match="No within-group"):
        make_hard_pairs(
            frame,
            group_col="hla_allele",
            negatives_per_positive=1,
            hard_negative_pool=10,
            random_state=17,
        )
