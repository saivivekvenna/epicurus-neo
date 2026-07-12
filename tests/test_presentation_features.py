"""Tests for genuine MHCflurry presentation-feature computation on peptide x HLA frames.

MHCflurry is a real, independent learned class-I presentation predictor (NOT the MixMHCpred backbone
that PRIME is built on, and NOT a neutral impute). We use it to give lossless-recovered candidates a
GENUINE presentation feature set (EL-percentile / affinity / processing) so the four-arm scorer stage
is attributed fairly rather than on a 0.5 imputation.

The heavy model load is guarded with importorskip so environments without the MHCflurry models skip.
"""

from __future__ import annotations

import pandas as pd
import pytest

# IMPORTANT: do NOT import mhcflurry (TensorFlow) at module/collection scope. Loading TensorFlow into the
# shared pytest process collides with the OpenMP runtime used by the xgboost-based tests and segfaults the
# interpreter. Deferring the import into the fixture keeps TensorFlow out of the process until this test
# actually runs (alphabetically after the xgboost tests), so the full suite stays green.
from benchmark.presentation_features import PRESENTATION_COLUMNS  # noqa: E402


@pytest.fixture(scope="module")
def scored():
    pytest.importorskip("mhcflurry")
    from benchmark.presentation_features import add_presentation_features

    frame = pd.DataFrame({
        "mutant_peptide": ["GILGFVFTL", "AAAAAAAAA", "", "NLVPMVATV"],
        "hla_allele": ["HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:01", "HLA-A*02:01"],
    })
    return add_presentation_features(frame)


def test_adds_the_three_presentation_columns(scored):
    for col in PRESENTATION_COLUMNS:
        assert col in scored.columns


def test_valid_peptides_get_numeric_features(scored):
    strong = scored.iloc[0]  # GILGFVFTL is a canonical A*02:01 binder
    assert pd.notna(strong["mhcflurry_el_percentile"])
    assert pd.notna(strong["mhcflurry_affinity_nm"])
    assert pd.notna(strong["mhcflurry_processing"])


def test_empty_peptide_is_nan_not_a_crash(scored):
    empty_row = scored[scored["mutant_peptide"] == ""].iloc[0]
    assert pd.isna(empty_row["mhcflurry_el_percentile"])


def test_known_binder_presents_better_than_random_nonamer(scored):
    # lower EL percentile = stronger presenter
    binder = scored[scored["mutant_peptide"] == "GILGFVFTL"].iloc[0]["mhcflurry_el_percentile"]
    poly_a = scored[scored["mutant_peptide"] == "AAAAAAAAA"].iloc[0]["mhcflurry_el_percentile"]
    assert binder < poly_a


def test_deterministic(scored):
    from benchmark.presentation_features import add_presentation_features

    frame = pd.DataFrame({
        "mutant_peptide": ["GILGFVFTL", "NLVPMVATV"],
        "hla_allele": ["HLA-A*02:01", "HLA-A*02:01"],
    })
    a = add_presentation_features(frame)["mhcflurry_el_percentile"].tolist()
    b = add_presentation_features(frame)["mhcflurry_el_percentile"].tolist()
    assert a == b


def test_input_rows_preserved_in_order(scored):
    assert list(scored["mutant_peptide"]) == ["GILGFVFTL", "AAAAAAAAA", "", "NLVPMVATV"]
