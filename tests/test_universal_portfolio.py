import numpy as np
import pandas as pd

from benchmark.universal_portfolio import (
    mutation_representatives,
    select_evidence_lane_portfolio,
    select_rank_fusion_cap1,
)


def _frame():
    return pd.DataFrame(
        {
            "mutation_id": ["m1", "m1", "m2", "m3", "m4"],
            "mutant_peptide": ["AAAAAAAA", "AAAAAAAK", "CCCCCCCC", "DDDDDDDD", "EEEEEEEE"],
            "hla_allele": ["HLA-A*01:01"] * 5,
            "prime": [0.9, 0.8, 0.7, 0.1, np.nan],
            "epicurus": [0.2, 0.1, 0.4, 0.9, 0.8],
            "presentation": [0.5, 0.4, 0.3, 0.2, 0.99],
        }
    )


def test_representatives_select_one_best_route_per_mutation():
    got = mutation_representatives(_frame(), "prime")
    assert got["mutation_id"].tolist() == ["m1", "m2", "m3"]
    assert got.iloc[0]["mutant_peptide"] == "AAAAAAAA"


def test_rank_fusion_is_unique_and_missing_evidence_has_no_advantage():
    got = select_rank_fusion_cap1(_frame(), ("prime", "epicurus"), k=4)
    assert got["mutation_id"].is_unique
    assert len(got) == 4
    assert got.iloc[-1]["mutation_id"] == "m4"


def test_lane_portfolio_exposes_distinct_experts_and_is_deterministic():
    a = select_evidence_lane_portfolio(_frame(), ("prime", "epicurus", "presentation"), k=4)
    b = select_evidence_lane_portfolio(_frame(), ("prime", "epicurus", "presentation"), k=4)
    assert a[["mutation_id", "portfolio_lane"]].equals(b[["mutation_id", "portfolio_lane"]])
    assert a["mutation_id"].tolist()[:3] == ["m1", "m3", "m4"]
    assert a["mutation_id"].is_unique


def test_lane_portfolio_skips_missing_lane_column():
    got = select_evidence_lane_portfolio(_frame(), ("absent", "prime"), k=2)
    assert got["mutation_id"].tolist() == ["m1", "m2"]
