import pandas as pd

from epicurus_neo.features import add_contrastive_features, anchor_positions, mutation_deltas
from epicurus_neo.portfolio import PortfolioConstraints, select_portfolio


def test_anchor_positions_for_class_i_peptide():
    assert anchor_positions(9) == {1, 8}


def test_mutation_deltas_split_anchor_and_tcr_face():
    deltas = mutation_deltas("SILNFEKLA", "SIINFEKLT")
    assert deltas["mutation_count"] == 2.0
    assert deltas["mutation_anchor_count"] == 1.0
    assert deltas["mutation_tcr_face_count"] == 1.0


def test_add_contrastive_features():
    frame = pd.DataFrame(
        {
            "mutant_peptide": ["SILNFEKLA"],
            "wildtype_peptide": ["SIINFEKLT"],
        }
    )
    out = add_contrastive_features(frame)
    assert out.loc[0, "mutation_count"] == 2.0
    assert "mutation_hydrophobicity_delta" in out.columns


def test_select_portfolio_respects_diversity_limits():
    frame = pd.DataFrame(
        {
            "candidate_id": ["a", "b", "c", "d"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:01", "HLA-B*07:02", "HLA-B*07:02"],
            "gene_symbol": ["KRAS", "KRAS", "TP53", "EGFR"],
            "epicurus_score": [0.9, 0.8, 0.7, 0.6],
        }
    )
    selected = select_portfolio(
        frame,
        constraints=PortfolioConstraints(k=3, max_per_hla=1, max_per_gene=1),
    )
    assert selected["candidate_id"].tolist() == ["a", "c"]

