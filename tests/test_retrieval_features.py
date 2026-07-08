import pandas as pd

from epicurus_neo.retrieval_features import add_retrieval_features, peptide_similarity


def test_peptide_similarity_handles_equal_and_shifted_lengths():
    assert peptide_similarity("AAAA", "AAAT") == 0.75
    assert peptide_similarity("AAAA", "XAAAAX") == 4 / 6


def test_add_retrieval_features_uses_labeled_same_hla_neighbors():
    frame = pd.DataFrame(
        {
            "candidate_id": ["query"],
            "hla_allele": ["HLA-A*02:01"],
            "mutant_peptide": ["AAAA"],
            "label": ["unknown"],
        }
    )
    reference = pd.DataFrame(
        {
            "candidate_id": ["pos", "neg", "other_hla"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:01", "HLA-B*07:02"],
            "mutant_peptide": ["AAAT", "CCCC", "AAAA"],
            "label": ["positive", "negative", "positive"],
        }
    )

    out = add_retrieval_features(frame, reference, top_k=2)

    assert out.loc[0, "retrieval_max_positive_similarity"] == 0.75
    assert out.loc[0, "retrieval_max_negative_similarity"] == 0.0
    assert out.loc[0, "retrieval_topk_positive_fraction"] == 0.5
    assert out.loc[0, "retrieval_reference_count"] == 2.0


def test_add_retrieval_features_excludes_self_candidate():
    frame = pd.DataFrame(
        {
            "candidate_id": ["same"],
            "hla_allele": ["HLA-A*02:01"],
            "mutant_peptide": ["AAAA"],
            "label": ["positive"],
        }
    )

    out = add_retrieval_features(frame, frame, top_k=1)

    assert pd.isna(out.loc[0, "retrieval_max_positive_similarity"])
    assert out.loc[0, "retrieval_reference_count"] == 0.0
