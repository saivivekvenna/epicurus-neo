import pandas as pd

from epicurus_neo.cli import build_parser, cmd_crossfit_retrieval_features, cmd_retrieval_features
from epicurus_neo.retrieval_features import (
    add_crossfit_retrieval_features,
    add_retrieval_features,
    embedding_cosine_similarity,
    peptide_biochemical_similarity,
    peptide_motif_embedding,
    peptide_similarity,
    residue_biochemical_similarity,
)


def test_peptide_similarity_handles_equal_and_shifted_lengths():
    assert peptide_similarity("AAAA", "AAAT") == 0.75
    assert peptide_similarity("AAAA", "XAAAAX") == 4 / 6


def test_biochemical_similarity_rewards_conservative_substitutions():
    assert residue_biochemical_similarity("A", "A") == 1.0
    assert residue_biochemical_similarity("K", "R") > residue_biochemical_similarity("K", "D")
    assert peptide_biochemical_similarity("AAAA", "AAAV") > peptide_similarity("AAAA", "AAAV")


def test_motif_embedding_produces_smooth_neighborhood():
    query = peptide_motif_embedding("AAAA")
    close = peptide_motif_embedding("AAAV")
    far = peptide_motif_embedding("WWWW")

    assert query.shape == close.shape
    assert embedding_cosine_similarity(query, close) > embedding_cosine_similarity(query, far)


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
    assert out.loc[0, "retrieval_biochemical_max_positive_similarity"] > 0.75
    assert out.loc[0, "retrieval_motif_max_positive_similarity"] > out.loc[
        0, "retrieval_motif_max_negative_similarity"
    ]
    assert "retrieval_motif_positive_prototype_similarity" in out.columns
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


def test_add_crossfit_retrieval_features_uses_other_fold_references_only():
    frame = pd.DataFrame(
        {
            "candidate_id": ["q", "same_fold_pos", "other_fold_pos", "other_fold_neg"],
            "hla_allele": ["HLA-A*02:01"] * 4,
            "mutant_peptide": ["AAAA", "AAAA", "AAAT", "CCCC"],
            "label": ["negative", "positive", "positive", "negative"],
            "retrieval_fold": [0, 0, 1, 1],
        }
    )

    out = add_crossfit_retrieval_features(
        frame,
        top_k=2,
        n_folds=2,
        fold_col="retrieval_fold",
    )

    query = out[out["candidate_id"] == "q"].iloc[0]
    assert query["retrieval_max_positive_similarity"] == 0.75
    assert query["retrieval_reference_count"] == 2.0


def test_retrieval_features_cli_writes_output(tmp_path):
    query = tmp_path / "query.csv"
    reference = tmp_path / "reference.csv"
    output = tmp_path / "out.csv"
    pd.DataFrame(
        {
            "candidate_id": ["query"],
            "hla_allele": ["HLA-A*02:01"],
            "mutant_peptide": ["AAAA"],
            "label": ["unknown"],
        }
    ).to_csv(query, index=False)
    pd.DataFrame(
        {
            "candidate_id": ["pos"],
            "hla_allele": ["HLA-A*02:01"],
            "mutant_peptide": ["AAAT"],
            "label": ["positive"],
        }
    ).to_csv(reference, index=False)

    parser = build_parser()
    args = parser.parse_args(
        [
            "add-retrieval-features",
            "--input",
            str(query),
            "--reference",
            str(reference),
            "--output",
            str(output),
        ]
    )
    assert cmd_retrieval_features(args) == 0
    assert "retrieval_max_positive_similarity" in pd.read_csv(output).columns


def test_crossfit_retrieval_features_cli_writes_output(tmp_path):
    input_path = tmp_path / "input.csv"
    output = tmp_path / "out.csv"
    pd.DataFrame(
        {
            "candidate_id": ["q", "pos"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:01"],
            "mutant_peptide": ["AAAA", "AAAT"],
            "label": ["negative", "positive"],
            "retrieval_fold": [0, 1],
        }
    ).to_csv(input_path, index=False)

    parser = build_parser()
    args = parser.parse_args(
        [
            "add-crossfit-retrieval-features",
            "--input",
            str(input_path),
            "--output",
            str(output),
            "--n-folds",
            "2",
        ]
    )
    assert cmd_crossfit_retrieval_features(args) == 0
    assert "retrieval_max_positive_similarity" in pd.read_csv(output).columns
