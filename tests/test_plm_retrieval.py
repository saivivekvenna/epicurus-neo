import numpy as np
import pandas as pd

from epicurus_neo.cli import build_parser
from epicurus_neo.plm_retrieval import add_embedding_retrieval_features


def test_embedding_retrieval_features_use_supplied_vectors():
    query = pd.DataFrame(
        {
            "candidate_id": ["q"],
            "hla_allele": ["HLA-A*02:01"],
            "mutant_peptide": ["AAAA"],
            "label": ["unknown"],
        }
    )
    reference = pd.DataFrame(
        {
            "candidate_id": ["pos", "neg"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:01"],
            "mutant_peptide": ["AAAV", "WWWW"],
            "label": ["positive", "negative"],
        }
    )
    embeddings = {
        "AAAA": np.array([1.0, 0.0]),
        "AAAV": np.array([0.9, 0.1]),
        "WWWW": np.array([0.0, 1.0]),
    }

    out = add_embedding_retrieval_features(query, reference, embeddings)

    assert out.loc[0, "retrieval_plm_max_positive_similarity"] > out.loc[
        0, "retrieval_plm_max_negative_similarity"
    ]
    assert out.loc[0, "retrieval_plm_reference_count"] == 2.0


def test_plm_retrieval_cli_parses_options():
    parser = build_parser()
    args = parser.parse_args(
        [
            "add-plm-retrieval-features",
            "--input",
            "query.csv",
            "--reference",
            "reference.csv",
            "--output",
            "out.csv",
            "--model-name",
            "facebook/esm2_t6_8M_UR50D",
            "--device",
            "cpu",
        ]
    )

    assert args.command == "add-plm-retrieval-features"
    assert args.device == "cpu"
