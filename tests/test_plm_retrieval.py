import numpy as np
import pandas as pd

from epicurus_neo.cli import build_parser
from epicurus_neo.plm_retrieval import (
    add_embedding_retrieval_features,
    load_plm_embedding_cache,
    save_plm_embedding_cache,
)


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


def test_plm_embedding_cache_round_trip(tmp_path):
    output = tmp_path / "embeddings.npz"
    embeddings = {
        "AAAA": np.array([1.0, 0.0], dtype=np.float32),
        "WWWW": np.array([0.0, 1.0], dtype=np.float32),
    }

    save_plm_embedding_cache(
        embeddings,
        output,
        model_name="test/model",
    )
    loaded, model_name = load_plm_embedding_cache(output)

    assert model_name == "test/model"
    assert set(loaded) == {"AAAA", "WWWW"}
    np.testing.assert_allclose(loaded["AAAA"], embeddings["AAAA"])


def test_plm_embedding_cache_cli_parses_multiple_inputs():
    parser = build_parser()
    args = parser.parse_args(
        [
            "build-plm-embedding-cache",
            "--input",
            "train.csv",
            "--input",
            "validation.csv",
            "--output",
            "embeddings.npz",
        ]
    )

    assert args.input == ["train.csv", "validation.csv"]
