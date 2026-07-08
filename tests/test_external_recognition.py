import numpy as np
import pandas as pd

from epicurus_neo.cli import build_parser
from epicurus_neo.external_recognition import (
    add_external_recognition_features,
    normalize_vdjdb,
)


def test_normalize_vdjdb_filters_and_aggregates_supported_class_i_epitopes():
    raw = pd.DataFrame(
        {
            "species": ["HomoSapiens", "HomoSapiens", "MusMusculus", "HomoSapiens"],
            "antigen.epitope": ["NLVPMVATV", "NLVPMVATV", "NLVPMVATV", "TOO_SHORT"],
            "antigen.species": ["CMV", "CMV", "CMV", "CMV"],
            "mhc.a": ["HLA-A*02:01"] * 4,
            "mhc.class": ["MHCI", "MHCI", "MHCI", "MHCII"],
            "vdjdb.score": [1, 3, 3, 3],
        }
    )

    out = normalize_vdjdb(raw)

    assert len(out) == 1
    assert out.loc[0, "mutant_peptide"] == "NLVPMVATV"
    assert out.loc[0, "recognition_support"] == 2
    assert out.loc[0, "recognition_max_evidence"] == 3
    assert out.loc[0, "recognition_origin"] == "pathogen"


def test_external_recognition_features_separate_hla_and_origin_neighbors():
    query = pd.DataFrame(
        {
            "mutant_peptide": ["NLVPMVATV"],
            "hla_allele": ["HLA-A*02:01"],
        }
    )
    reference = pd.DataFrame(
        {
            "mutant_peptide": ["NLVPMVATV", "AAAAAAAAL", "WWWWWWWWW"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*01:01", "HLA-A*02:01"],
            "recognition_origin": ["pathogen", "pathogen", "human"],
            "recognition_support": [10, 2, 1],
            "recognition_max_evidence": [3, 1, 1],
        }
    )
    embeddings = {
        "NLVPMVATV": np.array([1.0, 0.0], dtype=np.float32),
        "AAAAAAAAL": np.array([0.7, 0.3], dtype=np.float32),
        "WWWWWWWWW": np.array([0.0, 1.0], dtype=np.float32),
    }

    out = add_external_recognition_features(query, reference, embeddings)

    assert out.loc[0, "recognition_hla_max_similarity"] == 1.0
    assert out.loc[0, "recognition_pathogen_max_similarity"] == 1.0
    assert out.loc[0, "recognition_human_max_similarity"] == 0.0
    assert out.loc[0, "recognition_pathogen_minus_human_similarity"] == 1.0


def test_external_recognition_cli_parses_commands():
    parser = build_parser()
    normalize_args = parser.parse_args(
        ["normalize-vdjdb", "--input", "vdjdb.txt", "--output", "vdjdb.csv"]
    )
    feature_args = parser.parse_args(
        [
            "add-external-recognition-features",
            "--input",
            "query.csv",
            "--reference",
            "vdjdb.csv",
            "--embedding-cache",
            "embeddings.npz",
            "--output",
            "features.csv",
        ]
    )

    assert normalize_args.command == "normalize-vdjdb"
    assert feature_args.top_k == 5
