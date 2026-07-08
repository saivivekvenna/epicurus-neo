from pathlib import Path

import pandas as pd

from epicurus_neo.download import dataset_file_plans
from epicurus_neo.normalize import (
    normalize_gartner_table,
    normalize_neoranking_neopep,
    write_normalized,
)


def test_download_plan_includes_neoranking_files():
    plans = dataset_file_plans(
        "configs/datasets.yml",
        output_dir="data/raw",
        dataset_key="neoranking",
    )
    file_keys = {plan.file_key for plan in plans}
    assert "neopep_data_org" in file_keys
    assert "gartner_nmers_ranking" in file_keys


def test_normalize_neoranking_neopep_fixture(tmp_path: Path):
    source = tmp_path / "Neopep_data_org.txt"
    pd.DataFrame(
        {
            "patient": ["p1", "p1", "p2"],
            "dataset": ["NCI", "NCI", "TESLA"],
            "train_test": ["train", "train", "test"],
            "response_type": ["CD8", "negative", "not_tested"],
            "gene": ["KRAS", "TP53", "EGFR"],
            "mutant_seq": ["SILNFEKLA", "AAAAAAAAL", "BBBBBBBBL"],
            "wt_seq": ["SIINFEKLT", "AAAAAAAAT", "BBBBBBBBT"],
            "rnaseq_TPM": [50.0, 2.0, 3.0],
            "CCF": [0.9, 0.2, 0.1],
            "mutant_rank_netMHCpan": [0.1, 2.0, 4.0],
            "DAI_NetMHC": [5.0, 1.0, 0.5],
            "bestWTMatchScore_I": [0.1, 0.9, 0.8],
        }
    ).to_csv(source, sep="\t", index=False)

    normalized = normalize_neoranking_neopep(source)

    assert normalized["label"].tolist() == ["positive", "negative", "unknown"]
    assert normalized.loc[0, "expression_tpm"] == 50.0
    assert normalized.loc[0, "clonality_ccf"] == 0.9
    assert normalized.loc[0, "netmhcpan_mutant_rank"] == 0.1
    assert normalized.loc[0, "mutation_tcr_face_count"] == 1.0
    assert normalized.loc[0, "mutation_anchor_count"] == 1.0


def test_normalize_gartner_fixture(tmp_path: Path):
    source = tmp_path / "Gartner_nmers_ranking.txt"
    pd.DataFrame(
        {
            "Patient": ["3703", "3703"],
            "Cancer Type": ["Melanoma", "Melanoma"],
            "Gene Name": ["CCNE1", "NSDHL"],
            "Mutant Nmer": ["HEVLLPQYPQQILIQIAELLDLCVL", "FLSRILTGLNYEVPKYHIPYWVAYY"],
            "Rank NetMHC": [4, 24],
            "Rank Nmer Model": [2, 43],
        }
    ).to_csv(source, index=False)

    normalized = normalize_gartner_table(source)

    assert normalized["source_dataset"].unique().tolist() == ["gartner_nci"]
    assert normalized["patient_id"].tolist() == ["3703", "3703"]
    assert normalized["label"].tolist() == ["unknown", "unknown"]
    assert normalized.loc[0, "gene_symbol"] == "CCNE1"
    assert normalized.loc[0, "baseline_netmhc_rank"] == 4


def test_write_normalized(tmp_path: Path):
    out = tmp_path / "processed" / "table.csv"
    path = write_normalized(pd.DataFrame({"a": [1]}), out)
    assert path == out
    assert out.exists()

