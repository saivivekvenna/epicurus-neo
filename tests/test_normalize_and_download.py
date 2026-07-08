from pathlib import Path
from zipfile import ZipFile

import pandas as pd

from epicurus_neo.download import dataset_file_plans
from epicurus_neo.normalize import (
    normalize_bigmhc_table,
    normalize_gartner_table,
    normalize_neoranking_neopep,
    normalize_tesla_table,
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


def test_download_plan_includes_bigmhc_direct_files():
    plans = dataset_file_plans(
        "configs/datasets.yml",
        output_dir="data/raw",
        dataset_key="bigmhc",
    )
    file_keys = {plan.file_key for plan in plans}
    assert "datasets_zip" in file_keys
    assert "manafest" in file_keys


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


def test_normalize_gartner_tsv_fixture(tmp_path: Path):
    source = tmp_path / "NmersTestingSet.txt"
    pd.DataFrame(
        {
            "ID": ["3703", "3703"],
            "tumor type": ["MELANOMA", "MELANOMA"],
            "Gene Name": ["NSDHL", "FADS3"],
            "Wt Epitope": ["FLSRILTGLNYEAPKYHIPYWVAYY", "RHNYSRVAPLVKLLCAKHGLSYEVK"],
            "Mut Epitope": ["FLSRILTGLNYEVPKYHIPYWVAYY", "RHNYSRVAPLVKLLCAKHGLSYEVK"],
            "Screening Status": ["1", "unscreened"],
            "Gene Expression Decile for this sample(1=lowest expression-10=highest expression)": [6.0, 4.0],
        }
    ).to_csv(source, sep="\t", index=False)

    normalized = normalize_gartner_table(source)

    assert normalized["patient_id"].tolist() == ["3703", "3703"]
    assert normalized["label"].tolist() == ["positive", "unknown"]
    assert normalized.loc[0, "mutant_peptide"] == "FLSRILTGLNYEVPKYHIPYWVAYY"


def test_write_normalized(tmp_path: Path):
    out = tmp_path / "processed" / "table.csv"
    path = write_normalized(pd.DataFrame({"a": [1]}), out)
    assert path == out
    assert out.exists()


def test_normalize_tesla_fixture(tmp_path: Path):
    source = tmp_path / "TESLA_neoepitopes.xlsx"
    pd.DataFrame(
        {
            "peptide": ["NILGFTFDI", "ABCDEF"],
            "target_value": [0, 1],
            "allele": ["HLA-A02:01", "HLA-B44:02"],
        }
    ).to_excel(source, index=False)

    normalized = normalize_tesla_table(source)
    assert normalized["label"].tolist() == ["negative", "positive"]
    assert normalized.loc[0, "hla_allele_norm"] == "HLA-A*02:01"


def test_normalize_bigmhc_fixture(tmp_path: Path):
    source = tmp_path / "im_test.csv"
    pd.DataFrame(
        {
            "mhc": ["HLA-A*02:01", "HLA-B*07:02"],
            "pep": ["ALDKLSSQHLY", "ASIRNANLY"],
            "tgt": [1, 0],
            "BigMHC_EL": [0.9, 0.2],
            "BigMHC_IM": [0.7, 0.1],
            "BigMHC_ELIM": [0.8, 0.15],
        }
    ).to_csv(source, index=False)

    normalized = normalize_bigmhc_table(source)
    assert normalized["label"].tolist() == ["positive", "negative"]
    assert normalized.loc[0, "patient_id"] == "bigmhc_HLA-A*02:01"
    assert normalized.loc[0, "bigmhc_el_score"] == 0.9
    assert normalized.loc[0, "bigmhc_im_score"] == 0.7


def test_normalize_bigmhc_zip_member_fixture(tmp_path: Path):
    table = tmp_path / "im_test.csv"
    archive = tmp_path / "datasets.zip"
    pd.DataFrame(
        {
            "mhc": ["HLA-A*02:01"],
            "pep": ["ALDKLSSQHLY"],
            "tgt": [1],
        }
    ).to_csv(table, index=False)
    with ZipFile(archive, "w") as zip_file:
        zip_file.write(table, arcname="im_test.csv")

    normalized = normalize_bigmhc_table(archive, zip_member="im_test.csv")
    assert normalized.loc[0, "candidate_id"] == "bigmhc:im_test:0"


def test_normalize_reads_zipped_single_table(tmp_path: Path):
    table = tmp_path / "Neopep_data_org.txt"
    archive = tmp_path / "Neopep_data_org.txt.zip"
    pd.DataFrame(
        {
            "patient": ["p1"],
            "dataset": ["NCI"],
            "response_type": ["CD8"],
            "mutant_seq": ["SILNFEKLA"],
            "wt_seq": ["SIINFEKLT"],
        }
    ).to_csv(table, sep="\t", index=False)
    with ZipFile(archive, "w") as zip_file:
        zip_file.write(table, arcname="Neopep_data_org.txt")

    normalized = normalize_neoranking_neopep(archive)
    assert normalized.loc[0, "label"] == "positive"
