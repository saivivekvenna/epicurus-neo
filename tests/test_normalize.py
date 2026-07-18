from pathlib import Path
from zipfile import ZipFile

import pandas as pd

from epicurus_neo.normalize import (
    normalize_bigmhc_table,
    normalize_cd8_multimer_2025,
    normalize_gartner_table,
    normalize_improve_cv,
    normalize_neoranking_neopep,
    normalize_tesla_table,
    write_normalized,
)


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
            "ID": ["3703", "3703", "3703", "3703"],
            "tumor type": ["MELANOMA", "MELANOMA", "MELANOMA", "MELANOMA"],
            "Gene Name": ["NSDHL", "FADS3", "BRAF", "TP53"],
            "Wt Epitope": [
                "FLSRILTGLNYEAPKYHIPYWVAYY",
                "RHNYSRVAPLVKLLCAKHGLSYEVK",
                "AAAAAAAAL",
                "CCCCCCCCL",
            ],
            "Mut Epitope": [
                "FLSRILTGLNYEVPKYHIPYWVAYY",
                "RHNYSRVAPLVKLLCAKHGLSYEVK",
                "AAAAAAAAV",
                "CCCCCCCCV",
            ],
            "Screening Status": ["1", "0", "-", "unscreened"],
            "Gene Expression Decile for this sample(1=lowest expression-10=highest expression)": [
                6.0,
                4.0,
                2.0,
                1.0,
            ],
        }
    ).to_csv(source, sep="\t", index=False)

    normalized = normalize_gartner_table(source)

    assert normalized["patient_id"].tolist() == ["3703", "3703", "3703", "3703"]
    assert normalized["label"].tolist() == ["positive", "negative", "negative", "unknown"]
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


def test_normalize_cd8_multimer_2025_repairs_dimensions_and_deduplicates(tmp_path: Path):
    source = tmp_path / "mmc2.xlsx"
    pd.DataFrame(
        {
            "Patient ID": ["TIL#01", "TIL#01", "PBMC#01"],
            "Alt. ID": ["A", "A", "B"],
            "WT epitope": ["SIINFEKLT", "SIINFEKLT", "AAAAAAAAA"],
            "MUT epitope": ["SIINFEKLA", "SIINFEKLA", "AAAAAAAAV"],
            "HLA": ["HLA-A*02:01", "HLA-A*02:01", "HLA-B*07:02"],
            "Response": ["YES", "YES", "NO"],
            "Resp. Magnitude": [0.4, 0.4, None],
            "Binding affinity (%Rank score)": [0.2, 0.2, 5.0],
            "RNA expression (TPM)": [10.0, 10.0, 0.0],
            "Proteasomal processing score": [0.8, 0.8, 0.3],
            "EL (%Rank score)": [0.1, 0.1, 10.0],
            "RF classifier score": [0.9, 0.9, 0.1],
            "Agretopicity": [0.01, 0.01, 1.0],
            "Foreignness score": [0.7, 0.7, 0.0],
            "Dissimilarity": [1.0, 1.0, 0.0],
            "Dataset": ["TIL", "TIL", "PBMC"],
            "Tumor type": ["Melanoma_TIL", "Melanoma_TIL", "TNBC"],
            "SYMBOL": ["BRAF", "BRAF", "TP53"],
            "ENSG": ["ENSG1", "ENSG1", "ENSG2"],
            "Genome assembly": ["grch38", "grch38", "grch38"],
            "TMB": [5.0, 5.0, 2.0],
        }
    ).to_excel(source, index=False)

    normalized = normalize_cd8_multimer_2025(source)

    assert len(normalized) == 2
    assert normalized["label"].tolist() == ["positive", "negative"]
    assert normalized.loc[0, "patient_id"] == "cd8_multimer_2025:TIL#01"
    assert normalized.loc[0, "source_duplicate_count"] == "2"
    assert normalized.loc[0, "netmhcpan_binding_score"] == 0.998
    assert normalized.loc[0, "netmhcpan_el_score"] == 0.999
    assert normalized.loc[0, "agretopicity_score"] == 2.0
    assert "Resp. Magnitude" not in normalized.columns


def test_normalize_improve_cv_fixture(tmp_path: Path):
    source = tmp_path / "improve.tsv"
    pd.DataFrame(
        {
            "Patient": ["BC-1", "RH-1"],
            "HLA_allele": ["HLA-A02:01", "HLA-B07:02"],
            "Norm_peptide": ["SIINFEKLT", "AAAAAAAAA"],
            "Mut_peptide": ["SIINFEKLA", "AAAAAAAAV"],
            "response": [1, 0],
            "cohort": ["bladder", "melanoma"],
            "Partition": [0, 1],
            "Gene_Symbol": ["BRAF", "TP53"],
            "RankEL_4.1": [0.2, 10.0],
            "RankBA_4.1": [0.4, 20.0],
            "RankEL_wt_4.1": [5.0, 2.0],
            "Expression": [10.0, 0.0],
            "Stability": [2.0, 0.5],
            "Foreigness": [0.8, 0.0],
        }
    ).to_csv(source, sep="\t", index=False)

    normalized = normalize_improve_cv(source)

    assert normalized["label"].tolist() == ["positive", "negative"]
    assert normalized["patient_id"].tolist() == ["improve:BC-1", "improve:RH-1"]
    assert normalized["official_partition"].tolist() == ["0", "1"]
    assert normalized.loc[0, "netmhcpan_el_score"] == 0.998
    assert normalized.loc[1, "netmhcpan_binding_score"] == 0.8
