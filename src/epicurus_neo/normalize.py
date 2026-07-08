from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from epicurus_neo.features import add_contrastive_features
from epicurus_neo.schema import add_normalized_columns, validate_schema


PEPTIDE_COLS = [
    "mutant_seq",
    "mut_seq",
    "mutant_peptide",
    "Mutant Minimal Peptide",
    "Mutant Mmp",
    "Mutated Minimal Peptide",
    "Mutant Nmer",
    "mutant_nmer",
]

WILDTYPE_COLS = [
    "wt_seq",
    "wildtype_seq",
    "wildtype_peptide",
    "Wildtype Minimal Peptide",
    "Wildtype Mmp",
    "Wildtype Nmer",
    "wt_nmer",
]

HLA_COLS = [
    "hla_allele",
    "HLA",
    "Allele",
    "HLA Allele",
    "Best Allele",
    "best_allele",
    "mutant_best_allele",
]

PATIENT_COLS = ["patient", "Patient", "patient_id", "Sample", "sample"]
STUDY_COLS = ["dataset", "Dataset", "study_id", "Study"]
GENE_COLS = ["gene", "Gene Name", "Gene", "gene_symbol"]


def read_table(path: str | Path) -> pd.DataFrame:
    table_path = Path(path)
    suffix = table_path.suffix.lower()
    if suffix == ".csv":
        return pd.read_csv(table_path)
    if suffix in {".tsv", ".txt"}:
        return pd.read_csv(table_path, sep=None, engine="python")
    raise ValueError(f"Unsupported input table format: {table_path}")


def _first_existing(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    for column in candidates:
        if column in frame.columns:
            return column
    return None


def _series_or_default(frame: pd.DataFrame, candidates: list[str], default: Any = "") -> pd.Series:
    column = _first_existing(frame, candidates)
    if column is None:
        return pd.Series([default] * len(frame), index=frame.index)
    return frame[column]


def _label_from_response(value: Any) -> str:
    if value is None or pd.isna(value):
        return "unknown"
    normalized = str(value).strip().lower()
    if normalized in {"cd8", "positive", "pos", "reactive", "immunogenic", "1", "true", "yes"}:
        return "positive"
    if normalized in {"negative", "neg", "non-reactive", "nonreactive", "0", "false", "no"}:
        return "negative"
    if normalized in {"not_tested", "not tested", "unknown", "untested", "nan", ""}:
        return "unknown"
    return "unknown"


def _label_column(frame: pd.DataFrame) -> pd.Series:
    candidates = [
        "response_type",
        "Response Type",
        "response",
        "Response",
        "immunogenicity",
        "Immunogenicity",
        "reactivity",
        "Reactivity",
        "TIL Reactivity",
    ]
    source = _series_or_default(frame, candidates, "unknown")
    return source.map(_label_from_response)


def _candidate_ids(source_dataset: str, frame: pd.DataFrame) -> pd.Series:
    if "candidate_id" in frame.columns:
        return frame["candidate_id"].astype(str)
    return pd.Series(
        [f"{source_dataset}:{idx}" for idx in range(len(frame))],
        index=frame.index,
    )


def _copy_numeric_features(source: pd.DataFrame, target: pd.DataFrame) -> pd.DataFrame:
    out = target.copy()
    for column in source.columns:
        if column in out.columns:
            continue
        numeric = pd.to_numeric(source[column], errors="coerce")
        if numeric.notna().any():
            out[column] = numeric
    return out


def normalize_candidate_table(
    frame: pd.DataFrame,
    *,
    source_dataset: str,
    study_default: str | None = None,
) -> pd.DataFrame:
    labels = _label_column(frame)
    study = _series_or_default(frame, STUDY_COLS, study_default or source_dataset)
    patient = _series_or_default(frame, PATIENT_COLS, "")

    out = pd.DataFrame(
        {
            "candidate_id": _candidate_ids(source_dataset, frame),
            "source_dataset": source_dataset,
            "study_id": study.astype(str),
            "patient_id": patient.astype(str),
            "hla_allele": _series_or_default(frame, HLA_COLS, "").astype(str),
            "mutant_peptide": _series_or_default(frame, PEPTIDE_COLS, "").astype(str),
            "wildtype_peptide": _series_or_default(frame, WILDTYPE_COLS, "").astype(str),
            "label": labels,
            "label_weight": np.where(labels == "unknown", 0.0, 1.0),
            "assay_type": _series_or_default(frame, ["assay_type", "Assay", "Assay Type"], "unknown"),
        }
    )

    gene_col = _first_existing(frame, GENE_COLS)
    if gene_col is not None:
        out["gene_symbol"] = frame[gene_col].astype(str)

    out = _copy_numeric_features(frame, out)
    out = add_contrastive_features(out)
    out = add_normalized_columns(out)
    report = validate_schema(out)
    if not report.ok:
        raise ValueError(f"Normalized table failed canonical schema validation: {report}")
    return out


def normalize_neoranking_neopep(path: str | Path) -> pd.DataFrame:
    frame = read_table(path)
    out = normalize_candidate_table(
        frame,
        source_dataset="neoranking",
        study_default="neoranking",
    )

    rename = {
        "rnaseq_TPM": "expression_tpm",
        "CCF": "clonality_ccf",
        "mutant_rank_netMHCpan": "netmhcpan_mutant_rank",
        "mutant_rank_PRIME": "prime_mutant_rank",
        "mut_Rank_Stab": "netmhcstab_mutant_rank",
        "mut_binding_score": "binding_score",
        "DAI_NetMHC": "dai_netmhc",
        "DAI_MixMHC": "dai_mixmhc",
        "bestWTMatchScore_I": "self_similarity_score",
        "bestMutationScore_I": "foreignness_mutation_score",
    }
    for old, new in rename.items():
        if old in out.columns and new not in out.columns:
            out[new] = out[old]
    return out


def normalize_gartner_table(path: str | Path) -> pd.DataFrame:
    frame = read_table(path)
    out = normalize_candidate_table(
        frame,
        source_dataset="gartner_nci",
        study_default="gartner_nci",
    )
    rename = {
        "Rank NetMHC": "baseline_netmhc_rank",
        "Rank Nmer Model": "baseline_gartner_nmer_rank",
        "Rank Mmp Model": "baseline_gartner_mmp_rank",
        "AUC NetMHC": "source_auc_netmhc",
        "AUC NMER Model": "source_auc_nmer_model",
    }
    for old, new in rename.items():
        if old in out.columns and new not in out.columns:
            out[new] = out[old]
    return out


def write_normalized(frame: pd.DataFrame, output_path: str | Path) -> Path:
    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out, index=False)
    return out

