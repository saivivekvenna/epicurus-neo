from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import pandas as pd


VALID_LABELS = {"positive", "negative", "unknown"}

REQUIRED_COLUMNS = {
    "candidate_id",
    "source_dataset",
    "study_id",
    "patient_id",
    "hla_allele",
    "mutant_peptide",
    "wildtype_peptide",
    "label",
    "label_weight",
    "assay_type",
}


@dataclass(frozen=True)
class SchemaReport:
    row_count: int
    missing_columns: tuple[str, ...]
    invalid_labels: tuple[str, ...]
    duplicate_candidate_ids: int

    @property
    def ok(self) -> bool:
        return (
            not self.missing_columns
            and not self.invalid_labels
            and self.duplicate_candidate_ids == 0
        )


def normalize_peptide(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip().upper()


def normalize_hla(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    hla = str(value).strip().upper().replace("HLA-", "")
    if "*" not in hla and len(hla) >= 5 and hla[1].isdigit():
        hla = f"{hla[0]}*{hla[1:]}"
    if "*" in hla:
        locus, fields = hla.split("*", maxsplit=1)
        if ":" not in fields and len(fields) == 4 and fields.isdigit():
            fields = f"{fields[:2]}:{fields[2:]}"
        hla = f"{locus}*{fields}"
    return f"HLA-{hla}"


def add_normalized_columns(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["mutant_peptide_norm"] = out["mutant_peptide"].map(normalize_peptide)
    out["wildtype_peptide_norm"] = out["wildtype_peptide"].map(normalize_peptide)
    out["hla_allele_norm"] = out["hla_allele"].map(normalize_hla)
    out["mutant_hla_key"] = out["mutant_peptide_norm"] + "|" + out["hla_allele_norm"]
    out["wildtype_hla_key"] = out["wildtype_peptide_norm"] + "|" + out["hla_allele_norm"]
    return out


def validate_schema(frame: pd.DataFrame, required: Iterable[str] = REQUIRED_COLUMNS) -> SchemaReport:
    required_set = set(required)
    missing = tuple(sorted(required_set - set(frame.columns)))

    invalid_labels: tuple[str, ...] = ()
    if "label" in frame.columns:
        observed = {str(label) for label in frame["label"].dropna().unique()}
        invalid_labels = tuple(sorted(observed - VALID_LABELS))

    duplicate_candidate_ids = 0
    if "candidate_id" in frame.columns:
        duplicate_candidate_ids = int(frame["candidate_id"].duplicated().sum())

    return SchemaReport(
        row_count=len(frame),
        missing_columns=missing,
        invalid_labels=invalid_labels,
        duplicate_candidate_ids=duplicate_candidate_ids,
    )


def supervised_rows(frame: pd.DataFrame) -> pd.DataFrame:
    """Return experimentally labeled rows; unknown/unassayed rows are excluded."""
    return frame[frame["label"].isin(["positive", "negative"])].copy()
