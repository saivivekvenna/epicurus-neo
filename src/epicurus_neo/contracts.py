from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import pandas as pd


CANDIDATE_SCHEMA_VERSION = "epicurus-neo-candidate-1.0.0"
RANKED_SCHEMA_VERSION = "epicurus-neo-ranked-1.0.0"

CANDIDATE_REQUIRED_COLUMNS = {
    "candidate_id",
    "patient_id",
    "mutant_peptide",
}

RANKED_REQUIRED_COLUMNS = CANDIDATE_REQUIRED_COLUMNS | {
    "epicurus_neo_evidence_score",
    "epicurus_neo_lower_evidence_score",
    "evidence_tier",
    "selected",
    "rank",
}


@dataclass(frozen=True)
class ContractReport:
    schema_version: str
    row_count: int
    missing_columns: tuple[str, ...]
    empty_required_values: tuple[str, ...]
    duplicate_candidate_ids: int
    invalid_peptides: int

    @property
    def ok(self) -> bool:
        return not any(
            (
                self.missing_columns,
                self.empty_required_values,
                self.duplicate_candidate_ids,
                self.invalid_peptides,
            )
        )


def _empty_columns(frame: pd.DataFrame, columns: Iterable[str]) -> tuple[str, ...]:
    empty: list[str] = []
    for column in columns:
        if column not in frame.columns:
            continue
        values = frame[column]
        if values.isna().any() or values.astype(str).str.strip().eq("").any():
            empty.append(column)
    return tuple(sorted(empty))


def _invalid_peptide_count(frame: pd.DataFrame) -> int:
    if "mutant_peptide" not in frame.columns:
        return 0
    peptides = frame["mutant_peptide"].fillna("").astype(str).str.strip().str.upper()
    return int((~peptides.str.fullmatch(r"[ACDEFGHIKLMNPQRSTVWY]+", na=False)).sum())


def validate_candidate_contract(frame: pd.DataFrame) -> ContractReport:
    missing = tuple(sorted(CANDIDATE_REQUIRED_COLUMNS - set(frame.columns)))
    return ContractReport(
        schema_version=CANDIDATE_SCHEMA_VERSION,
        row_count=len(frame),
        missing_columns=missing,
        empty_required_values=_empty_columns(frame, CANDIDATE_REQUIRED_COLUMNS),
        duplicate_candidate_ids=(
            int(frame["candidate_id"].duplicated().sum())
            if "candidate_id" in frame.columns
            else 0
        ),
        invalid_peptides=_invalid_peptide_count(frame),
    )


def validate_ranked_contract(frame: pd.DataFrame) -> ContractReport:
    missing = tuple(sorted(RANKED_REQUIRED_COLUMNS - set(frame.columns)))
    return ContractReport(
        schema_version=RANKED_SCHEMA_VERSION,
        row_count=len(frame),
        missing_columns=missing,
        empty_required_values=_empty_columns(frame, CANDIDATE_REQUIRED_COLUMNS),
        duplicate_candidate_ids=(
            int(frame["candidate_id"].duplicated().sum())
            if "candidate_id" in frame.columns
            else 0
        ),
        invalid_peptides=_invalid_peptide_count(frame),
    )
