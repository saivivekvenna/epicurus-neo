from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from epicurus_neo.schema import add_normalized_columns


@dataclass(frozen=True)
class LeakageReport:
    shared_mutant_hla: tuple[str, ...]
    shared_wildtype_hla: tuple[str, ...]
    shared_patients: tuple[str, ...]
    shared_studies: tuple[str, ...]

    @property
    def has_leakage(self) -> bool:
        return bool(
            self.shared_mutant_hla
            or self.shared_wildtype_hla
            or self.shared_patients
            or self.shared_studies
        )

    def has_blocking_leakage(self, *, include_shared_studies: bool = True) -> bool:
        return bool(
            self.shared_mutant_hla
            or self.shared_wildtype_hla
            or self.shared_patients
            or (include_shared_studies and self.shared_studies)
        )


def _shared_values(train: pd.Series, test: pd.Series) -> tuple[str, ...]:
    shared = set(train.dropna().astype(str)) & set(test.dropna().astype(str))
    return tuple(sorted(value for value in shared if value))


def _peptide_hla_keys(frame: pd.DataFrame, peptide_kind: str) -> pd.Series:
    peptide_col = f"{peptide_kind}_peptide_norm"
    key_col = f"{peptide_kind}_hla_key"
    peptide = frame[peptide_col].fillna("").astype(str).str.strip()
    return frame.loc[peptide.ne(""), key_col]


def detect_exact_leakage(train: pd.DataFrame, test: pd.DataFrame) -> LeakageReport:
    train_norm = add_normalized_columns(train)
    test_norm = add_normalized_columns(test)
    return LeakageReport(
        shared_mutant_hla=_shared_values(
            _peptide_hla_keys(train_norm, "mutant"),
            _peptide_hla_keys(test_norm, "mutant"),
        ),
        shared_wildtype_hla=_shared_values(
            _peptide_hla_keys(train_norm, "wildtype"),
            _peptide_hla_keys(test_norm, "wildtype"),
        ),
        shared_patients=_shared_values(train_norm["patient_id"], test_norm["patient_id"]),
        shared_studies=_shared_values(train_norm["study_id"], test_norm["study_id"]),
    )


def purge_train_overlaps(train: pd.DataFrame, test: pd.DataFrame) -> pd.DataFrame:
    """Remove train rows that share exact peptide/HLA keys with the test fold."""
    train_norm = add_normalized_columns(train)
    test_norm = add_normalized_columns(test)
    blocked_mutant = set(_peptide_hla_keys(test_norm, "mutant").dropna().astype(str))
    blocked_wildtype = set(_peptide_hla_keys(test_norm, "wildtype").dropna().astype(str))
    train_mutant = _peptide_hla_keys(train_norm, "mutant")
    train_wildtype = _peptide_hla_keys(train_norm, "wildtype")
    keep = ~(
        train_norm.index.to_series().isin(train_mutant[train_mutant.isin(blocked_mutant)].index)
        | train_norm.index.to_series().isin(
            train_wildtype[train_wildtype.isin(blocked_wildtype)].index
        )
    )
    return train.loc[keep].copy()
