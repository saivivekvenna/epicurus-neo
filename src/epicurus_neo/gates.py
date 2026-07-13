from __future__ import annotations

import re
from dataclasses import dataclass

import numpy as np
import pandas as pd


STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")

# Class-I length validity: a deterministic gate must remove only the genuinely IMPOSSIBLE, not the
# merely atypical. Class I canonically presents 8-11mers, but 12-14mers are documented (bulged
# conformations, esp. HLA-B) and standard class-I predictors (NetMHCpan) score 8-14. A hard 8-11
# cutoff was found to delete a validated immunogenic 12mer on the Müller NCI cohort, so the
# impossible-length bound is 8-14 (see MILESTONE_7_SYNTHESIS / test_class_i_12mer_is_not_dropped).
CLASS_I_MIN_LEN = 8
CLASS_I_MAX_LEN = 14


@dataclass(frozen=True)
class GateSummary:
    input_count: int
    survivor_count: int
    removed_count: int
    removed_fraction: float
    reason_counts: dict[str, int]


def _upper(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame:
        return pd.Series("", index=frame.index, dtype="object")
    return frame[column].fillna("").astype(str).str.strip().str.upper()


def _missense_mutant_residue(value: object) -> str | None:
    match = re.fullmatch(r"p\.([A-Z])\d+([A-Z])", str(value).strip().upper())
    return match.group(2) if match else None


def apply_deterministic_gate(frame: pd.DataFrame) -> pd.DataFrame:
    """Remove only rule-verifiable invalid routes, preserving every row and reason.

    Rules are first-wins so one candidate has one stable primary removal reason. The two
    biological rules validated on the raw SHERPA pool are vendor/source calls: a lost HLA route
    and a gene explicitly called unexpressed. Sequence rules run only when their preconditions are
    source-resolved.
    """
    out = frame.copy()
    reason = np.full(len(out), "", dtype=object)
    peptide = _upper(out, "mutant_peptide")

    malformed = peptide.map(lambda value: not value or bool(set(value) - STANDARD_AMINO_ACIDS))
    reason[(reason == "") & malformed.to_numpy()] = "MALFORMED_AA"

    mhc_class = _upper(out, "mhc_class")
    class_i = mhc_class.isin({"I", "CLASS_I", "CLASS I"})
    bad_class_i_length = class_i & ~peptide.str.len().between(CLASS_I_MIN_LEN, CLASS_I_MAX_LEN)
    reason[(reason == "") & bad_class_i_length.to_numpy()] = "BAD_CLASS_I_LENGTH"

    identity_columns = [
        column
        for column in ("patient_id", "mutation_id", "mutant_peptide", "hla_allele")
        if column in out
    ]
    if len(identity_columns) >= 3:
        duplicate = out.duplicated(identity_columns, keep="first")
        reason[(reason == "") & duplicate.to_numpy()] = "DUP_CANDIDATE"

    variant_type = _upper(out, "source_variant_type")
    if "protein_variant" in out:
        mutant_residue = out["protein_variant"].map(_missense_mutant_residue)
        mutation_missing = pd.Series(
            [
                variant_type.iloc[position] in {"SNV", "SNP", "MISSENSE"}
                and isinstance(residue, str)
                and residue not in sequence
                for position, (residue, sequence) in enumerate(zip(mutant_residue, peptide))
            ],
            index=out.index,
        )
        reason[(reason == "") & mutation_missing.to_numpy()] = "MUT_NOT_IN_PEPTIDE"

    lost_hla = _upper(out, "hla_loh_call").isin({"Y", "YES", "TRUE", "1", "LOST"})
    reason[(reason == "") & lost_hla.to_numpy()] = "HLA_LOH_LOST_ALLELE"

    not_expressed = _upper(out, "expression_call").isin({"N", "NO", "FALSE", "0", "NOT_EXPRESSED"})
    reason[(reason == "") & not_expressed.to_numpy()] = "GENE_NOT_EXPRESSED"

    out["deterministic_gate_reason"] = reason
    out["deterministic_gate_pass"] = reason == ""
    out["deterministic_gate_policy"] = "epicurus-validity-gate-1.0.0"
    return out


def summarize_gate(frame: pd.DataFrame) -> GateSummary:
    if "deterministic_gate_pass" not in frame:
        raise ValueError("frame has not passed through apply_deterministic_gate")
    reason_counts = {
        str(reason): int(count)
        for reason, count in frame.loc[
            ~frame["deterministic_gate_pass"], "deterministic_gate_reason"
        ].value_counts().items()
    }
    input_count = len(frame)
    survivor_count = int(frame["deterministic_gate_pass"].sum())
    removed_count = input_count - survivor_count
    return GateSummary(
        input_count=input_count,
        survivor_count=survivor_count,
        removed_count=removed_count,
        removed_fraction=float(removed_count / input_count) if input_count else 0.0,
        reason_counts=reason_counts,
    )
