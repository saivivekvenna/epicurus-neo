from __future__ import annotations

import pandas as pd

from epicurus_neo.features import (
    add_contrastive_features,
    add_sequence_features,
    infer_numeric_feature_columns,
)

# Anchor/TCR-face positions use the class-I binding register; nulled elsewhere.
_CLASS_I_ONLY = ("mutation_anchor_count", "mutation_tcr_face_count")

# Belt-and-suspenders on top of features.NON_FEATURE_COLUMNS.
_EXTRA_BANNED = {
    "study_id",
    "patient_id",
    "candidate_id",
    "cancer_type",
    "vaccine_platform",
    "mhc_class",
    "hla_alleles",
    "hla_allele",
    "mutant_peptide",
    "wildtype_peptide",
    "response_label",
    "event_type",
    "sample_date",
    "timepoint",
    "sample_id",
    "genomic_variant",
    "gene",
    "transcript",
    "protein_change",
    "candidate_source",
    "vaccine_inclusion",
    "vaccine_inclusion_origin",
    "generation_provenance",
    "mutant_wildtype_verified",
    "provenance_id",
    "schema_version",
}


def _add_class_indicators(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    out["is_class_i"] = (out.mhc_class == "CLASS_I").astype(float)
    out["is_class_ii"] = (out.mhc_class == "CLASS_II").astype(float)
    return out


def build_feature_matrix(frame: pd.DataFrame, tier: str) -> pd.DataFrame:
    """Build the requested feature tier. Missingness is never a study label."""
    if tier not in {"core", "contrastive", "presentation"}:
        raise ValueError(f"unknown tier: {tier!r}")
    matrix = _add_class_indicators(add_sequence_features(frame))
    if tier == "contrastive":
        matrix = add_contrastive_features(matrix)
        non_class_i = matrix.mhc_class != "CLASS_I"
        for column in _CLASS_I_ONLY:
            if column in matrix.columns:
                matrix.loc[non_class_i, column] = float("nan")
    # "presentation" tier consumes presentation_score/mhcflurry_* already merged
    # onto ``frame`` by Task 5; if absent the tier degenerates to core.
    return matrix


def feature_columns(matrix: pd.DataFrame) -> list[str]:
    numeric = infer_numeric_feature_columns(matrix)
    return [column for column in numeric if column not in _EXTRA_BANNED]


def assert_no_banned_features(columns: list[str]) -> None:
    leaked = sorted(set(columns) & _EXTRA_BANNED)
    if leaked:
        raise AssertionError(f"banned columns leaked into feature set: {leaked}")
