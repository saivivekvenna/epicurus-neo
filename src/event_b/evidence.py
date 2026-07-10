"""Validation for separate, extensible recognition-evidence channels."""

from __future__ import annotations

import pandas as pd

from event_b.models import (
    AvailabilityStatus,
    EvidenceFamily,
    InformationTiming,
    SCHEMAS,
)


RELIABILITY_DIMENSIONS = (
    "patient_specificity",
    "functional_relevance",
    "vaccine_relevance",
    "candidate_specificity",
    "assay_directness",
    "temporal_clarity",
    "source_completeness",
)


def validate_evidence(frame: pd.DataFrame) -> pd.DataFrame:
    out = SCHEMAS["recognition_evidence"].normalize(frame)
    allowed_families = {item.value for item in EvidenceFamily}
    allowed_availability = {item.value for item in AvailabilityStatus}
    allowed_timing = {item.value for item in InformationTiming}
    families = out.evidence_family.astype(str).str.upper()
    if not families.isin(allowed_families).all():
        raise ValueError(f"Unknown evidence families: {sorted(set(families) - allowed_families)}")
    availability = out.availability_status.astype(str).str.upper()
    if not availability.isin(allowed_availability).all():
        raise ValueError("Unknown evidence availability status")
    timing = out.information_timing.astype(str).str.upper()
    if not timing.isin(allowed_timing).all():
        raise ValueError("Unknown evidence information timing")
    for column in RELIABILITY_DIMENSIONS:
        numeric = pd.to_numeric(out[column], errors="coerce")
        present = out[column].notna()
        if (present & (~numeric.between(0, 1))).any():
            raise ValueError(f"{column} must be independently recorded between 0 and 1")
        out[column] = numeric
    return out


def evidence_availability_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    """Return channel availability only; never collapse evidence into a score."""
    evidence = validate_evidence(frame)
    if evidence.empty:
        return pd.DataFrame(columns=["candidate_id", "patient_id"])
    available = (
        evidence.availability_status.astype(str).str.upper().eq(AvailabilityStatus.AVAILABLE.value)
    )
    matrix = (
        evidence.assign(_available=available.astype(int))
        .pivot_table(
            index=["candidate_id", "patient_id"],
            columns="evidence_family",
            values="_available",
            aggfunc="max",
            fill_value=0,
        )
        .reset_index()
    )
    matrix.columns.name = None
    return matrix
