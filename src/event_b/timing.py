"""Prevent post-vaccine outcomes from entering pre-selection feature matrices."""

from __future__ import annotations

import pandas as pd

from event_b.evidence import validate_evidence
from event_b.models import InformationTiming


def preselection_evidence(frame: pd.DataFrame, *, reject_unknown: bool = True) -> pd.DataFrame:
    evidence = validate_evidence(frame)
    timing = evidence.information_timing.astype(str).str.upper()
    if reject_unknown and timing.eq(InformationTiming.UNKNOWN.value).any():
        ids = evidence.loc[timing.eq(InformationTiming.UNKNOWN.value), "evidence_id"].tolist()
        raise ValueError(f"Unknown-timing evidence cannot enter pre-selection features: {ids[:5]}")
    return evidence.loc[timing.eq(InformationTiming.PRE_SELECTION.value)].copy()


def assert_preselection_columns(
    requested_columns: list[str] | tuple[str, ...],
    timing_registry: dict[str, str],
) -> None:
    unknown = set(requested_columns).difference(timing_registry)
    if unknown:
        raise ValueError(f"Feature timing is undeclared: {sorted(unknown)}")
    leaked = [
        column
        for column in requested_columns
        if str(timing_registry[column]).upper() != InformationTiming.PRE_SELECTION.value
    ]
    if leaked:
        raise ValueError(f"Outcome-only features requested before vaccine design: {sorted(leaked)}")
