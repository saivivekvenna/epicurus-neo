"""Explicit biological-event and assay metadata for three-state labels."""

from __future__ import annotations

from enum import Enum

import pandas as pd

from benchmark.labels import Label, coerce_label


class EventType(str, Enum):
    PRE_EXISTING_REACTIVITY = "pre_existing_reactivity"
    VACCINE_INDUCED_RESPONSE = "vaccine_induced_response"
    POST_TREATMENT_REACTIVITY = "post_treatment_reactivity"
    UNKNOWN = "unknown"


class Assay(str, Enum):
    ELISPOT = "elispot"
    ICS = "ics"
    TETRAMER = "tetramer"
    MANAFEST = "manafest"
    CYTOKINE_RELEASE = "cytokine_release"
    OTHER = "other"
    UNKNOWN = "unknown"


class Timepoint(str, Enum):
    PRE_VACCINE = "pre_vaccine"
    POST_PRIME = "post_prime"
    POST_BOOST = "post_boost"
    POST_VACCINE = "post_vaccine"
    POST_TREATMENT = "post_treatment"
    UNKNOWN = "unknown"


REQUIRED_LABEL_METADATA = {
    "record_id",
    "dataset_id",
    "event_type",
    "assay",
    "label",
    "timepoint",
    "provenance",
}


def _coerce_enum(value, enum_type: type[Enum], field: str):
    if isinstance(value, enum_type):
        return value.value
    try:
        return enum_type(str(value).strip().lower()).value
    except ValueError as error:
        raise ValueError(f"Unknown {field}: {value!r}") from error


def validate_label_metadata(frame: pd.DataFrame) -> pd.DataFrame:
    """Validate Event A/B semantics without guessing missing biological context."""
    missing = REQUIRED_LABEL_METADATA.difference(frame.columns)
    if missing:
        raise ValueError(f"Label metadata missing columns: {sorted(missing)}")
    if frame["record_id"].isna().any() or frame["record_id"].duplicated().any():
        raise ValueError("record_id must be non-null and unique")
    if frame["dataset_id"].isna().any():
        raise ValueError("dataset_id must be non-null")
    if frame["provenance"].fillna("").astype(str).str.strip().eq("").any():
        raise ValueError("Every label requires non-empty provenance")

    out = frame.copy()
    out["event_type"] = out["event_type"].map(
        lambda value: _coerce_enum(value, EventType, "event_type")
    )
    out["assay"] = out["assay"].map(lambda value: _coerce_enum(value, Assay, "assay"))
    out["label"] = out["label"].map(lambda value: coerce_label(value).name)
    out["timepoint"] = out["timepoint"].map(
        lambda value: _coerce_enum(value, Timepoint, "timepoint")
    )

    tested = out["label"].isin([Label.POSITIVE.name, Label.TESTED_NEGATIVE.name])
    unknown_assay = out["assay"].eq(Assay.UNKNOWN.value)
    if (tested & unknown_assay).any():
        raise ValueError("Tested labels require an explicit assay")
    vaccine_event = out["event_type"].eq(EventType.VACCINE_INDUCED_RESPONSE.value)
    pre_vaccine = out["timepoint"].eq(Timepoint.PRE_VACCINE.value)
    if (vaccine_event & pre_vaccine & tested).any():
        raise ValueError("A tested vaccine-induced response cannot be measured pre-vaccine")
    return out
