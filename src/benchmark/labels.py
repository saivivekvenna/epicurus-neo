"""Three-state label schema and dataset-manifest validation."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


class Label(Enum):
    POSITIVE = 1
    TESTED_NEGATIVE = 0
    UNTESTED = -1


_TEXT_LABELS = {
    "POSITIVE": Label.POSITIVE,
    "TESTED_NEGATIVE": Label.TESTED_NEGATIVE,
    "NEGATIVE": Label.TESTED_NEGATIVE,
    "UNTESTED": Label.UNTESTED,
    "UNKNOWN": Label.UNTESTED,
    "1": Label.POSITIVE,
    "0": Label.TESTED_NEGATIVE,
    "-1": Label.UNTESTED,
}


def coerce_label(value: Any) -> Label:
    if isinstance(value, Label):
        return value
    key = str(getattr(value, "name", value)).strip().upper()
    try:
        return _TEXT_LABELS[key]
    except KeyError as error:
        raise ValueError(f"Unknown label value: {value!r}") from error


def validate_labels(
    values: pd.Series | list[Any], *, require_three_states: bool = True
) -> pd.Series:
    """Validate and return Label values; reject binary loader outputs by default."""
    series = pd.Series(values, copy=True)
    labels = series.map(coerce_label)
    observed = set(labels)
    if require_three_states and observed != set(Label):
        missing = sorted(label.name for label in set(Label).difference(observed))
        raise ValueError(
            "Dataset loader emitted a two-valued/incomplete label schema; "
            f"missing states: {missing}"
        )
    return labels


@dataclass(frozen=True)
class DatasetManifest:
    dataset: str
    negative_provenance: str
    label_states: tuple[Label, ...]


def load_manifest(path: str | Path) -> DatasetManifest:
    payload = yaml.safe_load(Path(path).read_text())
    missing = {"dataset", "negative_provenance", "label_states"}.difference(payload)
    if missing:
        raise ValueError(f"Manifest missing required keys: {sorted(missing)}")
    states = tuple(coerce_label(value) for value in payload["label_states"])
    if set(states) != set(Label):
        raise ValueError("Manifest must declare POSITIVE, TESTED_NEGATIVE, and UNTESTED")
    if not str(payload["negative_provenance"]).strip():
        raise ValueError("negative_provenance must be explicit")
    return DatasetManifest(str(payload["dataset"]), str(payload["negative_provenance"]), states)
