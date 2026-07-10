"""Versioned registry for public Event-B studies and their ingestion state."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Any

import yaml


REGISTRY_VERSION = "event-b-registry-1.0.0"


class StudyStatus(str, Enum):
    REGISTERED = "REGISTERED"
    SOURCES_PENDING = "SOURCES_PENDING"
    SOURCES_PINNED = "SOURCES_PINNED"
    ADAPTER_IMPLEMENTED = "ADAPTER_IMPLEMENTED"
    REVIEW_REQUIRED = "REVIEW_REQUIRED"
    ACCEPTED = "ACCEPTED"
    BLOCKED_SOURCE_UNAVAILABLE = "BLOCKED_SOURCE_UNAVAILABLE"
    BLOCKED_LABEL_SEMANTICS = "BLOCKED_LABEL_SEMANTICS"
    BLOCKED_PATIENT_MAPPING = "BLOCKED_PATIENT_MAPPING"
    REJECTED_NOT_EVENT_B = "REJECTED_NOT_EVENT_B"


@dataclass(frozen=True)
class StudyRegistryEntry:
    canonical_study_id: str
    cohort_id: str
    publication_ids: tuple[str, ...]
    trial_id: str | None
    cancer_type: str
    vaccine_platform: str
    antigen_design: str
    source_availability: str
    adapter_status: str
    ingestion_status: StudyStatus
    expected_data_structures: tuple[str, ...]
    known_publication_overlap: tuple[str, ...]
    known_limitations: tuple[str, ...]
    current_blocker: str | None

    @classmethod
    def from_dict(cls, row: dict[str, Any]) -> "StudyRegistryEntry":
        tuple_fields = (
            "publication_ids",
            "expected_data_structures",
            "known_publication_overlap",
            "known_limitations",
        )
        normalized = dict(row)
        for field in tuple_fields:
            normalized[field] = tuple(normalized.get(field) or ())
        normalized["ingestion_status"] = StudyStatus(normalized["ingestion_status"])
        entry = cls(**normalized)
        entry.validate()
        return entry

    def validate(self) -> None:
        if not self.canonical_study_id or not self.cohort_id:
            raise ValueError("registry entries require canonical study and cohort IDs")
        if not self.publication_ids:
            raise ValueError(f"{self.canonical_study_id} requires at least one publication ID")
        if self.antigen_design not in {"SHARED", "PERSONALIZED"}:
            raise ValueError(f"{self.canonical_study_id} has invalid antigen_design")
        if self.ingestion_status.value.startswith("BLOCKED_") and not self.current_blocker:
            raise ValueError(f"{self.canonical_study_id} is blocked without a blocker")
        if self.ingestion_status is StudyStatus.ACCEPTED and self.adapter_status != "IMPLEMENTED":
            raise ValueError(f"{self.canonical_study_id} is accepted without an adapter")


@dataclass(frozen=True)
class StudyRegistry:
    registry_version: str
    studies: tuple[StudyRegistryEntry, ...]

    @classmethod
    def read(cls, path: str | Path) -> "StudyRegistry":
        payload = yaml.safe_load(Path(path).read_text())
        if payload.get("registry_version") != REGISTRY_VERSION:
            raise ValueError(f"unsupported registry version {payload.get('registry_version')!r}")
        studies = tuple(StudyRegistryEntry.from_dict(row) for row in payload.get("studies", ()))
        identifiers = [entry.canonical_study_id for entry in studies]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("canonical study IDs must be unique")
        cohorts = [entry.cohort_id for entry in studies]
        if len(cohorts) != len(set(cohorts)):
            raise ValueError("cohort IDs must be unique; publication overlap belongs in metadata")
        return cls(payload["registry_version"], studies)

    def get(self, study_id: str) -> StudyRegistryEntry:
        for entry in self.studies:
            if entry.canonical_study_id == study_id:
                return entry
        raise KeyError(f"study is not registered: {study_id}")

    def backbone(self) -> tuple[StudyRegistryEntry, ...]:
        order = ("mkras_vax_2026", "pdac_neovax_2023", "nous_209_2025", "fukuoka_dc")
        return tuple(self.get(study_id) for study_id in order)
