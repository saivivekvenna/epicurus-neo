"""Explicit adapter contract and limitations declaration."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from event_b.corpus import EventBCorpus
from event_b.manifest import SourceManifest


@dataclass(frozen=True)
class AdapterDeclaration:
    source_name: str
    source_version: str
    adapter_name: str
    adapter_version: str
    supported_entities: tuple[str, ...]
    supported_event_types: tuple[str, ...]
    known_limitations: tuple[str, ...]
    mapping_assumptions: tuple[str, ...]
    missing_fields: tuple[str, ...]


class StudyAdapter(Protocol):
    declaration: AdapterDeclaration

    def extract(self, manifest: SourceManifest) -> dict[str, object]: ...

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus: ...
