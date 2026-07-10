"""Generic CSV/TSV and JSONL adapters with source-specific mappings kept outside schemas."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from event_b.adapters.base import AdapterDeclaration
from event_b.corpus import EventBCorpus
from event_b.manifest import SourceManifest
from event_b.models import BiologicalEvent, SCHEMAS


def _read(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    if path.suffix.lower() == ".jsonl":
        return pd.read_json(path, lines=True)
    raise ValueError(f"Unsupported generic adapter file: {path}")


class GenericTableAdapter:
    def __init__(
        self,
        *,
        source_name: str,
        source_version: str,
        entity_paths: dict[str, str | Path],
        column_maps: dict[str, dict[str, str]] | None = None,
        constants: dict[str, dict[str, object]] | None = None,
        supported_event_types: tuple[str, ...] = (BiologicalEvent.UNKNOWN_EVENT.value,),
        limitations: tuple[str, ...] = (),
        assumptions: tuple[str, ...] = (),
    ):
        unknown = set(entity_paths).difference(SCHEMAS)
        if unknown:
            raise ValueError(f"Unknown canonical entities: {sorted(unknown)}")
        self.entity_paths = {key: Path(value) for key, value in entity_paths.items()}
        self.column_maps = column_maps or {}
        self.constants = constants or {}
        self.declaration = AdapterDeclaration(
            source_name,
            source_version,
            "generic_table",
            "1.0.0",
            tuple(sorted(entity_paths)),
            supported_event_types,
            limitations,
            assumptions,
            (),
        )

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        return {entity: _read(path) for entity, path in self.entity_paths.items()}

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        tables = {}
        for entity, value in extracted.items():
            frame = pd.DataFrame(value).rename(columns=self.column_maps.get(entity, {})).copy()
            for column, constant in self.constants.get(entity, {}).items():
                frame[column] = constant
            if "source_manifest_id" in SCHEMAS[entity].columns:
                frame["source_manifest_id"] = manifest.manifest_id
            tables[entity] = SCHEMAS[entity].normalize(frame)
        corpus = EventBCorpus()
        for entity, frame in tables.items():
            setattr(corpus, entity, frame)
        return corpus


class JsonlExtractionAdapter(GenericTableAdapter):
    """Import reviewable extraction JSONL; records remain unaccepted by default."""

    def __init__(self, source_name: str, source_version: str, path: str | Path):
        self.path = Path(path)
        self.declaration = AdapterDeclaration(
            source_name,
            source_version,
            "jsonl_extraction",
            "1.0.0",
            tuple(SCHEMAS),
            tuple(item.value for item in BiologicalEvent),
            ("LLM-extracted values require deterministic validation and review",),
            ("Each JSONL row contains entity_type and record",),
            (),
        )

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        grouped: dict[str, list[dict]] = {}
        for line in self.path.read_text().splitlines():
            payload = json.loads(line)
            entity = payload["entity_type"]
            if entity not in SCHEMAS:
                raise ValueError(f"Unknown extracted entity: {entity}")
            grouped.setdefault(entity, []).append(payload["record"])
        return grouped

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        corpus = EventBCorpus()
        for entity, records in extracted.items():
            frame = pd.DataFrame(records)
            if "source_manifest_id" in SCHEMAS[entity].columns:
                frame["source_manifest_id"] = manifest.manifest_id
            setattr(corpus, entity, SCHEMAS[entity].normalize(frame))
        return corpus
