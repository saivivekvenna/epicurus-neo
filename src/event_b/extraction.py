"""Endpoint-optional, cached, schema-constrained LLM extraction tasks."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path
from typing import Protocol

from jsonschema import ValidationError, validate

from event_b.models import (
    AssayType,
    BiologicalEvent,
    ResponseLabel,
    ReviewStatus,
    SCHEMA_VERSION,
    ValueOrigin,
)


EXTRACTION_VERSION = "event-b-extraction-1.0.0"

ASSAY_EXTRACTION_SCHEMA = {
    "$schema": "https://json-schema.org/draft/2020-12/schema",
    "title": "Event-B candidate-assay extraction",
    "type": "object",
    "additionalProperties": False,
    "required": ["assays", "provenance"],
    "properties": {
        "assays": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": True,
                "required": [
                    "assay_id",
                    "patient_id",
                    "study_id",
                    "event_type",
                    "response_label",
                    "relative_to_vaccine",
                    "explicit_assay_inclusion",
                    "review_status",
                    "provenance_id",
                ],
                "properties": {
                    "assay_id": {"type": "string", "minLength": 1},
                    "patient_id": {"type": "string", "minLength": 1},
                    "study_id": {"type": "string", "minLength": 1},
                    "candidate_id": {"type": ["string", "null"]},
                    "assay_type": {"enum": [item.value for item in AssayType]},
                    "event_type": {"enum": [item.value for item in BiologicalEvent]},
                    "response_label": {"enum": [item.value for item in ResponseLabel]},
                    "relative_to_vaccine": {"type": "string"},
                    "explicit_assay_inclusion": {"type": ["boolean", "null"]},
                    "review_status": {"const": ReviewStatus.NEEDS_REVIEW.value},
                    "provenance_id": {"type": "string", "minLength": 1},
                },
            },
        },
        "provenance": {
            "type": "array",
            "items": {
                "type": "object",
                "additionalProperties": True,
                "required": [
                    "provenance_id",
                    "entity_type",
                    "entity_id",
                    "field_name",
                    "source_document",
                    "source_fragment",
                    "extraction_confidence",
                    "value_origin",
                    "review_status",
                ],
                "properties": {
                    "value_origin": {"const": ValueOrigin.LLM_EXTRACTED.value},
                    "review_status": {"const": ReviewStatus.NEEDS_REVIEW.value},
                    "extraction_confidence": {"type": "number", "minimum": 0, "maximum": 1},
                },
            },
        },
    },
}


@dataclass(frozen=True)
class ExtractionTask:
    task_id: str
    study_id: str
    source_document: str
    source_checksum: str
    source_text: str
    prompt: str
    output_schema: dict
    extraction_version: str = EXTRACTION_VERSION
    canonical_schema_version: str = SCHEMA_VERSION

    @classmethod
    def create(
        cls,
        *,
        study_id: str,
        source_document: str,
        source_checksum: str,
        source_text: str,
    ) -> "ExtractionTask":
        prompt = (
            "Extract only candidate-resolved assay observations explicitly supported by the source. "
            "Never infer TESTED_NEGATIVE from omission, vaccine inclusion from candidate generation, "
            "or Event B from clinical outcome. Use UNTESTED/UNKNOWN_EVENT when resolution is unsafe. "
            "All records must remain NEEDS_REVIEW and preserve source fragments."
        )
        identity = f"{EXTRACTION_VERSION}|{study_id}|{source_checksum}|{source_text}"
        return cls(
            "extract:" + sha256(identity.encode()).hexdigest()[:24],
            study_id,
            source_document,
            source_checksum,
            source_text,
            prompt,
            ASSAY_EXTRACTION_SCHEMA,
        )


class InferenceProvider(Protocol):
    provider_id: str

    def infer(self, task: ExtractionTask) -> dict: ...


class NoEndpointProvider:
    provider_id = "none"

    def infer(self, task: ExtractionTask) -> dict:
        del task
        raise RuntimeError("No inference endpoint is configured")


@dataclass(frozen=True)
class ExtractionResult:
    task_id: str
    status: str
    provider_id: str
    attempts: int
    cache_path: str | None
    error: str | None = None


class ExtractionCache:
    def __init__(self, directory: str | Path):
        self.directory = Path(directory)
        self.directory.mkdir(parents=True, exist_ok=True)

    def path(self, task: ExtractionTask, provider_id: str) -> Path:
        key = sha256(f"{task.task_id}|{provider_id}".encode()).hexdigest()
        return self.directory / f"{key}.json"

    def load(self, task: ExtractionTask, provider_id: str) -> dict | None:
        path = self.path(task, provider_id)
        return json.loads(path.read_text()) if path.exists() else None

    def store(self, task: ExtractionTask, provider_id: str, raw_output: dict) -> Path:
        path = self.path(task, provider_id)
        path.write_text(json.dumps(raw_output, indent=2, sort_keys=True) + "\n")
        return path


def emit_extraction_tasks(tasks: list[ExtractionTask], output_dir: str | Path) -> tuple[Path, Path]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    task_path = output / "extraction_tasks.jsonl"
    schema_path = output / "assay_extraction.schema.json"
    rows = [
        json.dumps(asdict(task), sort_keys=True)
        for task in sorted(tasks, key=lambda item: item.task_id)
    ]
    task_path.write_text("\n".join(rows) + ("\n" if rows else ""))
    schema_path.write_text(json.dumps(ASSAY_EXTRACTION_SCHEMA, indent=2, sort_keys=True) + "\n")
    return task_path, schema_path


def run_extraction(
    task: ExtractionTask,
    provider: InferenceProvider | None,
    cache: ExtractionCache,
    *,
    retries: int = 2,
) -> ExtractionResult:
    provider = provider or NoEndpointProvider()
    cached = cache.load(task, provider.provider_id)
    if cached is not None:
        validate(instance=cached, schema=task.output_schema)
        return ExtractionResult(
            task.task_id,
            "CACHED",
            provider.provider_id,
            0,
            str(cache.path(task, provider.provider_id)),
        )
    if isinstance(provider, NoEndpointProvider):
        return ExtractionResult(
            task.task_id, "PENDING", provider.provider_id, 0, None, "No inference endpoint"
        )
    last_error = None
    for attempt in range(1, retries + 2):
        try:
            raw = provider.infer(task)
            validate(instance=raw, schema=task.output_schema)
            path = cache.store(task, provider.provider_id, raw)
            return ExtractionResult(
                task.task_id, "EXTRACTED_NEEDS_REVIEW", provider.provider_id, attempt, str(path)
            )
        except (RuntimeError, ValidationError, ValueError) as error:
            last_error = str(error)
    return ExtractionResult(
        task.task_id, "FAILED", provider.provider_id, retries + 1, None, last_error
    )
