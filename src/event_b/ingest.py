"""Provenance-first staged ingestion orchestration."""

from __future__ import annotations

from dataclasses import dataclass

from event_b.adapters.base import StudyAdapter
from event_b.corpus import EventBCorpus
from event_b.manifest import SourceManifest
from event_b.review import ReviewIssue
from event_b.validation import validate_corpus


@dataclass(frozen=True)
class IngestionResult:
    stage: str
    manifest: SourceManifest
    extracted: dict[str, object]
    normalized_corpus: EventBCorpus
    accepted_corpus: EventBCorpus
    review_queue: tuple[ReviewIssue, ...]


def ingest_source(adapter: StudyAdapter, manifest: SourceManifest) -> IngestionResult:
    if manifest.adapter_name != adapter.declaration.adapter_name:
        raise ValueError("Source manifest and adapter declaration do not match")
    extracted = adapter.extract(manifest)
    normalized = adapter.normalize(extracted, manifest)
    validation = validate_corpus(normalized)
    return IngestionResult(
        "VALIDATED",
        manifest,
        extracted,
        validation.normalized_corpus,
        validation.accepted_corpus,
        validation.review_queue,
    )
