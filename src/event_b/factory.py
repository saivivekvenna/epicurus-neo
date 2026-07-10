"""Idempotent, interruption-safe Event-B study ingestion jobs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
from typing import Callable

from event_b.adapters.base import StudyAdapter
from event_b.export import export_corpus
from event_b.ingest import IngestionResult, ingest_source
from event_b.manifest import SourceManifest
from event_b.models import SCHEMA_VERSION
from event_b.registry import REGISTRY_VERSION, StudyRegistryEntry
from event_b.sources import validate_source_file


JOB_VERSION = "event-b-job-1.0.0"


@dataclass(frozen=True)
class JobCheckpoint:
    job_version: str
    study_id: str
    status: str
    stage: str
    input_fingerprint: str
    schema_version: str
    adapter_version: str
    registry_version: str
    output_dir: str
    stale_reasons: tuple[str, ...]
    updated_at: str
    error: str | None = None

    @classmethod
    def read(cls, path: str | Path) -> "JobCheckpoint":
        payload = json.loads(Path(path).read_text())
        payload["stale_reasons"] = tuple(payload.get("stale_reasons", ()))
        return cls(**payload)


def _utcnow() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def _atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def input_fingerprint(manifest: SourceManifest, adapter_version: str) -> str:
    payload = "|".join(
        (JOB_VERSION, REGISTRY_VERSION, SCHEMA_VERSION, adapter_version, manifest.fingerprint)
    )
    return sha256(payload.encode()).hexdigest()


class EventBJobRunner:
    def __init__(self, root: str | Path):
        self.root = Path(root)

    def checkpoint_path(self, study_id: str) -> Path:
        return self.root / "jobs" / f"{study_id}.json"

    def _write(self, checkpoint: JobCheckpoint) -> None:
        _atomic_json(self.checkpoint_path(checkpoint.study_id), asdict(checkpoint))

    def block(self, entry: StudyRegistryEntry) -> JobCheckpoint:
        checkpoint = JobCheckpoint(
            JOB_VERSION,
            entry.canonical_study_id,
            "BLOCKED",
            "REGISTRY_REVIEWED",
            "",
            SCHEMA_VERSION,
            "",
            REGISTRY_VERSION,
            str((self.root / "studies" / entry.canonical_study_id).resolve()),
            (),
            _utcnow(),
            entry.current_blocker,
        )
        self._write(checkpoint)
        return checkpoint

    def ingest(
        self,
        entry: StudyRegistryEntry,
        adapter: StudyAdapter,
        manifest: SourceManifest,
        *,
        rebuild: bool = False,
        after_stage: Callable[[str], None] | None = None,
    ) -> tuple[JobCheckpoint, IngestionResult | None]:
        if manifest.schema_version != SCHEMA_VERSION:
            raise ValueError("manifest schema version does not match the current corpus schema")
        if manifest.adapter_version != adapter.declaration.adapter_version:
            raise ValueError("manifest adapter version does not match the adapter declaration")
        for document in manifest.documents:
            if not document.local_path:
                raise ValueError(f"source document {document.document_id} has no local path")
            validation = validate_source_file(
                document.local_path, expected_media_type=document.media_type
            )
            if validation.checksum_sha256 != document.checksum_sha256:
                raise ValueError(f"checksum changed for {document.document_id}")

        fingerprint = input_fingerprint(manifest, adapter.declaration.adapter_version)
        checkpoint_path = self.checkpoint_path(entry.canonical_study_id)
        output = self.root / "studies" / entry.canonical_study_id
        existing = JobCheckpoint.read(checkpoint_path) if checkpoint_path.exists() else None
        marker = output / "_SUCCESS.json"
        if existing and existing.status == "COMPLETE" and not rebuild:
            stale = []
            if existing.input_fingerprint != fingerprint:
                stale.append("INPUT_FINGERPRINT_CHANGED")
            if existing.schema_version != SCHEMA_VERSION:
                stale.append("SCHEMA_VERSION_CHANGED")
            if existing.adapter_version != adapter.declaration.adapter_version:
                stale.append("ADAPTER_VERSION_CHANGED")
            if not marker.exists():
                stale.append("SUCCESS_MARKER_MISSING")
            if not stale:
                return existing, None
            stale_checkpoint = JobCheckpoint(
                **{
                    **asdict(existing),
                    "status": "STALE",
                    "stale_reasons": tuple(stale),
                    "updated_at": _utcnow(),
                }
            )
            self._write(stale_checkpoint)
            raise RuntimeError("study outputs are stale: " + ", ".join(stale))

        running = JobCheckpoint(
            JOB_VERSION,
            entry.canonical_study_id,
            "RUNNING",
            "SOURCES_VALIDATED",
            fingerprint,
            SCHEMA_VERSION,
            adapter.declaration.adapter_version,
            REGISTRY_VERSION,
            str(output.resolve()),
            (),
            _utcnow(),
        )
        self._write(running)
        if after_stage:
            after_stage(running.stage)
        try:
            result = ingest_source(adapter, manifest)
            review_queue = tuple(result.review_queue) + tuple(
                getattr(adapter, "review_issues", ())
            )
            normalized = JobCheckpoint(
                **{**asdict(running), "stage": "CORPUS_VALIDATED", "updated_at": _utcnow()}
            )
            self._write(normalized)
            if after_stage:
                after_stage(normalized.stage)
            written = export_corpus(
                result.accepted_corpus,
                output,
                review_queue=review_queue,
                source_manifests=[manifest],
                adapter_declarations=[adapter.declaration],
            )
            _atomic_json(
                marker,
                {
                    "input_fingerprint": fingerprint,
                    "written": sorted(written),
                    "review_queue_n": len(review_queue),
                },
            )
            complete = JobCheckpoint(
                **{
                    **asdict(normalized),
                    "status": "COMPLETE",
                    "stage": "EXPORTED",
                    "updated_at": _utcnow(),
                }
            )
            self._write(complete)
            return complete, result
        except BaseException as error:
            failed = JobCheckpoint(
                **{
                    **asdict(running),
                    "status": "INTERRUPTED",
                    "error": f"{type(error).__name__}: {error}",
                    "updated_at": _utcnow(),
                }
            )
            self._write(failed)
            raise
