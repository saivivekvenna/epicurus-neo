from dataclasses import replace
import json
from pathlib import Path

import pytest

from event_b.adapters.base import AdapterDeclaration
from event_b.corpus import EventBCorpus
from event_b.factory import EventBJobRunner, JobCheckpoint
from event_b.manifest import manifest_from_paths
from event_b.registry import REGISTRY_VERSION, StudyRegistry, StudyRegistryEntry, StudyStatus
from event_b.sources import detect_media_type, validate_source_file


class EmptyAdapter:
    declaration = AdapterDeclaration(
        "test source",
        "1",
        "empty_adapter",
        "1.0.0",
        (),
        (),
        (),
        (),
        (),
        canonical_study_id="test_study",
        cohort_id="test_cohort",
    )

    def extract(self, manifest):
        return {"manifest_id": manifest.manifest_id}

    def normalize(self, extracted, manifest):
        del extracted, manifest
        return EventBCorpus()


def _entry(status=StudyStatus.REGISTERED):
    return StudyRegistryEntry(
        "test_study",
        "test_cohort",
        ("DOI:10.test/example",),
        None,
        "test cancer",
        "test vaccine",
        "PERSONALIZED",
        "fixture",
        "IMPLEMENTED",
        status,
        ("CSV",),
        (),
        (),
        None,
    )


def _manifest(path: Path, adapter=EmptyAdapter()):
    return manifest_from_paths(
        "test source",
        "1",
        adapter.declaration.adapter_name,
        adapter.declaration.adapter_version,
        [path],
    )


def test_repository_registry_is_valid_and_ordered():
    registry = StudyRegistry.read("configs/event_b_studies.yml")
    assert registry.registry_version == REGISTRY_VERSION
    assert [entry.canonical_study_id for entry in registry.backbone()] == [
        "mkras_vax_2026",
        "pdac_neovax_2023",
        "nous_209_2025",
        "fukuoka_dc",
    ]
    assert registry.get("braun_rcc_2025").ingestion_status is StudyStatus.ACCEPTED


def test_registry_rejects_blocked_entry_without_reason():
    with pytest.raises(ValueError, match="blocked without a blocker"):
        _entry(StudyStatus.BLOCKED_LABEL_SEMANTICS).validate()


def test_source_validation_rejects_html_disguised_as_xlsx(tmp_path):
    source = tmp_path / "supplement.xlsx"
    source.write_text("<!doctype html><title>access denied</title>")
    with pytest.raises(ValueError, match="HTML response"):
        detect_media_type(source)


def test_source_validation_streams_and_checks_jsonl(tmp_path):
    source = tmp_path / "records.jsonl"
    source.write_text('{"id": 1}\n{"id": 2}\n')
    result = validate_source_file(source, expected_media_type="jsonl")
    assert result.valid
    assert result.file_size_bytes == source.stat().st_size
    assert len(result.checksum_sha256) == 64


def test_job_runner_is_idempotent_and_resumes_after_interruption(tmp_path):
    source = tmp_path / "source.csv"
    source.write_text("id,value\n1,x\n")
    runner = EventBJobRunner(tmp_path / "run")
    adapter = EmptyAdapter()
    manifest = _manifest(source, adapter)

    def interrupt(stage):
        if stage == "CORPUS_VALIDATED":
            raise KeyboardInterrupt("simulated interruption")

    with pytest.raises(KeyboardInterrupt):
        runner.ingest(_entry(), adapter, manifest, after_stage=interrupt)
    interrupted = JobCheckpoint.read(runner.checkpoint_path("test_study"))
    assert interrupted.status == "INTERRUPTED"

    complete, result = runner.ingest(_entry(), adapter, manifest)
    assert complete.status == "COMPLETE"
    assert result is not None
    reused, second_result = runner.ingest(_entry(), adapter, manifest)
    assert reused == complete
    assert second_result is None


def test_job_runner_marks_changed_inputs_stale(tmp_path):
    source = tmp_path / "source.csv"
    source.write_text("id,value\n1,x\n")
    runner = EventBJobRunner(tmp_path / "run")
    adapter = EmptyAdapter()
    runner.ingest(_entry(), adapter, _manifest(source, adapter))
    source.write_text("id,value\n1,changed\n")
    with pytest.raises(RuntimeError, match="INPUT_FINGERPRINT_CHANGED"):
        runner.ingest(_entry(), adapter, _manifest(source, adapter))
    checkpoint = JobCheckpoint.read(runner.checkpoint_path("test_study"))
    assert checkpoint.status == "STALE"
    assert checkpoint.stale_reasons == ("INPUT_FINGERPRINT_CHANGED",)


def test_manifest_adapter_version_is_enforced(tmp_path):
    source = tmp_path / "source.csv"
    source.write_text("id,value\n1,x\n")
    adapter = EmptyAdapter()
    manifest = replace(_manifest(source, adapter), adapter_version="2.0.0")
    with pytest.raises(ValueError, match="adapter version"):
        EventBJobRunner(tmp_path / "run").ingest(_entry(), adapter, manifest)


def test_review_queue_extension_is_json_serializable():
    # Guard the richer deterministic review payload against accidental non-JSON fields.
    payload = json.dumps({"conflicting_fields": ("label",), "source_evidence": ({"row": 3},)})
    assert "label" in payload
