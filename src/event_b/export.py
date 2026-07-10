"""Deterministic model-ready exports without fitting a model."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from event_b.audit import corpus_audit, render_audit_markdown
from event_b.corpus import EventBCorpus
from event_b.evidence import evidence_availability_matrix
from event_b.manifest import SourceManifest
from event_b.review import ReviewIssue, write_review_queue
from event_b.schema import write_schema_bundle
from event_b.splits import SplitManifest


PARQUET_TABLES = {
    "studies": "studies.parquet",
    "patients": "patients.parquet",
    "vaccines": "vaccines.parquet",
    "candidates": "candidates.parquet",
    "assays": "assays.parquet",
    "clinical_outcomes": "clinical_outcomes.parquet",
    "recognition_evidence": "recognition_evidence.parquet",
    "candidate_funnel_links": "candidate_funnel_links.parquet",
    "provenance": "provenance.parquet",
}


def _write_parquet(frame: pd.DataFrame, path: Path) -> None:
    id_columns = [column for column in frame.columns if column.endswith("_id")]
    sort_columns = id_columns[:1] or list(frame.columns[:1])
    ordered = frame.sort_values(sort_columns, kind="mergesort", na_position="last").reset_index(
        drop=True
    )
    ordered.to_parquet(path, index=False, engine="pyarrow", compression="zstd")


def _event_labels(corpus: EventBCorpus) -> pd.DataFrame:
    columns = [
        "assay_id",
        "candidate_id",
        "patient_id",
        "study_id",
        "vaccine_id",
        "event_type",
        "response_label",
        "assay_type",
        "timepoint",
        "relative_to_vaccine",
        "review_status",
        "provenance_id",
        "schema_version",
    ]
    return corpus.assays.loc[:, columns].copy()


def _model_ready(corpus: EventBCorpus) -> pd.DataFrame:
    candidate_columns = [
        "candidate_id",
        "patient_id",
        "study_id",
        "mutant_peptide",
        "wildtype_peptide",
        "hla_alleles",
        "mhc_class",
        "vaccine_inclusion",
    ]
    assay_columns = [
        "assay_id",
        "candidate_id",
        "event_type",
        "response_label",
        "assay_type",
        "timepoint",
        "relative_to_vaccine",
        "review_status",
        "provenance_id",
    ]
    # One row per assay is deliberate: repeated assays and timepoints are never flattened away.
    model_ready = corpus.assays.loc[:, assay_columns].merge(
        corpus.candidates.loc[:, candidate_columns],
        on="candidate_id",
        how="left",
        validate="many_to_one",
    )
    availability = evidence_availability_matrix(corpus.recognition_evidence)
    if not availability.empty:
        model_ready = model_ready.merge(
            availability,
            on=["candidate_id", "patient_id"],
            how="left",
            validate="many_to_one",
        )
    return model_ready


def export_corpus(
    corpus: EventBCorpus,
    output_dir: str | Path,
    *,
    review_queue: tuple[ReviewIssue, ...] | list[ReviewIssue] = (),
    source_manifests: tuple[SourceManifest, ...] | list[SourceManifest] = (),
    split_manifests: tuple[SplitManifest, ...] | list[SplitManifest] = (),
    adapter_declarations=(),
) -> dict[str, str]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    written = {}
    for entity, path in write_schema_bundle(output / "schemas").items():
        written[f"schema:{entity}"] = path
    for entity, filename in PARQUET_TABLES.items():
        path = output / filename
        _write_parquet(getattr(corpus, entity), path)
        written[entity] = str(path)
    event_path = output / "event_labels.parquet"
    model_path = output / "model_ready_recognition.parquet"
    _write_parquet(_event_labels(corpus), event_path)
    _write_parquet(_model_ready(corpus), model_path)
    written["event_labels"] = str(event_path)
    written["model_ready_recognition"] = str(model_path)

    review_path = output / "review_queue.jsonl"
    write_review_queue(list(review_queue), review_path)
    written["review_queue"] = str(review_path)
    manifests_path = output / "source_manifest.json"
    manifests_path.write_text(
        json.dumps([manifest.to_dict() for manifest in source_manifests], indent=2, sort_keys=True)
        + "\n"
    )
    written["source_manifest"] = str(manifests_path)
    split_dir = output / "split_manifests"
    for manifest in sorted(split_manifests, key=lambda item: item.split_type):
        path = split_dir / f"{manifest.split_type.lower()}.json"
        manifest.write(path)
        written[f"split:{manifest.split_type}"] = str(path)
    audit = corpus_audit(corpus, review_queue, adapter_declarations)
    audit_json = output / "corpus_audit.json"
    audit_md = output / "corpus_audit.md"
    audit_json.write_text(json.dumps(audit, indent=2, sort_keys=True, default=str) + "\n")
    audit_md.write_text(render_audit_markdown(audit))
    written["corpus_audit_json"] = str(audit_json)
    written["corpus_audit_markdown"] = str(audit_md)
    return written
