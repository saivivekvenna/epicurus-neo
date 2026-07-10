#!/usr/bin/env python3
"""Build/audit the Event-B substrate or emit pending extraction tasks."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

from event_b.adapters import ImproveEventAAdapter
from event_b.export import export_corpus
from event_b.extraction import ExtractionTask, emit_extraction_tasks
from event_b.ingest import ingest_source
from event_b.manifest import manifest_from_paths
from event_b.splits import SplitType, generate_split_manifest


def cmd_improve(args) -> int:
    data_zip = args.improve_repo / "data.zip"
    adapter = ImproveEventAAdapter(data_zip)
    manifest = manifest_from_paths(
        adapter.declaration.source_name,
        adapter.declaration.source_version,
        adapter.declaration.adapter_name,
        adapter.declaration.adapter_version,
        [data_zip],
    )
    result = ingest_source(adapter, manifest)
    split_input = result.accepted_corpus.candidates.loc[
        :, ["candidate_id", "patient_id", "study_id", "mutant_peptide", "hla_alleles"]
    ].copy()
    split_manifest = generate_split_manifest(split_input, SplitType.PATIENT_HOLDOUT)
    written = export_corpus(
        result.accepted_corpus,
        args.output,
        review_queue=result.review_queue,
        source_manifests=[manifest],
        split_manifests=[split_manifest],
        adapter_declarations=[adapter.declaration],
    )
    print(json.dumps({"written": written, "review_queue_n": len(result.review_queue)}, indent=2))
    return 0


def cmd_emit_tasks(args) -> int:
    tasks = []
    for source in sorted(args.source):
        raw = source.read_bytes()
        try:
            text = raw.decode("utf-8")
        except UnicodeDecodeError as error:
            raise ValueError(
                f"{source} is not UTF-8 text; extract text before task emission"
            ) from error
        tasks.append(
            ExtractionTask.create(
                study_id=args.study_id,
                source_document=str(source.resolve()),
                source_checksum=sha256(raw).hexdigest(),
                source_text=text,
            )
        )
    task_path, schema_path = emit_extraction_tasks(tasks, args.output)
    print(
        json.dumps(
            {"tasks": str(task_path), "schema": str(schema_path), "status": "PENDING"}, indent=2
        )
    )
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(prog="event-b-corpus")
    sub = parser.add_subparsers(dest="command", required=True)
    improve = sub.add_parser("import-improve-event-a")
    improve.add_argument("improve_repo", type=Path)
    improve.add_argument("output", type=Path)
    improve.set_defaults(func=cmd_improve)
    tasks = sub.add_parser("emit-extraction-tasks")
    tasks.add_argument("--study-id", required=True)
    tasks.add_argument("--source", action="append", type=Path, required=True)
    tasks.add_argument("--output", type=Path, required=True)
    tasks.set_defaults(func=cmd_emit_tasks)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
