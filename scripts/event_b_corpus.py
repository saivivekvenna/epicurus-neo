#!/usr/bin/env python3
"""Build/audit the Event-B substrate or emit pending extraction tasks."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
from pathlib import Path

from event_b.adapters import (
    BraunRCCAdapter,
    HuNeoVaxAdapter,
    ImproveEventAAdapter,
    MKRASVaxAdapter,
    Nous209Adapter,
    PDACNeoVaxAdapter,
)
from event_b.adapters import hu_neovax
from event_b.adapters import mkras_vax
from event_b.adapters import nous_209
from event_b.adapters import pdac_neovax
from event_b.adapters.braun_rcc import (
    DOI,
    EXPECTED_SHA256,
    NCT,
    PMCID,
    STUDY_ID,
    SUPPL_URL,
    braun_source_paths,
)
from event_b.audit import corpus_audit, render_audit_markdown
from event_b.braun_pipeline import (
    build_braun_corpus,
    combine_corpora,
    load_corpus_from_parquet,
    reconcile_braun,
    render_reconciliation_markdown,
)
from event_b.export import export_corpus
from event_b.extraction import ExtractionTask, emit_extraction_tasks
from event_b.factory import EventBJobRunner
from event_b.hu_pipeline import (
    assert_quality_gates as assert_hu_gates,
    build_hu_corpus,
    reconcile_hu,
    render_reconciliation_markdown as render_hu_reconciliation_markdown,
)
from event_b.ingest import ingest_source
from event_b.manifest import manifest_from_paths
from event_b.registry import StudyRegistry
from event_b.review import read_review_queue
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


def _fetch_record(raw_dir: Path) -> dict:
    return {
        "publication_id": f"DOI:{DOI}; {PMCID}",
        "trial_id": NCT,
        "source_url": SUPPL_URL,
        "retrieval_date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "files": [
            {
                "name": name,
                "sha256": digest,
                "media_type": Path(name).suffix.lstrip("."),
                "local_path": str((raw_dir / "extracted" / name).resolve()),
            }
            for name, digest in sorted(EXPECTED_SHA256.items())
        ],
    }


def cmd_braun(args) -> int:
    raw_dir = args.raw_dir
    build = build_braun_corpus(raw_dir)
    accepted = build.result.accepted_corpus

    split_input = accepted.candidates.loc[
        :, ["candidate_id", "patient_id", "study_id", "mutant_peptide", "hla_alleles"]
    ].copy()
    braun_split = generate_split_manifest(split_input, SplitType.PATIENT_HOLDOUT)

    braun_out = Path(args.output)
    written = export_corpus(
        accepted,
        braun_out,
        review_queue=build.review_queue,
        source_manifests=[build.manifest],
        split_manifests=[braun_split],
        adapter_declarations=[build.adapter.declaration],
    )
    recon = reconcile_braun(raw_dir, build)
    recon_md = render_reconciliation_markdown(recon)
    (braun_out / "braun_reconciliation.md").write_text(recon_md)
    (braun_out / "fetch_record.json").write_text(
        json.dumps(_fetch_record(Path(raw_dir)), indent=2, sort_keys=True) + "\n"
    )

    milestone = Path(args.milestone_dir)
    milestone.mkdir(parents=True, exist_ok=True)
    (milestone / "braun_reconciliation.md").write_text(recon_md)
    (milestone / "braun_reconciliation.json").write_text(
        json.dumps(recon, indent=2, sort_keys=True, default=str) + "\n"
    )

    combined_summary: dict[str, object] = {"built": False}
    improve_dir = Path(args.improve_corpus)
    if (improve_dir / "assays.parquet").exists():
        improve = load_corpus_from_parquet(improve_dir)
        combined = combine_corpora(improve, accepted)
        combined_input = (
            combined.candidates.loc[
                :, ["candidate_id", "patient_id", "study_id", "mutant_peptide", "hla_alleles"]
            ]
            .drop_duplicates(subset="candidate_id")
            .reset_index(drop=True)
        )
        split_manifests = [generate_split_manifest(combined_input, SplitType.PATIENT_HOLDOUT)]
        try:
            split_manifests.append(
                generate_split_manifest(combined_input, SplitType.STUDY_HOLDOUT)
            )
        except ValueError:
            pass
        combined_issues = tuple(read_review_queue(improve_dir / "review_queue.jsonl")) + tuple(
            build.review_queue
        )
        combined_written = export_corpus(
            combined,
            args.combined_output,
            review_queue=combined_issues,
            source_manifests=[build.manifest],
            split_manifests=split_manifests,
            adapter_declarations=[
                ImproveEventAAdapter.declaration,
                build.adapter.declaration,
            ],
        )
        audit = corpus_audit(
            combined,
            combined_issues,
            [ImproveEventAAdapter.declaration, build.adapter.declaration],
        )
        (milestone / "combined_corpus_audit.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True, default=str) + "\n"
        )
        (milestone / "combined_corpus_audit.md").write_text(render_audit_markdown(audit))
        combined_summary = {
            "built": True,
            "output": str(args.combined_output),
            "written_tables": len(combined_written),
            "sample_sizes": audit["sample_sizes"],
            "event_counts": audit["event_counts"],
            "model_readiness": audit["model_readiness"],
        }

    print(
        json.dumps(
            {
                "study_id": STUDY_ID,
                "braun_output": str(braun_out),
                "written_tables": len(written),
                "accepted": recon["accepted"],
                "review_queue": recon["review_queue"],
                "summary_reconciles": recon["summary_reconciles"],
                "combined": combined_summary,
            },
            indent=2,
            default=str,
        )
    )
    return 0


def _hu_fetch_record(raw_dir: Path) -> dict:
    return {
        "publication_id": f"DOI:{hu_neovax.DOI}; {hu_neovax.PMCID}",
        "trial_id": hu_neovax.NCT,
        "acquisition": "manual (PMC author-manuscript supplement behind a proof-of-work gate)",
        "source_url": hu_neovax.DOWNLOAD_URL,
        "cross_check_source": hu_neovax.OTT_PMCID,
        "retrieval_date": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "files": [
            {
                "name": name,
                "sha256": digest,
                "media_type": Path(name).suffix.lstrip("."),
                "local_path": str(hu_neovax.hu_source_path(raw_dir).resolve()),
            }
            for name, digest in sorted(hu_neovax.EXPECTED_SHA256.items())
        ],
    }


def _split_input(corpus) -> "object":
    return (
        corpus.candidates.loc[
            :, ["candidate_id", "patient_id", "study_id", "mutant_peptide", "hla_alleles"]
        ]
        .drop_duplicates(subset="candidate_id")
        .reset_index(drop=True)
    )


def cmd_hu(args) -> int:
    raw_dir = args.raw_dir
    build = build_hu_corpus(raw_dir)
    accepted = build.result.accepted_corpus
    assert_hu_gates(accepted)

    hu_split = generate_split_manifest(_split_input(accepted), SplitType.PATIENT_HOLDOUT)
    hu_out = Path(args.output)
    written = export_corpus(
        accepted,
        hu_out,
        review_queue=build.review_queue,
        source_manifests=[build.manifest],
        split_manifests=[hu_split],
        adapter_declarations=[build.adapter.declaration],
    )
    recon = reconcile_hu(build)
    recon_md = render_hu_reconciliation_markdown(recon)
    (hu_out / "hu_reconciliation.md").write_text(recon_md)
    (hu_out / "fetch_record.json").write_text(
        json.dumps(_hu_fetch_record(Path(raw_dir)), indent=2, sort_keys=True) + "\n"
    )

    milestone = Path(args.milestone_dir)
    milestone.mkdir(parents=True, exist_ok=True)
    (milestone / "hu_reconciliation.md").write_text(recon_md)
    (milestone / "hu_reconciliation.json").write_text(
        json.dumps(recon, indent=2, sort_keys=True, default=str) + "\n"
    )

    # Combine with whichever prior corpora are materialised: IMPROVE (Event-A) and Braun (Event-B).
    corpora = []
    declarations = []
    combined_issues: tuple = tuple(build.review_queue)
    improve_dir = Path(args.improve_corpus)
    if (improve_dir / "assays.parquet").exists():
        corpora.append(load_corpus_from_parquet(improve_dir))
        declarations.append(ImproveEventAAdapter.declaration)
        combined_issues = tuple(read_review_queue(improve_dir / "review_queue.jsonl")) + combined_issues
    braun_dir = Path(args.braun_corpus)
    if (braun_dir / "assays.parquet").exists():
        corpora.append(load_corpus_from_parquet(braun_dir))
        declarations.append(BraunRCCAdapter.declaration)
        combined_issues = tuple(read_review_queue(braun_dir / "review_queue.jsonl")) + combined_issues

    combined_summary: dict[str, object] = {"built": False}
    if corpora:
        combined = combine_corpora(*corpora, accepted)
        declarations.append(build.adapter.declaration)
        combined_input = _split_input(combined)
        split_manifests = [generate_split_manifest(combined_input, SplitType.PATIENT_HOLDOUT)]
        try:
            split_manifests.append(generate_split_manifest(combined_input, SplitType.STUDY_HOLDOUT))
        except ValueError:
            pass
        combined_written = export_corpus(
            combined,
            args.combined_output,
            review_queue=combined_issues,
            source_manifests=[build.manifest],
            split_manifests=split_manifests,
            adapter_declarations=declarations,
        )
        audit = corpus_audit(combined, combined_issues, declarations)
        (milestone / "combined_corpus_audit.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True, default=str) + "\n"
        )
        (milestone / "combined_corpus_audit.md").write_text(render_audit_markdown(audit))
        combined_summary = {
            "built": True,
            "output": str(args.combined_output),
            "written_tables": len(combined_written),
            "sample_sizes": audit["sample_sizes"],
            "event_counts": audit["event_counts"],
            "model_readiness": audit["model_readiness"],
        }

    print(
        json.dumps(
            {
                "study_id": hu_neovax.STUDY_ID,
                "hu_output": str(hu_out),
                "written_tables": len(written),
                "reconciliation": recon["reconciliation"],
                "accepted": recon["accepted"],
                "epitope_spreading": recon["epitope_spreading"],
                "review_queue": recon["review_queue"],
                "combined": combined_summary,
            },
            indent=2,
            default=str,
        )
    )
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


def _factory_adapter(study_id: str, raw_root: Path):
    if study_id == STUDY_ID:
        adapter = BraunRCCAdapter(raw_root / "braun_rcc_2025")
        paths = braun_source_paths(raw_root / "braun_rcc_2025")
    elif study_id == hu_neovax.STUDY_ID:
        adapter = HuNeoVaxAdapter(raw_root / "hu_melanoma_2021")
        paths = hu_neovax.hu_source_paths(raw_root / "hu_melanoma_2021")
    elif study_id == mkras_vax.STUDY_ID:
        adapter = MKRASVaxAdapter(raw_root / "mkras_vax_2026")
        paths = mkras_vax.stage_sources(raw_root / "mkras_vax_2026")
    elif study_id == pdac_neovax.STUDY_ID:
        adapter = PDACNeoVaxAdapter(raw_root / "pdac_neovax_2023")
        paths = [pdac_neovax.stage_source(raw_root / "pdac_neovax_2023")]
    elif study_id == nous_209.STUDY_ID:
        adapter = Nous209Adapter(raw_root / "nous_209_2025")
        paths = nous_209.stage_sources(raw_root / "nous_209_2025")
    else:
        return None
    return adapter, manifest_from_paths(
        adapter.declaration.source_name,
        adapter.declaration.source_version,
        adapter.declaration.adapter_name,
        adapter.declaration.adapter_version,
        paths,
    )


def _checkpoint_summary(checkpoint) -> dict:
    return {
        "study_id": checkpoint.study_id,
        "status": checkpoint.status,
        "stage": checkpoint.stage,
        "output_dir": checkpoint.output_dir,
        "stale_reasons": list(checkpoint.stale_reasons),
        "error": checkpoint.error,
    }


def cmd_ingest_study(args) -> int:
    registry = StudyRegistry.read(args.registry)
    entry = registry.get(args.study_id)
    runner = EventBJobRunner(args.output_root)
    resolved = _factory_adapter(entry.canonical_study_id, args.raw_root)
    if resolved is None:
        if entry.ingestion_status.value.startswith("BLOCKED_"):
            checkpoint = runner.block(entry)
            print(json.dumps(_checkpoint_summary(checkpoint), indent=2, sort_keys=True))
            return 0
        raise RuntimeError(
            f"{entry.canonical_study_id} has no implemented adapter; "
            f"registry status is {entry.ingestion_status.value}"
        )
    adapter, manifest = resolved
    checkpoint, result = runner.ingest(entry, adapter, manifest, rebuild=args.rebuild)
    summary = _checkpoint_summary(checkpoint)
    summary["reused"] = result is None
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


def cmd_ingest_backbone(args) -> int:
    registry = StudyRegistry.read(args.registry)
    summaries = []
    failures = []
    for entry in registry.backbone():
        namespace = argparse.Namespace(
            registry=args.registry,
            study_id=entry.canonical_study_id,
            output_root=args.output_root,
            raw_root=args.raw_root,
            rebuild=args.rebuild,
        )
        try:
            resolved = _factory_adapter(entry.canonical_study_id, args.raw_root)
            runner = EventBJobRunner(args.output_root)
            if resolved is None:
                if entry.ingestion_status.value.startswith("BLOCKED_"):
                    summaries.append(_checkpoint_summary(runner.block(entry)))
                else:
                    failures.append(
                        {
                            "study_id": entry.canonical_study_id,
                            "status": entry.ingestion_status.value,
                            "error": "adapter not implemented",
                        }
                    )
                continue
            adapter, manifest = resolved
            checkpoint, result = runner.ingest(
                entry, adapter, manifest, rebuild=namespace.rebuild
            )
            row = _checkpoint_summary(checkpoint)
            row["reused"] = result is None
            summaries.append(row)
        except Exception as error:  # continue to the next registered study
            failures.append(
                {
                    "study_id": entry.canonical_study_id,
                    "error": f"{type(error).__name__}: {error}",
                }
            )
    payload = {"studies": summaries, "failures": failures}
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 1 if failures else 0


def main() -> int:
    parser = argparse.ArgumentParser(prog="event-b-corpus")
    sub = parser.add_subparsers(dest="command", required=True)
    improve = sub.add_parser("import-improve-event-a")
    improve.add_argument("improve_repo", type=Path)
    improve.add_argument("output", type=Path)
    improve.set_defaults(func=cmd_improve)
    braun = sub.add_parser("import-braun-rcc")
    braun.add_argument("--raw-dir", type=Path, default=Path("data/raw/braun_rcc_2025"))
    braun.add_argument("--output", type=Path, default=Path("outputs/event_b_braun"))
    braun.add_argument("--improve-corpus", type=Path, default=Path("outputs/event_b_corpus"))
    braun.add_argument(
        "--combined-output", type=Path, default=Path("outputs/event_b_corpus_combined")
    )
    braun.add_argument("--milestone-dir", type=Path, default=Path("artifacts/milestone_5a"))
    braun.set_defaults(func=cmd_braun)
    hu = sub.add_parser("import-hu-neovax")
    hu.add_argument("--raw-dir", type=Path, default=Path("data/raw/hu_melanoma_2021"))
    hu.add_argument("--output", type=Path, default=Path("outputs/event_b_hu"))
    hu.add_argument("--improve-corpus", type=Path, default=Path("outputs/event_b_corpus"))
    hu.add_argument("--braun-corpus", type=Path, default=Path("outputs/event_b_braun"))
    hu.add_argument("--combined-output", type=Path, default=Path("outputs/event_b_corpus_combined"))
    hu.add_argument("--milestone-dir", type=Path, default=Path("artifacts/milestone_5b1"))
    hu.set_defaults(func=cmd_hu)
    tasks = sub.add_parser("emit-extraction-tasks")
    tasks.add_argument("--study-id", required=True)
    tasks.add_argument("--source", action="append", type=Path, required=True)
    tasks.add_argument("--output", type=Path, required=True)
    tasks.set_defaults(func=cmd_emit_tasks)
    factory = sub.add_parser("ingest-study")
    factory.add_argument("study_id")
    factory.add_argument("--registry", type=Path, default=Path("configs/event_b_studies.yml"))
    factory.add_argument("--raw-root", type=Path, default=Path("data/raw"))
    factory.add_argument("--output-root", type=Path, default=Path("outputs/event_b_backbone"))
    factory.add_argument("--rebuild", action="store_true")
    factory.set_defaults(func=cmd_ingest_study)
    backbone = sub.add_parser("ingest-backbone")
    backbone.add_argument("--registry", type=Path, default=Path("configs/event_b_studies.yml"))
    backbone.add_argument("--raw-root", type=Path, default=Path("data/raw"))
    backbone.add_argument("--output-root", type=Path, default=Path("outputs/event_b_backbone"))
    backbone.add_argument("--resume", action="store_true", help="reuse checksum-matching outputs")
    backbone.add_argument("--rebuild", action="store_true")
    backbone.set_defaults(func=cmd_ingest_backbone)
    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
