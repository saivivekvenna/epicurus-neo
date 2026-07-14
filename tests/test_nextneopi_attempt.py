from __future__ import annotations

import csv
import hashlib
import json
import shutil
from pathlib import Path

import pytest

from benchmark.nextneopi_attempt import (
    VcfAllele,
    finalize_non_success,
    finalize_success,
    normalize_variant_rows,
    verify_attempt_manifest,
)


@pytest.mark.parametrize(("chrom", "pos", "ref", "alt", "expected"), [
    ("chrX", 41344255, "G", "C", "chrX-41344254-41344255-G-C"),
    ("chr1", 227091906, "GC", "G", "chr1-227091906-227091907-GC-G"),
    ("chr15", 44767453, "T", "TC", "chr15-44767453-44767453-T-TC"),
])
def test_pvac_id_matches_official_pvactools_snv_deletion_insertion_examples(
    chrom, pos, ref, alt, expected
):
    assert VcfAllele(chrom, pos, ref, alt, 1).pvac_id == expected


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _evidence(tmp_path: Path, *, variants: int = 21) -> dict:
    tmp_path.mkdir(parents=True, exist_ok=True)
    source = tmp_path / "nextNEOpi.nf"
    source.write_text("workflow { nextNEOpi() }\n")
    patch = tmp_path / "instrumentation.patch"
    patch.write_text("publish native evidence\n")
    config = tmp_path / "config.json"
    frozen = {
        "policy_id": "nextneopi-track-a-native-aggregate-v1",
        "upstream": {"commit": "a" * 40, "nextneopi_nf_sha256": _sha(source)},
        "runtime": {"instrumentation_patch_sha256": _sha(patch)},
    }
    config.write_text(json.dumps(frozen, sort_keys=True))
    input_manifest = tmp_path / "INPUT_MANIFEST.json"
    input_manifest.write_text(json.dumps({
        "policy_id": frozen["policy_id"], "patient_id": "Hu_test", "labels_opened": False,
        "command": "nextflow run nextNEOpi.nf --locked-inputs",
        "config": {"sha256": _sha(config)},
        "upstream": {"nextneopi_nf_sha256": _sha(source)},
        "instrumentation_patch": {"sha256": _sha(patch)},
    }, sort_keys=True))
    log = tmp_path / "nextflow.log"
    log.write_text("execution complete\n")
    trace = tmp_path / "trace.txt"
    trace.write_text("task_id\tstatus\n1\tCOMPLETED\n")
    reference = tmp_path / "GRCh38.fa"
    reference.write_text(">chr1\n" + "A" * 1000 + "\n")
    vcf = tmp_path / "benchmark_pvac_input.vcf"
    with vcf.open("w") as handle:
        handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for index in range(variants):
            handle.write(f"chr1\t{index + 10}\t.\tA\tC\t.\tPASS\t.\n")
    index = tmp_path / "benchmark_pvac_input.vcf.tbi"
    index.write_bytes(b"index evidence")
    aggregate = tmp_path / "native.all_epitopes.aggregated.tsv"
    with aggregate.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t",
            fieldnames=("ID", "Best Peptide", "Allele", "Tier"),
        )
        writer.writeheader()
        for idx in range(variants):
            pos = idx + 10
            writer.writerow({
                "ID": f"chr1-{pos - 1}-{pos}-A-C", "Best Peptide": f"PEPTIDE{idx}",
                "Allele": "HLA-A*01:01", "Tier": "Pass",
            })
    return locals()


def _copy_normalizer(source: Path, reference: Path, destination: Path) -> None:
    assert reference.is_file()
    shutil.copyfile(source, destination)


def _success_args(e: dict, tmp_path: Path) -> dict:
    return {
        "patient_id": "Hu_test", "input_manifest": e["input_manifest"],
        "source": e["source"], "config": e["config"], "instrumentation": e["patch"],
        "aggregate": e["aggregate"], "pvac_input_vcf": e["vcf"],
        "pvac_input_vcf_index": e["index"], "reference": e["reference"],
        "execution_log": e["log"], "execution_trace": e["trace"],
        "execution_command": "nextflow run nextNEOpi.nf --locked-inputs",
        "started_at": "2026-07-13T10:00:00-06:00",
        "finished_at": "2026-07-13T11:30:00-06:00",
        "runtime": "nextflow-22.10.8/singularity/linux-amd64",
        "output_dir": tmp_path / "attempt", "normalizer": _copy_normalizer,
    }


def test_success_binds_provenance_maps_every_id_and_freezes_canonical_top20(tmp_path):
    evidence = _evidence(tmp_path)
    result = finalize_success(**_success_args(evidence, tmp_path))

    assert result["comparator_id"] == "nextneopi_track_a"
    assert result["scope"] == "FULL_PIPELINE_IDENTICAL_RAW_INPUT"
    assert result["patient_id"] == "Hu_test" and result["labels_opened"] is False
    assert result["status"] == "SUCCEEDED"
    assert result["execution"]["exit_code"] == 0
    assert result["input_manifest"]["sha256"] == _sha(evidence["input_manifest"])
    assert result["native_aggregate"]["sha256"] == _sha(evidence["aggregate"])
    assert result["published_pvac_input"]["vcf"]["sha256"] == _sha(evidence["vcf"])
    assert result["published_pvac_input"]["index"]["sha256"] == _sha(evidence["index"])
    assert result["grch38_reference"]["sha256"] == _sha(evidence["reference"])
    assert result["normalization"]["mapped_aggregate_ids"] == 21
    assert result["portfolio"]["size"] == 20
    assert (tmp_path / "attempt" / "ATTEMPT_MANIFEST.json").is_file()

    with (tmp_path / "attempt" / "NORMALIZATION_MAP.tsv").open() as handle:
        mapping = list(csv.DictReader(handle, delimiter="\t"))
    assert len(mapping) == 21
    assert mapping[0]["canonical_variant"] == "chr1:10:A:C"
    assert mapping[19]["selected_top20"] == "true"
    assert mapping[20]["selected_top20"] == "false"
    with (tmp_path / "attempt" / "NEXTNEOPI_MUTATION_TOP20.tsv").open() as handle:
        top = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["rank"] for row in top] == [str(i) for i in range(1, 21)]
    assert [row["canonical_mutation"] for row in top[:2]] == ["chr1:10:A:C", "chr1:11:A:C"]


def test_success_rejects_missing_or_ambiguous_aggregate_mapping(tmp_path):
    evidence = _evidence(tmp_path / "missing", variants=1)
    evidence["aggregate"].write_text(
        "ID\tBest Peptide\tAllele\tTier\nchr1-999-1000-A-C\tPEP\tHLA-A*01:01\tPass\n"
    )
    with pytest.raises(ValueError, match="one-to-one"):
        finalize_success(**_success_args(evidence, tmp_path / "missing"))

    evidence = _evidence(tmp_path / "ambiguous", variants=1)
    with evidence["vcf"].open("a") as handle:
        handle.write("chr1\t10\tduplicate\tA\tC\t.\tPASS\t.\n")
    with pytest.raises(ValueError, match="one-to-one"):
        finalize_success(**_success_args(evidence, tmp_path / "ambiguous"))


def test_success_rejects_canonical_duplicates_from_normalization(tmp_path):
    evidence = _evidence(tmp_path, variants=2)

    def duplicate_normalizer(source: Path, reference: Path, destination: Path) -> None:
        lines = source.read_text().splitlines()
        variants = [line.split("\t") for line in lines if not line.startswith("#")]
        variants[1][0:2] = variants[0][0:2]
        destination.write_text(
            "\n".join([line for line in lines if line.startswith("#")] +
                      ["\t".join(row) for row in variants]) + "\n"
        )

    args = _success_args(evidence, tmp_path)
    args["normalizer"] = duplicate_normalizer
    with pytest.raises(ValueError, match="canonical duplicate"):
        finalize_success(**args)


def test_success_rejects_tampered_pin_and_invalid_execution_claim(tmp_path):
    evidence = _evidence(tmp_path)
    evidence["source"].write_text("tampered\n")
    with pytest.raises(ValueError, match="source hash"):
        finalize_success(**_success_args(evidence, tmp_path))

    evidence = _evidence(tmp_path / "exit")
    with pytest.raises(ValueError, match="exit_code 0"):
        finalize_success(**_success_args(evidence, tmp_path / "exit"), exit_code=1)


@pytest.mark.parametrize(("status", "exit_code"), [("FAILED", 17), ("ABSTAINED", None)])
def test_verified_non_success_is_evaluable_zero_portfolio(tmp_path, status, exit_code):
    evidence = _evidence(tmp_path)
    result = finalize_non_success(
        patient_id="Hu_test", input_manifest=evidence["input_manifest"],
        source=evidence["source"], config=evidence["config"],
        instrumentation=evidence["patch"], execution_log=evidence["log"],
        execution_trace=evidence["trace"], output_dir=tmp_path / "attempt",
        execution_status=status, exit_code=exit_code, reason="verified runtime limitation",
        execution_command="nextflow run nextNEOpi.nf --locked-inputs",
        started_at="2026-07-13T10:00:00-06:00",
        finished_at="2026-07-13T10:00:01-06:00", runtime="singularity/linux-amd64",
    )
    assert result["status"] == status and result["evaluable"] is True
    assert result["portfolio"]["size"] == 0
    assert result["native_aggregate"] is None
    assert (tmp_path / "attempt" / "ATTEMPT_MANIFEST.json").is_file()


def test_non_success_still_rejects_missing_or_tampered_provenance(tmp_path):
    evidence = _evidence(tmp_path)
    evidence["trace"].unlink()
    with pytest.raises(FileNotFoundError):
        finalize_non_success(
            patient_id="Hu_test", input_manifest=evidence["input_manifest"],
            source=evidence["source"], config=evidence["config"],
            instrumentation=evidence["patch"], execution_log=evidence["log"],
            execution_trace=evidence["trace"], output_dir=tmp_path / "attempt",
            execution_status="FAILED", exit_code=1, reason="failed",
            execution_command="nextflow run nextNEOpi.nf --locked-inputs",
            started_at="2026-07-13T10:00:00-06:00",
            finished_at="2026-07-13T10:00:01-06:00", runtime="singularity/linux-amd64",
        )


def test_aggregate_mapping_uses_each_alt_of_multiallelic_input_exactly_once(tmp_path):
    evidence = _evidence(tmp_path, variants=1)
    evidence["vcf"].write_text(
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr1\t10\t.\tA\tC,G\t.\tPASS\t.\n"
    )
    evidence["aggregate"].write_text(
        "ID\tBest Peptide\tAllele\tTier\n"
        "chr1-9-10-A-G\tPEP\tHLA-A*01:01\tPass\n"
    )
    result = finalize_success(**_success_args(evidence, tmp_path))
    assert result["portfolio"]["size"] == 1
    with (tmp_path / "attempt" / "NEXTNEOPI_MUTATION_TOP20.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["canonical_mutation"] == "chr1:10:A:G"


def test_verify_success_returns_ordered_snapshot_and_detects_artifact_tampering(tmp_path):
    evidence = _evidence(tmp_path)
    finalize_success(**_success_args(evidence, tmp_path))
    manifest = tmp_path / "attempt" / "ATTEMPT_MANIFEST.json"
    ok, snapshot = verify_attempt_manifest(manifest, "Hu_test")
    assert ok is True
    assert snapshot["status"] == "SUCCEEDED"
    assert snapshot["mutation_ids"][:2] == ["chr1:10:A:C", "chr1:11:A:C"]
    assert len(snapshot["mutation_ids"]) == 20

    (tmp_path / "attempt" / "NEXTNEOPI_MUTATION_TOP20.tsv").write_text("tampered\n")
    ok, detail = verify_attempt_manifest(manifest, "Hu_test")
    assert ok is False and "mismatch" in detail["reason"]


def test_verify_non_success_returns_empty_ordered_snapshot(tmp_path):
    evidence = _evidence(tmp_path)
    finalize_non_success(
        patient_id="Hu_test", input_manifest=evidence["input_manifest"],
        source=evidence["source"], config=evidence["config"],
        instrumentation=evidence["patch"], execution_log=evidence["log"],
        execution_trace=evidence["trace"], output_dir=tmp_path / "attempt",
        execution_status="FAILED", exit_code=9, reason="verified failure",
        execution_command="nextflow run nextNEOpi.nf --locked-inputs",
        started_at="2026-07-13T10:00:00-06:00",
        finished_at="2026-07-13T10:00:01-06:00", runtime="singularity/linux-amd64",
    )
    ok, snapshot = verify_attempt_manifest(
        tmp_path / "attempt" / "ATTEMPT_MANIFEST.json", "Hu_test"
    )
    assert ok is True and snapshot["status"] == "FAILED"
    assert snapshot["mutation_ids"] == [] and snapshot["evaluable"] is True


def test_verify_rejects_wrong_identity_path_escape_and_symlink(tmp_path):
    evidence = _evidence(tmp_path)
    finalize_success(**_success_args(evidence, tmp_path))
    manifest_path = tmp_path / "attempt" / "ATTEMPT_MANIFEST.json"
    ok, detail = verify_attempt_manifest(manifest_path, "wrong_patient")
    assert ok is False and "patient_id" in detail["reason"]

    manifest = json.loads(manifest_path.read_text())
    manifest["portfolio"]["path"] = str(evidence["log"].resolve())
    manifest["portfolio"]["bytes"] = evidence["log"].stat().st_size
    manifest["portfolio"]["sha256"] = _sha(evidence["log"])
    manifest_path.chmod(0o644)  # simulate an adversary who deliberately breaks the immutable bit
    manifest_path.write_text(json.dumps(manifest))
    ok, detail = verify_attempt_manifest(manifest_path, "Hu_test")
    assert ok is False and "escapes attempt directory" in detail["reason"]

    link = tmp_path / "manifest-link.json"
    link.symlink_to(manifest_path)
    ok, detail = verify_attempt_manifest(link, "Hu_test")
    assert ok is False


def test_finalization_rejects_symlink_evidence_and_never_overwrites_attempt(tmp_path):
    evidence = _evidence(tmp_path)
    source_link = tmp_path / "nextNEOpi-link.nf"
    source_link.symlink_to(evidence["source"])
    args = _success_args(evidence, tmp_path)
    args["source"] = source_link
    with pytest.raises(ValueError, match="symlink"):
        finalize_success(**args)

    finalize_success(**_success_args(evidence, tmp_path))
    manifest = tmp_path / "attempt" / "ATTEMPT_MANIFEST.json"
    original = manifest.read_bytes()
    assert manifest.stat().st_mode & 0o222 == 0
    with pytest.raises(FileExistsError, match="nonempty"):
        finalize_success(**_success_args(evidence, tmp_path))
    assert manifest.read_bytes() == original


def test_normalize_variant_rows_preserves_ids_order_and_mnvs(tmp_path):
    reference = tmp_path / "GRCh38.fa"
    reference.write_text(">chr1\n" + "A" * 100 + "\n")
    rows = [
        {"row_id": "assay-2", "chrom": "chr1", "pos": 20, "ref": "AC", "alt": "GT"},
        {"row_id": "assay-1", "chrom": "chr1", "pos": 10, "ref": "A", "alt": "C"},
    ]
    normalized = normalize_variant_rows(rows, reference, _copy_normalizer)
    assert [row["row_id"] for row in normalized] == ["assay-2", "assay-1"]
    assert normalized[0]["canonical_mutation"] == "chr1:20:AC:GT"
    assert len(normalized) == 2  # the MNV was not decomposed/rescued by components
    assert "canonical_mutation" not in rows[0]


def test_normalize_variant_rows_fails_closed_on_duplicate_ids_or_canonical_collapse(tmp_path):
    reference = tmp_path / "GRCh38.fa"
    reference.write_text(">chr1\n" + "A" * 100 + "\n")
    duplicate_ids = [
        {"row_id": "same", "chrom": "chr1", "pos": 10, "ref": "A", "alt": "C"},
        {"row_id": "same", "chrom": "chr1", "pos": 11, "ref": "A", "alt": "G"},
    ]
    with pytest.raises(ValueError, match="IDs must be unique"):
        normalize_variant_rows(duplicate_ids, reference, _copy_normalizer)

    duplicate_variants = [
        {"row_id": "one", "chrom": "chr1", "pos": 10, "ref": "A", "alt": "C"},
        {"row_id": "two", "chrom": "chr1", "pos": 10, "ref": "A", "alt": "C"},
    ]
    with pytest.raises(ValueError, match="canonical duplicate"):
        normalize_variant_rows(duplicate_variants, reference, _copy_normalizer)
