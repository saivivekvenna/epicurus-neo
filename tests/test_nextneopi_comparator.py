from __future__ import annotations

import csv
from pathlib import Path
import hashlib
import json

import pytest

from benchmark.nextneopi_comparator import freeze_native_portfolio, prepare_batch


def _fastq(path: Path) -> Path:
    path.write_bytes(b"@r\nAC\n+\n!!\n")
    return path


def test_prepare_batch_hashes_three_raw_sample_types(tmp_path, monkeypatch):
    monkeypatch.chdir(Path(__file__).parents[1])
    fastqs = [_fastq(tmp_path / f"r{i}.fastq") for i in range(6)]
    provenance = tmp_path / "CONVERT_PROVENANCE.json"
    provenance.write_text(json.dumps([
        {
            "status": "OK",
            "raw_input_identity_ready": True,
            "sha256_per_file": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()},
            "bytes_per_file": {path.name: path.stat().st_size},
        }
        for path in fastqs
    ]))
    manifest = prepare_batch(
        patient_id="Hu_test",
        tumor_r1=fastqs[0], tumor_r2=fastqs[1],
        normal_r1=fastqs[2], normal_r2=fastqs[3],
        rna_r1=fastqs[4], rna_r2=fastqs[5],
        conversion_provenance=provenance,
        output_dir=tmp_path / "out",
    )
    assert list(manifest["raw_fastqs"]) == ["tumor_DNA", "normal_DNA", "tumor_RNA"]
    assert manifest["labels_opened"] is False
    assert manifest["hla_file_supplied"] is False
    assert manifest["epicurus_conversion_provenance"]["identity_match"] is True
    assert manifest["upstream"]["nextneopi_nf_sha256"] == (
        "bc55f844133084366db6f23584df55c9bdaaa6bc478fd160a825b54da65f8c66"
    )
    assert manifest["instrumentation_patch"]["sha256"] == (
        "3a01cb41f7400f70e2319ac9dab4fc45443f6bdc928dba5e3ebed5246ebc4456"
    )
    with (tmp_path / "out" / "nextneopi_batch.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert [row["sampleType"] for row in rows] == ["tumor_DNA", "normal_DNA", "tumor_RNA"]
    assert all(row["HLAfile"] == "" for row in rows)


def test_freeze_native_portfolio_preserves_first_twenty(tmp_path):
    aggregate = tmp_path / "native.tsv"
    with aggregate.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("ID", "Best Peptide", "Allele", "Tier", "IC50 MT"),
            delimiter="\t",
        )
        writer.writeheader()
        for index in range(25):
            writer.writerow(
                {"ID": f"v{index}", "Best Peptide": "AAAAAAAA", "Allele": "A*01:01",
                 "Tier": "Pass", "IC50 MT": index}
            )
    manifest = freeze_native_portfolio(aggregate, tmp_path / "out")
    assert manifest["native_aggregate"]["row_count"] == 25
    assert manifest["portfolio"]["size"] == 20
    with (tmp_path / "out" / "nextneopi_native_top20.csv").open() as handle:
        selected = list(csv.DictReader(handle))
    assert [row["variant_id"] for row in selected] == [f"v{i}" for i in range(20)]


def test_freeze_native_portfolio_rejects_duplicate_mutations(tmp_path):
    aggregate = tmp_path / "native.tsv"
    aggregate.write_text(
        "ID\tBest Peptide\tAllele\tTier\n"
        "v1\tAAAAAAAA\tA*01:01\tPass\n"
        "v1\tBBBBBBBB\tA*01:01\tPass\n"
    )
    with pytest.raises(ValueError, match="unique"):
        freeze_native_portfolio(aggregate, tmp_path / "out")


def test_prepare_batch_rejects_same_name_fastq_with_different_bytes(tmp_path, monkeypatch):
    monkeypatch.chdir(Path(__file__).parents[1])
    fastqs = [_fastq(tmp_path / f"r{i}.fastq") for i in range(6)]
    provenance = tmp_path / "CONVERT_PROVENANCE.json"
    provenance.write_text(json.dumps([
        {
            "status": "OK",
            "raw_input_identity_ready": True,
            "sha256_per_file": {path.name: hashlib.sha256(path.read_bytes()).hexdigest()},
            "bytes_per_file": {path.name: path.stat().st_size},
        }
        for path in fastqs
    ]))
    fastqs[0].write_bytes(b"@different\nAC\n+\n!!\n")
    with pytest.raises(ValueError, match="differs from Epicurus"):
        prepare_batch(
            patient_id="Hu_test",
            tumor_r1=fastqs[0], tumor_r2=fastqs[1],
            normal_r1=fastqs[2], normal_r2=fastqs[3],
            rna_r1=fastqs[4], rna_r2=fastqs[5],
            conversion_provenance=provenance,
            output_dir=tmp_path / "out",
        )
