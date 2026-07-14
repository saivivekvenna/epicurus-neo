from __future__ import annotations

import csv
from pathlib import Path

import pytest

from benchmark.nextneopi_comparator import freeze_native_portfolio, prepare_batch


def _fastq(path: Path) -> Path:
    path.write_bytes(b"@r\nAC\n+\n!!\n")
    return path


def test_prepare_batch_hashes_three_raw_sample_types(tmp_path, monkeypatch):
    monkeypatch.chdir(Path(__file__).parents[1])
    fastqs = [_fastq(tmp_path / f"r{i}.fastq") for i in range(6)]
    manifest = prepare_batch(
        patient_id="Hu_test",
        tumor_r1=fastqs[0], tumor_r2=fastqs[1],
        normal_r1=fastqs[2], normal_r2=fastqs[3],
        rna_r1=fastqs[4], rna_r2=fastqs[5],
        output_dir=tmp_path / "out",
    )
    assert list(manifest["raw_fastqs"]) == ["tumor_DNA", "normal_DNA", "tumor_RNA"]
    assert manifest["labels_opened"] is False
    assert manifest["hla_file_supplied"] is False
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
