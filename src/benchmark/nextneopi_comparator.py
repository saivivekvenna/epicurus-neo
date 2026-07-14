"""Label-blind preparation and portfolio freezing for the nextNEOpi Track-A comparator."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


CONFIG = Path("configs/frozen/nextneopi_track_a_v1.json")
REQUIRED_SAMPLE_TYPES = ("tumor_DNA", "normal_DNA", "tumor_RNA")
AGGREGATE_REQUIRED = ("ID", "Best Peptide", "Allele", "Tier")


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _checked_fastq(path: str | Path) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)
    if not resolved.name.endswith((".fastq.gz", ".fq.gz", ".fastq", ".fq")):
        raise ValueError(f"not a FASTQ path: {resolved}")
    return resolved


def prepare_batch(
    *,
    patient_id: str,
    tumor_r1: str | Path,
    tumor_r2: str | Path,
    normal_r1: str | Path,
    normal_r2: str | Path,
    rna_r1: str | Path,
    rna_r2: str | Path,
    output_dir: str | Path,
    sex: str = "NA",
) -> dict:
    """Create nextNEOpi's FASTQ batch CSV and a checksum-bound, label-free input manifest."""
    if not patient_id or any(char in patient_id for char in ",\n\r"):
        raise ValueError("patient_id must be a non-empty CSV-safe string")
    if sex not in {"female", "male", "NA"}:
        raise ValueError("sex must be female, male, or NA")

    inputs = {
        "tumor_DNA": (_checked_fastq(tumor_r1), _checked_fastq(tumor_r2)),
        "normal_DNA": (_checked_fastq(normal_r1), _checked_fastq(normal_r2)),
        "tumor_RNA": (_checked_fastq(rna_r1), _checked_fastq(rna_r2)),
    }
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    batch = destination / "nextneopi_batch.csv"
    with batch.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("sampleName", "reads1", "reads2", "sampleType", "HLAfile", "sex"),
        )
        writer.writeheader()
        for sample_type in REQUIRED_SAMPLE_TYPES:
            read1, read2 = inputs[sample_type]
            writer.writerow(
                {
                    "sampleName": patient_id,
                    "reads1": str(read1),
                    "reads2": str(read2),
                    "sampleType": sample_type,
                    "HLAfile": "",
                    "sex": sex,
                }
            )

    config = CONFIG.resolve()
    if not config.is_file():
        raise FileNotFoundError(config)
    manifest = {
        "policy_id": "nextneopi-track-a-native-aggregate-v1",
        "patient_id": patient_id,
        "labels_opened": False,
        "batch_csv": {"path": str(batch.resolve()), "sha256": sha256_file(batch)},
        "config": {"path": str(config), "sha256": sha256_file(config)},
        "raw_fastqs": {
            sample_type: [
                {"path": str(path), "size": path.stat().st_size, "sha256": sha256_file(path)}
                for path in inputs[sample_type]
            ]
            for sample_type in REQUIRED_SAMPLE_TYPES
        },
        "hla_file_supplied": False,
        "command": (
            "NXF_VER=22.10.8 nextflow run nextNEOpi.nf --batchFile "
            f"{batch.resolve()} --WES true --RNA_tag_seq false "
            "--pVACseq_filter_set standard -profile singularity -resume"
        ),
    }
    manifest_path = destination / "INPUT_MANIFEST.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return manifest


def freeze_native_portfolio(aggregate_tsv: str | Path, output_dir: str | Path) -> dict:
    """Freeze the first 20 rows of pVACtools' native, already-sorted aggregate report."""
    aggregate = Path(aggregate_tsv).resolve()
    with aggregate.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [column for column in AGGREGATE_REQUIRED if column not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"native aggregate missing columns: {missing}")
        rows = list(reader)
    if not rows:
        raise ValueError("native aggregate is empty")
    ids = [row["ID"] for row in rows]
    if any(not value for value in ids) or len(ids) != len(set(ids)):
        raise ValueError("native aggregate IDs must be non-empty and unique")

    selected = rows[:20]
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    portfolio = destination / "nextneopi_native_top20.csv"
    columns = (
        "rank", "variant_id", "peptide", "hla", "tier", "ic50_mt",
        "rna_expr", "rna_vaf", "dna_vaf",
    )
    with portfolio.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for rank, row in enumerate(selected, 1):
            writer.writerow(
                {
                    "rank": rank,
                    "variant_id": row["ID"],
                    "peptide": row["Best Peptide"],
                    "hla": row["Allele"],
                    "tier": row["Tier"],
                    "ic50_mt": row.get("IC50 MT", ""),
                    "rna_expr": row.get("RNA Expr", ""),
                    "rna_vaf": row.get("RNA VAF", ""),
                    "dna_vaf": row.get("DNA VAF", ""),
                }
            )

    manifest = {
        "policy_id": "nextneopi-track-a-native-aggregate-v1",
        "labels_opened": False,
        "native_aggregate": {
            "path": str(aggregate),
            "sha256": sha256_file(aggregate),
            "row_count": len(rows),
        },
        "portfolio": {
            "path": str(portfolio.resolve()),
            "sha256": sha256_file(portfolio),
            "size": len(selected),
            "selection": "first 20 rows in native pVACtools aggregate order",
        },
    }
    (destination / "PORTFOLIO_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return manifest
