"""Frozen, label-blind common-input contract for pVACtools and Vaxrank."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


CONFIG = Path("configs/frozen/controlled_comparators_v1.json")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file(path: str | Path) -> Path:
    resolved = Path(path).expanduser().resolve()
    if not resolved.is_file() or resolved.is_symlink():
        raise FileNotFoundError(resolved)
    return resolved


def read_hla(path: str | Path) -> list[str]:
    alleles = [line.strip() for line in _file(path).read_text().splitlines() if line.strip()]
    if not alleles or len(alleles) != len(set(alleles)):
        raise ValueError("HLA panel must contain unique, non-empty alleles")
    if any(not allele.startswith("HLA-") for allele in alleles):
        raise ValueError("HLA panel must use HLA-* allele names")
    return alleles


def prepare_common_bundle(
    *,
    patient_id: str,
    pass_vcf: str | Path,
    pass_vcf_index: str | Path,
    rna_bam: str | Path,
    rna_bam_index: str | Path,
    hla_panel: str | Path,
    pvac_ready_vcf: str | Path,
    pvac_ready_vcf_index: str | Path,
    output_dir: str | Path,
) -> dict:
    """Hash the common evidence and emit exact, outcome-free comparator commands."""
    if not patient_id or any(char in patient_id for char in ",\n\r"):
        raise ValueError("patient_id must be non-empty and CSV-safe")
    paths = {
        "pass_vcf": _file(pass_vcf),
        "pass_vcf_index": _file(pass_vcf_index),
        "rna_bam": _file(rna_bam),
        "rna_bam_index": _file(rna_bam_index),
        "hla_panel": _file(hla_panel),
        "pvac_ready_vcf": _file(pvac_ready_vcf),
        "pvac_ready_vcf_index": _file(pvac_ready_vcf_index),
    }
    alleles = read_hla(paths["hla_panel"])
    config = _file(CONFIG)
    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    pvac_out = destination / "pvacseq"
    vaxrank_out = destination / "vaxrank"
    algorithms = "NetMHCpan NetMHCpanEL MHCflurry MHCflurryEL MHCnuggetsI BigMHC_IM DeepImmuno"
    common_filters = (
        "--binding-threshold 500 --top-score-metric median --minimum-fold-change 0.0 "
        "--normal-cov 5 --tdna-cov 10 --trna-cov 10 --normal-vaf 0.02 "
        "--tdna-vaf 0.25 --trna-vaf 0.25 --expn-val 1 --maximum-transcript-support-level 1"
    )
    manifest = {
        "policy_id": "epicurus-controlled-comparators-v1",
        "patient_id": patient_id,
        "labels_opened": False,
        "config": {"path": str(config), "sha256": _sha256(config)},
        "inputs": {
            name: {"path": str(path), "size": path.stat().st_size, "sha256": _sha256(path)}
            for name, path in paths.items()
        },
        "hla_alleles": alleles,
        "commands": {
            "pvacseq": (
                f"pvacseq run --normal-sample-name {patient_id}_N -e1 8,9,10,11 "
                + common_filters + " "
                f"{paths['pvac_ready_vcf']} {patient_id}_T "
                f"{','.join(alleles)} {algorithms} {pvac_out}"
            ),
            "vaxrank": (
                f"vaxrank --vcf {paths['pass_vcf']} --bam {paths['rna_bam']} "
                f"--mhc-predictor mhcflurry --mhc-alleles {','.join(alleles)} "
                f"--output-dir {vaxrank_out}"
            ),
        },
        "allow_dna_only_fallback": False,
    }
    (destination / "CONTROLLED_INPUT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return manifest


def freeze_vaxrank_portfolio(ranked_csv: str | Path, output_dir: str | Path) -> dict:
    """Freeze one top peptide for each of Vaxrank's first 20 natively ranked variants."""
    source = _file(ranked_csv)
    with source.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"chr", "pos", "ref", "alt", "variant_rank", "peptide_rank", "amino_acids"}
        missing = sorted(required.difference(reader.fieldnames or []))
        if missing:
            raise ValueError(f"Vaxrank ranked CSV missing columns: {missing}")
        rows = list(reader)
    if not rows:
        raise ValueError("Vaxrank ranked CSV is empty")

    by_rank: dict[int, dict] = {}
    for row in rows:
        rank = int(row["variant_rank"])
        peptide_rank = int(row["peptide_rank"])
        if rank < 1 or peptide_rank < 1:
            raise ValueError("Vaxrank ranks must be positive")
        current = by_rank.get(rank)
        if current is None or peptide_rank < int(current["peptide_rank"]):
            by_rank[rank] = row
    ordered = sorted(by_rank)
    if ordered != list(range(1, len(ordered) + 1)):
        raise ValueError("Vaxrank variant ranks must be contiguous from 1")
    selected = [by_rank[rank] for rank in ordered[:20]]

    destination = Path(output_dir).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    portfolio = destination / "vaxrank_native_top20.csv"
    fields = ("rank", "variant_id", "peptide", "combined_score", "expression_score", "epitope_score")
    with portfolio.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in selected:
            writer.writerow(
                {
                    "rank": int(row["variant_rank"]),
                    "variant_id": f"{row['chr']}:{row['pos']}:{row['ref']}:{row['alt']}",
                    "peptide": row["amino_acids"],
                    "combined_score": row.get("combined_score", ""),
                    "expression_score": row.get("expression_score", ""),
                    "epitope_score": row.get("target_epitope_score", ""),
                }
            )
    result = {
        "policy_id": "vaxrank-native-top20-v1",
        "labels_opened": False,
        "source": {"path": str(source), "sha256": _sha256(source), "rows": len(rows)},
        "portfolio": {"path": str(portfolio), "sha256": _sha256(portfolio), "size": len(selected)},
    }
    (destination / "VAXRANK_PORTFOLIO_MANIFEST.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    return result
