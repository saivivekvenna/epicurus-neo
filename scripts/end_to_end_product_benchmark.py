"""Canonical full Epicurus-owned product-path audit on Hu_287 and Sid."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from benchmark.end_to_end_product import (  # noqa: E402
    evaluate_frozen_pipeline,
    freeze_product_pipeline,
)
from epicurus_neo.product import normalize_product_candidates  # noqa: E402
from event_b.sid_full_pipeline import prepare_sid_gate_frame  # noqa: E402


OUT = ROOT / "artifacts/milestone_7_decision/end_to_end_product"
HU_UNIVERSE = ROOT / "data/raw/miller_ipv/hu_287/freeze/universe.csv"
HU_LABELS = ROOT / "data/raw/miller_ipv/miller_recognition_labels.csv"
SID_UNIVERSE = ROOT / "artifacts/milestone_7_decision/sid_benchmark/scored_candidates.csv.gz"
SID_VAF = ROOT / "data/raw/osteosarc/site_cache/variant_vafs_long.tsv"
CODE = (
    ROOT / "src/epicurus_neo/product.py",
    ROOT / "src/epicurus_neo/gates.py",
    ROOT / "src/benchmark/end_to_end_product.py",
    ROOT / "scripts/end_to_end_product_benchmark.py",
)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _hu_frame() -> pd.DataFrame:
    raw = pd.read_csv(HU_UNIVERSE, low_memory=False)
    return normalize_product_candidates(raw, source_name="hu287_lossless_raw_reconstruction")


def _sid_frame() -> pd.DataFrame:
    return prepare_sid_gate_frame(SID_UNIVERSE, SID_VAF)


def _hu_positives() -> set[str]:
    labels = pd.read_csv(HU_LABELS)
    pos = labels[(labels["patient_id"] == "Hu_287") & (labels["label"] == "POSITIVE")]
    return {
        f"{str(row.chrom).removeprefix('chr')}:{int(row.pos)}:{row.ref}:{row.alt}"
        for row in pos.itertuples()
    }


def _sid_positives() -> set[str]:
    from event_b.sid_benchmark import hudson_positive_variant_ids

    return set(hudson_positive_variant_ids())


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    frames = {"Hu_287": _hu_frame(), "Sid": _sid_frame()}
    frozen = {
        "status": "FROZEN_BEFORE_LABEL_JOIN",
        "protocol": "PROTOCOL.md",
        "product_boundary": "raw-derived mutation/HLA/RNA candidate universe -> Epicurus top20",
        "patients": {},
        "input_sha256": {
            str(path.relative_to(ROOT)): _sha(path)
            for path in (HU_UNIVERSE, SID_UNIVERSE, SID_VAF)
        },
        "code_sha256": {str(path.relative_to(ROOT)): _sha(path) for path in CODE},
    }
    for patient, frame in frames.items():
        frozen["patients"][patient] = freeze_product_pipeline(frame)
        frozen["patients"][patient]["evidence_status"] = (
            "complete local tumor/normal WES + tumor RNA + inferred HLA + lossless class-I generation"
            if patient == "Hu_287"
            else "longitudinal WES/WGS + matched T2 RNA + inferred HLA; generation 130/147 mutations"
        )

    freeze_path = OUT / "FROZEN_PIPELINE.json"
    freeze_path.write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")

    # The only outcome imports occur after all product stages and selections are on disk.
    positives = {"Hu_287": _hu_positives(), "Sid": _sid_positives()}
    evaluation = {
        patient: evaluate_frozen_pipeline(frozen["patients"][patient], positives[patient])
        for patient in frames
    }
    result = {
        "status": "POST_HOC_PRODUCT_INTEGRATION_AUDIT_NOT_INDEPENDENT_VALIDATION",
        "labels_opened_after_freeze": True,
        "freeze_sha256": _sha(freeze_path),
        "patients": {
            patient: {
                "frozen_pipeline": frozen["patients"][patient],
                "evaluation_only": evaluation[patient],
            }
            for patient in frames
        },
        "headline_scope": "actual Epicurus-owned product path only",
    }
    (OUT / "RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (OUT / "REPORT.md").write_text(_markdown(result))
    print(
        json.dumps(
            {
                patient: value["evaluation_only"]["product_hits_at_20"]
                for patient, value in result["patients"].items()
            },
            indent=2,
        )
    )
    return 0


def _markdown(result: dict) -> str:
    lines = [
        "# Canonical Epicurus product-path end-to-end audit",
        "",
        "> Actual production normalization, deterministic gate, evidence score, eligibility policy, and capped portfolio. Both patients were previously inspected; this is an integration audit, not independent validation.",
        "",
        "| Patient | Generated positive mutations | Deterministic-valid | Product-eligible | Product top 20 | PRIME plain | PRIME cap-2 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for patient, value in result["patients"].items():
        ev = value["evaluation_only"]
        stages = ev["stage_reachability"]
        lines.append(
            f"| {patient} | {stages['generated']['n']}/{stages['generated']['of']} | "
            f"{stages['deterministic_valid']['n']}/{stages['deterministic_valid']['of']} | "
            f"{stages['product_eligible']['n']}/{stages['product_eligible']['of']} | "
            f"**{ev['product_hits_at_20']['n']}/{ev['product_hits_at_20']['of']}** | "
            f"{ev['prime_plain_hits_at_20']['n']}/{ev['prime_plain_hits_at_20']['of']} | "
            f"{ev['prime_cap2_hits_at_20']['n']}/{ev['prime_cap2_hits_at_20']['of']} |"
        )
    lines += ["", "## Patient funnels", ""]
    for patient, value in result["patients"].items():
        frozen = value["frozen_pipeline"]
        ev = value["evaluation_only"]
        lines += [
            f"### {patient}",
            "",
            f"- Candidate rows: {frozen['counts']['candidate_rows']:,}",
            f"- Generated mutations: {frozen['counts']['generated_mutations']}",
            f"- Deterministic-valid rows: {frozen['counts']['deterministic_valid_rows']:,}",
            f"- Product-eligible rows: {frozen['counts']['product_eligible_rows']:,}",
            f"- Selected routes / unique mutations: {frozen['counts']['selected_routes']} / {frozen['counts']['selected_unique_mutations']}",
            f"- Exclusions: {frozen['removal_reasons']}",
            f"- Positive last-reached stages: {ev['last_reached_stage_by_positive']}",
            "",
        ]
    lines += [
        "## Honest interpretation",
        "",
        "This report is the deliverable-level check: every headline number comes from the same shipped product logic, not from a mix-and-match research arm. Hu_287 tests a complete local raw-data reconstruction; Sid exposes both incomplete generation (130/147 mutations) and any downstream product losses. Neither patient is blind, so this establishes runnable behavior and patient-specific outcomes—not general superiority.",
    ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
