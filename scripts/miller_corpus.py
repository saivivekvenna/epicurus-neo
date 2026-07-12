"""Miller IPV corpus ingestion — read the S1/S2 supplement, validate, emit the label corpus + audit.

Reads the (now-obtained) supplement from the gitignored raw location, parses S2 (the canonical 754-peptide
/ 13-patient tested table) via the frozen ingestion contract, writes the unified recognition-label frame to
the gitignored raw dir, and emits a COMMITTED aggregate-only corpus audit (counts + per-patient + crosswalk
+ the exact remaining gaps). Peptide sequences are NOT committed (they stay in the gitignored raw frame);
only aggregate statistics are published.

LOCKED_TEST: Miller is never used for development/tuning. This runner only ingests + audits.

    .venv/bin/python -m scripts.miller_corpus
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from benchmark.miller_ingest import (  # noqa: E402
    mutation_recognition,
    parse_miller_labels,
    parse_sra_runinfo,
    patient_input_crosswalk,
    validate_recognition_labels,
    SRA_RUNINFO_FIXTURE,
)

RAW = ROOT / "data" / "raw" / "miller_ipv"
SUPP = RAW / "supplement"
S2 = SUPP / "scitranslmed.abj9905_data_files_s1_to_s4" / "abj9905_Data_file_S2.xlsx"
S1_V2 = SUPP / "scitranslmed.abj9905_data_files_s1_to_s4" / "scitranslmed.abj9905.data_file_s1.v2.xlsx"
S1_V1 = SUPP / "scitranslmed.abj9905.data_file_s1.v1.xlsx"
OUT = ROOT / "artifacts" / "milestone_7_decision" / "external_validation" / "miller_ipv"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run() -> dict:
    sra = parse_sra_runinfo(SRA_RUNINFO_FIXTURE)
    sra_patients = set(sra["patient_id"])

    raw_s2 = pd.read_excel(S2, sheet_name="Sheet 1")
    labels = parse_miller_labels(raw_s2, readout="IFN-g")            # PRIMARY (CD8/class-I biased)
    labels_any = parse_miller_labels(raw_s2, readout="any")          # sensitivity (paper's 199 set)

    contract = validate_recognition_labels(labels, sra_patients=sra_patients)
    mut = mutation_recognition(labels)
    mut_any = mutation_recognition(labels_any)

    # write the unified label frame to the gitignored raw dir (peptides not committed)
    RAW.mkdir(parents=True, exist_ok=True)
    labels.to_csv(RAW / "miller_recognition_labels.csv", index=False)

    label_patients = set(labels["patient_id"])
    per_patient = (
        labels.assign(pos=labels["label"].eq("POSITIVE"))
        .groupby("patient_id")
        .agg(n_peptides=("label", "size"), n_ifng_pos=("pos", "sum"),
             n_mutations=("mutation_id", "nunique"))
        .reset_index()
    )
    per_patient = per_patient.merge(
        mut.groupby("patient_id")["recognized"].sum().rename("n_recognized_mutations_ifng").reset_index(),
        on="patient_id", how="left")

    audit = {
        "cohort": "miller_ipv",
        "status": "LABELS_INGESTED_AND_VALIDATED",
        "intended_use": "LOCKED_TEST",
        "doi": "10.1126/scitranslmed.abj9905",
        "pmid": 38416845,
        "checksums_sha256": {
            "S2": _sha256(S2), "S1_v2": _sha256(S1_V2),
            "S1_v1": _sha256(S1_V1) if S1_V1.exists() else None,
        },
        "canonical_label_source": "S2 (Data_file_S2.xlsx, 'Sheet 1'): 754 tested 20-mers x 13 patients",
        "sheet_map": {
            "S1_v2": "single-patient (Hu_159) detail with 'mut position in peptide' + ref peptide",
            "S2": "CANONICAL per-peptide tested table (patient/variant/peptide + IFN-g/IL-5/both/any)",
            "S3": "figure source data (IFN-g / IL-5 SFC/10^6 PBMC) — not a per-candidate label table",
            "S4": "reagents / RRID list — irrelevant",
        },
        "counts": {
            "n_peptides_tested": int(len(labels)),
            "n_patients": int(labels["patient_id"].nunique()),
            "n_mutations": int(labels["mutation_id"].nunique()),
            "label_ifng": labels["label"].value_counts().to_dict(),
            "n_ifng_pos_peptides": int((labels["ifn_g"]).sum()),
            "n_il5_pos_peptides": int((labels["il_5"]).sum()),
            "n_any_pos_peptides": int((labels["any_cytokine"]).sum()),
            "n_recognized_mutations_ifng": int(mut["recognized"].sum()),
            "n_recognized_mutations_any": int(mut_any["recognized"].sum()),
            "peptide_length_dist": labels["mutant_peptide"].str.len().value_counts().to_dict(),
            "variant_type_dist": labels["source_variant_type"].value_counts().to_dict(),
        },
        "contract_validation": contract,
        "sra_crosswalk": {
            "sra_patients": sorted(sra_patients),
            "label_patients": sorted(label_patients),
            "all_label_patients_have_inputs": bool(label_patients <= sra_patients),
            "patients_matched": int(len(label_patients & sra_patients)),
            "input_completeness": bool(patient_input_crosswalk(sra)["complete"].all()),
        },
        "per_patient": per_patient.to_dict("records"),
        "reconciliation_notes": [
            f"Paper reports 349 tested variants; S2 has {int(labels['mutation_id'].nunique())} distinct "
            f"(gene,chr,pos) mutations — a {349 - int(labels['mutation_id'].nunique())}-variant gap (likely "
            "a few multi-transcript/indel rows collapsing on coordinates); documented, not resolved.",
            f"Paper's '199 (26 pct) induced T cell responses' = the 'any' (IFN-g OR IL-5) set (verified "
            f"199); the PRIMARY class-I label here is IFN-g ({int((labels['ifn_g']).sum())} pos peptides).",
        ],
        "remaining_gaps_for_the_four_arm_run": [
            "HLA per patient: NOT in S2 (20-mer ELISpot, no MHC restriction) -> type from WES (OptiType).",
            "RNA expression/TPM: NOT in S2 -> quantify from the RNA-seq runs.",
            "Full candidate universe: S2 has 343 IPV-PREFILTERED tested mutations -> re-enumerate the full "
            "class-I mutanome from WES for a fair denominator (both arms share it).",
            "20-mer -> class-I minimal epitope: recover mutant residue position (present: derivable from "
            "mut vs ref peptide diff) and enumerate 8-11mers spanning it for PRIME/Epicurus scoring.",
            "=> the LABEL half is DONE; the ranking arms are BLOCKED on the SRA download + bioinformatics "
            "(HLA + expression + mutanome), not on any missing file.",
        ],
    }
    _write(audit)
    return audit


def _write(audit: dict) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "CORPUS_AUDIT.json").write_text(json.dumps(audit, indent=2, default=str) + "\n")
    c = audit["counts"]
    x = audit["sra_crosswalk"]
    lines = [
        "# Miller IPV — recognition-label corpus audit",
        "",
        f"> Status **{audit['status']}** · LOCKED_TEST · DOI {audit['doi']} (PMID {audit['pmid']}).",
        f"> Canonical label source: {audit['canonical_label_source']}.",
        "",
        "## Counts",
        f"- peptides tested: **{c['n_peptides_tested']}** · patients: **{c['n_patients']}** · "
        f"mutations: **{c['n_mutations']}**",
        f"- IFN-g label (primary): {c['label_ifng']}",
        f"- positive peptides — IFN-g {c['n_ifng_pos_peptides']} · IL-5 {c['n_il5_pos_peptides']} · "
        f"any {c['n_any_pos_peptides']} (paper's 199)",
        f"- recognized mutations — IFN-g {c['n_recognized_mutations_ifng']} · any "
        f"{c['n_recognized_mutations_any']}",
        f"- peptide lengths: {c['peptide_length_dist']} · variant types: {c['variant_type_dist']}",
        "",
        "## Contract validation",
        f"- ok=**{audit['contract_validation']['ok']}** · rows={audit['contract_validation']['n_rows']} · "
        f"invalid_labels={audit['contract_validation']['invalid_labels']} · "
        f"invalid_peptides={audit['contract_validation']['n_invalid_peptides']} · "
        f"conflicting_keys={audit['contract_validation']['n_conflicting_keys']}",
        "",
        "## SRA crosswalk",
        f"- all {x['patients_matched']} label patients have public inputs: "
        f"**{x['all_label_patients_have_inputs']}** · input trios complete: {x['input_completeness']}",
        "",
        "## Per-patient (IFN-g)",
        "",
        "| patient | peptides | IFN-g+ peptides | mutations | recognized mutations |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in audit["per_patient"]:
        lines.append(f"| {r['patient_id']} | {r['n_peptides']} | {int(r['n_ifng_pos'])} | "
                     f"{r['n_mutations']} | {int(r['n_recognized_mutations_ifng'])} |")
    lines += ["", "## Reconciliation notes", ""]
    lines += [f"- {n}" for n in audit["reconciliation_notes"]]
    lines += ["", "## Remaining gaps for the four-arm run", ""]
    lines += [f"- {g}" for g in audit["remaining_gaps_for_the_four_arm_run"]]
    (OUT / "CORPUS_AUDIT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    audit = run()
    print(json.dumps({"status": audit["status"], "counts": audit["counts"],
                      "crosswalk": audit["sra_crosswalk"]["all_label_patients_have_inputs"],
                      "contract_ok": audit["contract_validation"]["ok"]}, indent=2, default=str))


if __name__ == "__main__":
    main()
