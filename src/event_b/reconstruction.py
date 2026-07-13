"""NeoVax melanoma (Ott 2017 + Hu 2021) decision-problem reconstruction scaffold.

This module turns the data-quality *measurement* into a data-acquisition *action*. It records,
per vaccinated patient, exactly which elements of a reconstructable decision problem are available
locally, which sit behind controlled access (dbGaP phs001451), and which must be regenerated — then
applies a strict completeness gate that **cannot** call any patient decision-ready until every
element is verified present.

It also extracts the one genuinely new local asset: the MAF-style somatic mutation identities for
Hu patients 11-12 (Suppl Dataset 2a/2b). That is a *mutation-level* universe fragment — deliberately
NOT a peptide candidate universe and NOT a denominator:

  * it has no tumor VAF, no read depth, no RNA expression, no clonality;
  * it needs the patient HLA genotype and a candidate-generation run to become a candidate universe;
  * the vaccine peptides remain the INCLUDED subset; untested candidates stay UNTESTED, never negative.

Nothing here fabricates evidence or downloads controlled data. It reuses the Hu streaming reader so
the 257 MB workbook is never fully materialised, and it runs (in a degraded, clearly-flagged mode)
even when that manually-placed workbook is absent.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import yaml

from event_b.adapters.hu_neovax import (
    STUDY_ID,
    hu_source_path,
    read_workbook_sheets,
    stage_hu_supplement,
)


RECONSTRUCTION_VERSION = "event-b-reconstruction-1.0.0"

DEFAULT_MANIFEST = "configs/source_manifests/neovax_reconstruction.yml"

# Worksheet numbers inside the Hu workbook carrying the somatic mutation identities.
SHEET_SNP = 3   # Suppl Dataset 2a. snp
SHEET_INDEL = 4  # Suppl Dataset 2b. indel

# Protein-altering variant classifications that can yield a neoepitope. Everything else
# (Intron, Silent, 3'/5'UTR, IGR, RNA, Flank, Non-coding_Transcript) is upstream noise for
# candidate generation and is excluded from the neoantigen-generating count.
NEOANTIGEN_GENERATING_CLASSES = frozenset(
    {
        "Missense_Mutation",
        "Nonsense_Mutation",
        "Nonstop_Mutation",
        "Splice_Site",
        "Frame_Shift_Del",
        "Frame_Shift_Ins",
        "In_Frame_Del",
        "In_Frame_Ins",
        "Translation_Start_Site",
        "De_novo_Start_InFrame",
        "De_novo_Start_OutOfFrame",
    }
)

# Element readiness statuses. Only PRESENT_LOCAL counts toward decision-readiness.
STATUS_PRESENT_LOCAL = "PRESENT_LOCAL"
STATUS_PARTIAL_LOCAL = "PARTIAL_LOCAL"
STATUS_BLOCKED = "BLOCKED_CONTROLLED_ACCESS"
STATUS_REGEN = "REQUIRES_REGENERATION"
STATUS_ABSENT = "ABSENT"

_COHORT_OTT = "ott_1_6"
_COHORT_HU = "hu_11_12"


def cohort_of(source_patient: str) -> str:
    """Map a source patient number to its manifest cohort key."""
    text = str(source_patient).strip()
    if text in {"11", "12"}:
        return _COHORT_HU
    if text in {"1", "2", "3", "4", "5", "6"}:
        return _COHORT_OTT
    return _COHORT_OTT


# --------------------------------------------------------------------------------------
# Manifest
# --------------------------------------------------------------------------------------
def load_reconstruction_manifest(path: str | Path = DEFAULT_MANIFEST) -> dict:
    with Path(path).open() as handle:
        return yaml.safe_load(handle)


def vaccinated_patients(manifest: dict) -> list[str]:
    """The vaccinated source patients, taken from the manifest (not hard-coded here)."""
    supplies = (
        manifest.get("local_sources", {})
        .get("hu_supplement_workbook", {})
        .get("supplies", {})
        .get("vaccine_components_and_assays", {})
    )
    patients = [str(p) for p in supplies.get("scope_patients", [])]
    return sorted(patients, key=lambda p: int(p))


# --------------------------------------------------------------------------------------
# Somatic mutation universe extraction (pure grid parsing + guarded file reader)
# --------------------------------------------------------------------------------------
def _header_map(row: list[str]) -> dict[str, int]:
    return {str(cell).strip().lower(): i for i, cell in enumerate(row) if str(cell).strip()}


def _pick(headers: dict[str, int], *names: str) -> int | None:
    for name in names:
        idx = headers.get(name.lower())
        if idx is not None:
            return idx
    return None


def _get(row: list[str], idx: int | None) -> str:
    if idx is None or idx < 0 or idx >= len(row):
        return ""
    return str(row[idx]).strip()


def parse_mutation_grid(grid: list[list[str]], *, variant_source: str) -> list[dict]:
    """Parse a Suppl Dataset 2a/2b grid into normalized mutation rows.

    Pure function over an already-extracted grid, so it is testable without the workbook.
    """
    header_idx = next(
        (i for i, row in enumerate(grid) if row and str(row[0]).strip().lower() == "patient id"),
        None,
    )
    if header_idx is None:
        return []
    headers = _header_map(grid[header_idx])
    col_patient = _pick(headers, "patient id")
    col_gene = _pick(headers, "hugo symbol", "gene")
    col_transcript = _pick(headers, "annotation transcript", "transcript")
    col_class = _pick(headers, "variant classification")
    col_type = _pick(headers, "variant type")
    col_chrom = _pick(headers, "chromosome")
    col_start = _pick(headers, "start position", "start")
    col_end = _pick(headers, "end position", "end")
    col_ref = _pick(headers, "reference allele")
    col_alt = _pick(headers, "tumor seq allele2", "tumor seq allele1")
    col_cdna = _pick(headers, "cdna change")
    col_codon = _pick(headers, "codon change")
    col_protein = _pick(headers, "protein change")

    rows: list[dict] = []
    for row in grid[header_idx + 1 :]:
        patient = _get(row, col_patient)
        if not patient.replace(".", "", 1).isdigit():
            continue
        source_patient = str(int(float(patient)))
        classification = _get(row, col_class)
        rows.append(
            {
                "patient_id": f"{STUDY_ID}:{source_patient}",
                "source_patient": source_patient,
                "cohort": cohort_of(source_patient),
                "gene": _get(row, col_gene),
                "transcript": _get(row, col_transcript),
                "variant_classification": classification,
                "variant_type": _get(row, col_type),
                "chromosome": _get(row, col_chrom),
                "start_position": _get(row, col_start),
                "end_position": _get(row, col_end),
                "reference_allele": _get(row, col_ref),
                "tumor_allele": _get(row, col_alt),
                "cdna_change": _get(row, col_cdna),
                "codon_change": _get(row, col_codon),
                "protein_change": _get(row, col_protein),
                "variant_source": variant_source,
                "is_neoantigen_generating": classification in NEOANTIGEN_GENERATING_CLASSES,
            }
        )
    return rows


def extract_hu_mutation_universe(raw_dir: str | Path) -> pd.DataFrame:
    """Extract the Hu (patients 11-12) somatic mutation identities from the workbook.

    Reuses the streaming reader; raises the adapter's actionable message if the checksum-pinned
    workbook is absent (no fabrication). Callers that must tolerate absence should guard with
    ``hu_workbook_available``.
    """
    path = stage_hu_supplement(raw_dir)
    sheets = read_workbook_sheets(path, (SHEET_SNP, SHEET_INDEL))
    rows = parse_mutation_grid(sheets[SHEET_SNP], variant_source="snp")
    rows += parse_mutation_grid(sheets[SHEET_INDEL], variant_source="indel")
    columns = [
        "patient_id", "source_patient", "cohort", "gene", "transcript",
        "variant_classification", "variant_type", "chromosome", "start_position",
        "end_position", "reference_allele", "tumor_allele", "cdna_change",
        "codon_change", "protein_change", "variant_source", "is_neoantigen_generating",
    ]
    return pd.DataFrame(rows, columns=columns)


def hu_workbook_available(raw_dir: str | Path) -> bool:
    try:
        stage_hu_supplement(raw_dir)
        return True
    except RuntimeError:
        return False


def mutation_universe_summary(frame: pd.DataFrame) -> dict[str, dict]:
    """Per-patient partial-universe sizes. Explicitly a mutation-level fragment, not a denominator."""
    summary: dict[str, dict] = {}
    if frame.empty:
        return summary
    for source_patient, group in frame.groupby(frame.source_patient.astype(str)):
        summary[str(source_patient)] = {
            "n_somatic_mutations": int(len(group)),
            "n_neoantigen_generating": int(group["is_neoantigen_generating"].sum()),
            "n_snp": int((group["variant_source"] == "snp").sum()),
            "n_indel": int((group["variant_source"] == "indel").sum()),
        }
    return summary


# --------------------------------------------------------------------------------------
# Evidence ledger + strict readiness gate
# --------------------------------------------------------------------------------------
@dataclass(frozen=True)
class ReadinessVerdict:
    verdict: str
    n_present_local: int
    n_required: int
    blocking_elements: tuple[tuple[str, str], ...]  # (element_key, status)


def _element_status(element: dict, cohort: str) -> str:
    by_cohort = element.get("status_by_cohort", {})
    return str(by_cohort.get(cohort, STATUS_ABSENT))


def evaluate_patient_readiness(manifest: dict, cohort: str) -> ReadinessVerdict:
    """Strict gate: decision-ready only if EVERY required element is PRESENT_LOCAL for the cohort."""
    rule = manifest.get("readiness_rule", {})
    ready_status = rule.get("ready_status", STATUS_PRESENT_LOCAL)
    ready_verdict = rule.get("verdict_when_all_ready", "DECISION_PROBLEM_READY_UNIVERSE")
    not_ready_verdict = rule.get("verdict_otherwise", "NOT_DECISION_READY")
    elements = manifest.get("required_elements", [])
    present = 0
    blockers: list[tuple[str, str]] = []
    for element in elements:
        status = _element_status(element, cohort)
        if status == ready_status:
            present += 1
        else:
            blockers.append((str(element["key"]), status))
    verdict = ready_verdict if not blockers else not_ready_verdict
    return ReadinessVerdict(verdict, present, len(elements), tuple(blockers))


def build_evidence_ledger(
    manifest: dict,
    mutation_summary: dict[str, dict] | None = None,
    workbook_present: bool = False,
) -> pd.DataFrame:
    """One row per (patient, required element): declared status, readiness, and verification detail."""
    mutation_summary = mutation_summary or {}
    elements = manifest.get("required_elements", [])
    rows: list[dict] = []
    for source_patient in vaccinated_patients(manifest):
        cohort = cohort_of(source_patient)
        for element in elements:
            key = str(element["key"])
            status = _element_status(element, cohort)
            verified = False
            detail = ""
            if key == "somatic_mutations" and cohort == _COHORT_HU:
                if workbook_present and source_patient in mutation_summary:
                    verified = True
                    info = mutation_summary[source_patient]
                    detail = (
                        f"{info['n_somatic_mutations']} mutations "
                        f"({info['n_neoantigen_generating']} neoantigen-generating); "
                        "identities only, no VAF/expression"
                    )
                else:
                    detail = "declared PARTIAL_LOCAL but workbook not verified in this run"
            rows.append(
                {
                    "patient_id": f"{STUDY_ID}:{source_patient}",
                    "source_patient": source_patient,
                    "cohort": cohort,
                    "element": key,
                    "status": status,
                    "is_ready_element": status == STATUS_PRESENT_LOCAL,
                    "verified_local": verified,
                    "detail": detail,
                }
            )
    return pd.DataFrame(rows)


def patient_readiness_frame(manifest: dict, mutation_summary: dict[str, dict] | None = None,
                            workbook_present: bool = False) -> pd.DataFrame:
    """Per-patient roll-up: verdict, element-status counts, and the partial-universe size."""
    mutation_summary = mutation_summary or {}
    elements = manifest.get("required_elements", [])
    rows: list[dict] = []
    for source_patient in vaccinated_patients(manifest):
        cohort = cohort_of(source_patient)
        verdict = evaluate_patient_readiness(manifest, cohort)
        counts = {STATUS_PRESENT_LOCAL: 0, STATUS_PARTIAL_LOCAL: 0, STATUS_BLOCKED: 0,
                  STATUS_REGEN: 0, STATUS_ABSENT: 0}
        for element in elements:
            counts[_element_status(element, cohort)] = counts.get(_element_status(element, cohort), 0) + 1
        info = mutation_summary.get(source_patient, {})
        rows.append(
            {
                "patient_id": f"{STUDY_ID}:{source_patient}",
                "source_patient": source_patient,
                "cohort": cohort,
                "verdict": verdict.verdict,
                "n_present_local": verdict.n_present_local,
                "n_required": verdict.n_required,
                "n_partial_local": counts[STATUS_PARTIAL_LOCAL],
                "n_blocked_controlled": counts[STATUS_BLOCKED],
                "n_requires_regen": counts[STATUS_REGEN],
                "n_absent": counts[STATUS_ABSENT],
                "partial_universe_mutations": int(info.get("n_somatic_mutations", 0)),
                "partial_universe_neoantigen_generating": int(info.get("n_neoantigen_generating", 0)),
                "blocking_elements": ", ".join(f"{k}:{s}" for k, s in verdict.blocking_elements),
            }
        )
    frame = pd.DataFrame(rows)
    return frame.sort_values(
        ["n_present_local", "n_partial_local"], ascending=False, kind="mergesort"
    ).reset_index(drop=True)


# --------------------------------------------------------------------------------------
# Orchestration + rendering
# --------------------------------------------------------------------------------------
def reconstruction_audit(raw_dir: str | Path = "data/raw/hu_melanoma_2021",
                         manifest_path: str | Path = DEFAULT_MANIFEST) -> dict:
    """Run the full reconstruction audit; degrades cleanly if the workbook is absent."""
    manifest = load_reconstruction_manifest(manifest_path)
    present = hu_workbook_available(raw_dir)
    if present:
        mutation_frame = extract_hu_mutation_universe(raw_dir)
        mut_summary = mutation_universe_summary(mutation_frame)
    else:
        mut_summary = {}
    readiness = patient_readiness_frame(manifest, mut_summary, workbook_present=present)
    ledger = build_evidence_ledger(manifest, mut_summary, workbook_present=present)
    any_ready = bool((readiness["verdict"] == manifest["readiness_rule"]["verdict_when_all_ready"]).any()) if not readiness.empty else False
    accessions = manifest.get("sequencing_accessions", {})
    return {
        "version": RECONSTRUCTION_VERSION,
        "workbook_present": present,
        "dbgap_study": accessions.get("dbgap_study"),
        "access_level": accessions.get("access_level"),
        "redistribution": accessions.get("redistribution"),
        "n_vaccinated_patients": len(vaccinated_patients(manifest)),
        "any_decision_ready_patient": any_ready,
        "total_partial_universe_mutations": int(readiness["partial_universe_mutations"].sum()) if not readiness.empty else 0,
        "patients_with_local_mutation_universe": int((readiness["partial_universe_mutations"] > 0).sum()) if not readiness.empty else 0,
        "mutation_universe_summary": mut_summary,
        "readiness": readiness,
        "ledger": ledger,
        "manifest": manifest,
    }


def render_reconstruction_markdown(audit: dict) -> str:
    manifest = audit["manifest"]
    accessions = manifest.get("sequencing_accessions", {})
    readiness: pd.DataFrame = audit["readiness"]
    lines = [
        "# NeoVax reconstruction & acquisition audit",
        "",
        f"Version `{audit['version']}`. No controlled data is downloaded and no evidence is fabricated.",
        "",
        "## Controlled-access sequencing",
        "",
        f"- dbGaP study: **{accessions.get('dbgap_study')}** "
        f"({', '.join(v['id'] for v in accessions.get('versions', []))})",
        f"- Access level: **{accessions.get('access_level')}**",
        f"- Redistribution: **{accessions.get('redistribution')}**",
        f"- Workbook present this run: **{audit['workbook_present']}**",
        "",
        "## Decision-readiness (strict gate)",
        "",
        f"- Vaccinated patients: **{audit['n_vaccinated_patients']}**",
        f"- Any decision-problem-ready patient: **{audit['any_decision_ready_patient']}**",
        f"- Patients with a local mutation-universe fragment: "
        f"**{audit['patients_with_local_mutation_universe']}**",
        f"- Total local somatic mutations (partial universe, identities only): "
        f"**{audit['total_partial_universe_mutations']}**",
        "",
        "| Patient | Cohort | Verdict | Present | Partial | Blocked | Regen | Absent | Local muts (neoAg) |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for _, r in readiness.iterrows():
        muts = (
            f"{r['partial_universe_mutations']} ({r['partial_universe_neoantigen_generating']})"
            if r["partial_universe_mutations"] else "—"
        )
        lines.append(
            f"| {r['source_patient']} | {r['cohort']} | {r['verdict']} | "
            f"{r['n_present_local']}/{r['n_required']} | {r['n_partial_local']} | "
            f"{r['n_blocked_controlled']} | {r['n_requires_regen']} | {r['n_absent']} | {muts} |"
        )
    lines += [
        "",
        "## Shortest route to the first decision-ready patient",
        "",
        "Patients 11-12 are closest: their somatic mutation identities are already local. The "
        "remaining blockers are all in controlled-access dbGaP or require regeneration:",
        "",
    ]
    for step in manifest.get("acquisition_checklist", []):
        lines.append(f"1. {step}")
    lines.append("")
    return "\n".join(lines)
