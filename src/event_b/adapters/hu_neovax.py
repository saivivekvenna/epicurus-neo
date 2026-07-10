"""Hu 2021 melanoma NeoVax (NCT01970358) Event-B adapter.

Long-term follow-up of the Ott 2017 personal neoantigen vaccine trial. The source
supplement consolidates per-peptide CD8 (class-I predicted minimal epitopes) and CD4
(class-II overlapping assay peptides) IFN-gamma ELISpot immunogenicity for the eight
vaccinated patients (Ott's original 1-6 plus new 11-12), at week 16 and at 3-4.5 years,
together with epitope-spreading responses to neoantigens that were NOT in the vaccine.
Epitope spreading is kept strictly separate and never counted as vaccine-candidate
recognition.

Source (author manuscript; not in the open-access bulk subset; supplied manually and
pinned by sha256):
    Hu Z, Leet DE, Allesoe RL, et al. Personal neoantigen vaccines induce persistent
    memory T cell responses and epitope spreading in patients with melanoma.
    Nat Med 2021;27:515-525. doi:10.1038/s41591-020-01206-4 ; PMC8273876 ;
    ClinicalTrials NCT01970358.

Positivity: the source's own 0/1 ELISpot calls are ingested as reported (Hu 2021 scored
a response positive at >=2.5x the DMSO control). Unlike Braun, the table carries no raw
replicates to recompute from, so the call is author-reported and reconciled against Ott
2017's published totals (CD8 15/97, CD4 58/97 neoantigens across patients 1-6).

De-novo basis: both Ott 2017 and Hu 2021 state neoantigen-reactive T cells were absent
before vaccination. Unlike Braun, this table carries no pre-vaccine (week-0) column, so
de-novo here is author-asserted, not re-verified from a baseline; recognition-evidence
temporal_clarity is set lower accordingly (evidence strength is not flattened).

The workbook is ~2.2 GB uncompressed (mostly bulk TCR-repertoire sheets that are not
recognition labels); only the small recognition sheets are read, via a streaming reader
that never materialises the giant shared-string table or the TCR sheets.
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from xml.etree import ElementTree as ET
import zipfile

import pandas as pd

from benchmark.funnel import ReachabilityStatus
from event_b.adapters.base import AdapterDeclaration
from event_b.corpus import EventBCorpus
from event_b.manifest import SourceManifest
from event_b.models import (
    AssayType,
    AvailabilityStatus,
    BiologicalEvent,
    EvidenceFamily,
    InformationTiming,
    MHCClass,
    ResponseLabel,
    ReviewStatus,
    SCHEMAS,
    ValueOrigin,
    VaccineInclusion,
    stable_candidate_id,
)
from event_b.review import ReviewIssue


STUDY_ID = "hu_neovax_2021"
NCT = "NCT01970358"
DOI = "10.1038/s41591-020-01206-4"
PMCID = "PMC8273876"
OTT_PMCID = "PMC5577644"  # cross-check: original per-peptide screen for patients 1-6

DATASET_FILE = "NIHMS1707651-supplement-Suppl_DataSet.xlsx"
EXPECTED_SHA256 = {
    DATASET_FILE: "9ed88622f0b8dcc63ef7c345dd2d2d86bc086b600e222ae628a648ce1eae1034",
}
# Manual placement (author-manuscript supplement behind a proof-of-work download gate).
MANUAL_SUBDIR = "manual"
SOURCE_PAGE = f"https://pmc.ncbi.nlm.nih.gov/articles/{PMCID}/"
DOWNLOAD_URL = (
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/8273876/bin/"
    "NIHMS1707651-supplement-Suppl_DataSet.xlsx"
)

# Worksheet numbers (xl/worksheets/sheetN.xml) inside the workbook.
SHEET_CD8 = 6          # Suppl Dataset 4a — CD8 vaccine-peptide reactivity
SHEET_CD4 = 7          # Suppl Dataset 4b — CD4 vaccine-peptide reactivity
SHEET_SPREADING = (41, 42, 43)  # Suppl Dataset 11a/b/c — class-I epitope spreading
# 11d (44) and 11e (45) are MSEC-nominated, wild-type-ligand only (11e has #REF!-corrupted
# cells); deferred rather than mis-parsed. Dataset 12 is class-II spreading *predictions*
# (no measured reactivity) and is not ingested as labels. Both are declared limitations.

# Ott 2017 reconciliation anchors (patients 1-6, week-16 readout).
OTT_CD8_NEOANTIGENS = (15, 97)
OTT_CD4_NEOANTIGENS = (58, 97)
OTT_COHORT = {"1", "2", "3", "4", "5", "6"}

_REACHED = ReachabilityStatus.REACHED.value
_NOT_ASSESSED = ReachabilityStatus.NOT_ASSESSED.value

# Differentiated reliability vectors. Distinct channels get distinct strengths so a CD8
# minimal-epitope call, a pre-stimulation-amplified CD4 call, and an epitope-spreading
# response are never flattened into identical-strength evidence.
_RELIABILITY = {
    "CD8_VACCINE": {
        "patient_specificity": 1.0,
        "functional_relevance": 1.0,
        "vaccine_relevance": 1.0,
        "candidate_specificity": 1.0,   # specific minimal epitope + restricting HLA
        "assay_directness": 0.9,        # peptide-pulsed autologous APC
        "temporal_clarity": 0.7,        # de-novo author-asserted; no week-0 baseline in table
        "source_completeness": 0.9,
    },
    "CD4_VACCINE_EXVIVO": {
        "patient_specificity": 1.0,
        "functional_relevance": 1.0,
        "vaccine_relevance": 1.0,
        "candidate_specificity": 0.6,   # overlapping ~15mer; no single minimal epitope/HLA
        "assay_directness": 0.9,        # ex-vivo (no in-vitro expansion)
        "temporal_clarity": 0.7,
        "source_completeness": 0.9,
    },
    "CD4_VACCINE_PRESTIM": {
        "patient_specificity": 1.0,
        "functional_relevance": 1.0,
        "vaccine_relevance": 1.0,
        "candidate_specificity": 0.6,
        "assay_directness": 0.6,        # after in-vitro pre-stimulation (amplified)
        "temporal_clarity": 0.7,
        "source_completeness": 0.9,
    },
    "SPREADING": {
        "patient_specificity": 1.0,
        "functional_relevance": 1.0,
        "vaccine_relevance": 0.0,       # explicitly NOT a vaccine candidate
        "candidate_specificity": 1.0,
        "assay_directness": 0.9,
        "temporal_clarity": 0.6,        # spreading onset timing less precisely bracketed
        "source_completeness": 0.8,
    },
    "PERSISTENCE": {
        "patient_specificity": 1.0,
        "functional_relevance": 1.0,
        "vaccine_relevance": 1.0,
        "candidate_specificity": 0.8,
        "assay_directness": 0.7,
        "temporal_clarity": 0.8,        # clear late (3-4.5 yr) measurement
        "source_completeness": 0.8,
    },
}


# --------------------------------------------------------------------------------------
# Source staging + streaming reader
# --------------------------------------------------------------------------------------
def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def hu_source_path(raw_dir: str | Path) -> Path:
    return Path(raw_dir) / MANUAL_SUBDIR / DATASET_FILE


def _manual_placement_message(path: Path) -> str:
    return "\n".join(
        [
            "Hu 2021 NeoVax supplement is missing or its checksum does not match.",
            "It is a PMC author manuscript behind a JavaScript proof-of-work download gate,",
            "so it cannot be fetched programmatically. Download it in a browser from:",
            f"    {DOWNLOAD_URL}",
            f"(or via the article page {SOURCE_PAGE}), then place it EXACTLY at:",
            f"    {path}",
            f"    expected sha256={EXPECTED_SHA256[DATASET_FILE]}",
            "No records are fabricated; ingestion refuses to proceed without the verified source.",
        ]
    )


def stage_hu_supplement(raw_dir: str | Path) -> Path:
    """Return the verified path to the manually-placed workbook or fail with guidance."""
    path = hu_source_path(raw_dir)
    if not path.exists() or sha256_file(path) != EXPECTED_SHA256[DATASET_FILE]:
        raise RuntimeError(_manual_placement_message(path))
    return path


def hu_source_paths(raw_dir: str | Path) -> list[Path]:
    return [stage_hu_supplement(raw_dir)]


_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_ROW, _C, _V, _T, _IS, _SI = _NS + "row", _NS + "c", _NS + "v", _NS + "t", _NS + "is", _NS + "si"


def _col_idx(ref: str) -> int:
    letters = "".join(ch for ch in ref if ch.isalpha())
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch.upper()) - 64)
    return n - 1


def read_workbook_sheets(path: str | Path, sheet_numbers) -> dict[int, list[list[str]]]:
    """Read specific worksheets from a huge xlsx without loading the whole file.

    One streaming pass parses the (small) target sheets into sparse cells and collects the
    shared-string indices they reference; a single streaming pass over sharedStrings.xml
    resolves only those indices (stopping at the max needed), so the multi-hundred-MB
    shared-string table and the giant bulk sheets are never fully materialised.
    """
    archive = zipfile.ZipFile(path)
    raw_sheets: dict[int, list[list[tuple[int, str, object]]]] = {}
    wanted: set[int] = set()
    for number in sheet_numbers:
        part = f"xl/worksheets/sheet{number}.xml"
        rows: list[list[tuple[int, str, object]]] = []
        with archive.open(part) as handle:
            for _, element in ET.iterparse(handle, events=("end",)):
                if element.tag != _ROW:
                    continue
                cells: list[tuple[int, str, object]] = []
                for i, cell in enumerate(element.iter(_C)):
                    ref = cell.get("r", "")
                    column = _col_idx(ref) if ref else i
                    kind = cell.get("t")
                    if kind == "s":
                        value = cell.find(_V)
                        index = int(value.text) if (value is not None and value.text) else None
                        if index is not None:
                            wanted.add(index)
                        cells.append((column, "s", index))
                    elif kind == "inlineStr":
                        inline = cell.find(_IS)
                        text = (
                            "".join(t.text or "" for t in inline.iter(_T))
                            if inline is not None
                            else ""
                        )
                        cells.append((column, "x", text))
                    else:
                        value = cell.find(_V)
                        cells.append((column, "x", value.text if value is not None else ""))
                rows.append(cells)
                element.clear()
        raw_sheets[number] = rows

    strings: dict[int, str] = {}
    if wanted:
        maximum, index = max(wanted), -1
        with archive.open("xl/sharedStrings.xml") as handle:
            for _, element in ET.iterparse(handle, events=("end",)):
                if element.tag != _SI:
                    continue
                index += 1
                if index in wanted:
                    strings[index] = "".join(t.text or "" for t in element.iter(_T))
                element.clear()
                if index >= maximum:
                    break

    materialised: dict[int, list[list[str]]] = {}
    for number, rows in raw_sheets.items():
        grid: list[list[str]] = []
        for cells in rows:
            if not cells:
                grid.append([])
                continue
            width = max(column for column, _, _ in cells) + 1
            row = [""] * width
            for column, kind, payload in cells:
                row[column] = strings.get(payload, "") if kind == "s" else (str(payload) if payload is not None else "")
            grid.append(row)
        materialised[number] = grid
    return materialised


# --------------------------------------------------------------------------------------
# Cell / header helpers
# --------------------------------------------------------------------------------------
def _cell(row: list[str], column: int) -> str:
    return row[column].strip() if 0 <= column < len(row) else ""


def _is_patient(value: str) -> bool:
    text = str(value).strip()
    return bool(text) and text.replace(".", "", 1).isdigit()


def _patient_number(value: str) -> str:
    return str(int(float(value)))


def _peptide(value: str) -> str:
    text = str(value).strip().upper()
    return text if text and text.isalpha() else ""


def _normalize_react(value: str) -> str:
    """Map an ELISpot cell to POS / NEG / ND (not-done / unscored)."""
    text = str(value).strip().lower()
    if text in {"1", "1.0"}:
        return "POS"
    if text in {"0", "0.0"}:
        return "NEG"
    return "ND"


def _split_header_data(rows: list[list[str]]) -> tuple[list[list[str]], list[list[str]]]:
    for i, row in enumerate(rows):
        if row and _is_patient(row[0]):
            return rows[:i], rows[i:]
    return rows, []


def _header_strings(header_rows: list[list[str]]) -> dict[int, str]:
    """Concatenate the header text stacked above each column."""
    combined: dict[int, list[str]] = {}
    for row in header_rows:
        for column, value in enumerate(row):
            text = str(value).strip()
            if text:
                combined.setdefault(column, []).append(text)
    return {column: " | ".join(parts) for column, parts in combined.items()}


def _find_column(headers: dict[int, str], *keywords: str) -> int | None:
    for column in sorted(headers):
        text = headers[column].lower()
        if all(keyword.lower() in text for keyword in keywords):
            return column
    return None


def _hla_list(value: str) -> str:
    text = str(value).strip()
    # Drop trailing artefacts like "B56:01 Rank" seen in a couple of spreading rows.
    token = text.split()[0] if text else ""
    return json.dumps([token]) if token else json.dumps([])


def _prov(
    entity: str,
    entity_id: str,
    *,
    document: str,
    table: str,
    row: int,
    fragment: str,
    origin: str = ValueOrigin.SOURCE_REPORTED.value,
) -> dict:
    provenance_id = "prov:" + sha256(f"{entity}|{entity_id}".encode()).hexdigest()[:20]
    return {
        "provenance_id": provenance_id,
        "entity_type": entity,
        "entity_id": entity_id,
        "field_name": "*",
        "source_document": document,
        "table": table,
        "row": row,
        "source_fragment": fragment,
        "extraction_method": "streaming_xlsx_adapter",
        "extraction_confidence": 1.0,
        "value_origin": origin,
        "review_status": ReviewStatus.ACCEPTED.value,
    }


class HuNeoVaxAdapter:
    """Hu 2021 melanoma NeoVax adapter (CD8/CD4 vaccine recognition + epitope spreading)."""

    declaration = AdapterDeclaration(
        "Hu 2021 melanoma NeoVax personal neoantigen vaccine",
        f"{PMCID}/NatMed-2021",
        "hu_neovax_event_b",
        "1.0.0",
        (
            "studies",
            "patients",
            "vaccines",
            "candidates",
            "assays",
            "recognition_evidence",
            "candidate_funnel_links",
            "provenance",
        ),
        (
            BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
            BiologicalEvent.EPITOPE_SPREADING.value,
        ),
        (
            "The source is a ~2.2 GB author-manuscript workbook; only the small recognition sheets "
            "(Datasets 4a/4b vaccine, 11a-c epitope spreading) are read via a streaming reader. The "
            "bulk TCR-repertoire sheets (Dataset 8) are not recognition labels and are not ingested.",
            "ELISpot calls are the source's own 0/1 (>=2.5x DMSO, Hu 2021); the table carries no raw "
            "replicates, so calls are ingested as reported and reconciled against Ott 2017's reported "
            "neoantigen totals rather than recomputed.",
            "The table has no pre-vaccine (week-0) column; de-novo status is author-asserted (Ott/Hu: "
            "no pre-vaccination reactivity), not re-verified from a baseline, reflected in a lower "
            "recognition_evidence temporal_clarity than a baseline-verified de-novo claim.",
            "Per-peptide class-II restriction for CD4 assay peptides is not resolved in the source; "
            "CD4 candidate hla_alleles are left empty and mhc_class is CLASS_II.",
            "Patient HLA genotypes are not in this table; patient hla_alleles left empty (the HLA-subset "
            "check is skipped, not fabricated). Class-I restricting alleles are stored per candidate.",
            "Epitope-spreading Datasets 11d/11e (MSEC-nominated, wild-type-ligand only; 11e has "
            "#REF!-corrupted cells) are deferred, not mis-parsed; Dataset 12 (class-II spreading "
            "predictions, no measured reactivity) is not ingested as labels.",
        ),
        (
            "CD8 rows (Dataset 4a) are class-I: the tested entity is the predicted minimal epitope "
            "(mutant_peptide) with its restricting HLA; the parent long immunizing peptide is kept in "
            "candidate_source, preserving the long-peptide/minimal-epitope relationship.",
            "CD4 rows (Dataset 4b) are class-II: the tested entity is the overlapping assay peptide; the "
            "week-16 label is positive if reactive ex-vivo or after pre-stimulation, with the assay "
            "condition recorded so ex-vivo and pre-stimulation evidence keep distinct assay_directness.",
            "Vaccine-included peptides (4a/4b) are EVENT_B; responses to non-vaccine neoantigens (11a-c) "
            "are EPITOPE_SPREADING with vaccine_inclusion NOT_INCLUDED, never a vaccine-candidate label.",
            "Patients 1-6 are the Ott 2017 cohort and 11-12 are new; a single consolidated source avoids "
            "cross-source patient double-counting.",
            "Week-16 is the primary recognition timepoint; a scored 3-4.5-year readout is recorded as "
            "LONGITUDINAL_PERSISTENCE recognition evidence, not a separate label.",
        ),
        (
            "patient_hla_genotype",
            "pre_vaccine_baseline",
            "cd4_class_ii_restriction",
            "raw_elispot_replicates",
            "per_peptide_clinical_outcome",
        ),
    )

    def __init__(self, raw_dir: str | Path):
        self.raw_dir = Path(raw_dir)
        self.review_issues: list[ReviewIssue] = []
        # Accumulators, reset per normalize().
        self._studies: dict[str, dict] = {}
        self._patients: dict[str, dict] = {}
        self._vaccines: dict[str, dict] = {}
        self._candidates: list[dict] = []
        self._assays: list[dict] = []
        self._evidence: list[dict] = []
        self._funnel: list[dict] = []
        self._provenance: list[dict] = []
        self._seen_candidates: set[str] = set()
        self._document = f"{PMCID}:{DATASET_FILE}"

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        path = stage_hu_supplement(self.raw_dir)
        sheets = read_workbook_sheets(path, (SHEET_CD8, SHEET_CD4, *SHEET_SPREADING))
        return {
            "cd8": sheets[SHEET_CD8],
            "cd4": sheets[SHEET_CD4],
            "spreading": {n: sheets[n] for n in SHEET_SPREADING},
        }

    # -- entity ensures -----------------------------------------------------------------
    def _ensure_patient(self, source_patient: str, manifest: SourceManifest) -> tuple[str, str]:
        patient_id = f"{STUDY_ID}:{source_patient}"
        vaccine_id = f"{STUDY_ID}:{source_patient}:pcv"
        if patient_id in self._patients:
            return patient_id, vaccine_id
        cohort = "Ott-2017 cohort" if source_patient in OTT_COHORT else "Hu-2021 additional"
        patient_prov = _prov(
            "patients",
            patient_id,
            document=self._document,
            table="Dataset 3/4",
            row=0,
            fragment=f"source_patient={source_patient}; {cohort}; melanoma NeoVax",
        )
        self._provenance.append(patient_prov)
        self._patients[patient_id] = {
            "patient_id": patient_id,
            "source_patient_id": source_patient,
            "study_id": STUDY_ID,
            "cancer_type": "melanoma",
            "treatment_context": cohort,
            "provenance_id": patient_prov["provenance_id"],
        }
        vaccine_prov = _prov(
            "vaccines",
            vaccine_id,
            document=self._document,
            table="-",
            row=0,
            fragment=f"personalized long-peptide NeoVax + poly-ICLC; patient {source_patient}",
        )
        self._provenance.append(vaccine_prov)
        self._vaccines[vaccine_id] = {
            "vaccine_id": vaccine_id,
            "patient_id": patient_id,
            "study_id": STUDY_ID,
            "vaccine_platform": "synthetic long peptide (SLP) + poly-ICLC",
            "formulation": "up to 20 long peptides in 4 pools admixed with poly-ICLC",
            "mhc_class_intent": MHCClass.BOTH.value,
            "concurrent_therapy": "none",
            "provenance_id": vaccine_prov["provenance_id"],
        }
        return patient_id, vaccine_id

    def _ensure_study(self, manifest: SourceManifest) -> None:
        if STUDY_ID in self._studies:
            return
        study_prov = _prov(
            "studies",
            STUDY_ID,
            document=self._document,
            table="-",
            row=0,
            fragment=(
                f"{NCT}; melanoma personal neoantigen vaccine (NeoVax); follow-up of Ott 2017; "
                "de-novo Event-B is author-asserted (no pre-vaccine reactivity reported)"
            ),
        )
        self._provenance.append(study_prov)
        self._studies[STUDY_ID] = {
            "study_id": STUDY_ID,
            "title": (
                "Personal neoantigen vaccines induce persistent memory T cell responses and "
                "epitope spreading in patients with melanoma"
            ),
            "publication_ids": f"DOI:{DOI}; {PMCID}; cross-check {OTT_PMCID}",
            "trial_id": NCT,
            "cancer_type": "melanoma",
            "vaccine_platform": "synthetic long peptide (SLP) + poly-ICLC",
            "adjuvant": "poly-ICLC (Hiltonol)",
            "vaccination_schedule": "prime series then boosts (Ott 2017 NeoVax schedule)",
            "source_urls": SOURCE_PAGE,
            "source_paths": str(hu_source_path(self.raw_dir)),
            "source_checksums": EXPECTED_SHA256[DATASET_FILE],
            "source_manifest_id": manifest.manifest_id,
            "provenance_id": study_prov["provenance_id"],
        }

    # -- emission -----------------------------------------------------------------------
    def _emit(
        self,
        *,
        patient_id: str,
        vaccine_id: str,
        source_patient: str,
        gene: str,
        protein_change: str,
        peptide: str,
        wildtype: str,
        hla_json: str,
        mhc_class: str,
        parent_imp: str,
        vaccine_inclusion: str,
        event_type: str,
        label: str,
        reliability_key: str,
        source_dataset: str,
        table: str,
        source_row: int,
        detail: str,
        persistence: str | None = None,
    ) -> None:
        """Append candidate + week-16 assay + recognition evidence + funnel + provenance."""
        identity = {
            "study_id": STUDY_ID,
            "patient_id": patient_id,
            "sample_id": "",
            "timepoint": "",
            "genomic_variant": f"{gene}|{protein_change}",
            "transcript": "",
            "mutant_peptide": peptide,
            "hla_alleles": hla_json,
        }
        candidate_id = stable_candidate_id(identity)
        if candidate_id in self._seen_candidates:
            self.review_issues.append(
                ReviewIssue.create(
                    "candidates",
                    candidate_id,
                    "DUPLICATE_IDENTITY",
                    f"{source_dataset}: {gene}|{protein_change}|{peptide} repeats an existing identity",
                )
            )
            return
        self._seen_candidates.add(candidate_id)

        candidate_prov = _prov(
            "candidates",
            candidate_id,
            document=self._document,
            table=table,
            row=source_row,
            fragment=f"{gene}|{protein_change}|{peptide}|hla={hla_json}|parent={parent_imp}",
        )
        self._provenance.append(candidate_prov)
        self._candidates.append(
            {
                "candidate_id": candidate_id,
                "patient_id": patient_id,
                "study_id": STUDY_ID,
                "genomic_variant": f"{gene}|{protein_change}" if gene or protein_change else pd.NA,
                "gene": gene or pd.NA,
                "protein_change": protein_change or pd.NA,
                "mutant_peptide": peptide,
                "wildtype_peptide": wildtype or pd.NA,
                "peptide_length": len(peptide),
                "hla_alleles": hla_json,
                "mhc_class": mhc_class,
                "candidate_source": (
                    f"{source_dataset}; parent immunizing peptide {parent_imp or 'n/a'}"
                    if vaccine_inclusion == VaccineInclusion.INCLUDED.value
                    else f"{source_dataset}; non-vaccine neoantigen (epitope spreading)"
                ),
                "vaccine_inclusion": vaccine_inclusion,
                "vaccine_inclusion_origin": ValueOrigin.SOURCE_REPORTED.value,
                "generation_provenance": f"parent_imp={parent_imp}" if parent_imp else pd.NA,
                "mutant_wildtype_verified": True,
                "provenance_id": candidate_prov["provenance_id"],
            }
        )

        funnel_prov = _prov(
            "candidate_funnel_links",
            candidate_id,
            document=self._document,
            table=table,
            row=source_row,
            fragment="assayed for T-cell recognition; upstream reachability not reported",
        )
        self._provenance.append(funnel_prov)
        self._funnel.append(
            {
                "funnel_link_id": "funnel:" + sha256(candidate_id.encode()).hexdigest()[:20],
                "candidate_id": candidate_id,
                "patient_id": patient_id,
                "study_id": STUDY_ID,
                "mutation_called": _REACHED,
                "transcript_represented": _NOT_ASSESSED,
                "peptide_generated": _REACHED,
                "survives_gating": _NOT_ASSESSED,
                "hla_included": _REACHED if mhc_class == MHCClass.CLASS_I.value else _NOT_ASSESSED,
                "presentation_candidate": _NOT_ASSESSED,
                "ranking_stage": _NOT_ASSESSED,
                "top_k": _NOT_ASSESSED,
                "recognition_scored": _REACHED,
                "vaccine_inclusion": (
                    _REACHED if vaccine_inclusion == VaccineInclusion.INCLUDED.value else _NOT_ASSESSED
                ),
                "functional_assay": _REACHED,
                "provenance_id": funnel_prov["provenance_id"],
            }
        )

        qualitative = {"POS": "POSITIVE", "NEG": "NEGATIVE", "ND": "UNSCORED"}[label]
        response_label = {
            "POS": ResponseLabel.POSITIVE.value,
            "NEG": ResponseLabel.TESTED_NEGATIVE.value,
            "ND": ResponseLabel.UNTESTED.value,
        }[label]
        review_status = (
            ReviewStatus.NEEDS_REVIEW.value if label == "ND" else ReviewStatus.ACCEPTED.value
        )
        if label == "ND":
            self.review_issues.append(
                ReviewIssue.create(
                    "candidates",
                    candidate_id,
                    "UNSCORABLE_ASSAY",
                    f"{source_dataset}: {gene}|{protein_change} recorded n.d. at week 16",
                )
            )

        assay_id = "assay:" + sha256(f"hu|{candidate_id}|wk16".encode()).hexdigest()[:20]
        assay_prov = _prov(
            "assays",
            assay_id,
            document=self._document,
            table=table,
            row=source_row,
            fragment=f"{detail}; week16 -> {response_label}",
            origin=ValueOrigin.SOURCE_REPORTED.value,
        )
        self._provenance.append(assay_prov)
        self._assays.append(
            {
                "assay_id": assay_id,
                "patient_id": patient_id,
                "study_id": STUDY_ID,
                "candidate_id": candidate_id,
                "vaccine_id": vaccine_id,
                "assay_type": AssayType.ELISPOT.value,
                "sample_type": "PBMC, IFN-gamma ELISpot (peptide-pulsed autologous APC)",
                "timepoint": "week16_post_vaccine",
                "relative_to_vaccine": "POST_VACCINE",
                "stimulation_protocol": detail,
                "positivity_threshold": (
                    "author-reported ELISpot call (>=2.5x DMSO control, Hu 2021); ingested as reported"
                ),
                "qualitative_result": qualitative,
                "source_interpretation": (
                    "recognized" if label == "POS" else ("assayed, not recognized" if label == "NEG" else "not done")
                ),
                "event_type": event_type,
                "response_label": response_label,
                "explicit_assay_inclusion": True,
                "review_status": review_status,
                "provenance_id": assay_prov["provenance_id"],
            }
        )

        # Week-16 recognition evidence with a differentiated reliability vector. An n.d.
        # (not-done) row carries no evidence, so none is recorded for it.
        if label in {"POS", "NEG"}:
            self._append_evidence(
                candidate_id,
                patient_id,
                evidence_family=(
                    EvidenceFamily.VACCINE_EVENT_B.value
                    if event_type == BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
                    else EvidenceFamily.FUNCTIONAL_T_CELL_ASSAY.value
                ),
                source_dataset=source_dataset,
                value=qualitative,
                reliability=_RELIABILITY[reliability_key],
                tag="wk16",
            )
        # Optional 3-4.5-year persistence, only when actually scored.
        if persistence in {"POS", "NEG"}:
            self._append_evidence(
                candidate_id,
                patient_id,
                evidence_family=EvidenceFamily.LONGITUDINAL_PERSISTENCE.value,
                source_dataset=source_dataset,
                value="POSITIVE" if persistence == "POS" else "NEGATIVE",
                reliability=_RELIABILITY["PERSISTENCE"],
                tag="yr3_4.5",
            )

    def _append_evidence(
        self,
        candidate_id: str,
        patient_id: str,
        *,
        evidence_family: str,
        source_dataset: str,
        value: str,
        reliability: dict,
        tag: str,
    ) -> None:
        evidence_id = "evid:" + sha256(f"{candidate_id}|{tag}".encode()).hexdigest()[:20]
        evidence_prov = _prov(
            "recognition_evidence",
            evidence_id,
            document=self._document,
            table=source_dataset,
            row=0,
            fragment=f"{source_dataset}:{tag}={value}",
            origin=ValueOrigin.SOURCE_REPORTED.value,
        )
        self._provenance.append(evidence_prov)
        self._evidence.append(
            {
                "evidence_id": evidence_id,
                "candidate_id": candidate_id,
                "patient_id": patient_id,
                "evidence_family": evidence_family,
                "source_dataset": f"{PMCID}:{source_dataset}",
                "measured_or_predicted": "MEASURED",
                "value": value,
                "units": "IFN-gamma ELISpot qualitative call",
                "directionality": "positive_indicates_recognition",
                "availability_status": AvailabilityStatus.AVAILABLE.value,
                "information_timing": InformationTiming.OUTCOME_ONLY.value,
                "evidence_quality": "author_reported_binary_elispot",
                **reliability,
                "provenance_id": evidence_prov["provenance_id"],
            }
        )

    # -- channel parsers ----------------------------------------------------------------
    def _parse_cd8(self, rows: list[list[str]], manifest: SourceManifest) -> None:
        header_rows, data = _split_header_data(rows)
        headers = _header_strings(header_rows)
        col_hla = _find_column(headers, "hla allele")
        col_mut = _find_column(headers, "mutated peptide")
        col_wt = _find_column(headers, "wild type peptide")
        col_imp = _find_column(headers, "immunizing peptide")
        col_pool = _find_column(headers, "immunizing pool")
        col_r16 = _find_column(headers, "16 week")
        col_ryr = _find_column(headers, "3-4.5") or _find_column(headers, "4-4.5")
        if None in (col_hla, col_mut, col_wt, col_r16):
            raise ValueError("Hu Dataset 4a (CD8) header layout not recognised")
        for offset, row in enumerate(data):
            if not (row and _is_patient(row[0])):
                continue
            peptide = _peptide(_cell(row, col_mut))
            if not peptide:
                continue
            source_patient = _patient_number(row[0])
            patient_id, vaccine_id = self._ensure_patient(source_patient, manifest)
            gene = _cell(row, _find_column(headers, "gene") or 2)
            protein_change = _cell(row, _find_column(headers, "protein change") or 3)
            hla_json = _hla_list(_cell(row, col_hla))
            wildtype = _peptide(_cell(row, col_wt))
            imp = _cell(row, col_imp) if col_imp is not None else ""
            pool = _cell(row, col_pool) if col_pool is not None else ""
            mhc_class = MHCClass.CLASS_I.value if 8 <= len(peptide) <= 14 else MHCClass.UNKNOWN.value
            self._emit(
                patient_id=patient_id,
                vaccine_id=vaccine_id,
                source_patient=source_patient,
                gene=gene,
                protein_change=protein_change,
                peptide=peptide,
                wildtype=wildtype if wildtype != peptide else "",
                hla_json=hla_json,
                mhc_class=mhc_class,
                parent_imp=imp,
                vaccine_inclusion=VaccineInclusion.INCLUDED.value,
                event_type=BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                label=_normalize_react(_cell(row, col_r16)),
                reliability_key="CD8_VACCINE",
                source_dataset="Dataset_4a_CD8",
                table="Suppl Dataset 4a",
                source_row=offset,
                detail=f"CD8 minimal epitope, HLA {json.loads(hla_json)}, pool {pool}",
                persistence=(
                    _normalize_react(_cell(row, col_ryr)) if col_ryr is not None else None
                ),
            )

    def _parse_cd4(self, rows: list[list[str]], manifest: SourceManifest) -> None:
        header_rows, data = _split_header_data(rows)
        headers = _header_strings(header_rows)
        col_gene = _find_column(headers, "gene") or 1
        col_change = _find_column(headers, "protein change") or 2
        col_asp = _find_column(headers, "assay peptide")
        col_imp = _find_column(headers, "immunizing peptide")
        # Week-16 conditions are at fixed offsets under the "Week 16" header; validate then use.
        col_ex16, col_ps16, col_exyr, col_psyr = 9, 10, 14, 15
        if col_asp is None:
            raise ValueError("Hu Dataset 4b (CD4) header layout not recognised")
        for offset, row in enumerate(data):
            if not (row and _is_patient(row[0])):
                continue
            # The assay-peptide sequence is the column after its ID; find by keyword then +1.
            peptide = _peptide(_cell(row, col_asp + 1)) or _peptide(_cell(row, 8))
            if not peptide:
                continue
            source_patient = _patient_number(row[0])
            patient_id, vaccine_id = self._ensure_patient(source_patient, manifest)
            gene = _cell(row, col_gene)
            protein_change = _cell(row, col_change)
            imp = _cell(row, col_imp + 1) if col_imp is not None else ""
            ex16 = _normalize_react(_cell(row, col_ex16))
            ps16 = _normalize_react(_cell(row, col_ps16))
            # CD4 week-16 label: positive if reactive ex-vivo or after pre-stimulation.
            if ex16 == "POS" or ps16 == "POS":
                label = "POS"
                reliability_key = "CD4_VACCINE_EXVIVO" if ex16 == "POS" else "CD4_VACCINE_PRESTIM"
            elif ex16 == "NEG" or ps16 == "NEG":
                label = "NEG"
                reliability_key = "CD4_VACCINE_PRESTIM"
            else:
                label = "ND"
                reliability_key = "CD4_VACCINE_PRESTIM"
            mhc_class = MHCClass.CLASS_II.value if 12 <= len(peptide) <= 30 else MHCClass.UNKNOWN.value
            exyr = _normalize_react(_cell(row, col_exyr))
            psyr = _normalize_react(_cell(row, col_psyr))
            persistence = "POS" if "POS" in (exyr, psyr) else ("NEG" if "NEG" in (exyr, psyr) else None)
            self._emit(
                patient_id=patient_id,
                vaccine_id=vaccine_id,
                source_patient=source_patient,
                gene=gene,
                protein_change=protein_change,
                peptide=peptide,
                wildtype="",
                hla_json=json.dumps([]),
                mhc_class=mhc_class,
                parent_imp=imp,
                vaccine_inclusion=VaccineInclusion.INCLUDED.value,
                event_type=BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                label=label,
                reliability_key=reliability_key,
                source_dataset="Dataset_4b_CD4",
                table="Suppl Dataset 4b",
                source_row=offset,
                detail=f"CD4 assay peptide; wk16 ex-vivo={ex16}, pre-stim={ps16}",
                persistence=persistence,
            )

    def _parse_spreading(self, sheets: dict[int, list[list[str]]], manifest: SourceManifest) -> None:
        for number, rows in sheets.items():
            header_rows, data = _split_header_data(rows)
            headers = _header_strings(header_rows)
            col_hla = _find_column(headers, "hla allele")
            col_mut = _find_column(headers, "mutated peptide")
            col_wt = _find_column(headers, "wild type peptide")
            col_gene = _find_column(headers, "gene")
            col_change = _find_column(headers, "protein change")
            col_r16 = _find_column(headers, "16 week")
            if None in (col_hla, col_mut, col_r16):
                # 11d/11e (MSEC, wild-type-ligand only) do not expose a mutant epitope column;
                # they are declared out of scope rather than mis-parsed.
                self.review_issues.append(
                    ReviewIssue.create(
                        "recognition_evidence",
                        f"sheet{number}",
                        "SPREADING_LAYOUT_DEFERRED",
                        f"epitope-spreading sheet{number} lacks a mutant-epitope/reactivity layout; deferred",
                    )
                )
                continue
            for offset, row in enumerate(data):
                if not (row and _is_patient(row[0])):
                    continue
                peptide = _peptide(_cell(row, col_mut))
                if not peptide:
                    continue
                source_patient = _patient_number(row[0])
                patient_id, vaccine_id = self._ensure_patient(source_patient, manifest)
                gene = _cell(row, col_gene) if col_gene is not None else ""
                protein_change = _cell(row, col_change) if col_change is not None else ""
                hla_json = _hla_list(_cell(row, col_hla))
                wildtype = _peptide(_cell(row, col_wt)) if col_wt is not None else ""
                mhc_class = MHCClass.CLASS_I.value if 8 <= len(peptide) <= 14 else MHCClass.UNKNOWN.value
                self._emit(
                    patient_id=patient_id,
                    vaccine_id=vaccine_id,
                    source_patient=source_patient,
                    gene=gene,
                    protein_change=protein_change,
                    peptide=peptide,
                    wildtype=wildtype if wildtype != peptide else "",
                    hla_json=hla_json,
                    mhc_class=mhc_class,
                    parent_imp="",
                    vaccine_inclusion=VaccineInclusion.NOT_INCLUDED.value,
                    event_type=BiologicalEvent.EPITOPE_SPREADING.value,
                    label=_normalize_react(_cell(row, col_r16)),
                    reliability_key="SPREADING",
                    source_dataset=f"Dataset_11_sheet{number}",
                    table=f"Suppl Dataset 11 (sheet{number})",
                    source_row=offset,
                    detail=f"epitope spreading (non-vaccine neoantigen), HLA {json.loads(hla_json)}",
                    persistence=None,
                )

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        self.review_issues = []
        self._studies, self._patients, self._vaccines = {}, {}, {}
        self._candidates, self._assays, self._evidence = [], [], []
        self._funnel, self._provenance = [], []
        self._seen_candidates = set()

        self._ensure_study(manifest)
        self._parse_cd8(list(extracted["cd8"]), manifest)
        self._parse_cd4(list(extracted["cd4"]), manifest)
        self._parse_spreading(dict(extracted["spreading"]), manifest)

        # Record candidate counts per vaccine now that all peptides are seen.
        counts: dict[str, int] = {}
        for candidate in self._candidates:
            if candidate["vaccine_inclusion"] == VaccineInclusion.INCLUDED.value:
                counts[candidate["patient_id"]] = counts.get(candidate["patient_id"], 0) + 1
        for vaccine in self._vaccines.values():
            vaccine["candidate_count"] = counts.get(vaccine["patient_id"], 0)

        return EventBCorpus(
            studies=SCHEMAS["studies"].normalize(pd.DataFrame(list(self._studies.values()))),
            patients=SCHEMAS["patients"].normalize(pd.DataFrame(list(self._patients.values()))),
            vaccines=SCHEMAS["vaccines"].normalize(pd.DataFrame(list(self._vaccines.values()))),
            candidates=SCHEMAS["candidates"].normalize(pd.DataFrame(self._candidates)),
            assays=SCHEMAS["assays"].normalize(pd.DataFrame(self._assays)),
            recognition_evidence=SCHEMAS["recognition_evidence"].normalize(
                pd.DataFrame(self._evidence)
            ),
            candidate_funnel_links=SCHEMAS["candidate_funnel_links"].normalize(
                pd.DataFrame(self._funnel)
            ),
            provenance=SCHEMAS["provenance"].normalize(pd.DataFrame(self._provenance)),
        )
