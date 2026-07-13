"""Zhao 2026 Fukuoka DC-vaccine Event-B adapter.

Ingests 2,317 short (8-11mer, class I) SNV-derived neoantigen peptides administered
to 352 cancer patients, each scored by a single post-vaccination IFN-gamma ELISpot
fold-increase. The immunogenicity call is recomputed from the source ``ELSPOT ratio``
column using the paper's own stated rule; no counts are hard-coded, no negative is
inferred from omission, and repeated assays of an identical peptide/HLA within a
patient are preserved rather than collapsed.

Source (CC BY, open access):
    Zhao P et al. (incl. Morisaki T). Profiling immunogenic neoantigen peptides elicited
    by personalized neoantigen vaccine in cancer patients. Front Immunol 2026.
    doi:10.3389/fimmu.2026.1829509 ; PMID 42344930 ; PMC13286890.

Positivity rule (verbatim, Zhao 2026):
    "Immunogenic neoantigen peptides were defined as those inducing a >=2.0-fold increase
     in IFN-gamma ELISPOT responses after vaccination."

Event-B justification:
    Every peptide is a personalized vaccine component measured only after vaccination by
    post-vaccine IFN-gamma ELISpot; the fold-increase is post-vaccine over pre-vaccine
    baseline, so the readout is vaccine-induced (Event B) by construction.

Critical limitation (encoded in metadata, not silently dropped):
    Peptides were pre-selected for vaccination by an HLA-binding workflow. This is
    vaccinated-subset recognition discrimination, NOT full-denominator top-K ranking.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import shutil

import pandas as pd

from benchmark.funnel import ReachabilityStatus
from event_b.adapters.base import AdapterDeclaration
from event_b.corpus import EventBCorpus
from event_b.manifest import SourceManifest
from event_b.models import (
    AssayType,
    BiologicalEvent,
    MHCClass,
    ResponseLabel,
    ReviewStatus,
    SCHEMAS,
    ValueOrigin,
    VaccineInclusion,
    stable_candidate_id,
)
from event_b.review import ReviewIssue


STUDY_ID = "zhao_dc_2026"
DOI = "10.3389/fimmu.2026.1829509"
PMID = "42344930"
PMCID = "PMC13286890"
SUPPL_URL = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"
ARTICLE_URL = (
    "https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2026.1829509/full"
)

TABLE_FILE = "Table1.xlsx"
DATASHEET_FILE = "DataSheet1.pdf"
PEPTIDE_SHEET = "S1"

# Frozen content identity: a silently changed remote file can never be ingested undetected.
EXPECTED_SHA256 = {
    TABLE_FILE: "1fa76cf45435c39dc28e9d52e584d56938cee531dda953ad11cfbcc9c2617aea",
    DATASHEET_FILE: "55030d203fc284fdb505399e789cb7a8a832d1769029f820f78282ef733eaa09",
}

# Paper's stated positivity rule.
RATIO_POSITIVE_MIN = 2.0
# The source right-censors very large fold-increases as the literal string ">=5.0"; those
# are unambiguously >= the 2.0 positivity threshold, so they are POSITIVE with a recorded
# censor floor (never dropped, never treated as missing).
CENSORED_TOKEN = "≥5.0"  # '≥5.0'
CENSORED_TOKEN_ASCII = ">=5.0"
CENSORED_VALUE = 5.0

POSITIVITY_RULE = (
    ">=2.0-fold increase in IFN-gamma ELISPOT responses after vaccination "
    "(source 'ELSPOT ratio' >= 2.0; the right-censored '>=5.0' cell counts as >=2.0)"
)

_REACHED = ReachabilityStatus.REACHED.value
_NOT_ASSESSED = ReachabilityStatus.NOT_ASSESSED.value


def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _checksums_ok(paths: dict[str, Path]) -> bool:
    return all(
        path.exists() and sha256_file(path) == EXPECTED_SHA256[name]
        for name, path in paths.items()
    )


def _manual_placement_message(manual_dir: Path) -> str:
    lines = [
        "Could not fetch Zhao 2026 supplements automatically (no network?).",
        f"Download the open-access (CC BY) supplement bundle from Europe PMC:\n    {SUPPL_URL}",
        f"unzip it and place these files into:\n    {manual_dir}",
    ]
    lines += [f"    {name}  sha256={digest}" for name, digest in sorted(EXPECTED_SHA256.items())]
    lines.append("No records were fabricated; ingestion refuses to proceed without the source.")
    return "\n".join(lines)


def stage_zhao_supplements(raw_dir: str | Path) -> dict[str, Path]:
    """Return verified paths to the source files, fetching on cache miss.

    Local ``manual/`` cache first; download the Europe PMC supplementaryFiles zip on
    miss; actionable failure when offline. Every returned file is checksum-verified.
    """
    manual = Path(raw_dir) / "manual"
    paths = {name: manual / name for name in EXPECTED_SHA256}
    if _checksums_ok(paths):
        return paths
    manual.mkdir(parents=True, exist_ok=True)
    zip_path = Path(raw_dir) / f"{PMCID}_supplementaryFiles.zip"
    if not zip_path.exists():
        try:
            import urllib.request

            request = urllib.request.Request(SUPPL_URL, headers={"User-Agent": "epicurus-neo/1.0"})
            with urllib.request.urlopen(request, timeout=180) as response, zip_path.open("wb") as fh:
                shutil.copyfileobj(response, fh)
        except Exception as error:  # noqa: BLE001 - network failure surfaces as guidance
            raise RuntimeError(_manual_placement_message(manual)) from error
    from zipfile import ZipFile

    with ZipFile(zip_path) as archive:
        members = {Path(name).name: name for name in archive.namelist()}
        for name in EXPECTED_SHA256:
            if name in members:
                with archive.open(members[name]) as src, (manual / name).open("wb") as dst:
                    shutil.copyfileobj(src, dst)
    if not _checksums_ok(paths):
        observed = {name: (sha256_file(p) if p.exists() else "MISSING") for name, p in paths.items()}
        raise RuntimeError(
            "Zhao supplement checksums did not match after staging; refusing to ingest. "
            f"Observed: {observed}"
        )
    return paths


def zhao_source_paths(raw_dir: str | Path) -> list[Path]:
    """Ordered list of source files used to build the source manifest."""
    paths = stage_zhao_supplements(raw_dir)
    return [paths[name] for name in sorted(EXPECTED_SHA256)]


def read_peptide_sheet(path: str | Path) -> pd.DataFrame:
    """Read the S1 per-peptide sheet, forward-filling the merged patient/cancer cells."""
    frame = pd.read_excel(path, sheet_name=PEPTIDE_SHEET, header=1, engine="openpyxl")
    frame.columns = [str(column).strip() for column in frame.columns]
    # Patient ID and Cancer type are merged cells: only the first peptide row of each
    # patient carries the value. Forward-fill reconstructs the per-row patient identity.
    frame["Patient ID"] = frame["Patient ID"].ffill()
    frame["Cancer type"] = frame["Cancer type"].ffill()
    return frame


def _txt(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)) or value is pd.NA:
        return ""
    return str(value).strip()


def _clean_int_text(value: object) -> str:
    try:
        number = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return _txt(value)
    if pd.isna(number):
        return ""
    return str(int(number)) if float(number).is_integer() else _txt(value)


def parse_elspot_ratio(value: object) -> tuple[float | None, bool]:
    """Parse the source 'ELSPOT ratio' cell.

    Returns ``(ratio, censored)``. The literal '>=5.0' cell is right-censored: it maps
    to ``(5.0, True)``. A numeric cell maps to ``(float, False)``. An unparseable/blank
    cell maps to ``(None, False)`` and must be routed to review, never assumed negative.
    """
    text = _txt(value)
    if text in {CENSORED_TOKEN, CENSORED_TOKEN_ASCII}:
        return CENSORED_VALUE, True
    try:
        number = float(text)
    except ValueError:
        return None, False
    if pd.isna(number):
        return None, False
    return number, False


def immunogenic_label(value: object) -> tuple[str, float | None, bool]:
    """Apply the paper's rule. Returns ``(response_label, ratio, censored)``.

    POSITIVE iff ELSpot ratio >= 2.0 (including censored '>=5.0'); TESTED_NEGATIVE iff
    ratio < 2.0. A cell that cannot be scored yields UNTESTED for the caller to route to
    review rather than being coerced to a negative.
    """
    ratio, censored = parse_elspot_ratio(value)
    if ratio is None:
        return ResponseLabel.UNTESTED.value, None, False
    if ratio >= RATIO_POSITIVE_MIN:
        return ResponseLabel.POSITIVE.value, ratio, censored
    return ResponseLabel.TESTED_NEGATIVE.value, ratio, censored


def _prov(
    entity: str,
    entity_id: str,
    *,
    row: int | str,
    fragment: str,
    origin: str = ValueOrigin.SOURCE_REPORTED.value,
) -> dict:
    provenance_id = "prov:" + sha256(f"{entity}|{entity_id}".encode()).hexdigest()[:20]
    return {
        "provenance_id": provenance_id,
        "entity_type": entity,
        "entity_id": entity_id,
        "field_name": "*",
        "source_document": f"{PMCID}:{TABLE_FILE}",
        "table": PEPTIDE_SHEET,
        "row": str(row),
        "source_fragment": fragment,
        "extraction_method": "deterministic_xlsx_adapter",
        "extraction_confidence": 1.0,
        "value_origin": origin,
        "review_status": ReviewStatus.ACCEPTED.value,
    }


class ZhaoDCAdapter:
    """Study-specific adapter; keeps source column names out of the canonical schema."""

    declaration = AdapterDeclaration(
        "Zhao 2026 Fukuoka DC personalized neoantigen peptide vaccine",
        f"{PMCID}/FrontImmunol-2026",
        "zhao_dc_event_b",
        "1.0.0",
        ("studies", "patients", "vaccines", "candidates", "assays", "provenance"),
        (BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,),
        (
            "Vaccinated-subset discrimination only: peptides were pre-selected for vaccination by "
            "an HLA-binding workflow, so this asset cannot support full-denominator top-K claims.",
            "One post-vaccine ELISpot fold-increase per peptide; no pre-vaccine per-peptide "
            "magnitude is published, so the pre/post baseline is not stored separately.",
            "HLA per peptide is the predicted restricting allele, not the patient's experimentally "
            "typed genotype; patient hla_alleles is left empty so the HLA-subset check is skipped "
            "rather than fabricated.",
            "No genomic coordinates, gene symbol, or transcript are published; only peptide-level "
            "mutant/wild-type sequences and the 1-based mutated position within the peptide.",
            "Single clinic / single platform (Fukuoka General Cancer Clinic peptide-pulsed DC "
            "vaccine); likely platform-specific bias.",
        ),
        (
            "A peptide is POSITIVE iff the source 'ELSPOT ratio' >= 2.0 (the paper's rule); every "
            "other tested peptide is TESTED_NEGATIVE because all 2,317 were in the post-vaccine "
            "ELISpot denominator. The right-censored '>=5.0' cell counts as POSITIVE (>= 2.0).",
            "All peptides are EVENT_B: each is a personalized vaccine component measured only after "
            "vaccination by post-vaccine IFN-gamma ELISpot fold-increase.",
            "Only SNV-derived short (8-11mer) class I peptides are in scope; the paper excluded 12 "
            "short indel and one short multi-nucleotide peptide from this 2,317 set.",
        ),
        ("gene", "protein_change", "genomic_variant", "transcript", "patient_hla_genotype"),
        canonical_study_id=STUDY_ID,
        cohort_id="fukuoka_general_cancer_clinic_dc_2018_2023",
        source_files=(TABLE_FILE, DATASHEET_FILE),
        supported_timepoints=("POST_VACCINE",),
        positivity_rules=(POSITIVITY_RULE,),
        baseline_semantics=(
            "Post-vaccine over pre-vaccine IFN-gamma ELISpot fold-increase; per-peptide baseline "
            "magnitudes are not published (only the ratio)."
        ),
        vaccine_component_structure="Personalized short SNV neoantigen peptides (HLA-binding pre-selected).",
        assay_target_structure="One post-vaccine IFN-gamma ELISpot fold-increase per vaccine peptide.",
        candidate_identity_completeness="PATIENT_PEPTIDE_AND_RESTRICTING_HLA_RESOLVED",
        unresolved_ambiguities=(
            "Patient overlap with the blocked fukuoka_dc (Morisaki 2021) cohort is plausible; "
            "flagged by the overlap audit, not silently assumed absent.",
        ),
    )

    def __init__(self, raw_dir: str | Path):
        self.raw_dir = Path(raw_dir)
        self.review_issues: list[ReviewIssue] = []

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        paths = stage_zhao_supplements(self.raw_dir)
        return {"peptides": read_peptide_sheet(paths[TABLE_FILE])}

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        self.review_issues = []
        peptides = pd.DataFrame(extracted["peptides"]).reset_index(drop=True)

        peptide_counts = (
            peptides["Patient ID"].map(_clean_int_text).replace("", pd.NA).dropna().value_counts()
        ).to_dict()

        studies: dict[str, dict] = {}
        patients: dict[str, dict] = {}
        vaccines: dict[str, dict] = {}
        candidates: list[dict] = []
        assays: list[dict] = []
        funnel_links: list[dict] = []
        provenance: list[dict] = []
        seen_candidates: set[str] = set()

        study_prov = _prov(
            "studies",
            STUDY_ID,
            row=0,
            fragment=(
                f"{DOI}; Fukuoka General Cancer Clinic personalized neoantigen peptide-pulsed DC "
                "vaccine; post-vaccine IFN-gamma ELISpot; 352 patients / 2,317 short SNV peptides"
            ),
        )
        provenance.append(study_prov)
        studies[STUDY_ID] = {
            "study_id": STUDY_ID,
            "title": (
                "Profiling immunogenic neoantigen peptides elicited by personalized neoantigen "
                "vaccine in cancer patients"
            ),
            "publication_ids": f"DOI:{DOI}; PMID:{PMID}; {PMCID}",
            "trial_id": pd.NA,
            "cancer_type": "mixed advanced solid tumors",
            "vaccine_platform": "personalized neoantigen peptide-pulsed dendritic-cell vaccine",
            "adjuvant": pd.NA,
            "vaccination_schedule": pd.NA,
            "source_urls": f"{ARTICLE_URL} ; {SUPPL_URL}",
            "source_paths": str((self.raw_dir / "manual" / TABLE_FILE)),
            "source_checksums": EXPECTED_SHA256[TABLE_FILE],
            "source_manifest_id": manifest.manifest_id,
            "provenance_id": study_prov["provenance_id"],
        }

        for source_row, row in peptides.iterrows():
            source_patient = _clean_int_text(row.get("Patient ID"))
            if not source_patient:
                continue
            spreadsheet_row = int(source_row) + 3  # header on sheet rows 1-2; data from row 3
            patient_id = f"{STUDY_ID}:{source_patient}"
            vaccine_id = f"{STUDY_ID}:{source_patient}:pcv"
            cancer_type = _txt(row.get("Cancer type")) or "unspecified"

            if patient_id not in patients:
                patient_prov = _prov(
                    "patients",
                    patient_id,
                    row=spreadsheet_row,
                    fragment=f"source_patient={source_patient}; cancer_type={cancer_type}",
                )
                provenance.append(patient_prov)
                patients[patient_id] = {
                    "patient_id": patient_id,
                    "source_patient_id": source_patient,
                    "study_id": STUDY_ID,
                    "cancer_type": cancer_type,
                    # Full HLA genotype is not published; left empty so no fabricated genotype
                    # gates the per-candidate restricting allele.
                    "hla_alleles": pd.NA,
                    "provenance_id": patient_prov["provenance_id"],
                }
                vaccine_prov = _prov(
                    "vaccines",
                    vaccine_id,
                    row=spreadsheet_row,
                    fragment=(
                        f"personalized peptide-pulsed DC vaccine; assayed_peptides="
                        f"{peptide_counts.get(source_patient, 0)}"
                    ),
                )
                provenance.append(vaccine_prov)
                vaccines[vaccine_id] = {
                    "vaccine_id": vaccine_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "vaccine_platform": "personalized neoantigen peptide-pulsed dendritic-cell vaccine",
                    "mhc_class_intent": MHCClass.CLASS_I.value,
                    "candidate_count": int(peptide_counts.get(source_patient, 0)),
                    "provenance_id": vaccine_prov["provenance_id"],
                }

            peptide = _txt(row.get("peptide(mut)")).upper()
            wildtype = _txt(row.get("peptide(wt)")).upper()
            hla = _txt(row.get("HLA type")).upper()
            mutation_type = _txt(row.get("mutation type")) or "SNV"
            position = _clean_int_text(row.get("position"))
            peptide_source_id = _clean_int_text(row.get("Peptide ID")) or f"row{spreadsheet_row}"
            try:
                length = int(float(row.get("Length")))
            except (TypeError, ValueError):
                length = len(peptide)

            identity = {
                "study_id": STUDY_ID,
                "patient_id": patient_id,
                "sample_id": "",
                "timepoint": "",
                "genomic_variant": "",
                "transcript": "",
                "mutant_peptide": peptide,
                "hla_alleles": hla,
            }
            candidate_id = stable_candidate_id(identity)

            # Repeated assays of an identical peptide/HLA within a patient share one candidate
            # identity but are kept as distinct assay rows (never collapsed). Only the candidate
            # row is de-duplicated.
            if candidate_id not in seen_candidates:
                seen_candidates.add(candidate_id)
                candidate_prov = _prov(
                    "candidates",
                    candidate_id,
                    row=spreadsheet_row,
                    fragment=(
                        f"peptide_id={peptide_source_id}; mut={peptide}; wt={wildtype}; "
                        f"hla={hla}; mutation_type={mutation_type}; mut_position={position}"
                    ),
                )
                provenance.append(candidate_prov)
                candidates.append(
                    {
                        "candidate_id": candidate_id,
                        "patient_id": patient_id,
                        "study_id": STUDY_ID,
                        "sample_id": pd.NA,
                        "timepoint": pd.NA,
                        "genomic_variant": pd.NA,
                        "gene": pd.NA,
                        "transcript": pd.NA,
                        "protein_change": pd.NA,
                        "mutant_peptide": peptide,
                        "wildtype_peptide": wildtype or pd.NA,
                        "peptide_length": length,
                        "hla_alleles": hla or pd.NA,
                        "mhc_class": MHCClass.CLASS_I.value,
                        "candidate_source": (
                            "vaccine neoantigen peptide (HLA-binding pre-selected); post-vaccine "
                            "ELISpot"
                        ),
                        "vaccine_inclusion": VaccineInclusion.INCLUDED.value,
                        "vaccine_inclusion_origin": ValueOrigin.SOURCE_REPORTED.value,
                        "mutant_wildtype_verified": bool(
                            peptide and wildtype and peptide != wildtype
                        ),
                        "provenance_id": candidate_prov["provenance_id"],
                    }
                )
                funnel_prov = _prov(
                    "candidate_funnel_links",
                    candidate_id,
                    row=spreadsheet_row,
                    fragment="vaccine-included and post-vaccine ELISpot-assayed; upstream funnel not reported",
                )
                provenance.append(funnel_prov)
                funnel_links.append(
                    {
                        "funnel_link_id": "funnel:" + sha256(candidate_id.encode()).hexdigest()[:20],
                        "candidate_id": candidate_id,
                        "patient_id": patient_id,
                        "study_id": STUDY_ID,
                        "mutation_called": _REACHED,
                        "transcript_represented": _NOT_ASSESSED,
                        "peptide_generated": _REACHED,
                        "survives_gating": _NOT_ASSESSED,
                        "hla_included": _REACHED,
                        "presentation_candidate": _NOT_ASSESSED,
                        "ranking_stage": _NOT_ASSESSED,
                        "top_k": _NOT_ASSESSED,
                        "recognition_scored": _REACHED,
                        "vaccine_inclusion": _REACHED,
                        "functional_assay": _REACHED,
                        "provenance_id": funnel_prov["provenance_id"],
                    }
                )

            label, ratio, censored = immunogenic_label(row.get("ELSPOT ratio"))
            assay_id = "assay:" + sha256(
                f"zhao|{patient_id}|{peptide_source_id}".encode()
            ).hexdigest()[:20]
            if label == ResponseLabel.UNTESTED.value:
                review_status = ReviewStatus.NEEDS_REVIEW.value
                self.review_issues.append(
                    ReviewIssue.create(
                        "assays",
                        assay_id,
                        "UNSCORABLE_ASSAY",
                        f"peptide_id={peptide_source_id}: ELSpot ratio "
                        f"{row.get('ELSPOT ratio')!r} could not be scored; held for review, not "
                        "assumed negative",
                    )
                )
            else:
                review_status = ReviewStatus.ACCEPTED.value

            qualitative = "POSITIVE" if label == ResponseLabel.POSITIVE.value else (
                "NEGATIVE" if label == ResponseLabel.TESTED_NEGATIVE.value else "UNSCORED"
            )
            ratio_display = ">=5.0(censored)" if censored else (
                "N/A" if ratio is None else f"{ratio:g}"
            )
            assay_prov = _prov(
                "assays",
                assay_id,
                row=spreadsheet_row,
                fragment=(
                    f"peptide_id={peptide_source_id}; ELSpot_ratio={row.get('ELSPOT ratio')!r}; "
                    f"rule=[ratio>={RATIO_POSITIVE_MIN}] -> {label}"
                ),
                origin=ValueOrigin.DETERMINISTICALLY_DERIVED.value,
            )
            provenance.append(assay_prov)
            assays.append(
                {
                    "assay_id": assay_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "candidate_id": candidate_id,
                    "vaccine_id": vaccine_id,
                    "assay_type": AssayType.ELISPOT.value,
                    "sample_type": "PBMC (post-vaccination)",
                    "timepoint": "post_vaccine",
                    "relative_to_vaccine": "POST_VACCINE",
                    "stimulation_protocol": (
                        "IFN-gamma ELISpot (Mabtech kit); post-vaccine over pre-vaccine "
                        "fold-increase in spot counts"
                    ),
                    "positivity_threshold": POSITIVITY_RULE,
                    "quantitative_result": ratio,
                    "result_units": (
                        "fold increase (post/pre IFN-gamma ELISpot spot-count ratio; "
                        ">=5.0 right-censored)"
                    ),
                    "qualitative_result": qualitative,
                    "source_interpretation": (
                        "immunogenic" if label == ResponseLabel.POSITIVE.value else (
                            "assayed, non-immunogenic"
                            if label == ResponseLabel.TESTED_NEGATIVE.value
                            else "unscorable"
                        )
                    ),
                    "event_type": BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                    "response_label": label,
                    "explicit_assay_inclusion": True,
                    "review_status": review_status,
                    "provenance_id": assay_prov["provenance_id"],
                }
            )

        return EventBCorpus(
            studies=SCHEMAS["studies"].normalize(pd.DataFrame(list(studies.values()))),
            patients=SCHEMAS["patients"].normalize(pd.DataFrame(list(patients.values()))),
            vaccines=SCHEMAS["vaccines"].normalize(pd.DataFrame(list(vaccines.values()))),
            candidates=SCHEMAS["candidates"].normalize(pd.DataFrame(candidates)),
            assays=SCHEMAS["assays"].normalize(pd.DataFrame(assays)),
            candidate_funnel_links=SCHEMAS["candidate_funnel_links"].normalize(
                pd.DataFrame(funnel_links)
            ),
            provenance=SCHEMAS["provenance"].normalize(pd.DataFrame(provenance)),
        )
