"""Braun 2025 clear-cell RCC NeoVax (NCT02950766) Event-B adapter.

This ingests peptide-resolved, de-novo vaccine-induced immunogenicity from an
open-access personalized neoantigen vaccine trial. The immunogenicity call is
recomputed from the raw IFN-gamma ELISpot replicates using the paper's own
stated positivity rule; no counts are hard-coded, no negatives are inferred
from omission, and no pooled response is decomposed to the peptide level.

Source (CC BY-NC-ND, open access):
    Braun DA, Moranzoni G, Chea V, et al. A neoantigen vaccine generates
    antitumour immunity in renal cell carcinoma. Nature 2025;639:474-482.
    doi:10.1038/s41586-024-08507-5 ; PMC11903305 ; ClinicalTrials NCT02950766.

Positivity rule (Braun 2025, Methods / Extended Data Fig. 3 legend, verbatim):
    "P < 0.05 by two-sided t-test and mean spot count at least three-fold
     higher than DMSO [no-stim] control."

Event-B justification (Braun 2025, Results, verbatim):
    "For all peptides, no pre-existing immune responses were detected"
plus ex-vivo week-0 (pre-vaccine) pool baselines that are background
(background-subtracted magnitudes well below the ELISpot positivity floor).
"""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
import shutil
from statistics import mean

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


STUDY_ID = "braun_rcc_2025"
NCT = "NCT02950766"
DOI = "10.1038/s41586-024-08507-5"
PMCID = "PMC11903305"

# Nature MOESM supplementary files. The frozen sha256 pins content identity so a
# silently changed remote file can never be ingested without detection.
INVITRO_FILE = "41586_2024_8507_MOESM4_ESM.xlsx"  # 'In Vitro' per-peptide; 'Ex Vivo' pool baseline
METADATA_FILE = "41586_2024_8507_MOESM9_ESM.xlsx"  # sheet '1c' patient cohort/stage
SUMMARY_FILE = "41586_2024_8507_MOESM10_ESM.xlsx"  # sheet '2e' driver/passenger reconciliation
EXPECTED_SHA256 = {
    INVITRO_FILE: "c113c42b0773049fe7e3f6b983485d15cd00cb847c6fd5de532cea4c9715d0c1",
    METADATA_FILE: "f67adad4fa4dc61ffcd50ee60c2fd2a63ea8d36b57a844a036e59953f0ff95c8",
    SUMMARY_FILE: "3030a7dfa494a0e7068ae22789e64857d9d3113fadd248bfe3204bcc8e3ce8a7",
}
SUPPL_URL = f"https://www.ebi.ac.uk/europepmc/webservices/rest/{PMCID}/supplementaryFiles"

# Paper's stated positivity rule constants.
P_MAX = 0.05
FOLD_MIN = 3.0
STIM_COLS = [f"InVitro_PeptideStim_Replicate0{i}" for i in (1, 2, 3)]
NOSTIM_COLS = [f"InVitro_NoStim_Replicate0{i}" for i in (1, 2, 3)]
# ELISpot positivity floor (net background-subtracted SFC per 1e6 PBMC) above which a
# pre-vaccine (week-0) pool response would count as pre-existing (Event-A) and block a
# clean de-novo Event-B assignment. Braun's week-0 pool magnitudes sit well below this.
PRE_VACCINE_SFC_FLOOR = 50.0

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


def _manual_placement_message(extracted: Path) -> str:
    lines = [
        "Could not fetch Braun RCC 2025 supplements automatically (no network?).",
        f"Manually download the open-access supplement bundle from:\n    {SUPPL_URL}",
        f"then place these files into:\n    {extracted}",
    ]
    lines += [f"    {name}  sha256={digest}" for name, digest in sorted(EXPECTED_SHA256.items())]
    lines.append("No records were fabricated; ingestion refuses to proceed without the source.")
    return "\n".join(lines)


def stage_braun_supplements(raw_dir: str | Path) -> dict[str, Path]:
    """Return verified paths to the three source files, fetching on cache miss.

    Local cache first; download from Europe PMC on miss; actionable failure when
    offline. Every returned file is checksum-verified against ``EXPECTED_SHA256``.
    """
    extracted = Path(raw_dir) / "extracted"
    paths = {name: extracted / name for name in EXPECTED_SHA256}
    if _checksums_ok(paths):
        return paths
    extracted.mkdir(parents=True, exist_ok=True)
    zip_path = Path(raw_dir) / "suppl.zip"
    try:
        import urllib.request

        request = urllib.request.Request(SUPPL_URL, headers={"User-Agent": "epicurus-neo/1.0"})
        with urllib.request.urlopen(request, timeout=180) as response, zip_path.open("wb") as fh:
            shutil.copyfileobj(response, fh)
    except Exception as error:  # noqa: BLE001 - network failure surfaces as actionable guidance
        raise RuntimeError(_manual_placement_message(extracted)) from error
    from zipfile import ZipFile

    with ZipFile(zip_path) as archive:
        for name in EXPECTED_SHA256:
            archive.extract(name, extracted)
    if not _checksums_ok(paths):
        observed = {name: sha256_file(path) for name, path in paths.items()}
        raise RuntimeError(
            "Braun supplement checksums did not match after download; refusing to ingest. "
            f"Observed: {observed}"
        )
    return paths


def braun_source_paths(raw_dir: str | Path) -> list[Path]:
    """Ordered list of the source files used to build the source manifest."""
    return [stage_braun_supplements(raw_dir)[name] for name in sorted(EXPECTED_SHA256)]


def _read_sheet(path: Path, sheet: str) -> pd.DataFrame:
    frame = pd.read_excel(path, sheet_name=sheet, engine="openpyxl")
    frame.columns = [str(column).strip() for column in frame.columns]
    return frame


def _num(value: object) -> float | None:
    try:
        result = float(value)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None
    return None if pd.isna(result) else result


def _txt(value: object) -> str:
    """Blank-safe text: pandas reads empty cells as NaN, which must not become 'nan'."""
    if value is None or (isinstance(value, float) and pd.isna(value)) or value is pd.NA:
        return ""
    return str(value).strip()


def _clean_int(value: object) -> str:
    number = _num(value)
    if number is not None and float(number).is_integer():
        return str(int(number))
    return _txt(value)


def _means(row: pd.Series) -> tuple[float | None, float | None]:
    stim = [v for v in (_num(row.get(c)) for c in STIM_COLS) if v is not None]
    nostim = [v for v in (_num(row.get(c)) for c in NOSTIM_COLS) if v is not None]
    return (mean(stim) if stim else None, mean(nostim) if nostim else None)


def immunogenic_call(row: pd.Series) -> bool | None:
    """Apply the paper's stated rule. Returns None when the row cannot be scored."""
    pvalue = _num(row.get("Ttest_pvalue_InVitroStim"))
    stim_mean, nostim_mean = _means(row)
    if pvalue is None or stim_mean is None or nostim_mean is None:
        return None
    return (pvalue < P_MAX) and (stim_mean >= FOLD_MIN * nostim_mean)


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
        "extraction_method": "deterministic_xlsx_adapter",
        "extraction_confidence": 1.0,
        "value_origin": origin,
        "review_status": ReviewStatus.ACCEPTED.value,
    }


class BraunRCCAdapter:
    """Study-specific adapter; keeps source column names out of the canonical schema."""

    declaration = AdapterDeclaration(
        "Braun 2025 RCC NeoVax personalized neoantigen vaccine",
        f"{PMCID}/Nature-2025",
        "braun_rcc_event_b",
        "1.0.0",
        ("studies", "patients", "vaccines", "candidates", "assays", "provenance"),
        (BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,),
        (
            "Per-peptide HLA restriction is a prediction (best short epitope), not experimentally "
            "deconvolved; it is not stored as the assay restriction.",
            "Immunogenicity is measured after in vitro peptide stimulation and pool deconvolution "
            "of the post-vaccine (week-16) leukapheresis: one post-vaccine readout per peptide.",
            "Patient HLA genotypes are absent from the public supplement; patient hla_alleles left "
            "empty (the HLA-subset check is thereby skipped rather than fabricated).",
            "CD4 vs CD8 attribution is pool-level (Extended Data Fig. 4), not per-peptide; the long "
            "vaccine peptide's mhc_class is left UNKNOWN.",
            "Absolute vaccination/sample calendar dates are de-identified (relative months only).",
        ),
        (
            "A peptide is POSITIVE iff P<0.05 two-sided t-test AND mean IFN-gamma spot count is at "
            "least 3x the no-stim DMSO control (the paper's rule); every other assayed peptide is "
            "TESTED_NEGATIVE because it was explicitly in the deconvolution denominator.",
            "All assayed vaccine peptides are EVENT_B de-novo: the paper states 'no pre-existing "
            "immune responses were detected' and ex-vivo week-0 baselines are background.",
            "The long synthetic vaccine peptide is the tested entity (mutant_peptide); mhc_class "
            "UNKNOWN avoids imposing a spurious class-I length constraint on a >=20mer.",
        ),
        ("transcript", "wildtype_peptide", "patient_hla_genotype", "per_peptide_clinical_outcome"),
    )

    def __init__(self, raw_dir: str | Path):
        self.raw_dir = Path(raw_dir)
        self.review_issues: list[ReviewIssue] = []

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        paths = stage_braun_supplements(self.raw_dir)
        return {
            "invitro": _read_sheet(paths[INVITRO_FILE], "In Vitro"),
            "exvivo": _read_sheet(paths[INVITRO_FILE], "Ex Vivo"),
            "meta": _read_sheet(paths[METADATA_FILE], "1c"),
        }

    def _baseline_max(self, exvivo: pd.DataFrame) -> float:
        week0 = [c for c in exvivo.columns if c.startswith("Week 0")]
        if not week0:
            return 0.0
        values = [abs(v) for v in (_num(x) for x in exvivo[week0[0]]) if v is not None]
        return max(values) if values else 0.0

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        self.review_issues = []
        invitro = pd.DataFrame(extracted["invitro"])
        exvivo = pd.DataFrame(extracted["exvivo"])
        meta = pd.DataFrame(extracted["meta"])
        document = f"{PMCID}:{INVITRO_FILE}"

        baseline_max = self._baseline_max(exvivo)
        # Corroborates the paper's per-peptide "no pre-existing responses" statement. If a week-0
        # pool ever crossed the positivity floor we would have to demote those peptides; it does not.
        de_novo_supported = baseline_max < PRE_VACCINE_SFC_FLOOR

        patient_meta: dict[str, dict] = {}
        for _, row in meta.iterrows():
            if pd.isna(row.get("ID")):
                continue
            source_id = _txt(row["ID"]).split("-")[-1]
            patient_meta[source_id] = {
                "cohort": _txt(row.get("Cohort")),
                "stage": _txt(row.get("Stage")),
            }

        peptide_counts = invitro["Patient_ID"].astype(str).str.strip().value_counts().to_dict()

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
            document=document,
            table="-",
            row=0,
            fragment=(
                f"{NCT}; phase I high-risk resected ccRCC; personalized SLP + poly-ICLC; "
                f"de-novo Event-B (paper: no pre-existing responses; ex-vivo wk0 |max|="
                f"{baseline_max:.1f} SFC < {PRE_VACCINE_SFC_FLOOR:.0f} floor)"
            ),
        )
        provenance.append(study_prov)
        studies[STUDY_ID] = {
            "study_id": STUDY_ID,
            "title": "A neoantigen vaccine generates antitumour immunity in renal cell carcinoma",
            "publication_ids": f"DOI:{DOI}; {PMCID}",
            "trial_id": NCT,
            "cancer_type": "clear cell renal cell carcinoma",
            "vaccine_platform": "synthetic long peptide (SLP) + poly-ICLC",
            "adjuvant": "poly-ICLC (Hiltonol, Oncovir)",
            "vaccination_schedule": (
                "prime series then 2 boosts; 4 pools x up to 5 SLPs (300 ug/peptide) admixed with "
                "0.5 mg poly-ICLC per pool"
            ),
            "source_urls": f"https://pmc.ncbi.nlm.nih.gov/articles/{PMCID}/",
            "source_paths": str((self.raw_dir / "extracted" / INVITRO_FILE)),
            "source_checksums": EXPECTED_SHA256[INVITRO_FILE],
            "source_manifest_id": manifest.manifest_id,
            "provenance_id": study_prov["provenance_id"],
        }

        for source_row, row in invitro.reset_index(drop=True).iterrows():
            source_patient = str(row["Patient_ID"]).strip()
            patient_id = f"{STUDY_ID}:{source_patient}"
            vaccine_id = f"{STUDY_ID}:{source_patient}:pcv"
            info = patient_meta.get(source_patient, {"cohort": "", "stage": ""})
            cohort = info["cohort"]
            concurrent = "ipilimumab" if "ipilimumab" in cohort.lower() else "none"

            if patient_id not in patients:
                patient_prov = _prov(
                    "patients",
                    patient_id,
                    document=f"{PMCID}:{METADATA_FILE}",
                    table="1c",
                    row=0,
                    fragment=f"source_patient={source_patient}; cohort={cohort}; stage={info['stage']}",
                )
                provenance.append(patient_prov)
                patients[patient_id] = {
                    "patient_id": patient_id,
                    "source_patient_id": source_patient,
                    "study_id": STUDY_ID,
                    "cancer_type": "clear cell renal cell carcinoma",
                    "disease_stage": info["stage"] or pd.NA,
                    "treatment_context": cohort or pd.NA,
                    "provenance_id": patient_prov["provenance_id"],
                }
                vaccine_prov = _prov(
                    "vaccines",
                    vaccine_id,
                    document=document,
                    table="-",
                    row=0,
                    fragment=f"personalized SLP; concurrent={concurrent}; assayed_peptides="
                    f"{peptide_counts.get(source_patient, 0)}",
                )
                provenance.append(vaccine_prov)
                vaccines[vaccine_id] = {
                    "vaccine_id": vaccine_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "vaccine_platform": "synthetic long peptide (SLP) + poly-ICLC",
                    "formulation": (
                        "4 pools x up to 5 SLPs (300 ug/peptide) admixed with 0.5 mg poly-ICLC per "
                        "pool"
                    ),
                    "mhc_class_intent": MHCClass.BOTH.value,
                    "candidate_count": int(peptide_counts.get(source_patient, 0)),
                    "concurrent_therapy": concurrent,
                    "provenance_id": vaccine_prov["provenance_id"],
                }

            peptide = _txt(row["Vaccine_Peptide"]).upper()
            gene = _txt(row.get("Hugo_Symbol"))
            change_full = _txt(row.get("Gene_and_Protein_Change"))
            protein_change = change_full.split("|", 1)[1] if "|" in change_full else change_full
            variant_type = _txt(row.get("Variant_Type"))
            chromosome = _clean_int(row.get("Chromosome"))
            position = _clean_int(row.get("Start_position"))
            genomic_variant = f"chr{chromosome}:{position}:{variant_type}".strip(":")
            peptide_id = _txt(row.get("Peptide_ID")) or f"row{source_row}"
            mutation_type = _txt(row.get("Mutation_type"))

            identity = {
                "study_id": STUDY_ID,
                "patient_id": patient_id,
                "sample_id": "",
                "timepoint": "",
                "genomic_variant": genomic_variant,
                "transcript": "",
                "mutant_peptide": peptide,
                "hla_alleles": "",
            }
            candidate_id = stable_candidate_id(identity)
            if candidate_id in seen_candidates:
                self.review_issues.append(
                    ReviewIssue.create(
                        "candidates",
                        candidate_id,
                        "DUPLICATE_IDENTITY",
                        f"peptide {peptide_id} collides with an existing candidate identity",
                    )
                )
                continue
            seen_candidates.add(candidate_id)

            call = immunogenic_call(row)
            stim_mean, nostim_mean = _means(row)
            pvalue = _num(row.get("Ttest_pvalue_InVitroStim"))
            stim_raw = [row.get(c) for c in STIM_COLS]
            nostim_raw = [row.get(c) for c in NOSTIM_COLS]

            candidate_prov = _prov(
                "candidates",
                candidate_id,
                document=document,
                table="In Vitro",
                row=int(source_row) + 2,
                fragment=f"{gene}|{protein_change}|{variant_type}|peptide_id={peptide_id}",
            )
            provenance.append(candidate_prov)
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "sample_id": pd.NA,
                    "timepoint": pd.NA,
                    "genomic_variant": genomic_variant,
                    "gene": gene or pd.NA,
                    "transcript": pd.NA,
                    "protein_change": protein_change or pd.NA,
                    "mutant_peptide": peptide,
                    "wildtype_peptide": pd.NA,
                    "peptide_length": len(peptide),
                    "hla_alleles": pd.NA,
                    "mhc_class": MHCClass.UNKNOWN.value,
                    "candidate_source": "vaccine SLP, in-vitro ELISpot deconvolution",
                    "vaccine_inclusion": VaccineInclusion.INCLUDED.value,
                    "vaccine_inclusion_origin": ValueOrigin.SOURCE_REPORTED.value,
                    "mutant_wildtype_verified": True,
                    "provenance_id": candidate_prov["provenance_id"],
                }
            )

            # Only stages with direct source evidence are 'reached'; the Braun paper does not report
            # the upstream pVACtools-style reachability funnel (gating/ranking/top-k), so those stay
            # 'not_assessed' rather than being inferred as reached from vaccine inclusion.
            funnel_prov = _prov(
                "candidate_funnel_links",
                candidate_id,
                document=document,
                table="In Vitro",
                row=int(source_row) + 2,
                fragment="vaccine-included and ELISpot-assayed; upstream reachability not reported",
            )
            provenance.append(funnel_prov)
            funnel_links.append(
                {
                    "funnel_link_id": "funnel:"
                    + sha256(candidate_id.encode()).hexdigest()[:20],
                    "candidate_id": candidate_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "mutation_called": _REACHED,
                    "transcript_represented": _NOT_ASSESSED,
                    "peptide_generated": _REACHED,
                    "survives_gating": _NOT_ASSESSED,
                    "hla_included": _NOT_ASSESSED,
                    "presentation_candidate": _NOT_ASSESSED,
                    "ranking_stage": _NOT_ASSESSED,
                    "top_k": _NOT_ASSESSED,
                    "recognition_scored": _REACHED,
                    "vaccine_inclusion": _REACHED,
                    "functional_assay": _REACHED,
                    "provenance_id": funnel_prov["provenance_id"],
                }
            )

            assay_id = "assay:" + sha256(f"braun|{patient_id}|{peptide_id}".encode()).hexdigest()[
                :20
            ]
            if call is None:
                label = ResponseLabel.UNTESTED.value
                review_status = ReviewStatus.NEEDS_REVIEW.value
                self.review_issues.append(
                    ReviewIssue.create(
                        "assays",
                        assay_id,
                        "UNSCORABLE_ASSAY",
                        f"peptide {peptide_id}: missing replicates or p-value; cannot apply rule",
                    )
                )
            elif not mutation_type:
                # Immunogenic-or-not by the rule, but the source left it unclassified (blank
                # Mutation_type) and excluded it from the driver/passenger summary (Fig. 2e). Do
                # not silently inflate the accepted positive count: hold it for review.
                label = ResponseLabel.POSITIVE.value if call else ResponseLabel.TESTED_NEGATIVE.value
                review_status = ReviewStatus.NEEDS_REVIEW.value
                self.review_issues.append(
                    ReviewIssue.create(
                        "assays",
                        assay_id,
                        "UNCLASSIFIED_MUTATION",
                        f"{gene}|{protein_change} scored {'immunogenic' if call else 'negative'} by "
                        "the paper rule but has a blank Mutation_type and is excluded from the "
                        "paper's Fig-2e driver/passenger summary; held for review",
                    )
                )
            else:
                label = ResponseLabel.POSITIVE.value if call else ResponseLabel.TESTED_NEGATIVE.value
                review_status = ReviewStatus.ACCEPTED.value

            event_type = (
                BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
                if de_novo_supported
                else BiologicalEvent.UNKNOWN_EVENT.value
            )
            # An unscorable row (UNTESTED) asserts no net magnitude; raw replicates stay in
            # provenance for the reviewer, but we do not report a result we could not score.
            net = (
                None
                if call is None or stim_mean is None or nostim_mean is None
                else round(stim_mean - nostim_mean, 2)
            )
            qualitative = "POSITIVE" if call is True else ("NEGATIVE" if call is False else "UNSCORED")
            assay_prov = _prov(
                "assays",
                assay_id,
                document=document,
                table="In Vitro",
                row=int(source_row) + 2,
                fragment=(
                    f"peptide_id={peptide_id}; stim={stim_raw}; nostim={nostim_raw}; p={pvalue}; "
                    f"rule=[p<{P_MAX} AND mean_stim>={FOLD_MIN:.0f}x mean_nostim] -> {label}; "
                    f"ex_vivo_wk0_pool_baseline_max={baseline_max:.1f} SFC"
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
                    "sample_type": "PBMC (post-vaccine leukapheresis), in-vitro peptide-stimulated",
                    "timepoint": "week16_post_vaccine_in_vitro_deconvolution",
                    "relative_to_vaccine": "POST_VACCINE",
                    "stimulation_protocol": (
                        "in vitro peptide stimulation; IFN-gamma ELISpot triplicate vs no-stim DMSO"
                    ),
                    "replicate_information": json.dumps(
                        {"peptide_stim": stim_raw, "no_stim": nostim_raw, "ttest_p": pvalue},
                        default=str,
                    ),
                    "positivity_threshold": (
                        "P<0.05 two-sided t-test AND mean spot count >=3x no-stim DMSO control "
                        "(Braun 2025 Methods / Extended Data Fig. 3)"
                    ),
                    "quantitative_result": net,
                    "result_units": "net SFC per 1e6 PBMC (mean triplicate minus no-stim)",
                    "qualitative_result": qualitative,
                    "source_interpretation": (
                        "immunogenic" if call else "assayed, non-immunogenic"
                    ),
                    "event_type": event_type,
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
