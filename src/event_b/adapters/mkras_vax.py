"""mKRAS-VAX shared-antigen Event-B adapter (NCT04117087, PDAC cohort)."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import urllib.request

import pandas as pd

from event_b.adapters.base import AdapterDeclaration
from event_b.adapters.common import entity_frame, provenance_record, stable_record_id
from event_b.corpus import EventBCorpus
from event_b.manifest import SourceManifest, manifest_from_paths, sha256_file
from event_b.models import (
    AssayType,
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


STUDY_ID = "mkras_vax_2026"
COHORT_ID = "NCT04117087_resected_pdac"
DOI = "10.1038/s41467-026-68324-4"
NCT = "NCT04117087"
PDF_FILE = "41467_2026_68324_MOESM1_ESM.pdf"
DATA_FILE = "41467_2026_68324_MOESM4_ESM.xlsx"
BASE_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41467-026-68324-4/MediaObjects/"
)
SOURCE_URLS = {name: BASE_URL + name for name in (PDF_FILE, DATA_FILE)}
EXPECTED_SHA256 = {
    PDF_FILE: "8bc2ad373da7dabcd6efb2b25926f2eacc04456333e9474f06a4fb9ebd6ff61b",
    DATA_FILE: "b7ac94e8c1427b63eb35b503c60ea12bea863165402b7a11cc2d83d927e09a7a",
}

ANTIGENS = {
    "G12V": ("p.Gly12Val", "YKLVVVGAVGVGKSALTIQLI"),
    "G12A": ("p.Gly12Ala", "YKLVVVGAAGVGKSALTIQLI"),
    "G12R": ("p.Gly12Arg", "YKLVVVGARGVGKSALTIQLI"),
    "G12C": ("p.Gly12Cys", "YKLVVVGACGVGKSALTIQLI"),
    "G12D": ("p.Gly12Asp", "YKLVVVGADGVGKSALTIQLI"),
    "G13D": ("p.Gly13Asp", "YKLVVVGAGDVGKSALTIQLI"),
}
WILDTYPE_SEQUENCE = "YKLVVVGAGGVGKSALTIQLI"
TUMOR_MUTATIONS = {
    "J1994_1": "G12R",
    "J1994_2": "G12V",
    "J1994_3": "G12D",
    "J1994_5": "G12V",
    "J1994_6": "G12V",
    "J1994_7": "G12V",
    "J1994_9": "G12V",
    "J1994_10": "G12D",
    "J1994_12": "G12R",
    "J1994_13": "G12D",
    "J1994_14": "G12V",
    "J1994_18": "G12D",
}


def stage_sources(raw_dir: str | Path) -> list[Path]:
    root = Path(raw_dir)
    root.mkdir(parents=True, exist_ok=True)
    paths = []
    for name in (PDF_FILE, DATA_FILE):
        path = root / name
        if not path.exists() or sha256_file(path) != EXPECTED_SHA256[name]:
            request = urllib.request.Request(
                SOURCE_URLS[name], headers={"User-Agent": "epicurus-neo/1.0"}
            )
            try:
                with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as out:
                    shutil.copyfileobj(response, out)
            except Exception as error:  # noqa: BLE001 - converted to a manual-source contract
                raise RuntimeError(
                    f"Download {SOURCE_URLS[name]} manually to {path}; "
                    f"expected sha256={EXPECTED_SHA256[name]}"
                ) from error
        observed = sha256_file(path)
        if observed != EXPECTED_SHA256[name]:
            raise RuntimeError(f"checksum mismatch for {path}: observed {observed}")
        paths.append(path)
    return paths


def source_manifest(raw_dir: str | Path) -> SourceManifest:
    adapter = MKRASVaxAdapter(raw_dir)
    return manifest_from_paths(
        adapter.declaration.source_name,
        adapter.declaration.source_version,
        adapter.declaration.adapter_name,
        adapter.declaration.adapter_version,
        stage_sources(raw_dir),
    )


_id = stable_record_id
_prov = provenance_record
_frame = entity_frame


class MKRASVaxAdapter:
    declaration = AdapterDeclaration(
        "mKRAS-VAX resected PDAC phase I",
        "Nature-Communications-2026",
        "mkras_vax_event_b",
        "1.0.0",
        (
            "antigens",
            "studies",
            "patients",
            "vaccines",
            "candidates",
            "assays",
            "recognition_evidence",
            "candidate_funnel_links",
            "entity_relationships",
            "provenance",
        ),
        (
            BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value,
            BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
        ),
        (
            "The primary 1/0 component calls are author-reported; raw time-series values are "
            "preserved in the source workbook but not re-thresholded.",
            "Patient HLA genotypes are reported only in a PDF table and are not transcribed by "
            "this adapter; no HLA restriction is inferred for a 21-mer response.",
            "All participants also received nivolumab and ipilimumab.",
        ),
        (
            "Figure 2c is an explicit patient-by-antigen denominator: 1 is POSITIVE and 0 is "
            "TESTED_NEGATIVE for response within 17 weeks.",
            "Baseline individual-antigen assays are explicit and reported as lacking "
            "pre-vaccination responses; they remain Event-A observations.",
            "Tumor KRAS mutation is distinct from response to the other five shared vaccine "
            "components and is recorded only on the matching patient observation.",
        ),
        ("patient_hla_genotype", "experimentally_resolved_hla_restriction"),
        canonical_study_id=STUDY_ID,
        cohort_id=COHORT_ID,
        source_files=(PDF_FILE, DATA_FILE),
        supported_timepoints=("BASELINE", "MAXIMUM_WITHIN_17_WEEKS"),
        positivity_rules=("Figure 2c author call: 1 positive, 0 negative",),
        baseline_semantics="All six antigens tested at baseline; no pre-vaccine response detected.",
        vaccine_component_structure="Pool of six shared mutant-KRAS 21-mer peptides.",
        assay_target_structure="Each 21-mer individually restimulated in ex-vivo IFNg ELISpot.",
        candidate_identity_completeness="PATIENT_AND_SHARED_ANTIGEN_RESOLVED",
    )

    def __init__(self, raw_dir: str | Path):
        self.raw_dir = Path(raw_dir)
        self.review_issues = []

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        paths = {path.name: path for path in stage_sources(self.raw_dir)}
        calls = pd.read_excel(
            paths[DATA_FILE], sheet_name="Figure 2c", header=4, usecols="B:H", engine="openpyxl"
        )
        calls = calls.rename(columns={calls.columns[0]: "patient"}).dropna(subset=["patient"])
        calls["patient"] = calls["patient"].astype(str).str.strip()
        return {"calls": calls}

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        calls = pd.DataFrame(extracted["calls"])
        expected_columns = {"patient", *ANTIGENS}
        if set(calls.columns) != expected_columns:
            raise ValueError(f"unexpected Figure 2c columns: {list(calls.columns)}")
        if set(calls.patient) != set(TUMOR_MUTATIONS):
            raise ValueError("Figure 2c patient denominator changed")

        rows: dict[str, list[dict]] = {entity: [] for entity in SCHEMAS}

        def add(entity: str, row: dict, provenance: dict) -> None:
            row["provenance_id"] = provenance["provenance_id"]
            rows[entity].append(row)
            rows["provenance"].append(provenance)

        study_prov = _prov(
            "studies",
            STUDY_ID,
            document=f"DOI:{DOI}",
            table="main article",
            row="Methods",
            column="study identity",
            fragment="NCT04117087; 12 treated participants; pooled six-peptide mKRAS-VAX",
            method="manual_primary_source_curation",
            origin=ValueOrigin.MANUALLY_CURATED.value,
        )
        add(
            "studies",
            {
                "study_id": STUDY_ID,
                "title": "Mutant KRAS vaccine with dual checkpoint blockade in resected PDAC",
                "publication_ids": json.dumps([f"DOI:{DOI}", "PMID:42336869"]),
                "trial_id": NCT,
                "cancer_type": "pancreatic ductal adenocarcinoma",
                "vaccine_platform": "pooled synthetic long-peptide vaccine",
                "adjuvant": "poly-ICLC",
                "vaccination_schedule": "prime/boost with concurrent nivolumab and ipilimumab",
                "source_urls": json.dumps([SOURCE_URLS[PDF_FILE], SOURCE_URLS[DATA_FILE]]),
                "source_paths": json.dumps([doc.local_path for doc in manifest.documents]),
                "source_checksums": json.dumps(EXPECTED_SHA256, sort_keys=True),
                "source_manifest_id": manifest.manifest_id,
            },
            study_prov,
        )

        antigen_ids = {}
        for mutation, (protein_change, sequence) in ANTIGENS.items():
            antigen_id = f"antigen:mkras:{mutation.lower()}"
            antigen_ids[mutation] = antigen_id
            prov = _prov(
                "antigens",
                antigen_id,
                document=f"DOI:{DOI}",
                table="Methods: Vaccine formulation and administration",
                row=1637,
                column=mutation,
                fragment=f"{mutation}: {sequence}",
                method="manual_primary_source_curation",
                origin=ValueOrigin.MANUALLY_CURATED.value,
            )
            add(
                "antigens",
                {
                    "antigen_id": antigen_id,
                    "study_id": STUDY_ID,
                    "gene": "KRAS",
                    "protein_change": protein_change,
                    "mutant_sequence": sequence,
                    "wildtype_sequence": WILDTYPE_SEQUENCE,
                    "component_type": "IMMUNIZING_LONG_PEPTIDE",
                    "peptide_length": len(sequence),
                    "hla_alleles": json.dumps([]),
                    "hla_evidence_type": "NOT_ASSESSED_FOR_LONG_PEPTIDE",
                },
                prov,
            )

        for source_row, call_row in calls.reset_index(drop=True).iterrows():
            source_patient = str(call_row.patient)
            patient_id = f"{STUDY_ID}:{source_patient}"
            patient_prov = _prov(
                "patients",
                patient_id,
                document=PDF_FILE,
                table="Supplementary Table 5 and 6",
                row=source_patient,
                column="Patient ID and tumor mutation",
                fragment=f"{source_patient}; tumor KRAS {TUMOR_MUTATIONS[source_patient]}",
                method="manual_primary_source_curation",
                origin=ValueOrigin.MANUALLY_CURATED.value,
            )
            add(
                "patients",
                {
                    "patient_id": patient_id,
                    "source_patient_id": source_patient,
                    "study_id": STUDY_ID,
                    "cancer_type": "pancreatic ductal adenocarcinoma",
                    "disease_stage": "resected",
                    "treatment_context": "adjuvant after surgery and standard therapy",
                    "prior_therapies": "perioperative chemotherapy for 11/12 cohort participants",
                    "hla_alleles": json.dumps([]),
                    "tumor_context": f"directly confirmed KRAS {TUMOR_MUTATIONS[source_patient]}",
                },
                patient_prov,
            )
            vaccine_id = f"vaccine:{STUDY_ID}:{source_patient.lower()}"
            vaccine_prov = _prov(
                "vaccines",
                vaccine_id,
                document=f"DOI:{DOI}",
                table="Methods: Vaccine formulation and administration",
                row=1637,
                column="formulation",
                fragment="six 21-mer peptides, 0.3 mg each, admixed with poly-ICLC",
                method="manual_primary_source_curation",
                origin=ValueOrigin.MANUALLY_CURATED.value,
            )
            add(
                "vaccines",
                {
                    "vaccine_id": vaccine_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "vaccine_platform": "pooled synthetic long-peptide vaccine",
                    "formulation": "six mKRAS 21-mers plus poly-ICLC",
                    "dose": "0.3 mg/peptide; 1.8 mg total",
                    "vaccination_dates": json.dumps([]),
                    "relative_schedule": "prime/boost; primary response within 17 weeks",
                    "candidate_count": 6,
                    "mhc_class_intent": MHCClass.BOTH.value,
                    "concurrent_therapy": "nivolumab and ipilimumab",
                },
                vaccine_prov,
            )

            for mutation, (protein_change, sequence) in ANTIGENS.items():
                raw_call = call_row[mutation]
                if raw_call not in (0, 1, 0.0, 1.0):
                    raise ValueError(f"non-binary Figure 2c call for {source_patient}/{mutation}")
                record = {
                    "study_id": STUDY_ID,
                    "patient_id": patient_id,
                    "sample_id": "PBMC",
                    "timepoint": "VACCINE_COMPONENT",
                    "genomic_variant": (
                        f"KRAS {mutation}" if TUMOR_MUTATIONS[source_patient] == mutation else pd.NA
                    ),
                    "gene": "KRAS",
                    "protein_change": protein_change,
                    "mutant_peptide": sequence,
                    "wildtype_peptide": WILDTYPE_SEQUENCE,
                    "peptide_length": len(sequence),
                    "hla_alleles": json.dumps([]),
                    "mhc_class": MHCClass.UNKNOWN.value,
                    "candidate_source": (
                        "shared vaccine component; tumor mutation directly confirmed"
                        if TUMOR_MUTATIONS[source_patient] == mutation
                        else "shared vaccine component; not the patient's confirmed tumor mutation"
                    ),
                    "vaccine_inclusion": VaccineInclusion.INCLUDED.value,
                    "vaccine_inclusion_origin": ValueOrigin.SOURCE_REPORTED.value,
                    "generation_provenance": "fixed six-component mKRAS-VAX formulation",
                    "mutant_wildtype_verified": True,
                }
                candidate_id = stable_candidate_id(record)
                candidate_prov = _prov(
                    "candidates",
                    candidate_id,
                    document=DATA_FILE,
                    table="Figure 2c",
                    row=source_row + 6,
                    column=mutation,
                    fragment=f"{source_patient}; {mutation}; call={int(raw_call)}",
                )
                add("candidates", {"candidate_id": candidate_id, **record}, candidate_prov)

                for relative, event, label, timepoint, value in (
                    (
                        "PRE_VACCINE",
                        BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value,
                        ResponseLabel.TESTED_NEGATIVE.value,
                        "BASELINE",
                        0,
                    ),
                    (
                        "POST_VACCINE",
                        BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                        ResponseLabel.POSITIVE.value
                        if int(raw_call) == 1
                        else ResponseLabel.TESTED_NEGATIVE.value,
                        "MAXIMUM_WITHIN_17_WEEKS",
                        int(raw_call),
                    ),
                ):
                    assay_id = _id("assay", candidate_id, event, timepoint)
                    assay_prov = _prov(
                        "assays",
                        assay_id,
                        document=DATA_FILE if relative == "POST_VACCINE" else PDF_FILE,
                        table="Figure 2c" if relative == "POST_VACCINE" else "Supplementary Table 5",
                        row=source_row + 6 if relative == "POST_VACCINE" else source_patient,
                        column=mutation,
                        fragment=(
                            f"author call={value}; 1=positive, 0=negative"
                            if relative == "POST_VACCINE"
                            else "baseline comparison explicitly assayed; no pre-vaccine response"
                        ),
                    )
                    add(
                        "assays",
                        {
                            "assay_id": assay_id,
                            "patient_id": patient_id,
                            "study_id": STUDY_ID,
                            "candidate_id": candidate_id,
                            "vaccine_id": vaccine_id,
                            "assay_type": AssayType.ELISPOT.value,
                            "sample_type": "PBMC",
                            "timepoint": timepoint,
                            "relative_to_vaccine": relative,
                            "stimulation_protocol": "individual 21-mer; 24 h ex-vivo IFNg ELISpot",
                            "replicate_information": "source statistics use two-way ANOVA/Dunnett",
                            "positivity_threshold": "source author call",
                            "quantitative_result": value,
                            "result_units": "binary author call",
                            "qualitative_result": label,
                            "source_interpretation": "explicit individual-antigen test",
                            "event_type": event,
                            "response_label": label,
                            "explicit_assay_inclusion": True,
                            "review_status": ReviewStatus.ACCEPTED.value,
                        },
                        assay_prov,
                    )

                    relationship_id = _id("rel", assay_id, "TESTS_RESPONSE_TO", antigen_ids[mutation])
                    rel_prov = _prov(
                        "entity_relationships",
                        relationship_id,
                        document=DATA_FILE,
                        table="Figure 2c",
                        row=source_row + 6,
                        column=mutation,
                        fragment="individual antigen stimulation",
                    )
                    add(
                        "entity_relationships",
                        {
                            "relationship_id": relationship_id,
                            "study_id": STUDY_ID,
                            "source_entity_type": "assays",
                            "source_entity_id": assay_id,
                            "target_entity_type": "antigens",
                            "target_entity_id": antigen_ids[mutation],
                            "relationship_type": "TESTS_RESPONSE_TO",
                        },
                        rel_prov,
                    )

                evidence_id = _id("evidence", candidate_id, "event_b")
                evidence_prov = _prov(
                    "recognition_evidence",
                    evidence_id,
                    document=DATA_FILE,
                    table="Figure 2c",
                    row=source_row + 6,
                    column=mutation,
                    fragment=f"primary Event-B call={int(raw_call)}",
                )
                add(
                    "recognition_evidence",
                    {
                        "evidence_id": evidence_id,
                        "candidate_id": candidate_id,
                        "patient_id": patient_id,
                        "evidence_family": EvidenceFamily.VACCINE_EVENT_B.value,
                        "source_dataset": f"{STUDY_ID}:Figure_2c",
                        "measured_or_predicted": "MEASURED_AUTHOR_CALL",
                        "value": int(raw_call),
                        "units": "binary",
                        "directionality": "higher_is_response",
                        "uncertainty": "binary source call; no calibrated probability",
                        "assay_or_model_version": "IFNg ELISpot primary endpoint",
                        "evidence_quality": "candidate-resolved explicit denominator",
                        "availability_status": "AVAILABLE",
                        "information_timing": InformationTiming.OUTCOME_ONLY.value,
                        "patient_specificity": 1.0,
                        "functional_relevance": 0.9,
                        "vaccine_relevance": 1.0,
                        "candidate_specificity": 1.0,
                        "assay_directness": 1.0,
                        "temporal_clarity": 1.0,
                        "source_completeness": 0.9,
                        "replication_status": "TECHNICAL_REPLICATES",
                    },
                    evidence_prov,
                )

                funnel_id = _id("funnel", candidate_id)
                funnel_prov = _prov(
                    "candidate_funnel_links",
                    funnel_id,
                    document=DATA_FILE,
                    table="Figure 2c",
                    row=source_row + 6,
                    column=mutation,
                    fragment="fixed vaccine component and explicit functional assay",
                )
                add(
                    "candidate_funnel_links",
                    {
                        "funnel_link_id": funnel_id,
                        "candidate_id": candidate_id,
                        "patient_id": patient_id,
                        "study_id": STUDY_ID,
                        "mutation_called": (
                            "REACHED" if TUMOR_MUTATIONS[source_patient] == mutation else "NOT_APPLICABLE"
                        ),
                        "transcript_represented": "NOT_ASSESSED",
                        "peptide_generated": "REACHED",
                        "survives_gating": "REACHED",
                        "hla_included": "NOT_ASSESSED",
                        "presentation_candidate": "NOT_ASSESSED",
                        "ranking_stage": "NOT_APPLICABLE_SHARED_FORMULATION",
                        "top_k": "NOT_APPLICABLE_SHARED_FORMULATION",
                        "recognition_scored": "REACHED",
                        "vaccine_inclusion": "REACHED",
                        "functional_assay": "REACHED",
                    },
                    funnel_prov,
                )

                for relationship_type, source_type, source_id, target_type, target_id in (
                    ("DERIVED_FROM", "candidates", candidate_id, "antigens", antigen_ids[mutation]),
                    (
                        "COMPONENT_OF_VACCINE",
                        "antigens",
                        antigen_ids[mutation],
                        "vaccines",
                        vaccine_id,
                    ),
                ):
                    relationship_id = _id(
                        "rel", source_id, relationship_type, target_id
                    )
                    rel_prov = _prov(
                        "entity_relationships",
                        relationship_id,
                        document=f"DOI:{DOI}",
                        table="Methods and Figure 2c",
                        row=source_patient,
                        column=mutation,
                        fragment=relationship_type,
                        method="deterministic_relationship",
                        origin=ValueOrigin.DETERMINISTICALLY_DERIVED.value,
                    )
                    add(
                        "entity_relationships",
                        {
                            "relationship_id": relationship_id,
                            "study_id": STUDY_ID,
                            "source_entity_type": source_type,
                            "source_entity_id": source_id,
                            "target_entity_type": target_type,
                            "target_entity_id": target_id,
                            "relationship_type": relationship_type,
                        },
                        rel_prov,
                    )

        return EventBCorpus(**{entity: _frame(entity, rows[entity]) for entity in SCHEMAS})


def reconcile(corpus: EventBCorpus, review_queue=()) -> dict:
    event_b = corpus.assays.event_type.astype(str).eq(
        BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
    )
    labels = corpus.assays.loc[event_b, "response_label"].astype(str)
    return {
        "source_reported": {
            "patients": 12,
            "primary_candidate_labels": 72,
            "positive_primary_labels": 60,
            "tested_negative_primary_labels": 12,
        },
        "extracted": {
            "patients": int(corpus.patients.patient_id.nunique()),
            "global_antigens": int(corpus.antigens.antigen_id.nunique()),
            "patient_candidates": int(corpus.candidates.candidate_id.nunique()),
            "assay_observations": int(corpus.assays.assay_id.nunique()),
            "event_a_observations": int((~event_b).sum()),
            "event_b_observations": int(event_b.sum()),
            "positive_primary_labels": int(labels.eq(ResponseLabel.POSITIVE.value).sum()),
            "tested_negative_primary_labels": int(
                labels.eq(ResponseLabel.TESTED_NEGATIVE.value).sum()
            ),
            "untested_candidates": int(labels.eq(ResponseLabel.UNTESTED.value).sum()),
            "review_queue": len(tuple(review_queue)),
        },
        "reconciles": bool(
            corpus.patients.patient_id.nunique() == 12
            and event_b.sum() == 72
            and labels.eq(ResponseLabel.POSITIVE.value).sum() == 60
            and labels.eq(ResponseLabel.TESTED_NEGATIVE.value).sum() == 12
        ),
        "material_discrepancies": [],
    }
