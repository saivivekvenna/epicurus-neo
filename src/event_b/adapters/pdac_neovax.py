"""Personalized PDAC mRNA NeoVax Event-B adapter (NCT04161755)."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import urllib.request

import pandas as pd

from event_b.adapters.base import AdapterDeclaration
from event_b.adapters.common import (
    entity_frame as _frame,
    provenance_record as _prov,
    stable_record_id as _id,
)
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


STUDY_ID = "pdac_neovax_2023"
COHORT_ID = "NCT04161755_phase1"
DOI = "10.1038/s41586-023-06063-y"
PMCID = "PMC10171177"
NCT = "NCT04161755"
DATA_FILE = "41586_2023_6063_MOESM4_ESM.xlsx"
SOURCE_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41586-023-06063-y/MediaObjects/" + DATA_FILE
)
EXPECTED_SHA256 = "c9942caa0de461e87c3725ae42377fea2a283cd170074fb2e36820303b429595"

SOURCE_COLUMNS = {
    "Patient number",
    "Neoantigen number",
    "Gene",
    "RefSeq transcript",
    "Substitution",
    "Mutant Neoantigen Sequence",
    "WT Neoantigen Sequence",
    "mRNA (+-13 AA (SNV); -15 AA to STOP (Indels))",
    "MHC-I  Allele (Best Prediction)",
    "MHC-I Mutant Epitope (Best Prediction)",
    "MHC-I WT Epitope",
    "MHC-II  Allele (Best Prediction)",
    "MHC-II Mutant Epitope (Best Prediction)",
    "MHC-II WT Epitope",
    "ELISpot Response",
}
RESPONSE_MAP = {
    "De novo response": ResponseLabel.POSITIVE.value,
    "No response": ResponseLabel.TESTED_NEGATIVE.value,
    "De novo response in pool": ResponseLabel.UNTESTED.value,
    "No data": ResponseLabel.UNTESTED.value,
}


def stage_source(raw_dir: str | Path) -> Path:
    root = Path(raw_dir)
    root.mkdir(parents=True, exist_ok=True)
    path = root / DATA_FILE
    if not path.exists() or sha256_file(path) != EXPECTED_SHA256:
        request = urllib.request.Request(SOURCE_URL, headers={"User-Agent": "epicurus-neo/1.0"})
        try:
            with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as out:
                shutil.copyfileobj(response, out)
        except Exception as error:  # noqa: BLE001 - converted to a manual-source contract
            raise RuntimeError(
                f"Download {SOURCE_URL} manually to {path}; expected sha256={EXPECTED_SHA256}"
            ) from error
    observed = sha256_file(path)
    if observed != EXPECTED_SHA256:
        raise RuntimeError(f"checksum mismatch for {path}: observed {observed}")
    return path


def source_manifest(raw_dir: str | Path) -> SourceManifest:
    adapter = PDACNeoVaxAdapter(raw_dir)
    return manifest_from_paths(
        adapter.declaration.source_name,
        adapter.declaration.source_version,
        adapter.declaration.adapter_name,
        adapter.declaration.adapter_version,
        [stage_source(raw_dir)],
    )


def _text(value: object) -> str:
    return "" if value is None or pd.isna(value) else str(value).strip()


def _int_text(value: object) -> str:
    return str(int(float(value)))


def _source_identifier(value: object) -> str:
    text = _text(value)
    try:
        number = float(text)
    except ValueError:
        return text
    return str(int(number)) if number.is_integer() else text


class PDACNeoVaxAdapter:
    declaration = AdapterDeclaration(
        "Rojas 2023 personalized PDAC mRNA NeoVax",
        f"{PMCID}/Nature-2023",
        "pdac_neovax_event_b",
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
            "Seven patient-25 candidates occur only in two positive multi-neoantigen pools; their "
            "individual labels remain UNTESTED.",
            "Two vaccine targets have no ELISpot data and remain UNTESTED.",
            "The follow-up Nature 2025 publication reuses this cohort and is not an independent study.",
            "All participants received atezolizumab and most subsequently received mFOLFIRINOX.",
        ),
        (
            "The source's 'De novo response' and 'No response' candidate rows are accepted as the "
            "explicit single-target author calls.",
            "'De novo response in pool' is never decomposed into positive candidate labels.",
            "The encoded neoantigen, overlapping 15-mer assay pool, and predicted minimal class-I/II "
            "epitopes are distinct antigen entities.",
        ),
        ("patient_hla_genotype", "raw_elispot_replicates_for_nonresponding_targets"),
        canonical_study_id=STUDY_ID,
        cohort_id=COHORT_ID,
        source_files=(DATA_FILE,),
        supported_timepoints=("POST_ATEZOLIZUMAB_PRE_VACCINE", "POST_VACCINE_PRIMING"),
        positivity_rules=(
            "at least 7 spots/300,000 PBMC and significant increase versus medium control",
        ),
        baseline_semantics="No vaccine-neoantigen responses detected before vaccination.",
        vaccine_component_structure="Up to 20 personalized mRNA-encoded long neoantigens.",
        assay_target_structure="One overlapping 15-mer pool per target, except two combined pools.",
        candidate_identity_completeness="PATIENT_AND_ENCODED_NEOANTIGEN_RESOLVED",
        unresolved_ambiguities=(
            "The public target table does not partition the seven patient-25 targets between its "
            "two positive combined pools.",
        ),
    )

    def __init__(self, raw_dir: str | Path):
        self.raw_dir = Path(raw_dir)
        self.review_issues = []

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        targets = pd.read_excel(
            stage_source(self.raw_dir),
            sheet_name="ml41081_targets_with_elispot",
            engine="openpyxl",
        )
        return {"targets": targets}

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        targets = pd.DataFrame(extracted["targets"])
        if set(targets.columns) != SOURCE_COLUMNS:
            raise ValueError("PDAC NeoVax target-table columns changed")
        responses = set(targets["ELISpot Response"].astype(str))
        if responses != set(RESPONSE_MAP):
            raise ValueError(f"unrecognized ELISpot response values: {sorted(responses)}")
        if len(targets) != 232 or targets["Patient number"].nunique() != 16:
            raise ValueError("PDAC NeoVax public denominator changed")

        rows: dict[str, list[dict]] = {entity: [] for entity in SCHEMAS}

        def add(entity: str, row: dict, provenance: dict) -> None:
            row["provenance_id"] = provenance["provenance_id"]
            rows[entity].append(row)
            rows["provenance"].append(provenance)

        add(
            "studies",
            {
                "study_id": STUDY_ID,
                "title": "Personalized RNA neoantigen vaccines stimulate T cells in PDAC",
                "publication_ids": json.dumps(
                    [f"DOI:{DOI}", f"{PMCID}", "PMID:37165196", "DOI:10.1038/s41586-024-08508-4"]
                ),
                "trial_id": NCT,
                "cancer_type": "pancreatic ductal adenocarcinoma",
                "vaccine_platform": "individualized mRNA-lipoplex neoantigen vaccine",
                "adjuvant": "RNA-lipoplex innate stimulation",
                "vaccination_schedule": "eight prime doses and one boost",
                "source_urls": json.dumps([SOURCE_URL]),
                "source_paths": json.dumps([doc.local_path for doc in manifest.documents]),
                "source_checksums": json.dumps({DATA_FILE: EXPECTED_SHA256}),
                "source_manifest_id": manifest.manifest_id,
            },
            _prov(
                "studies",
                STUDY_ID,
                document=f"DOI:{DOI}",
                table="main article",
                row="trial identity",
                column="cohort",
                fragment="NCT04161755; 16 vaccinated biomarker-evaluable participants",
                method="manual_primary_source_curation",
                origin=ValueOrigin.MANUALLY_CURATED.value,
            ),
        )

        patient_numbers = sorted({_int_text(value) for value in targets["Patient number"]})
        for source_patient in patient_numbers:
            patient_id = f"{STUDY_ID}:patient_{source_patient}"
            count = int(
                targets["Patient number"].map(_int_text).eq(source_patient).sum()
            )
            add(
                "patients",
                {
                    "patient_id": patient_id,
                    "source_patient_id": source_patient,
                    "study_id": STUDY_ID,
                    "cancer_type": "pancreatic ductal adenocarcinoma",
                    "disease_stage": "surgically resected",
                    "treatment_context": "sequential atezolizumab, vaccine, mFOLFIRINOX",
                    "prior_therapies": "surgery; no neoadjuvant chemotherapy",
                    "hla_alleles": json.dumps([]),
                    "tumor_context": "personalized tumor-derived somatic mutations",
                },
                _prov(
                    "patients",
                    patient_id,
                    document=DATA_FILE,
                    table="Supplementary Table 5",
                    row=source_patient,
                    column="Patient number",
                    fragment=f"patient {source_patient}; {count} vaccine targets",
                ),
            )
            vaccine_id = f"vaccine:{STUDY_ID}:patient_{source_patient}"
            add(
                "vaccines",
                {
                    "vaccine_id": vaccine_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "vaccine_platform": "individualized uridine mRNA-lipoplex vaccine",
                    "formulation": "autogene cevumeran personalized target set",
                    "vaccination_dates": json.dumps([]),
                    "relative_schedule": "eight priming doses followed by booster",
                    "candidate_count": count,
                    "mhc_class_intent": MHCClass.BOTH.value,
                    "concurrent_therapy": "atezolizumab before vaccine; mFOLFIRINOX after priming",
                },
                _prov(
                    "vaccines",
                    vaccine_id,
                    document=DATA_FILE,
                    table="Supplementary Table 5",
                    row=source_patient,
                    column="Neoantigen number",
                    fragment=f"{count} manufactured personalized targets",
                ),
            )

        for source_index, source in targets.reset_index(drop=True).iterrows():
            patient_number = _int_text(source["Patient number"])
            neo_number = _source_identifier(source["Neoantigen number"])
            patient_id = f"{STUDY_ID}:patient_{patient_number}"
            vaccine_id = f"vaccine:{STUDY_ID}:patient_{patient_number}"
            base = f"antigen:pdac:{patient_number}:{neo_number}"
            mutant = _text(source["Mutant Neoantigen Sequence"])
            wildtype = _text(source["WT Neoantigen Sequence"])
            hla_i = _text(source["MHC-I  Allele (Best Prediction)"])
            hla_ii = _text(source["MHC-II  Allele (Best Prediction)"])
            antigen_specs = (
                (base + ":encoded", mutant, wildtype, "ENCODED_MRNA_NEOANTIGEN", [], "NOT_APPLICABLE"),
                (
                    base + ":assay_pool",
                    mutant,
                    wildtype,
                    "OVERLAPPING_15MER_ASSAY_POOL_COVERAGE",
                    [],
                    "NOT_ASSESSED",
                ),
                (
                    base + ":mhci_pred",
                    _text(source["MHC-I Mutant Epitope (Best Prediction)"]),
                    _text(source["MHC-I WT Epitope"]),
                    "PREDICTED_MINIMAL_EPITOPE_CLASS_I",
                    [hla_i],
                    "PREDICTED_BEST_BINDER",
                ),
                (
                    base + ":mhcii_pred",
                    _text(source["MHC-II Mutant Epitope (Best Prediction)"]),
                    _text(source["MHC-II WT Epitope"]),
                    "PREDICTED_MINIMAL_EPITOPE_CLASS_II",
                    [hla_ii],
                    "PREDICTED_BEST_BINDER",
                ),
            )
            antigen_ids = {}
            for antigen_id, sequence, wt_sequence, component_type, alleles, evidence_type in antigen_specs:
                antigen_ids[component_type] = antigen_id
                add(
                    "antigens",
                    {
                        "antigen_id": antigen_id,
                        "study_id": STUDY_ID,
                        "gene": _text(source["Gene"]),
                        "protein_change": _text(source["Substitution"]),
                        "mutant_sequence": sequence,
                        "wildtype_sequence": wt_sequence,
                        "component_type": component_type,
                        "peptide_length": len(sequence),
                        "hla_alleles": json.dumps(alleles),
                        "hla_evidence_type": evidence_type,
                    },
                    _prov(
                        "antigens",
                        antigen_id,
                        document=DATA_FILE,
                        table="Supplementary Table 5",
                        row=source_index + 2,
                        column=component_type,
                        fragment=f"patient {patient_number}; target {neo_number}; {sequence}",
                    ),
                )

            encoded_id = antigen_ids["ENCODED_MRNA_NEOANTIGEN"]
            assay_target_id = antigen_ids["OVERLAPPING_15MER_ASSAY_POOL_COVERAGE"]
            record = {
                "study_id": STUDY_ID,
                "patient_id": patient_id,
                "sample_id": "resected_tumor",
                "timepoint": "MANUFACTURED_VACCINE_TARGET",
                "genomic_variant": _text(source["Substitution"]),
                "gene": _text(source["Gene"]),
                "transcript": _text(source["RefSeq transcript"]),
                "protein_change": _text(source["Substitution"]),
                "mutant_peptide": mutant,
                "wildtype_peptide": wildtype,
                "peptide_length": len(mutant),
                "hla_alleles": json.dumps([]),
                "mhc_class": MHCClass.UNKNOWN.value,
                "candidate_source": "personalized encoded target; HLA fields are predictions only",
                "vaccine_inclusion": VaccineInclusion.INCLUDED.value,
                "vaccine_inclusion_origin": ValueOrigin.SOURCE_REPORTED.value,
                "generation_provenance": _text(
                    source["mRNA (+-13 AA (SNV); -15 AA to STOP (Indels))"]
                ),
                "mutant_wildtype_verified": True,
            }
            candidate_id = stable_candidate_id(record)
            response_text = _text(source["ELISpot Response"])
            add(
                "candidates",
                {"candidate_id": candidate_id, **record},
                _prov(
                    "candidates",
                    candidate_id,
                    document=DATA_FILE,
                    table="Supplementary Table 5",
                    row=source_index + 2,
                    column="Mutant Neoantigen Sequence",
                    fragment=f"patient {patient_number}; target {neo_number}; {response_text}",
                ),
            )

            candidate_resolved = response_text in {"De novo response", "No response"}
            post_label = RESPONSE_MAP[response_text]
            for relative, event, timepoint, label in (
                (
                    "PRE_VACCINE",
                    BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value,
                    "POST_ATEZOLIZUMAB_PRE_VACCINE",
                    ResponseLabel.TESTED_NEGATIVE.value
                    if candidate_resolved
                    else ResponseLabel.UNTESTED.value,
                ),
                (
                    "POST_PRIME",
                    BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                    "POST_VACCINE_PRIMING",
                    post_label,
                ),
            ):
                assay_id = _id("assay", candidate_id, event)
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
                        "stimulation_protocol": "overlapping 15-mer pool covering one encoded target",
                        "replicate_information": "duplicate wells",
                        "positivity_threshold": "source-defined >=7 spots and significant increase",
                        "qualitative_result": label,
                        "source_interpretation": response_text,
                        "event_type": event,
                        "response_label": label,
                        "explicit_assay_inclusion": candidate_resolved,
                        "review_status": ReviewStatus.ACCEPTED.value,
                    },
                    _prov(
                        "assays",
                        assay_id,
                        document=DATA_FILE,
                        table="Supplementary Table 5",
                        row=source_index + 2,
                        column="ELISpot Response",
                        fragment=response_text,
                    ),
                )
                relationship_id = _id("rel", assay_id, "TESTS_RESPONSE_TO", assay_target_id)
                add(
                    "entity_relationships",
                    {
                        "relationship_id": relationship_id,
                        "study_id": STUDY_ID,
                        "source_entity_type": "assays",
                        "source_entity_id": assay_id,
                        "target_entity_type": "antigens",
                        "target_entity_id": assay_target_id,
                        "relationship_type": "TESTS_RESPONSE_TO",
                    },
                    _prov(
                        "entity_relationships",
                        relationship_id,
                        document=DATA_FILE,
                        table="Supplementary Table 5",
                        row=source_index + 2,
                        column="ELISpot Response",
                        fragment="assay target is overlapping 15-mer coverage pool",
                    ),
                )

            evidence_id = _id("evidence", candidate_id, "event_b")
            add(
                "recognition_evidence",
                {
                    "evidence_id": evidence_id,
                    "candidate_id": candidate_id,
                    "patient_id": patient_id,
                    "evidence_family": EvidenceFamily.VACCINE_EVENT_B.value,
                    "source_dataset": f"{STUDY_ID}:Supplementary_Table_5",
                    "measured_or_predicted": "MEASURED_AUTHOR_CALL",
                    "value": post_label,
                    "units": "three-state label",
                    "directionality": "positive_is_response",
                    "uncertainty": "pool-only and no-data rows remain UNTESTED",
                    "assay_or_model_version": "ex-vivo IFNg ELISpot",
                    "evidence_quality": "candidate-resolved" if candidate_resolved else "pool-only",
                    "availability_status": "AVAILABLE" if candidate_resolved else "NOT_ASSESSED",
                    "information_timing": InformationTiming.OUTCOME_ONLY.value,
                    "patient_specificity": 1.0,
                    "functional_relevance": 0.9,
                    "vaccine_relevance": 1.0,
                    "candidate_specificity": 1.0 if candidate_resolved else 0.0,
                    "assay_directness": 1.0,
                    "temporal_clarity": 1.0,
                    "source_completeness": 0.9 if candidate_resolved else 0.5,
                    "replication_status": "DUPLICATE_WELLS",
                },
                _prov(
                    "recognition_evidence",
                    evidence_id,
                    document=DATA_FILE,
                    table="Supplementary Table 5",
                    row=source_index + 2,
                    column="ELISpot Response",
                    fragment=response_text,
                ),
            )

            funnel_id = _id("funnel", candidate_id)
            add(
                "candidate_funnel_links",
                {
                    "funnel_link_id": funnel_id,
                    "candidate_id": candidate_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "mutation_called": "REACHED",
                    "transcript_represented": "REACHED",
                    "peptide_generated": "REACHED",
                    "survives_gating": "REACHED",
                    "hla_included": "REACHED_PREDICTED",
                    "presentation_candidate": "REACHED_PREDICTED",
                    "ranking_stage": "REACHED",
                    "top_k": "REACHED",
                    "recognition_scored": "REACHED" if candidate_resolved else "NOT_ASSESSED",
                    "vaccine_inclusion": "REACHED",
                    "functional_assay": "REACHED" if candidate_resolved else "POOL_ONLY",
                },
                _prov(
                    "candidate_funnel_links",
                    funnel_id,
                    document=DATA_FILE,
                    table="Supplementary Table 5",
                    row=source_index + 2,
                    column="ELISpot Response",
                    fragment=response_text,
                ),
            )

            relationships = [
                ("candidates", candidate_id, "antigens", encoded_id, "DERIVED_FROM"),
                ("antigens", encoded_id, "vaccines", vaccine_id, "COMPONENT_OF_VACCINE"),
                (
                    "antigens",
                    assay_target_id,
                    "antigens",
                    encoded_id,
                    "DERIVED_FROM",
                ),
            ]
            for component_type in (
                "PREDICTED_MINIMAL_EPITOPE_CLASS_I",
                "PREDICTED_MINIMAL_EPITOPE_CLASS_II",
            ):
                relationships.append(
                    ("antigens", antigen_ids[component_type], "antigens", encoded_id, "CONTAINED_WITHIN")
                )
            for source_type, source_id, target_type, target_id, relationship_type in relationships:
                relationship_id = _id("rel", source_id, relationship_type, target_id)
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
                    _prov(
                        "entity_relationships",
                        relationship_id,
                        document=DATA_FILE,
                        table="Supplementary Table 5",
                        row=source_index + 2,
                        column="entity relationship",
                        fragment=relationship_type,
                        method="deterministic_relationship",
                        origin=ValueOrigin.DETERMINISTICALLY_DERIVED.value,
                    ),
                )

        patient25 = f"{STUDY_ID}:patient_25"
        pool_assay_id = "assay:pdac:patient25:two_positive_combined_pools"
        add(
            "assays",
            {
                "assay_id": pool_assay_id,
                "patient_id": patient25,
                "study_id": STUDY_ID,
                "vaccine_id": "vaccine:pdac_neovax_2023:patient_25",
                "assay_type": AssayType.ELISPOT.value,
                "sample_type": "PBMC",
                "timepoint": "POST_VACCINE_PRIMING",
                "relative_to_vaccine": "POST_PRIME",
                "stimulation_protocol": "two positive combined pools spanning seven neoantigens",
                "positivity_threshold": "source-defined >=7 spots and significant increase",
                "quantitative_result": 2,
                "result_units": "positive combined pools",
                "qualitative_result": ResponseLabel.POSITIVE.value,
                "source_interpretation": "pool-level response only; no candidate decomposition",
                "event_type": BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                "response_label": ResponseLabel.POSITIVE.value,
                "explicit_assay_inclusion": True,
                "review_status": ReviewStatus.ACCEPTED.value,
            },
            _prov(
                "assays",
                pool_assay_id,
                document=f"DOI:{DOI}",
                table="Figure 1e caption",
                row="patient 25",
                column="pool response",
                fragment="2 positive pools containing 2 and 5 neoantigens",
                method="manual_primary_source_curation",
                origin=ValueOrigin.MANUALLY_CURATED.value,
            ),
        )
        return EventBCorpus(**{entity: _frame(entity, rows[entity]) for entity in SCHEMAS})


def reconcile(corpus: EventBCorpus, review_queue=()) -> dict:
    assays = corpus.assays
    event_b = assays.event_type.astype(str).eq(
        BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
    )
    primary = assays[event_b & assays.candidate_id.notna()]
    counts = primary.response_label.astype(str).value_counts().to_dict()
    pool_only = assays[event_b & assays.candidate_id.isna()]
    return {
        "source_reported": {
            "patients": 16,
            "vaccine_targets": 232,
            "single_target_positive": 23,
            "pool_only_candidates": 7,
            "no_data_candidates": 2,
            "positive_combined_pools": 2,
        },
        "extracted": {
            "patients": int(corpus.patients.patient_id.nunique()),
            "vaccine_targets": int(corpus.candidates.candidate_id.nunique()),
            "primary_candidate_labels": int(len(primary)),
            "positive": int(counts.get(ResponseLabel.POSITIVE.value, 0)),
            "tested_negative": int(counts.get(ResponseLabel.TESTED_NEGATIVE.value, 0)),
            "untested": int(counts.get(ResponseLabel.UNTESTED.value, 0)),
            "pool_level_observations": int(len(pool_only)),
            "assay_observations": int(len(assays)),
            "review_queue": len(tuple(review_queue)),
        },
        "reconciles": bool(
            corpus.patients.patient_id.nunique() == 16
            and len(primary) == 232
            and counts.get(ResponseLabel.POSITIVE.value, 0) == 23
            and counts.get(ResponseLabel.TESTED_NEGATIVE.value, 0) == 200
            and counts.get(ResponseLabel.UNTESTED.value, 0) == 9
            and len(pool_only) == 1
        ),
        "material_discrepancies": [
            "The paper's 25 response entities comprise 23 single-target calls plus 2 positive "
            "combined pools; the seven pool members remain UNTESTED individually."
        ],
    }
