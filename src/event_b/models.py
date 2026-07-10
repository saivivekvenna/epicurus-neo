"""Versioned enums and canonical table schemas for vaccine-response evidence."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from hashlib import sha256
from typing import Any

import pandas as pd


SCHEMA_VERSION = "event-b-1.0.0"


class BiologicalEvent(str, Enum):
    EVENT_A_PREEXISTING_REACTIVITY = "EVENT_A_PREEXISTING_REACTIVITY"
    EVENT_B_VACCINE_INDUCED_RESPONSE = "EVENT_B_VACCINE_INDUCED_RESPONSE"
    # Post-vaccine T-cell response to a neoantigen that was NOT in the vaccine (epitope
    # spreading). Vaccine-induced in the broad sense, but not recognition of a vaccine
    # candidate, so it is kept distinct and never counted as an Event-B training label.
    EPITOPE_SPREADING = "EPITOPE_SPREADING"
    EVENT_C_CLINICAL_OUTCOME = "EVENT_C_CLINICAL_OUTCOME"
    PRESENTATION_ONLY = "PRESENTATION_ONLY"
    UNKNOWN_EVENT = "UNKNOWN_EVENT"


class ResponseLabel(str, Enum):
    POSITIVE = "POSITIVE"
    TESTED_NEGATIVE = "TESTED_NEGATIVE"
    UNTESTED = "UNTESTED"


class ReviewStatus(str, Enum):
    ACCEPTED = "ACCEPTED"
    NEEDS_REVIEW = "NEEDS_REVIEW"
    CONFLICTED = "CONFLICTED"
    REJECTED = "REJECTED"


class ValueOrigin(str, Enum):
    SOURCE_REPORTED = "SOURCE_REPORTED"
    DETERMINISTICALLY_DERIVED = "DETERMINISTICALLY_DERIVED"
    LLM_EXTRACTED = "LLM_EXTRACTED"
    MANUALLY_CURATED = "MANUALLY_CURATED"
    UNKNOWN = "UNKNOWN"


class AssayType(str, Enum):
    ELISPOT = "ELISPOT"
    ICS = "ICS"
    TETRAMER = "TETRAMER"
    PMHC_MULTIMER = "PMHC_MULTIMER"
    MANAFEST = "MANAFEST"
    TCR_EXPANSION = "TCR_EXPANSION"
    CYTOKINE_RELEASE = "CYTOKINE_RELEASE"
    KILLING_ASSAY = "KILLING_ASSAY"
    HEALTHY_DONOR_PRIMING = "HEALTHY_DONOR_PRIMING"
    SINGLE_CELL_TCR_PMHC = "SINGLE_CELL_TCR_PMHC"
    OTHER = "OTHER"


class EvidenceFamily(str, Enum):
    PRESENTATION_PREDICTION = "PRESENTATION_PREDICTION"
    IMMUNOPEPTIDOMICS = "IMMUNOPEPTIDOMICS"
    MUTANT_WILDTYPE_DIFFERENCE = "MUTANT_WILDTYPE_DIFFERENCE"
    SELF_PEPTIDOME_SIMILARITY = "SELF_PEPTIDOME_SIMILARITY"
    THYMIC_TOLERANCE_PROXY = "THYMIC_TOLERANCE_PROXY"
    KNOWN_IMMUNOGENICITY = "KNOWN_IMMUNOGENICITY"
    HEALTHY_DONOR_PRIMABILITY = "HEALTHY_DONOR_PRIMABILITY"
    TCR_BINDING = "TCR_BINDING"
    TCR_EXPANSION = "TCR_EXPANSION"
    FUNCTIONAL_T_CELL_ASSAY = "FUNCTIONAL_T_CELL_ASSAY"
    VACCINE_EVENT_B = "VACCINE_EVENT_B"
    CLINICAL_CONTEXT = "CLINICAL_CONTEXT"
    TUMOR_CLONALITY = "TUMOR_CLONALITY"
    SPATIAL_CONTEXT = "SPATIAL_CONTEXT"
    LONGITUDINAL_PERSISTENCE = "LONGITUDINAL_PERSISTENCE"


class InformationTiming(str, Enum):
    PRE_SELECTION = "PRE_SELECTION"
    OUTCOME_ONLY = "OUTCOME_ONLY"
    UNKNOWN = "UNKNOWN"


class AvailabilityStatus(str, Enum):
    AVAILABLE = "AVAILABLE"
    NOT_AVAILABLE = "NOT_AVAILABLE"
    NOT_ASSESSED = "NOT_ASSESSED"


class VaccineInclusion(str, Enum):
    INCLUDED = "INCLUDED"
    NOT_INCLUDED = "NOT_INCLUDED"
    UNKNOWN = "UNKNOWN"


class MHCClass(str, Enum):
    CLASS_I = "CLASS_I"
    CLASS_II = "CLASS_II"
    BOTH = "BOTH"
    UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class EntitySchema:
    name: str
    id_column: str
    required: tuple[str, ...]
    columns: tuple[str, ...]
    version: str = SCHEMA_VERSION

    def normalize(self, frame: pd.DataFrame) -> pd.DataFrame:
        missing = set(self.required).difference(frame.columns)
        if missing:
            raise ValueError(f"{self.name} missing required columns: {sorted(missing)}")
        unexpected = set(frame.columns).difference(self.columns)
        if unexpected:
            raise ValueError(
                f"{self.name} contains source-specific fields outside the canonical schema: "
                f"{sorted(unexpected)}"
            )
        out = frame.copy()
        for column in self.columns:
            if column not in out:
                out[column] = pd.NA
        if out[self.id_column].isna().any() or out[self.id_column].duplicated().any():
            raise ValueError(f"{self.name}.{self.id_column} must be non-null and unique")
        out["schema_version"] = self.version
        return out.loc[:, list(self.columns)]


SCHEMAS = {
    "studies": EntitySchema(
        "studies",
        "study_id",
        ("study_id", "title", "source_manifest_id", "provenance_id"),
        (
            "study_id",
            "title",
            "publication_ids",
            "trial_id",
            "cancer_type",
            "vaccine_platform",
            "adjuvant",
            "vaccination_schedule",
            "source_urls",
            "source_paths",
            "source_checksums",
            "source_manifest_id",
            "provenance_id",
            "schema_version",
        ),
    ),
    "patients": EntitySchema(
        "patients",
        "patient_id",
        ("patient_id", "source_patient_id", "study_id", "provenance_id"),
        (
            "patient_id",
            "source_patient_id",
            "study_id",
            "cancer_type",
            "disease_stage",
            "treatment_context",
            "prior_therapies",
            "hla_alleles",
            "tumor_context",
            "provenance_id",
            "schema_version",
        ),
    ),
    "vaccines": EntitySchema(
        "vaccines",
        "vaccine_id",
        ("vaccine_id", "patient_id", "study_id", "provenance_id"),
        (
            "vaccine_id",
            "patient_id",
            "study_id",
            "vaccine_platform",
            "formulation",
            "dose",
            "vaccination_dates",
            "relative_schedule",
            "candidate_count",
            "mhc_class_intent",
            "concurrent_therapy",
            "provenance_id",
            "schema_version",
        ),
    ),
    "candidates": EntitySchema(
        "candidates",
        "candidate_id",
        ("candidate_id", "patient_id", "study_id", "mutant_peptide", "provenance_id"),
        (
            "candidate_id",
            "patient_id",
            "study_id",
            "sample_id",
            "sample_date",
            "timepoint",
            "genomic_variant",
            "gene",
            "transcript",
            "protein_change",
            "mutant_peptide",
            "wildtype_peptide",
            "peptide_length",
            "hla_alleles",
            "mhc_class",
            "candidate_source",
            "vaccine_inclusion",
            "vaccine_inclusion_origin",
            "generation_provenance",
            "mutant_wildtype_verified",
            "provenance_id",
            "schema_version",
        ),
    ),
    "assays": EntitySchema(
        "assays",
        "assay_id",
        ("assay_id", "patient_id", "study_id", "event_type", "response_label", "provenance_id"),
        (
            "assay_id",
            "patient_id",
            "study_id",
            "candidate_id",
            "vaccine_id",
            "assay_type",
            "sample_type",
            "sample_date",
            "timepoint",
            "relative_to_vaccine",
            "stimulation_protocol",
            "replicate_information",
            "positivity_threshold",
            "quantitative_result",
            "result_units",
            "qualitative_result",
            "source_interpretation",
            "event_type",
            "response_label",
            "explicit_assay_inclusion",
            "review_status",
            "provenance_id",
            "schema_version",
        ),
    ),
    "clinical_outcomes": EntitySchema(
        "clinical_outcomes",
        "outcome_id",
        ("outcome_id", "patient_id", "study_id", "outcome_type", "provenance_id"),
        (
            "outcome_id",
            "patient_id",
            "study_id",
            "outcome_type",
            "outcome_value",
            "assessment_date",
            "response_criteria",
            "progression_free_survival",
            "overall_survival",
            "recurrence",
            "provenance_id",
            "schema_version",
        ),
    ),
    "recognition_evidence": EntitySchema(
        "recognition_evidence",
        "evidence_id",
        ("evidence_id", "candidate_id", "evidence_family", "source_dataset", "provenance_id"),
        (
            "evidence_id",
            "candidate_id",
            "patient_id",
            "evidence_family",
            "source_dataset",
            "measured_or_predicted",
            "value",
            "units",
            "directionality",
            "uncertainty",
            "assay_or_model_version",
            "evidence_quality",
            "availability_status",
            "information_timing",
            "patient_specificity",
            "functional_relevance",
            "vaccine_relevance",
            "candidate_specificity",
            "assay_directness",
            "temporal_clarity",
            "source_completeness",
            "replication_status",
            "provenance_id",
            "schema_version",
        ),
    ),
    "candidate_funnel_links": EntitySchema(
        "candidate_funnel_links",
        "funnel_link_id",
        ("funnel_link_id", "candidate_id", "patient_id", "provenance_id"),
        (
            "funnel_link_id",
            "candidate_id",
            "patient_id",
            "study_id",
            "mutation_called",
            "transcript_represented",
            "peptide_generated",
            "survives_gating",
            "hla_included",
            "presentation_candidate",
            "ranking_stage",
            "top_k",
            "recognition_scored",
            "vaccine_inclusion",
            "functional_assay",
            "provenance_id",
            "schema_version",
        ),
    ),
    "provenance": EntitySchema(
        "provenance",
        "provenance_id",
        (
            "provenance_id",
            "entity_type",
            "entity_id",
            "field_name",
            "value_origin",
            "review_status",
        ),
        (
            "provenance_id",
            "entity_type",
            "entity_id",
            "field_name",
            "source_document",
            "page",
            "table",
            "figure",
            "supplementary_file",
            "row",
            "column",
            "source_fragment",
            "extraction_method",
            "extraction_confidence",
            "value_origin",
            "review_status",
            "schema_version",
        ),
    ),
}


def stable_candidate_id(record: dict[str, Any]) -> str:
    """Hash patient-specific biological identity; peptide sequence alone is insufficient."""
    fields = (
        "study_id",
        "patient_id",
        "sample_id",
        "timepoint",
        "genomic_variant",
        "transcript",
        "mutant_peptide",
        "hla_alleles",
    )
    values = []
    for field in fields:
        value = record.get(field, "")
        if field == "hla_alleles" and isinstance(value, (list, tuple, set)):
            value = ",".join(sorted(str(item).upper() for item in value))
        values.append(str(value).strip().upper())
    if not values[0] or not values[1] or not values[6]:
        raise ValueError("stable candidate identity requires study, patient, and mutant peptide")
    return "eventb:" + sha256("|".join(values).encode()).hexdigest()[:24]
