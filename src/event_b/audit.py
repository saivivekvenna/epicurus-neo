"""Machine-readable and Markdown Event-B corpus audits."""

from __future__ import annotations

from dataclasses import asdict
import json
from typing import Iterable

import pandas as pd

from event_b.adapters.base import AdapterDeclaration
from event_b.corpus import EventBCorpus
from event_b.models import BiologicalEvent, ResponseLabel
from event_b.review import ReviewIssue


def _counts(series: pd.Series) -> dict[str, int]:
    return {
        str(key): int(value)
        for key, value in series.fillna("<MISSING>").astype(str).value_counts().sort_index().items()
    }


def _multi_counts(series: pd.Series) -> dict[str, int]:
    values = []
    for value in series.dropna():
        text = str(value).strip()
        if text.startswith("["):
            try:
                values.extend(str(item) for item in json.loads(text))
                continue
            except json.JSONDecodeError:
                pass
        values.extend(item.strip() for item in text.replace(";", ",").split(",") if item.strip())
    return _counts(pd.Series(values, dtype="object"))


def _missing(corpus: EventBCorpus) -> dict[str, dict[str, float]]:
    result = {}
    for entity, frame in corpus.tables().items():
        result[entity] = {
            column: float(frame[column].isna().mean()) if len(frame) else 1.0
            for column in sorted(frame.columns)
        }
    return result


def corpus_audit(
    corpus: EventBCorpus,
    issues: Iterable[ReviewIssue] = (),
    adapter_declarations: Iterable[AdapterDeclaration] = (),
    *,
    minimum_event_b_patients: int = 100,
    minimum_event_b_studies: int = 2,
    minimum_positive_patients: int = 30,
) -> dict:
    assays = corpus.assays.copy()
    candidates = corpus.candidates.copy()
    patients = corpus.patients.copy()
    event_b = (
        assays.event_type.astype(str)
        .str.upper()
        .eq(BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value)
    )
    positive = assays.response_label.astype(str).str.upper().eq(ResponseLabel.POSITIVE.value)
    tested_negative = (
        assays.response_label.astype(str).str.upper().eq(ResponseLabel.TESTED_NEGATIVE.value)
    )
    candidate_resolved = assays.candidate_id.notna()
    event_b_patients = set(
        assays.loc[event_b & (positive | tested_negative), "patient_id"].astype(str)
    )
    positive_patients = set(assays.loc[event_b & positive, "patient_id"].astype(str))
    event_b_studies = set(
        assays.loc[event_b & (positive | tested_negative), "study_id"].astype(str)
    )
    tested_patients = set(
        assays.loc[event_b & (positive | tested_negative), "patient_id"].astype(str)
    )
    negative_patients = set(assays.loc[event_b & tested_negative, "patient_id"].astype(str))
    resolved_patients = set(assays.loc[candidate_resolved, "patient_id"].astype(str))
    issue_rows = [asdict(issue) for issue in issues]

    candidate_groups = assays.dropna(subset=["candidate_id"]).groupby("candidate_id")
    multiple_assays = int((candidate_groups.assay_id.nunique() > 1).sum()) if len(assays) else 0
    conflicting = 0
    a_b_differ = 0
    for _, group in candidate_groups:
        labels = set(group.response_label.astype(str).str.upper())
        if {ResponseLabel.POSITIVE.value, ResponseLabel.TESTED_NEGATIVE.value}.issubset(labels):
            conflicting += 1
        by_event = group.groupby(group.event_type.astype(str).str.upper()).response_label.apply(
            lambda values: set(values.astype(str).str.upper())
        )
        if (
            BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value in by_event
            and BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value in by_event
            and by_event[BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value]
            != by_event[BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value]
        ):
            a_b_differ += 1

    patient_label_counts = (
        assays.groupby(["patient_id", "response_label"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .to_dict("records")
        if len(assays)
        else []
    )
    study_label_counts = (
        assays.groupby(["study_id", "response_label"], dropna=False)
        .size()
        .unstack(fill_value=0)
        .reset_index()
        .to_dict("records")
        if len(assays)
        else []
    )
    sufficient = (
        len(event_b_patients) >= minimum_event_b_patients
        and len(event_b_studies) >= minimum_event_b_studies
        and len(positive_patients) >= minimum_positive_patients
    )
    return {
        "sample_sizes": {
            "peptide_n": int(candidates.mutant_peptide.nunique()),
            "patient_n": int(patients.patient_id.nunique()),
            "study_n": int(corpus.studies.study_id.nunique()),
            "vaccine_n": int(corpus.vaccines.vaccine_id.nunique()),
            "candidate_n": int(candidates.candidate_id.nunique()),
            "assay_n": int(assays.assay_id.nunique()),
        },
        "event_counts": _counts(assays.event_type),
        "response_counts": _counts(assays.response_label),
        "review_status_counts": _counts(assays.review_status),
        "mhc_class_breakdown": _counts(candidates.mhc_class),
        "assay_type_breakdown": _counts(assays.assay_type),
        "hla_coverage": _multi_counts(candidates.hla_alleles),
        "cancer_type_coverage": _counts(patients.cancer_type),
        "vaccine_timepoint_coverage": _counts(assays.relative_to_vaccine),
        "patient_label_counts": patient_label_counts,
        "study_label_counts": study_label_counts,
        "candidate_linkage_rate": float(candidate_resolved.mean()) if len(assays) else 0.0,
        "funnel_linkage_rate": (
            float(
                corpus.candidate_funnel_links.candidate_id.nunique()
                / max(candidates.candidate_id.nunique(), 1)
            )
        ),
        "unresolved_contradictions": sum(
            row["code"] == "CONTRADICTORY_LABELS" for row in issue_rows
        ),
        "review_queue_n": len(issue_rows),
        "missing_data_rates": _missing(corpus),
        "evidence_source_coverage": _counts(corpus.recognition_evidence.source_dataset),
        "adapter_coverage": [asdict(declaration) for declaration in adapter_declarations],
        "patients_with_event_b_positive": len(positive_patients),
        "patients_with_only_tested_negatives": len(
            (negative_patients & tested_patients) - positive_patients
        ),
        "patients_without_candidate_resolved_assays": len(
            set(patients.patient_id.astype(str)) - resolved_patients
        ),
        "candidates_tested_by_multiple_assays": multiple_assays,
        "candidates_with_conflicting_assays": conflicting,
        "candidates_where_event_a_and_b_differ": a_b_differ,
        "model_readiness": {
            "sufficient_for_recognition_model_development": sufficient,
            "event_b_patient_n": len(event_b_patients),
            "event_b_study_n": len(event_b_studies),
            "event_b_positive_patient_n": len(positive_patients),
            "registered_minimums": {
                "event_b_patients": minimum_event_b_patients,
                "event_b_studies": minimum_event_b_studies,
                "positive_patients": minimum_positive_patients,
            },
            "decision": (
                "SUFFICIENT_FOR_BASELINE_DIAGNOSTICS"
                if sufficient
                else "INSUFFICIENT_DATA_DO_NOT_FIT_RECOGNITION_MODEL"
            ),
        },
    }


def render_audit_markdown(audit: dict) -> str:
    sizes = audit["sample_sizes"]
    readiness = audit["model_readiness"]
    lines = [
        "# Event-B corpus audit",
        "",
        "This corpus is a recognition-evidence substrate. It is not proof of clinical benefit.",
        "",
        "## Independent sample sizes",
        "",
        f"- Peptides: {sizes['peptide_n']}",
        f"- Patients: {sizes['patient_n']}",
        f"- Studies: {sizes['study_n']}",
        f"- Vaccines: {sizes['vaccine_n']}",
        f"- Candidate-assay observations: {sizes['assay_n']}",
        "",
        "## Event and response counts",
        "",
        "```json",
        json.dumps(
            {"events": audit["event_counts"], "responses": audit["response_counts"]},
            indent=2,
            sort_keys=True,
        ),
        "```",
        "",
        "## Model-readiness decision",
        "",
        f"**{readiness['decision']}**",
        "",
        f"Event-B patients: {readiness['event_b_patient_n']}; studies: {readiness['event_b_study_n']}; "
        f"positive patients: {readiness['event_b_positive_patient_n']}.",
        "",
        f"Review queue: {audit['review_queue_n']} records; unresolved contradictions: "
        f"{audit['unresolved_contradictions']}.",
        "",
    ]
    return "\n".join(lines)
