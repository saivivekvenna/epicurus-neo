"""Deterministic corpus validation and contradiction detection."""

from __future__ import annotations

from dataclasses import dataclass
import json
from typing import Any

import pandas as pd

from event_b.corpus import EventBCorpus
from event_b.evidence import validate_evidence
from event_b.models import (
    BiologicalEvent,
    MHCClass,
    ResponseLabel,
    ReviewStatus,
    ValueOrigin,
    VaccineInclusion,
)
from event_b.review import ReviewIssue


@dataclass(frozen=True)
class ValidationResult:
    normalized_corpus: EventBCorpus
    accepted_corpus: EventBCorpus
    review_queue: tuple[ReviewIssue, ...]

    @property
    def ok(self) -> bool:
        return not self.review_queue


def _values(value: Any) -> list[str]:
    if value is None or value is pd.NA or (isinstance(value, float) and pd.isna(value)):
        return []
    if isinstance(value, (list, tuple, set)):
        return [str(item).strip().upper() for item in value if str(item).strip()]
    text = str(value).strip()
    if not text:
        return []
    if text.startswith("["):
        try:
            return [str(item).strip().upper() for item in json.loads(text)]
        except json.JSONDecodeError:
            pass
    return [item.strip().upper() for item in text.replace(";", ",").split(",") if item.strip()]


def _enum(value: Any, enum_type, default=None):
    if isinstance(value, enum_type):
        return value
    text = str(value).strip().upper()
    try:
        return enum_type(text)
    except ValueError:
        return default


def _issue(issues, entity, entity_id, code, message):
    issues.append(ReviewIssue.create(entity, str(entity_id), code, message))


def _validate_links(corpus: EventBCorpus, issues: list[ReviewIssue]) -> None:
    patients = corpus.patients.set_index("patient_id", drop=False)
    studies = set(corpus.studies["study_id"].astype(str))
    candidates = corpus.candidates.set_index("candidate_id", drop=False)
    vaccines = corpus.vaccines.set_index("vaccine_id", drop=False)
    for entity in ("patients", "vaccines", "candidates", "assays", "clinical_outcomes"):
        frame = getattr(corpus, entity)
        id_col = next(column for column in frame.columns if column.endswith("_id"))
        for _, row in frame.iterrows():
            entity_id = row[id_col]
            study_id = str(row.get("study_id", ""))
            if study_id and study_id not in studies:
                _issue(issues, entity, entity_id, "UNKNOWN_STUDY", f"study {study_id} is absent")
            patient_id = row.get("patient_id")
            if patient_id is not None and not pd.isna(patient_id):
                if patient_id not in patients.index:
                    _issue(
                        issues,
                        entity,
                        entity_id,
                        "UNKNOWN_PATIENT",
                        f"patient {patient_id} is absent",
                    )
                elif study_id and str(patients.loc[patient_id, "study_id"]) != study_id:
                    _issue(
                        issues,
                        entity,
                        entity_id,
                        "STUDY_MISMATCH",
                        "patient and record studies differ",
                    )
    for _, assay in corpus.assays.iterrows():
        candidate_id = assay.get("candidate_id")
        if candidate_id is not None and not pd.isna(candidate_id):
            if candidate_id not in candidates.index:
                _issue(issues, "assays", assay.assay_id, "UNKNOWN_CANDIDATE", str(candidate_id))
            elif str(candidates.loc[candidate_id, "study_id"]) != str(assay.study_id):
                _issue(
                    issues,
                    "assays",
                    assay.assay_id,
                    "STUDY_MISMATCH",
                    "candidate and assay studies differ",
                )
        vaccine_id = assay.get("vaccine_id")
        if vaccine_id is not None and not pd.isna(vaccine_id):
            if vaccine_id not in vaccines.index:
                _issue(issues, "assays", assay.assay_id, "UNKNOWN_VACCINE", str(vaccine_id))
            elif str(vaccines.loc[vaccine_id, "patient_id"]) != str(assay.patient_id):
                _issue(
                    issues,
                    "assays",
                    assay.assay_id,
                    "PATIENT_MISMATCH",
                    "vaccine and assay patients differ",
                )
    for _, evidence in corpus.recognition_evidence.iterrows():
        if str(evidence.candidate_id) not in candidates.index:
            _issue(
                issues,
                "recognition_evidence",
                evidence.evidence_id,
                "UNKNOWN_CANDIDATE",
                str(evidence.candidate_id),
            )
        elif evidence.patient_id is not None and not pd.isna(evidence.patient_id):
            if str(candidates.loc[str(evidence.candidate_id), "patient_id"]) != str(
                evidence.patient_id
            ):
                _issue(
                    issues,
                    "recognition_evidence",
                    evidence.evidence_id,
                    "PATIENT_MISMATCH",
                    "candidate and evidence patients differ",
                )


def _validate_candidates(corpus: EventBCorpus, issues: list[ReviewIssue]) -> None:
    patients = corpus.patients.set_index("patient_id", drop=False)
    for _, row in corpus.candidates.iterrows():
        candidate_id = row.candidate_id
        peptide = str(row.mutant_peptide).strip().upper()
        try:
            length = int(row.peptide_length)
            if length != len(peptide):
                _issue(
                    issues,
                    "candidates",
                    candidate_id,
                    "PEPTIDE_LENGTH",
                    f"stored {length}, actual {len(peptide)}",
                )
        except (TypeError, ValueError):
            _issue(
                issues,
                "candidates",
                candidate_id,
                "PEPTIDE_LENGTH",
                "missing or invalid peptide length",
            )
        mhc_class = _enum(row.mhc_class, MHCClass, MHCClass.UNKNOWN)
        if mhc_class is MHCClass.CLASS_I and not 8 <= len(peptide) <= 14:
            _issue(
                issues,
                "candidates",
                candidate_id,
                "MHC_LENGTH",
                "class I peptide length outside 8-14",
            )
        if mhc_class is MHCClass.CLASS_II and not 12 <= len(peptide) <= 30:
            _issue(
                issues,
                "candidates",
                candidate_id,
                "MHC_LENGTH",
                "class II peptide length outside 12-30",
            )
        wildtype = str(row.wildtype_peptide).strip().upper()
        if wildtype and wildtype not in {"<NA>", "NAN"} and wildtype == peptide:
            _issue(
                issues,
                "candidates",
                candidate_id,
                "PEPTIDE_ROLE",
                "mutant and wild-type peptides are identical",
            )
        if row.get("mutant_wildtype_verified") is False:
            _issue(
                issues,
                "candidates",
                candidate_id,
                "PEPTIDE_ROLE",
                "mutant/wild-type roles are unverified",
            )
        inclusion = _enum(row.vaccine_inclusion, VaccineInclusion, VaccineInclusion.UNKNOWN)
        origin = _enum(row.vaccine_inclusion_origin, ValueOrigin, ValueOrigin.UNKNOWN)
        if inclusion is VaccineInclusion.INCLUDED and origin is ValueOrigin.UNKNOWN:
            _issue(
                issues,
                "candidates",
                candidate_id,
                "VACCINE_INCLUSION",
                "included status has no value origin",
            )
        if row.patient_id in patients.index:
            patient_hla = set(_values(patients.loc[row.patient_id, "hla_alleles"]))
            candidate_hla = set(_values(row.hla_alleles))
            if patient_hla and candidate_hla and not candidate_hla.issubset(patient_hla):
                _issue(
                    issues,
                    "candidates",
                    candidate_id,
                    "HLA_MISMATCH",
                    f"{sorted(candidate_hla - patient_hla)} not in genotype",
                )


def _validate_assays(corpus: EventBCorpus, issues: list[ReviewIssue]) -> None:
    first_vaccine_dates: dict[str, pd.Timestamp] = {}
    for _, vaccine in corpus.vaccines.iterrows():
        dates = _values(vaccine.vaccination_dates)
        parsed = [pd.to_datetime(value, errors="coerce", utc=True) for value in dates]
        parsed = [value for value in parsed if not pd.isna(value)]
        if parsed:
            first_vaccine_dates[str(vaccine.vaccine_id)] = min(parsed)

    for _, row in corpus.assays.iterrows():
        assay_id = row.assay_id
        event = _enum(row.event_type, BiologicalEvent, BiologicalEvent.UNKNOWN_EVENT)
        label = _enum(row.response_label, ResponseLabel)
        relative = str(row.relative_to_vaccine).strip().upper()
        if label is None:
            _issue(
                issues,
                "assays",
                assay_id,
                "LABEL",
                f"unknown response label {row.response_label!r}",
            )
            continue
        if event is BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE and relative not in {
            "POST_VACCINE",
            "POST_PRIME",
            "POST_BOOST",
        }:
            _issue(
                issues,
                "assays",
                assay_id,
                "EVENT_B_TIMEPOINT",
                "Event B requires explicit post-vaccine evidence",
            )
        if relative == "PRE_VACCINE" and event is BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE:
            _issue(
                issues,
                "assays",
                assay_id,
                "EVENT_B_PRE",
                "pre-vaccine-only evidence cannot be Event B",
            )
        explicit = row.explicit_assay_inclusion
        explicit_inclusion = explicit is True or str(explicit).strip().upper() in {
            "TRUE",
            "1",
            "YES",
        }
        if label is ResponseLabel.TESTED_NEGATIVE and not explicit_inclusion:
            _issue(
                issues,
                "assays",
                assay_id,
                "NEGATIVE_DENOMINATOR",
                "tested negative requires explicit assay inclusion",
            )
        quantitative = pd.to_numeric(pd.Series([row.quantitative_result]), errors="coerce").iloc[0]
        positive_qualitative = str(row.qualitative_result).strip().upper() in {
            "POSITIVE",
            "REACTIVE",
            "YES",
        }
        if label is ResponseLabel.UNTESTED and (
            positive_qualitative or (pd.notna(quantitative) and quantitative > 0)
        ):
            _issue(
                issues,
                "assays",
                assay_id,
                "UNTESTED_RESULT",
                "untested record has a positive assay result",
            )
        if event is BiologicalEvent.EVENT_C_CLINICAL_OUTCOME:
            _issue(
                issues,
                "assays",
                assay_id,
                "CLINICAL_AS_ASSAY",
                "clinical outcomes must remain separate",
            )
        vaccine_id = str(row.vaccine_id)
        sample_date = pd.to_datetime(row.sample_date, errors="coerce", utc=True)
        first_date = first_vaccine_dates.get(vaccine_id)
        if first_date is not None and pd.notna(sample_date):
            if relative.startswith("POST") and sample_date < first_date:
                _issue(
                    issues,
                    "assays",
                    assay_id,
                    "CHRONOLOGY",
                    "post-vaccine sample predates vaccination",
                )
            if relative == "PRE_VACCINE" and sample_date > first_date:
                _issue(
                    issues,
                    "assays",
                    assay_id,
                    "CHRONOLOGY",
                    "pre-vaccine sample postdates vaccination",
                )

    accepted = corpus.assays.loc[
        corpus.assays["review_status"].astype(str).str.upper().eq(ReviewStatus.ACCEPTED.value)
    ].copy()
    accepted["_label"] = accepted["response_label"].astype(str).str.upper()
    group_cols = ["candidate_id", "assay_type", "timepoint"]
    for keys, group in accepted.dropna(subset=["candidate_id"]).groupby(group_cols, dropna=False):
        labels = set(group["_label"])
        if {ResponseLabel.POSITIVE.value, ResponseLabel.TESTED_NEGATIVE.value}.issubset(labels):
            for assay_id in group.assay_id:
                _issue(
                    issues,
                    "assays",
                    assay_id,
                    "CONTRADICTORY_LABELS",
                    f"conflicting accepted labels for {keys}",
                )


def _validate_provenance(corpus: EventBCorpus, issues: list[ReviewIssue]) -> None:
    available = set(corpus.provenance["provenance_id"].astype(str))
    for entity, frame in corpus.tables().items():
        if entity == "provenance" or "provenance_id" not in frame:
            continue
        for entity_id, provenance_id in zip(frame.iloc[:, 0], frame.provenance_id, strict=True):
            if str(provenance_id) not in available:
                _issue(issues, entity, entity_id, "MISSING_PROVENANCE", str(provenance_id))
    for _, row in corpus.provenance.iterrows():
        if _enum(row.value_origin, ValueOrigin) is None:
            _issue(issues, "provenance", row.provenance_id, "VALUE_ORIGIN", str(row.value_origin))


def validate_corpus(corpus: EventBCorpus) -> ValidationResult:
    normalized = corpus.normalized()
    issues: list[ReviewIssue] = []
    _validate_links(normalized, issues)
    _validate_candidates(normalized, issues)
    _validate_assays(normalized, issues)
    _validate_provenance(normalized, issues)
    try:
        normalized.recognition_evidence = validate_evidence(normalized.recognition_evidence)
    except ValueError as error:
        _issue(issues, "recognition_evidence", "corpus", "EVIDENCE_SCHEMA", str(error))
    affected_assays = {issue.entity_id for issue in issues if issue.entity_type == "assays"}
    affected_candidates = {issue.entity_id for issue in issues if issue.entity_type == "candidates"}
    affected_assays.update(
        normalized.assays.loc[
            normalized.assays.candidate_id.astype(str).isin(affected_candidates), "assay_id"
        ].astype(str)
    )
    accepted = normalized.assays.loc[
        normalized.assays.review_status.astype(str).str.upper().eq(ReviewStatus.ACCEPTED.value)
        & ~normalized.assays.assay_id.astype(str).isin(affected_assays)
    ].copy()
    accepted_corpus = EventBCorpus(**normalized.tables())
    accepted_corpus.assays = accepted
    if any(
        issue.entity_type == "recognition_evidence" and issue.entity_id == "corpus"
        for issue in issues
    ):
        accepted_corpus.recognition_evidence = normalized.recognition_evidence.iloc[0:0].copy()
    return ValidationResult(
        normalized, accepted_corpus, tuple(sorted(issues, key=lambda x: x.issue_id))
    )
