import pandas as pd
import pytest

from event_b.adapters.mkras_vax import (
    ANTIGENS,
    COHORT_ID,
    MKRASVaxAdapter,
    STUDY_ID,
    TUMOR_MUTATIONS,
    reconcile,
)
from event_b.ingest import ingest_source
from event_b.manifest import SourceManifest
from event_b.models import BiologicalEvent, ResponseLabel, SCHEMA_VERSION


CALLS = {
    "J1994_12": [1, 1, 1, 1, 1, 1],
    "J1994_5": [1, 1, 1, 1, 1, 0],
    "J1994_10": [1, 0, 0, 0, 0, 0],
    "J1994_9": [1, 1, 1, 0, 0, 1],
    "J1994_1": [1, 1, 1, 1, 1, 1],
    "J1994_2": [1, 1, 1, 1, 1, 1],
    "J1994_14": [1, 1, 1, 1, 1, 1],
    "J1994_3": [1, 1, 1, 1, 1, 1],
    "J1994_13": [1, 1, 1, 0, 0, 0],
    "J1994_6": [1, 1, 1, 1, 0, 1],
    "J1994_7": [1, 1, 1, 1, 1, 1],
    "J1994_18": [1, 1, 1, 1, 1, 1],
}


def _manifest():
    declaration = MKRASVaxAdapter.declaration
    return SourceManifest(
        "manifest:test-mkras",
        declaration.source_name,
        declaration.source_version,
        declaration.adapter_name,
        declaration.adapter_version,
        (),
        SCHEMA_VERSION,
    )


def _calls():
    return pd.DataFrame(
        [
            {"patient": patient, **dict(zip(ANTIGENS, values, strict=True))}
            for patient, values in CALLS.items()
        ]
    )


def _result():
    adapter = MKRASVaxAdapter("unused")
    adapter.extract = lambda manifest: {"calls": _calls()}
    return ingest_source(adapter, _manifest())


def test_mkras_reconciles_exact_explicit_denominator():
    result = _result()
    assert result.review_queue == ()
    report = reconcile(result.accepted_corpus)
    assert report["reconciles"]
    assert report["extracted"] == {
        "patients": 12,
        "global_antigens": 6,
        "patient_candidates": 72,
        "assay_observations": 144,
        "event_a_observations": 72,
        "event_b_observations": 72,
        "positive_primary_labels": 60,
        "tested_negative_primary_labels": 12,
        "untested_candidates": 0,
        "review_queue": 0,
    }


def test_shared_antigen_identity_does_not_collapse_patients():
    corpus = _result().accepted_corpus
    assert corpus.antigens.antigen_id.nunique() == 6
    assert corpus.candidates.candidate_id.nunique() == 72
    assert corpus.candidates.groupby("mutant_peptide").patient_id.nunique().eq(12).all()
    links = corpus.entity_relationships
    derived = links[links.relationship_type.eq("DERIVED_FROM")]
    assert len(derived) == 72
    assert derived.source_entity_id.nunique() == 72
    assert derived.target_entity_id.nunique() == 6


def test_tumor_mutation_is_patient_specific_not_cohort_inferred():
    corpus = _result().accepted_corpus
    candidates = corpus.candidates.copy()
    candidates["source_patient"] = candidates.patient_id.str.rsplit(":", n=1).str[-1]
    confirmed = candidates[candidates.genomic_variant.notna()]
    assert len(confirmed) == 12
    assert {
        row.source_patient: str(row.genomic_variant).replace("KRAS ", "")
        for row in confirmed.itertuples()
    } == TUMOR_MUTATIONS


def test_baseline_and_post_vaccine_events_remain_separate():
    assays = _result().accepted_corpus.assays
    baseline = assays[assays.relative_to_vaccine.eq("PRE_VACCINE")]
    post = assays[assays.relative_to_vaccine.eq("POST_VACCINE")]
    assert baseline.event_type.eq(BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value).all()
    assert baseline.response_label.eq(ResponseLabel.TESTED_NEGATIVE.value).all()
    assert post.event_type.eq(BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value).all()
    assert post.response_label.value_counts().to_dict() == {"POSITIVE": 60, "TESTED_NEGATIVE": 12}


def test_patient_level_response_is_not_inflated_over_component_calls():
    post = _result().accepted_corpus.assays
    post = post[post.event_type.eq(BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value)]
    patient = f"{STUDY_ID}:J1994_10"
    labels = post.loc[post.patient_id.eq(patient), "response_label"].value_counts().to_dict()
    assert labels == {"TESTED_NEGATIVE": 5, "POSITIVE": 1}


def test_invalid_or_missing_component_call_is_rejected():
    calls = _calls()
    calls.loc[0, "G12V"] = 2
    adapter = MKRASVaxAdapter("unused")
    with pytest.raises(ValueError, match="non-binary"):
        adapter.normalize({"calls": calls}, _manifest())


def test_declaration_locks_shared_structure_and_cohort():
    declaration = MKRASVaxAdapter.declaration
    assert declaration.canonical_study_id == STUDY_ID
    assert declaration.cohort_id == COHORT_ID
    assert declaration.candidate_identity_completeness == "PATIENT_AND_SHARED_ANTIGEN_RESOLVED"
