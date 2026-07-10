import pandas as pd

from event_b.adapters.pdac_neovax import (
    COHORT_ID,
    PDACNeoVaxAdapter,
    SOURCE_COLUMNS,
    STUDY_ID,
    reconcile,
)
from event_b.ingest import ingest_source
from event_b.manifest import SourceManifest
from event_b.models import BiologicalEvent, ResponseLabel, SCHEMA_VERSION


PATIENT_COUNTS = {
    1: 7,
    3: 20,
    4: 19,
    5: 16,
    6: 10,
    9: 17,
    10: 14,
    11: 14,
    14: 10,
    18: 10,
    19: 20,
    20: 19,
    23: 10,
    25: 17,
    28: 9,
    29: 20,
}


def _manifest():
    declaration = PDACNeoVaxAdapter.declaration
    return SourceManifest(
        "manifest:test-pdac",
        declaration.source_name,
        declaration.source_version,
        declaration.adapter_name,
        declaration.adapter_version,
        (),
        SCHEMA_VERSION,
    )


def _targets():
    rows = []
    index = 0
    for patient, count in PATIENT_COUNTS.items():
        for neo in range(1, count + 1):
            if patient == 25 and neo <= 7:
                response = "De novo response in pool"
            elif index < 23:
                response = "De novo response"
            elif 30 <= index < 32:
                response = "No data"
            else:
                response = "No response"
            mutant = "ACDEFGHIKLMNPQ" + ("R" if index % 2 else "S")
            wildtype = "ACDEFGHIKLMNPQ" + ("T" if index % 2 else "V")
            rows.append(
                {
                    "Patient number": patient,
                    "Neoantigen number": neo,
                    "Gene": f"GENE{index}",
                    "RefSeq transcript": f"NM_{index:06d}",
                    "Substitution": f"A{index + 1}V",
                    "Mutant Neoantigen Sequence": mutant,
                    "WT Neoantigen Sequence": wildtype,
                    "mRNA (+-13 AA (SNV); -15 AA to STOP (Indels))": "ACGT",
                    "MHC-I  Allele (Best Prediction)": "HLA-A*02:01",
                    "MHC-I Mutant Epitope (Best Prediction)": "ACDEFGHIK",
                    "MHC-I WT Epitope": "ACDEFGHIL",
                    "MHC-II  Allele (Best Prediction)": "HLA-DRB1*03:01",
                    "MHC-II Mutant Epitope (Best Prediction)": "ACDEFGHIKLMNPQ",
                    "MHC-II WT Epitope": "ACDEFGHILMNPQ",
                    "ELISpot Response": response,
                }
            )
            index += 1
    frame = pd.DataFrame(rows)
    assert set(frame.columns) == SOURCE_COLUMNS
    assert len(frame) == 232
    return frame


def _result():
    adapter = PDACNeoVaxAdapter("unused")
    adapter.extract = lambda manifest: {"targets": _targets()}
    return ingest_source(adapter, _manifest())


def test_pdac_reconciliation_preserves_three_state_candidate_labels():
    result = _result()
    assert result.review_queue == ()
    report = reconcile(result.accepted_corpus)
    assert report["reconciles"]
    assert report["extracted"] == {
        "patients": 16,
        "vaccine_targets": 232,
        "primary_candidate_labels": 232,
        "positive": 23,
        "tested_negative": 200,
        "untested": 9,
        "pool_level_observations": 1,
        "assay_observations": 465,
        "review_queue": 0,
    }


def test_positive_pool_is_not_decomposed_to_candidate_positives():
    assays = _result().accepted_corpus.assays
    post = assays[
        assays.event_type.eq(BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value)
    ]
    patient25 = post[post.patient_id.eq(f"{STUDY_ID}:patient_25")]
    candidate_rows = patient25[patient25.candidate_id.notna()]
    assert candidate_rows.response_label.eq(ResponseLabel.UNTESTED.value).sum() == 7
    pool_rows = patient25[patient25.candidate_id.isna()]
    assert len(pool_rows) == 1
    assert pool_rows.iloc[0].response_label == ResponseLabel.POSITIVE.value
    assert pool_rows.iloc[0].quantitative_result == 2


def test_encoded_assay_and_predicted_minimal_entities_are_distinct():
    corpus = _result().accepted_corpus
    component_counts = corpus.antigens.component_type.value_counts().to_dict()
    assert component_counts == {
        "ENCODED_MRNA_NEOANTIGEN": 232,
        "OVERLAPPING_15MER_ASSAY_POOL_COVERAGE": 232,
        "PREDICTED_MINIMAL_EPITOPE_CLASS_I": 232,
        "PREDICTED_MINIMAL_EPITOPE_CLASS_II": 232,
    }
    predicted = corpus.antigens[
        corpus.antigens.component_type.str.startswith("PREDICTED_MINIMAL")
    ]
    assert predicted.hla_evidence_type.eq("PREDICTED_BEST_BINDER").all()
    assert not predicted.hla_evidence_type.str.contains("EXPERIMENTAL").any()


def test_followup_publication_is_same_cohort_not_new_patients():
    corpus = _result().accepted_corpus
    assert corpus.studies.study_id.tolist() == [STUDY_ID]
    publication_ids = corpus.studies.iloc[0].publication_ids
    assert "10.1038/s41586-024-08508-4" in publication_ids
    assert corpus.patients.patient_id.nunique() == 16


def test_declaration_locks_cohort_and_candidate_granularity():
    declaration = PDACNeoVaxAdapter.declaration
    assert declaration.canonical_study_id == STUDY_ID
    assert declaration.cohort_id == COHORT_ID
    assert declaration.candidate_identity_completeness == "PATIENT_AND_ENCODED_NEOANTIGEN_RESOLVED"
