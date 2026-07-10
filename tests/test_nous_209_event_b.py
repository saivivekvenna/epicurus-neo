import pandas as pd

from event_b.adapters.nous_209 import COHORT_ID, Nous209Adapter, STUDY_ID, reconcile
from event_b.ingest import ingest_source
from event_b.manifest import SourceManifest
from event_b.models import ResponseLabel, SCHEMA_VERSION


def _manifest():
    declaration = Nous209Adapter.declaration
    return SourceManifest(
        "manifest:test-nous",
        declaration.source_name,
        declaration.source_version,
        declaration.adapter_name,
        declaration.adapter_version,
        (),
        SCHEMA_VERSION,
    )


def _fsp():
    rows = []
    pools = [f"Pool{index}" for index in range(1, 17) if index != 15]
    for index in range(115):
        rows.append(
            {
                "Immunogenic FSPS": f"FSP_GENE{index}_{index % 2 + 1}",
                "Peptide pools": pools[index % len(pools)],
            }
        )
    return pd.DataFrame(rows)


def _patients():
    return pd.DataFrame(
        {
            "patient": [f"Pt {index}" for index in range(1, 38)],
            "gender": ["F" if index % 2 else "M" for index in range(1, 38)],
            "reactive_pools": [index % 16 + 1 for index in range(37)],
        }
    )


def _result():
    adapter = Nous209Adapter("unused")
    adapter.extract = lambda manifest: {"fsp": _fsp(), "patients": _patients()}
    return ingest_source(adapter, _manifest())


def test_nous_reconciles_patient_level_event_b_without_candidate_labels():
    result = _result()
    assert result.review_queue == ()
    report = reconcile(result.accepted_corpus)
    assert report["reconciles"]
    assert report["extracted"] == {
        "patients": 37,
        "patient_level_event_b_observations": 37,
        "primary_candidate_labels": 0,
        "pool_entities": 16,
        "fsp_identifiers": 115,
        "review_queue": 0,
    }


def test_pool_response_is_never_inflated_to_fsp_labels():
    corpus = _result().accepted_corpus
    assert corpus.candidates.empty
    assert corpus.assays.candidate_id.isna().all()
    assert corpus.assays.response_label.eq(ResponseLabel.POSITIVE.value).all()
    assert corpus.assays.quantitative_result.between(1, 16).all()
    assert not corpus.entity_relationships.relationship_type.eq("TESTS_RESPONSE_TO").any()


def test_global_fsp_to_pool_structure_is_preserved_without_patient_mapping():
    corpus = _result().accepted_corpus
    relationships = corpus.entity_relationships
    assert len(relationships) == 115
    assert relationships.relationship_type.eq("CONTAINED_WITHIN").all()
    assert relationships.source_entity_id.nunique() == 115
    assert relationships.target_entity_id.nunique() == 15


def test_no_event_a_identity_is_fabricated_from_cohort_level_baseline_statement():
    corpus = _result().accepted_corpus
    assert corpus.assays.event_type.nunique() == 1
    assert corpus.assays.event_type.iloc[0] == "EVENT_B_VACCINE_INDUCED_RESPONSE"
    assert corpus.assays.relative_to_vaccine.eq("POST_BOOST").all()


def test_declaration_locks_patient_level_granularity():
    declaration = Nous209Adapter.declaration
    assert declaration.canonical_study_id == STUDY_ID
    assert declaration.cohort_id == COHORT_ID
    assert declaration.candidate_identity_completeness == "PATIENT_LEVEL_ONLY_NO_CANDIDATE_MAPPING"
