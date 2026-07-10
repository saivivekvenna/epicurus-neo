import json
from pathlib import Path

import pandas as pd
import pytest

from event_b.adapters import GenericTableAdapter, ImproveEventAAdapter, OsteosarcCaseStudyAdapter
from event_b.audit import corpus_audit, render_audit_markdown
from event_b.corpus import EventBCorpus
from event_b.evidence import evidence_availability_matrix, validate_evidence
from event_b.export import PARQUET_TABLES, export_corpus
from event_b.extraction import (
    ExtractionCache,
    ExtractionTask,
    emit_extraction_tasks,
    run_extraction,
)
from event_b.funnel import link_event_b_to_funnel
from event_b.manifest import manifest_from_paths
from event_b.models import (
    AssayType,
    AvailabilityStatus,
    BiologicalEvent,
    EvidenceFamily,
    InformationTiming,
    MHCClass,
    ResponseLabel,
    ReviewStatus,
    ValueOrigin,
    VaccineInclusion,
    stable_candidate_id,
)
from event_b.schema import entity_json_schema
from event_b.splits import SplitType, generate_split_manifest, peptide_cluster_ids
from event_b.timing import assert_preselection_columns, preselection_evidence
from event_b.validation import validate_corpus


def _provenance(entity, entity_id, provenance_id):
    return {
        "provenance_id": provenance_id,
        "entity_type": entity,
        "entity_id": entity_id,
        "field_name": "*",
        "source_document": "paper.pdf",
        "page": 4,
        "source_fragment": f"reported {entity_id}",
        "extraction_method": "manual_fixture",
        "extraction_confidence": 1.0,
        "value_origin": ValueOrigin.SOURCE_REPORTED.value,
        "review_status": ReviewStatus.ACCEPTED.value,
    }


def _corpus() -> EventBCorpus:
    study = {
        "study_id": "s1",
        "title": "Vaccine trial",
        "publication_ids": "doi:1",
        "cancer_type": "melanoma",
        "source_manifest_id": "m1",
        "provenance_id": "prov-study",
    }
    patient = {
        "patient_id": "p1",
        "source_patient_id": "P-01",
        "study_id": "s1",
        "cancer_type": "melanoma",
        "hla_alleles": json.dumps(["HLA-A*02:01", "HLA-B*07:02"]),
        "provenance_id": "prov-patient",
    }
    vaccine = {
        "vaccine_id": "v1",
        "patient_id": "p1",
        "study_id": "s1",
        "vaccine_platform": "peptide",
        "vaccination_dates": json.dumps(["2025-01-10"]),
        "candidate_count": 2,
        "mhc_class_intent": MHCClass.CLASS_I.value,
        "provenance_id": "prov-vaccine",
    }
    candidates = []
    for index, peptide in enumerate(["SLYNTVATL", "GILGFVFTL"], start=1):
        candidates.append(
            {
                "candidate_id": f"c{index}",
                "patient_id": "p1",
                "study_id": "s1",
                "sample_id": "tumor-0",
                "sample_date": "2024-12-01",
                "timepoint": "PRE_VACCINE",
                "genomic_variant": f"chr1:{index}:A:T",
                "gene": f"GENE{index}",
                "transcript": f"TX{index}",
                "protein_change": f"p.A{index}T",
                "mutant_peptide": peptide,
                "wildtype_peptide": "A" + peptide[1:],
                "peptide_length": len(peptide),
                "hla_alleles": json.dumps(["HLA-A*02:01"]),
                "mhc_class": MHCClass.CLASS_I.value,
                "candidate_source": "source table",
                "vaccine_inclusion": VaccineInclusion.INCLUDED.value,
                "vaccine_inclusion_origin": ValueOrigin.SOURCE_REPORTED.value,
                "generation_provenance": "source supplement",
                "mutant_wildtype_verified": True,
                "provenance_id": f"prov-c{index}",
            }
        )
    assays = [
        {
            "assay_id": "a-pre",
            "patient_id": "p1",
            "study_id": "s1",
            "candidate_id": "c1",
            "vaccine_id": "v1",
            "assay_type": AssayType.TETRAMER.value,
            "sample_type": "PBMC",
            "sample_date": "2025-01-01",
            "timepoint": "PRE_VACCINE",
            "relative_to_vaccine": "PRE_VACCINE",
            "qualitative_result": "POSITIVE",
            "event_type": BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value,
            "response_label": ResponseLabel.POSITIVE.value,
            "explicit_assay_inclusion": True,
            "review_status": ReviewStatus.ACCEPTED.value,
            "provenance_id": "prov-a-pre",
        },
        {
            "assay_id": "a-post",
            "patient_id": "p1",
            "study_id": "s1",
            "candidate_id": "c1",
            "vaccine_id": "v1",
            "assay_type": AssayType.ELISPOT.value,
            "sample_type": "PBMC",
            "sample_date": "2025-02-01",
            "timepoint": "POST_BOOST",
            "relative_to_vaccine": "POST_BOOST",
            "quantitative_result": 150,
            "result_units": "SFC/1e6",
            "qualitative_result": "POSITIVE",
            "event_type": BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
            "response_label": ResponseLabel.POSITIVE.value,
            "explicit_assay_inclusion": True,
            "review_status": ReviewStatus.ACCEPTED.value,
            "provenance_id": "prov-a-post",
        },
        {
            "assay_id": "a-neg",
            "patient_id": "p1",
            "study_id": "s1",
            "candidate_id": "c2",
            "vaccine_id": "v1",
            "assay_type": AssayType.ELISPOT.value,
            "sample_type": "PBMC",
            "sample_date": "2025-02-01",
            "timepoint": "POST_BOOST",
            "relative_to_vaccine": "POST_BOOST",
            "quantitative_result": 0,
            "qualitative_result": "NEGATIVE",
            "event_type": BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
            "response_label": ResponseLabel.TESTED_NEGATIVE.value,
            "explicit_assay_inclusion": True,
            "review_status": ReviewStatus.ACCEPTED.value,
            "provenance_id": "prov-a-neg",
        },
    ]
    outcome = {
        "outcome_id": "o1",
        "patient_id": "p1",
        "study_id": "s1",
        "outcome_type": "RECIST",
        "outcome_value": "PR",
        "assessment_date": "2025-03-01",
        "provenance_id": "prov-outcome",
    }
    evidence = [
        {
            "evidence_id": "e-pre",
            "candidate_id": "c1",
            "patient_id": "p1",
            "evidence_family": EvidenceFamily.PRESENTATION_PREDICTION.value,
            "source_dataset": "predictor",
            "measured_or_predicted": "PREDICTED",
            "value": 0.8,
            "units": "probability",
            "directionality": "higher_better",
            "uncertainty": 0.1,
            "assay_or_model_version": "frozen-1",
            "evidence_quality": "DIRECT",
            "availability_status": AvailabilityStatus.AVAILABLE.value,
            "information_timing": InformationTiming.PRE_SELECTION.value,
            "patient_specificity": 1.0,
            "functional_relevance": 0.3,
            "vaccine_relevance": 0.2,
            "candidate_specificity": 1.0,
            "assay_directness": 0.2,
            "temporal_clarity": 1.0,
            "source_completeness": 1.0,
            "replication_status": "SINGLE",
            "provenance_id": "prov-e-pre",
        },
        {
            "evidence_id": "e-outcome",
            "candidate_id": "c1",
            "patient_id": "p1",
            "evidence_family": EvidenceFamily.TCR_EXPANSION.value,
            "source_dataset": "trial",
            "measured_or_predicted": "MEASURED",
            "value": 10,
            "units": "fold_change",
            "availability_status": AvailabilityStatus.AVAILABLE.value,
            "information_timing": InformationTiming.OUTCOME_ONLY.value,
            "provenance_id": "prov-e-outcome",
        },
    ]
    entities = [
        ("studies", "s1", "prov-study"),
        ("patients", "p1", "prov-patient"),
        ("vaccines", "v1", "prov-vaccine"),
        ("candidates", "c1", "prov-c1"),
        ("candidates", "c2", "prov-c2"),
        ("assays", "a-pre", "prov-a-pre"),
        ("assays", "a-post", "prov-a-post"),
        ("assays", "a-neg", "prov-a-neg"),
        ("clinical_outcomes", "o1", "prov-outcome"),
        ("recognition_evidence", "e-pre", "prov-e-pre"),
        ("recognition_evidence", "e-outcome", "prov-e-outcome"),
    ]
    return EventBCorpus(
        studies=pd.DataFrame([study]),
        patients=pd.DataFrame([patient]),
        vaccines=pd.DataFrame([vaccine]),
        candidates=pd.DataFrame(candidates),
        assays=pd.DataFrame(assays),
        clinical_outcomes=pd.DataFrame([outcome]),
        recognition_evidence=pd.DataFrame(evidence),
        provenance=pd.DataFrame([_provenance(*entity) for entity in entities]),
    )


def test_stable_candidate_identity_is_patient_and_hla_specific():
    base = {
        "study_id": "s",
        "patient_id": "p1",
        "mutant_peptide": "AAAAAAAAA",
        "hla_alleles": ["A02"],
    }
    first = stable_candidate_id(base)
    assert first == stable_candidate_id(base)
    assert first != stable_candidate_id(base | {"patient_id": "p2"})
    assert first != stable_candidate_id(base | {"hla_alleles": ["B07"]})
    schema = entity_json_schema("candidates")
    assert schema["properties"]["schema_version"]["const"] == "event-b-1.0.0"
    assert schema["additionalProperties"] is False


def test_valid_corpus_separates_events_and_clinical_outcomes():
    result = validate_corpus(_corpus())
    assert result.ok
    assert len(result.accepted_corpus.assays) == 3
    assert set(result.accepted_corpus.assays.event_type) == {
        BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value,
        BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
    }
    assert len(result.accepted_corpus.clinical_outcomes) == 1


@pytest.mark.parametrize(
    ("column", "value", "code"),
    [
        ("explicit_assay_inclusion", False, "NEGATIVE_DENOMINATOR"),
        ("relative_to_vaccine", "PRE_VACCINE", "EVENT_B_TIMEPOINT"),
    ],
)
def test_invalid_event_b_records_enter_review_queue(column, value, code):
    corpus = _corpus()
    corpus.assays.loc[corpus.assays.assay_id.eq("a-neg"), column] = value
    result = validate_corpus(corpus)
    assert code in {issue.code for issue in result.review_queue}
    assert "a-neg" not in set(result.accepted_corpus.assays.assay_id)


def test_untested_positive_result_and_clinical_as_assay_are_rejected():
    corpus = _corpus()
    corpus.assays.loc[0, "response_label"] = ResponseLabel.UNTESTED.value
    corpus.assays.loc[0, "event_type"] = BiologicalEvent.EVENT_C_CLINICAL_OUTCOME.value
    result = validate_corpus(corpus)
    codes = {issue.code for issue in result.review_queue}
    assert {"UNTESTED_RESULT", "CLINICAL_AS_ASSAY"}.issubset(codes)


def test_candidate_linkage_hla_length_and_vaccine_inclusion_are_reviewed():
    corpus = _corpus()
    corpus.candidates.loc[0, "study_id"] = "wrong-study"
    corpus.candidates.loc[0, "hla_alleles"] = json.dumps(["HLA-C*99:99"])
    corpus.candidates.loc[0, "peptide_length"] = 8
    corpus.candidates.loc[0, "vaccine_inclusion_origin"] = ValueOrigin.UNKNOWN.value
    result = validate_corpus(corpus)
    codes = {issue.code for issue in result.review_queue}
    assert {
        "UNKNOWN_STUDY",
        "STUDY_MISMATCH",
        "HLA_MISMATCH",
        "PEPTIDE_LENGTH",
        "VACCINE_INCLUSION",
    }.issubset(codes)
    assert "a-post" not in set(result.accepted_corpus.assays.assay_id)


def test_contradictory_accepted_labels_are_not_silently_resolved():
    corpus = _corpus()
    conflict = corpus.assays.loc[corpus.assays.assay_id.eq("a-post")].copy()
    conflict["assay_id"] = "a-conflict"
    conflict["response_label"] = ResponseLabel.TESTED_NEGATIVE.value
    conflict["provenance_id"] = "prov-a-conflict"
    corpus.assays = pd.concat([corpus.assays, conflict], ignore_index=True)
    corpus.provenance = pd.concat(
        [corpus.provenance, pd.DataFrame([_provenance("assays", "a-conflict", "prov-a-conflict")])],
        ignore_index=True,
    )
    result = validate_corpus(corpus)
    assert any(issue.code == "CONTRADICTORY_LABELS" for issue in result.review_queue)
    assert "a-conflict" not in set(result.accepted_corpus.assays.assay_id)


def test_evidence_channels_and_timing_remain_separate():
    evidence = validate_evidence(_corpus().recognition_evidence)
    assert set(evidence.evidence_family) == {
        EvidenceFamily.PRESENTATION_PREDICTION.value,
        EvidenceFamily.TCR_EXPANSION.value,
    }
    assert preselection_evidence(evidence).evidence_id.tolist() == ["e-pre"]
    matrix = evidence_availability_matrix(evidence)
    assert EvidenceFamily.PRESENTATION_PREDICTION.value in matrix
    with pytest.raises(ValueError, match="Outcome-only"):
        assert_preselection_columns(
            ["post_vaccine_tcr_expansion"],
            {"post_vaccine_tcr_expansion": InformationTiming.OUTCOME_ONLY.value},
        )


def test_llm_extraction_without_endpoint_emits_tasks_and_stays_pending(tmp_path):
    task = ExtractionTask.create(
        study_id="s1",
        source_document="paper.pdf",
        source_checksum="abc",
        source_text="Post-vaccine ELISPOT results table",
    )
    task_path, schema_path = emit_extraction_tasks([task], tmp_path / "tasks")
    assert task_path.exists() and schema_path.exists()
    result = run_extraction(task, None, ExtractionCache(tmp_path / "cache"))
    assert result.status == "PENDING"
    assert not list((tmp_path / "cache").glob("*.json"))


def test_llm_schema_forces_review_and_preserves_raw_output(tmp_path):
    task = ExtractionTask.create(
        study_id="s1",
        source_document="paper.pdf",
        source_checksum="abc",
        source_text="candidate-resolved post-vaccine assay",
    )

    class Provider:
        provider_id = "fake-v1"

        def infer(self, unused_task):
            del unused_task
            return {
                "assays": [
                    {
                        "assay_id": "a1",
                        "patient_id": "p1",
                        "study_id": "s1",
                        "candidate_id": "c1",
                        "assay_type": AssayType.ELISPOT.value,
                        "event_type": BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                        "response_label": ResponseLabel.POSITIVE.value,
                        "relative_to_vaccine": "POST_BOOST",
                        "explicit_assay_inclusion": True,
                        "review_status": ReviewStatus.NEEDS_REVIEW.value,
                        "provenance_id": "prov1",
                    }
                ],
                "provenance": [
                    {
                        "provenance_id": "prov1",
                        "entity_type": "assays",
                        "entity_id": "a1",
                        "field_name": "response_label",
                        "source_document": "paper.pdf",
                        "source_fragment": "positive after boost",
                        "extraction_confidence": 0.8,
                        "value_origin": ValueOrigin.LLM_EXTRACTED.value,
                        "review_status": ReviewStatus.NEEDS_REVIEW.value,
                    }
                ],
            }

    cache = ExtractionCache(tmp_path / "cache")
    result = run_extraction(task, Provider(), cache)
    assert result.status == "EXTRACTED_NEEDS_REVIEW"
    raw = json.loads(Path(result.cache_path).read_text())
    assert raw["assays"][0]["review_status"] == ReviewStatus.NEEDS_REVIEW.value

    class InvalidProvider(Provider):
        provider_id = "bad-v1"

        def infer(self, unused_task):
            payload = super().infer(unused_task)
            payload["assays"][0]["review_status"] = ReviewStatus.ACCEPTED.value
            return payload

    invalid = run_extraction(task, InvalidProvider(), cache, retries=0)
    assert invalid.status == "FAILED"


def test_adapter_declarations_are_explicit_and_osteosarc_does_not_fabricate(tmp_path):
    source = tmp_path / "studies.csv"
    pd.DataFrame(
        [{"study_id": "s", "title": "x", "source_manifest_id": "m", "provenance_id": "p"}]
    ).to_csv(source, index=False)
    adapter = GenericTableAdapter(
        source_name="toy", source_version="1", entity_paths={"studies": source}
    )
    assert adapter.declaration.supported_entities == ("studies",)
    manifest = manifest_from_paths("toy", "1", "generic_table", "1.0.0", [source])
    normalized = adapter.normalize(adapter.extract(manifest), manifest)
    assert normalized.studies.loc[0, "study_id"] == "s"
    assert ImproveEventAAdapter.declaration.supported_event_types == (
        BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value,
    )
    with pytest.raises(RuntimeError, match="not supplied"):
        OsteosarcCaseStudyAdapter().extract(None)


def _split_frame():
    return pd.DataFrame(
        {
            "candidate_id": [f"c{i}" for i in range(6)],
            "patient_id": [f"p{i}" for i in range(6)],
            "study_id": ["s1", "s1", "s2", "s2", "s3", "s3"],
            "mutant_peptide": [
                "AAAAAAAAA",
                "AAAAAAAAT",
                "CCCCCCCCC",
                "DDDDDDDDD",
                "EEEEEEEEE",
                "FFFFFFFFF",
            ],
            "hla_alleles": [f"HLA-{i}" for i in range(6)],
            "cancer_type": ["a", "a", "b", "b", "c", "c"],
            "sample_date": pd.date_range("2024-01-01", periods=6, freq="180D").astype(str),
        }
    )


def test_all_split_types_are_deterministic_and_patient_safe():
    frame = _split_frame()
    for split_type in SplitType:
        kwargs = (
            {"temporal_cutoff": "2025-01-01"} if split_type is SplitType.TEMPORAL_HOLDOUT else {}
        )
        first = generate_split_manifest(frame, split_type, **kwargs)
        second = generate_split_manifest(frame, split_type, **kwargs)
        assert first == second
        assignments = pd.DataFrame(first.assignments)
        assert assignments.groupby("patient_id").split.nunique().max() == 1
    clusters = peptide_cluster_ids(frame.mutant_peptide, threshold=0.8)
    assert clusters[0] == clusters[1]


def test_hla_holdout_connects_overlapping_allele_sets():
    frame = _split_frame()
    frame.loc[0, "hla_alleles"] = json.dumps(["HLA-A*02:01", "HLA-B*07:02"])
    frame.loc[1, "hla_alleles"] = json.dumps(["HLA-C*07:01", "HLA-A*02:01"])
    manifest = generate_split_manifest(frame, SplitType.HLA_HOLDOUT)
    assignments = pd.DataFrame(manifest.assignments).set_index("candidate_id")
    assert assignments.loc["c0", "split"] == assignments.loc["c1", "split"]


def test_event_b_links_to_funnel_without_inference():
    corpus = _corpus().normalized()
    stages = {
        "mutation_called": "reached",
        "transcript_represented": "reached",
        "peptide_generated": "reached",
        "survives_gating": "not_assessed",
        "hla_included": "not_assessed",
        "presentation_candidate": "not_assessed",
        "ranking_stage": "not_assessed",
        "top_k": "not_assessed",
    }
    ledger = pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "patient_id": "p1",
                "study_id": "s1",
                "provenance_id": "prov-c1",
                **stages,
            }
        ]
    )
    links = link_event_b_to_funnel(corpus, ledger)
    assert links.loc[0, "survives_gating"] == "not_assessed"
    assert links.loc[0, "vaccine_inclusion"] == "reached"
    assert links.loc[0, "functional_assay"] == "reached"


def test_audit_reports_peptide_patient_study_counts_and_insufficient_data():
    corpus = _corpus().normalized()
    audit = corpus_audit(corpus)
    assert audit["sample_sizes"]["peptide_n"] == 2
    assert audit["sample_sizes"]["patient_n"] == 1
    assert audit["sample_sizes"]["study_n"] == 1
    assert audit["patients_with_event_b_positive"] == 1
    assert audit["candidates_where_event_a_and_b_differ"] == 0
    # A nonzero-but-sub-threshold Event-B corpus resolves to the intermediate verdict, not the
    # zero-Event-B "insufficient" state, and is still not sufficient for a general model.
    assert (
        audit["model_readiness"]["decision"]
        == "EVENT_B_VERTICAL_SLICE_VALIDATED_NOT_YET_SUFFICIENT_FOR_GENERAL_MODEL"
    )
    assert audit["model_readiness"]["sufficient_for_recognition_model_development"] is False
    assert "not proof of clinical benefit" in render_audit_markdown(audit)


def test_deterministic_model_ready_exports_preserve_repeated_assays(tmp_path):
    corpus = _corpus().normalized()
    manifest_source = tmp_path / "source.txt"
    manifest_source.write_text("source")
    manifest = manifest_from_paths("toy", "1", "generic_table", "1.0.0", [manifest_source])
    split = generate_split_manifest(_split_frame(), SplitType.PATIENT_HOLDOUT)
    first = export_corpus(
        corpus, tmp_path / "first", source_manifests=[manifest], split_manifests=[split]
    )
    second = export_corpus(
        corpus, tmp_path / "second", source_manifests=[manifest], split_manifests=[split]
    )
    for entity, filename in PARQUET_TABLES.items():
        pd.testing.assert_frame_equal(
            pd.read_parquet(tmp_path / "first" / filename),
            pd.read_parquet(tmp_path / "second" / filename),
        )
        assert entity in first and entity in second
    model_ready = pd.read_parquet(tmp_path / "first/model_ready_recognition.parquet")
    assert len(model_ready) == len(corpus.assays)
    assert model_ready.assay_id.nunique() == len(corpus.assays)
    assert EvidenceFamily.PRESENTATION_PREDICTION.value in model_ready.columns
    assert (tmp_path / "first/corpus_audit.md").exists()
    assert (tmp_path / "first/schemas/assays.schema.json").exists()
