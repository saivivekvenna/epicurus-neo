import json

import pandas as pd

from event_b.corpus import EventBCorpus
from event_b.models import BiologicalEvent, ResponseLabel
from event_b.registry import REGISTRY_VERSION, StudyRegistry, StudyRegistryEntry, StudyStatus
from event_b.sufficiency import render_sufficiency_markdown, sufficiency_audit


def _entry(study_id, design):
    return StudyRegistryEntry(
        study_id,
        f"{study_id}_cohort",
        (f"DOI:10.test/{study_id}",),
        None,
        f"cancer_{study_id}",
        f"platform_{study_id}",
        design,
        "fixture",
        "IMPLEMENTED",
        StudyStatus.ACCEPTED,
        ("CSV",),
        (),
        (),
        None,
    )


def _fixture():
    registry = StudyRegistry(
        REGISTRY_VERSION,
        (_entry("s1", "SHARED"), _entry("s2", "PERSONALIZED"), _entry("s3", "PERSONALIZED")),
    )
    studies = []
    patients = []
    candidates = []
    assays = []
    for study_index, entry in enumerate(registry.studies, start=1):
        studies.append(
            {
                "study_id": entry.canonical_study_id,
                "publication_ids": json.dumps(entry.publication_ids),
                "vaccine_platform": entry.vaccine_platform,
            }
        )
        for patient_index in range(2):
            patient_id = f"{entry.canonical_study_id}:p{patient_index}"
            patients.append(
                {
                    "patient_id": patient_id,
                    "study_id": entry.canonical_study_id,
                    "cancer_type": entry.cancer_type,
                }
            )
            for candidate_index in range(2):
                candidate_id = f"c:{entry.canonical_study_id}:{patient_index}:{candidate_index}"
                if entry.antigen_design == "SHARED":
                    peptide = f"SHAREDPEPTIDE{candidate_index}"
                else:
                    peptide = f"PEPTIDE{study_index}{patient_index}{candidate_index}"
                candidates.append(
                    {
                        "candidate_id": candidate_id,
                        "patient_id": patient_id,
                        "study_id": entry.canonical_study_id,
                        "mutant_peptide": peptide,
                        "hla_alleles": json.dumps([f"HLA-A*0{study_index}:01"]),
                        "mhc_class": "CLASS_I" if study_index != 2 else "CLASS_II",
                        "sample_date": f"202{study_index}-0{patient_index + 1}-01",
                    }
                )
                label = (
                    ResponseLabel.POSITIVE.value
                    if patient_index == candidate_index
                    else ResponseLabel.TESTED_NEGATIVE.value
                )
                assays.append(
                    {
                        "assay_id": f"a:{candidate_id}",
                        "candidate_id": candidate_id,
                        "patient_id": patient_id,
                        "study_id": entry.canonical_study_id,
                        "event_type": BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                        "response_label": label,
                        "assay_type": "ELISPOT",
                        "source_interpretation": "candidate-level",
                    }
                )
    corpus = EventBCorpus(
        studies=pd.DataFrame(studies),
        patients=pd.DataFrame(patients),
        candidates=pd.DataFrame(candidates),
        assays=pd.DataFrame(assays),
    )
    return corpus, registry


def test_sufficiency_keeps_patient_study_and_label_counts_distinct():
    corpus, registry = _fixture()
    audit = sufficiency_audit(corpus, registry)
    counts = audit["required_counts"]
    assert counts["event_b_patients"] == 6
    assert counts["event_b_studies"] == 3
    # With no patient-level-only evidence, the candidate-resolved tier equals the headline.
    assert counts["candidate_resolved_patient_n"] == 6
    assert counts["candidate_resolved_study_n"] == 3
    assert counts["patient_level_only_patient_n"] == 0
    assert counts["primary_candidate_labels"] == 12
    assert counts["assay_observations"] == 12
    assert counts["event_b_positives"] == 6
    assert counts["event_b_tested_negatives"] == 6
    assert counts["class_i_observations"] == 8
    assert counts["class_ii_observations"] == 4


def test_registered_minimum_and_verdict_are_conservative():
    corpus, registry = _fixture()
    audit = sufficiency_audit(corpus, registry)
    assert not audit["registered_minimum"]["met"]
    assert audit["verdict"] == "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"
    assert audit["label_quality"]["inferred_negatives"] == 0


def test_split_feasibility_requires_both_labels_on_both_sides():
    corpus, registry = _fixture()
    splits = sufficiency_audit(corpus, registry)["split_feasibility"]
    assert splits["PATIENT_HOLDOUT"]["feasible"]
    assert splits["STUDY_HOLDOUT"]["feasible"]
    assert splits["HLA_HOLDOUT"]["feasible"]
    assert splits["PEPTIDE_CLUSTER_HOLDOUT"]["feasible"]
    assert splits["CANCER_TYPE_HOLDOUT"]["feasible"]
    assert splits["SHARED_ANTIGEN_GROUP_HOLDOUT"]["feasible"]
    for result in splits.values():
        if not result["feasible"]:
            continue
        assert result["label_counts"]["evaluation"]["POSITIVE"] > 0
        assert result["label_counts"]["evaluation"]["TESTED_NEGATIVE"] > 0


def test_sufficiency_markdown_contains_verdict_and_study_status():
    corpus, registry = _fixture()
    markdown = render_sufficiency_markdown(sufficiency_audit(corpus, registry))
    assert "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA" in markdown
    assert "Evidence tiers" in markdown
    assert "candidate_resolved_patient_n" in markdown
    assert "`s1`" in markdown
    assert "Split feasibility" in markdown


def test_patient_level_only_evidence_never_satisfies_candidate_resolved_gate():
    """A patient-level-only cohort must lift the headline count but not the peptide-ranking tier."""
    corpus, registry = _fixture()
    extra_patient = pd.DataFrame(
        [{"patient_id": "s1:pl", "study_id": "s1", "cancer_type": "cancer_s1"}]
    )
    extra_assay = pd.DataFrame(
        [
            {
                "assay_id": "a:s1:pl",
                "candidate_id": pd.NA,
                "patient_id": "s1:pl",
                "study_id": "s1",
                "event_type": BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                "response_label": ResponseLabel.POSITIVE.value,
                "assay_type": "ELISPOT",
                "source_interpretation": "patient-level only",
            }
        ]
    )
    corpus = EventBCorpus(
        studies=corpus.studies,
        patients=pd.concat([corpus.patients, extra_patient], ignore_index=True),
        candidates=corpus.candidates,
        assays=pd.concat([corpus.assays, extra_assay], ignore_index=True),
    )
    counts = sufficiency_audit(corpus, registry)["required_counts"]
    # Headline superset absorbs the patient-level-only positive ...
    assert counts["event_b_patients"] == 7
    assert counts["event_b_positive_patients"] == 7
    # ... but the candidate-resolved tier that the gate depends on is untouched.
    assert counts["candidate_resolved_patient_n"] == 6
    assert counts["candidate_resolved_positive_patient_n"] == 6
    assert counts["patient_level_only_patient_n"] == 1
