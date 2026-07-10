"""Hu 2021 melanoma NeoVax Event-B vertical-slice tests.

Core logic runs against synthetic in-memory sheets (no network, no 2.2 GB file) so it
runs in the standard suite; the streaming reader is exercised on a tiny real xlsx built
with openpyxl. The real-source reconciliation, determinism, and multi-study combination
tests skip when the gitignored supplement is absent; the milestone-5b1 CI job supplies it.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from event_b.adapters.hu_neovax import (
    EXPECTED_SHA256,
    HuNeoVaxAdapter,
    hu_source_path,
    read_workbook_sheets,
)
from event_b.audit import corpus_audit
from event_b.hu_pipeline import (
    assert_quality_gates,
    build_hu_corpus,
    reconcile_hu,
)
from event_b.manifest import manifest_from_paths
from event_b.models import BiologicalEvent, MHCClass, ResponseLabel, VaccineInclusion
from event_b.validation import validate_corpus


RAW_DIR = Path("data/raw/hu_melanoma_2021")
HAS_SOURCE = hu_source_path(RAW_DIR).exists()
skip_no_source = pytest.mark.skipif(
    not HAS_SOURCE, reason="Hu supplement absent (gitignored); supplied in the milestone-5b1 CI job"
)


# --------------------------------------------------------------------------------------
# Synthetic sheets that mimic the real Dataset 4a/4b/11a header geometry
# --------------------------------------------------------------------------------------
def _cd8_sheet() -> list[list[str]]:
    r16 = "Peptide pulsed autologous APC (16 weeks)"
    ryr = "Peptide pulsed autologous APC (3-4.5 years)"
    return [
        ["Supplementary Dataset 4a."],
        ["Patient ID", "Immunizing pool", "Gene", "Protein change", "Peptide length",
         "HLA allele", "Mutated peptide", "", "Wild type peptide", "", "Immunizing peptide",
         "", "EPT peptide ID", "Gene expression (TPM)", "CD8+ T cell reactivity by IFN-g ELISPOT", ""],
        ["", "", "", "", "", "", "Sequence", "Affinity (nM)", "Sequence", "Affinity  (nM)",
         "Sequence", "ID", "", "", r16, ryr],
        ["1", "A", "ACPP", "p.E34K", "9", "A24:02", "KLKFVTLVF", "142", "ELKFVTLVF", "1140",
         "DRSVLAKKLKFVTLVFRHGDRSPID", "1-IMP04", "1-EPT4A", "1.3", "1", "0"],
        ["1", "A", "PRTG", "p.F1055L", "9", "A02:01", "FLFQDSKKI", "83", "FFFQDSKKI", "16441",
         "NNSKKKWFLFQDSKKIQVEQPQ", "1-IMP03", "1-EPT3A", "0.6", "0", "n.d."],
        ["2", "B", "ARHGEF", "p.V651A", "9", "B15:01", "ALFASRPRF", "63", "VLFASRPRF", "76",
         "RRGGALFASRPRFTPL", "2-IMP12", "2-EPT12A", "2.8", "n.d.", "n.d."],
    ]


def _cd4_sheet() -> list[list[str]]:
    return [
        ["Supplementary Dataset 4b."],
        ["", "", "", "Predicted class I peptide", "", "Immunizing peptide", "", "Assay peptide",
         "", "CD4+ T cell reactivity by IFN-g ELISPOT", "", "", "", "", ""],
        ["Patient ID", "Gene", "Protein change", "Mutant", "WT", "ID", "Sequence", "ID",
         "Sequence", "Week 16", "", "", "", "", "Year 3-4.5", ""],
        ["", "", "", "", "", "", "", "", "", "Ex vivo peptide pulsed",
         "After pre-stimulation Peptide pulsed", "Minigene", "Minigene blocked", "Autologus tumor",
         "Ex vivo peptide pulsed", "After pre-stimulation"],
        # ex vivo positive -> EXVIVO reliability
        ["1", "ZBED4", "p.S218F", "VQKVASKIPF", "VQKVASKIPS", "1-IMP02", "SPIKLVQKVASKIPFPDRITEESV",
         "1-ASP4", "SPIKLVQKVASKIPF", "1", "1", "n.d.", "n.d.", "n.d.", "1", "1"],
        # pre-stim only positive -> PRESTIM reliability
        ["1", "PRTG", "p.F1055L", "FLFQDSKKI", "FFFQDSKKI", "1-IMP03", "NNSKKKWFLFQDSKKIQVEQPQ",
         "1-ASP8", "KKWFLFQDSKKIQVEQ", "0", "1", "n.d.", "n.d.", "n.d.", "0", "0"],
        # negative
        ["2", "GENEX", "p.A1T", "AAAAAAAAA", "TAAAAAAAA", "2-IMP01", "AAAAAAAAAAAAAAAAAAAA",
         "2-ASP1", "AAAAAAAAAAAAAAA", "0", "0", "n.d.", "n.d.", "n.d.", "n.d.", "n.d."],
    ]


def _spreading_sheet() -> list[list[str]]:
    r16 = "Peptide pulsed autologous APC (16 weeks)"
    return [
        ["Supplementary Dataset 11a."],
        ["Patient ID", "Gene", "Protein change", "Peptide length", "HLA allele",
         "Mutated peptide", "", "Wild type peptide", "", "CD8+ T cell reactivity by IFN-g ELISPOT", ""],
        ["", "", "", "", "", "Sequence", "Affinity (nM)", "Sequence", "Affinity  (nM)", r16,
         "Peptide pulsed autologous APC (post-pembro)"],
        # a spreading POSITIVE: must never become an Event-B vaccine label
        ["2", "ASCC3", "p.S1373L", "10", "A03:01", "RVFNKYPTLK", "5", "RVFNKYPTSK", "7", "1", "0"],
        ["6", "SCAMP4", "p.S159F", "10", "B56:01", "LPAIMFFVSA", "5", "LPAIMFSVSA", "5", "0", "0"],
    ]


def _normalize(tmp_path: Path):
    source = tmp_path / "source.txt"
    source.write_text("synthetic")
    manifest = manifest_from_paths("Hu", "test", "hu_neovax_event_b", "1.0.0", [source])
    adapter = HuNeoVaxAdapter(tmp_path)
    extracted = {"cd8": _cd8_sheet(), "cd4": _cd4_sheet(), "spreading": {41: _spreading_sheet()}}
    corpus = adapter.normalize(extracted, manifest)
    return adapter, corpus


def test_streaming_reader_reads_targeted_sheets(tmp_path):
    openpyxl = pytest.importorskip("openpyxl")
    path = tmp_path / "mini.xlsx"
    workbook = openpyxl.Workbook()
    first = workbook.active
    first.title = "first"
    first.append(["Patient ID", "HLA allele", "Mutated peptide"])
    first.append(["1", "A02:01", "KLKFVTLVF"])
    second = workbook.create_sheet("second")
    second.append(["only", "sheet", "two"])
    workbook.save(path)

    sheets = read_workbook_sheets(path, (1, 2))
    assert sheets[1][0] == ["Patient ID", "HLA allele", "Mutated peptide"]
    assert sheets[1][1] == ["1", "A02:01", "KLKFVTLVF"]
    assert sheets[2][0] == ["only", "sheet", "two"]


def test_normalize_three_state_labels_and_mhc_class(tmp_path):
    adapter, corpus = _normalize(tmp_path)
    accepted = validate_corpus(corpus).accepted_corpus

    # Distinct patients came through (1-2 from vaccine sheets, 6 from the spreading sheet);
    # no validator issues on the clean synthetic rows.
    assert set(corpus.patients.source_patient_id.astype(str)) == {"1", "2", "6"}
    assert validate_corpus(corpus).review_queue == ()

    # CD8 minimal epitopes are class I; CD4 assay peptides are class II.
    by_source = corpus.candidates.set_index("candidate_id")
    cd8 = corpus.candidates[corpus.candidates.candidate_source.str.contains("4a")]
    cd4 = corpus.candidates[corpus.candidates.candidate_source.str.contains("4b")]
    assert set(cd8.mhc_class) == {MHCClass.CLASS_I.value}
    assert set(cd4.mhc_class) == {MHCClass.CLASS_II.value}

    # Three-state labels present; the n.d. CD8 rows are UNTESTED and held for review.
    labels = accepted.assays.response_label.astype(str).value_counts().to_dict()
    assert labels.get(ResponseLabel.POSITIVE.value, 0) >= 2
    assert labels.get(ResponseLabel.TESTED_NEGATIVE.value, 0) >= 2
    assert ResponseLabel.UNTESTED.value not in accepted.assays.response_label.astype(str).tolist()
    assert any(issue.code == "UNSCORABLE_ASSAY" for issue in adapter.review_issues)
    del by_source


def test_epitope_spreading_is_separated_and_never_vaccine(tmp_path):
    _, corpus = _normalize(tmp_path)
    accepted = validate_corpus(corpus).accepted_corpus

    spreading = accepted.assays[
        accepted.assays.event_type.astype(str).eq(BiologicalEvent.EPITOPE_SPREADING.value)
    ]
    assert len(spreading) == 2  # one positive, one negative
    # The spreading POSITIVE must NOT be counted as an Event-B vaccine label.
    spreading_candidates = corpus.candidates[
        corpus.candidates.candidate_id.isin(spreading.candidate_id)
    ]
    assert set(spreading_candidates.vaccine_inclusion) == {VaccineInclusion.NOT_INCLUDED.value}

    audit = corpus_audit(accepted)
    # Event-B patient count excludes the spreading-only signal.
    assert BiologicalEvent.EPITOPE_SPREADING.value in audit["event_counts"]
    event_b_positive_patients = audit["model_readiness"]["event_b_positive_patient_n"]
    # Patient 6 appears only via a negative spreading row -> not an Event-B positive patient.
    assert "6" not in {
        pid.split(":")[-1]
        for pid in accepted.assays.loc[
            accepted.assays.event_type.astype(str).eq(
                BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
            )
            & accepted.assays.response_label.astype(str).eq(ResponseLabel.POSITIVE.value),
            "patient_id",
        ].astype(str)
    }
    assert event_b_positive_patients >= 1
    assert_quality_gates(accepted)


def test_reliability_is_not_flattened(tmp_path):
    _, corpus = _normalize(tmp_path)
    evidence = corpus.recognition_evidence.merge(
        corpus.candidates[["candidate_id", "candidate_source"]], on="candidate_id", how="left"
    )
    # Compare like-for-like: the week-16 recognition channel (VACCINE_EVENT_B), not the
    # separate LONGITUDINAL_PERSISTENCE evidence which carries its own reliability vector.
    vaccine = evidence[evidence.evidence_family.eq("VACCINE_EVENT_B")]
    cd8 = vaccine[vaccine.candidate_source.str.contains("4a", na=False)]
    cd4 = vaccine[vaccine.candidate_source.str.contains("4b", na=False)]
    spread = evidence[evidence.evidence_family.eq("FUNCTIONAL_T_CELL_ASSAY")]

    # CD8 minimal epitope is more candidate-specific than a CD4 overlapping peptide.
    assert set(cd8.candidate_specificity) == {1.0}
    assert set(cd4.candidate_specificity) == {0.6}
    # Ex-vivo and pre-stimulation CD4 keep distinct assay directness (not flattened).
    assert set(cd4.assay_directness) == {0.9, 0.6}
    # Epitope spreading carries zero vaccine relevance.
    assert set(spread.vaccine_relevance) == {0.0}


def test_all_vaccine_recognition_is_post_vaccine_event_b(tmp_path):
    _, corpus = _normalize(tmp_path)
    vaccine = corpus.assays[
        corpus.assays.event_type.astype(str).eq(
            BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
        )
    ]
    assert (vaccine.relative_to_vaccine.astype(str).str.upper() == "POST_VACCINE").all()
    # IMPROVE's Event-A is never introduced by this adapter.
    assert not corpus.assays.event_type.astype(str).eq(
        BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value
    ).any()


def test_quality_gates_catch_spreading_labelled_as_vaccine(tmp_path):
    _, corpus = _normalize(tmp_path)
    accepted = validate_corpus(corpus).accepted_corpus
    assert_quality_gates(accepted)  # clean slice passes

    tampered = accepted.candidates.copy()
    spreading_ids = accepted.assays.loc[
        accepted.assays.event_type.astype(str).eq(BiologicalEvent.EPITOPE_SPREADING.value),
        "candidate_id",
    ]
    tampered.loc[
        tampered.candidate_id.isin(spreading_ids), "vaccine_inclusion"
    ] = VaccineInclusion.INCLUDED.value
    broken = validate_corpus(corpus).accepted_corpus
    broken.candidates = tampered
    with pytest.raises(AssertionError, match="epitope-spreading"):
        assert_quality_gates(broken)


def test_multi_study_verdict_from_two_event_b_studies(tmp_path):
    """A two-Event-B-study corpus (patients still sub-threshold) yields the multi-study verdict."""
    _, corpus = _normalize(tmp_path)
    accepted = validate_corpus(corpus).accepted_corpus
    assays = accepted.assays.copy()
    # Clone the accepted assays under a second study id to simulate a second Event-B study.
    second = assays.copy()
    second["study_id"] = "second_study"
    second["patient_id"] = second["patient_id"].astype(str) + ":s2"
    second["assay_id"] = second["assay_id"].astype(str) + ":s2"
    accepted.assays = pd.concat([assays, second], ignore_index=True)

    audit = corpus_audit(accepted)
    assert audit["model_readiness"]["event_b_study_n"] == 2
    assert (
        audit["model_readiness"]["decision"]
        == "EVENT_B_MULTI_STUDY_CORPUS_VALIDATED_INSUFFICIENT_PATIENTS_FOR_GENERAL_MODEL"
    )


@skip_no_source
def test_real_source_reconciles_to_ott(tmp_path):
    build = build_hu_corpus(RAW_DIR)
    assert_quality_gates(build.result.accepted_corpus)
    recon = reconcile_hu(build)

    assert recon["source_observed"]["patient_n"] == 8
    assert set(recon["source_observed"]["patients"]) == {"1", "2", "3", "4", "5", "6", "11", "12"}
    cd8 = recon["reconciliation"]["cd8"]
    cd4 = recon["reconciliation"]["cd4"]
    assert cd8["observed"]["positive_neoantigens"] == 15  # matches Ott 2017 exactly
    assert cd8["reconciles"] is True
    assert cd4["reconciles_within_tolerance"] is True
    assert recon["epitope_spreading"]["all_non_vaccine"] is True
    assert (
        recon["model_readiness"]["decision"]
        == "EVENT_B_VERTICAL_SLICE_VALIDATED_NOT_YET_SUFFICIENT_FOR_GENERAL_MODEL"
    )


@skip_no_source
def test_real_source_checksum_pinned_and_deterministic():
    path = hu_source_path(RAW_DIR)
    from event_b.adapters.hu_neovax import sha256_file

    assert sha256_file(path) == EXPECTED_SHA256[path.name]

    build_a = build_hu_corpus(RAW_DIR)
    build_b = build_hu_corpus(RAW_DIR)
    key = ["assay_id", "response_label", "event_type", "provenance_id"]
    left = build_a.result.accepted_corpus.assays.sort_values("assay_id")[key].reset_index(drop=True)
    right = build_b.result.accepted_corpus.assays.sort_values("assay_id")[key].reset_index(drop=True)
    pd.testing.assert_frame_equal(left, right)


@skip_no_source
def test_combined_three_study_audit_multi_study_verdict():
    from event_b.braun_pipeline import build_braun_corpus, combine_corpora, load_corpus_from_parquet

    braun_raw = Path("data/raw/braun_rcc_2025")
    if not (braun_raw / "extracted").exists():
        pytest.skip("Braun source absent; combined multi-study audit needs both Event-B studies")
    hu = build_hu_corpus(RAW_DIR)
    braun = build_braun_corpus(braun_raw)
    corpora = [braun.result.accepted_corpus, hu.result.accepted_corpus]
    improve_dir = Path("outputs/event_b_corpus")
    if (improve_dir / "assays.parquet").exists():
        corpora.insert(0, load_corpus_from_parquet(improve_dir))
    combined = combine_corpora(*corpora)
    audit = corpus_audit(combined, hu.review_queue)
    assert audit["event_counts"][BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value] == 670
    assert audit["event_counts"][BiologicalEvent.EPITOPE_SPREADING.value] == 82
    assert audit["model_readiness"]["event_b_study_n"] == 2
    assert audit["model_readiness"]["event_b_patient_n"] == 17
    assert (
        audit["model_readiness"]["decision"]
        == "EVENT_B_MULTI_STUDY_CORPUS_VALIDATED_INSUFFICIENT_PATIENTS_FOR_GENERAL_MODEL"
    )
