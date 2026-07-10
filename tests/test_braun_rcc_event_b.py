"""Braun RCC 2025 Event-B vertical-slice tests.

Core logic is exercised with synthetic in-memory sheets (no network, no large data,
so it runs in the standard suite). The real-source reconciliation and determinism
tests skip when the gitignored supplements are absent; the milestone-5a CI job
fetches them from Europe PMC and runs this file end-to-end.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from event_b.adapters.braun_rcc import (
    EXPECTED_SHA256,
    BraunRCCAdapter,
    braun_source_paths,
    immunogenic_call,
)
from event_b.audit import corpus_audit
from event_b.braun_pipeline import (
    assert_quality_gates,
    build_braun_corpus,
    combine_corpora,
    reconcile_braun,
)
from event_b.manifest import manifest_from_paths
from event_b.models import BiologicalEvent, ResponseLabel
from event_b.validation import validate_corpus


RAW_DIR = Path("data/raw/braun_rcc_2025")
HAS_SOURCE = all((RAW_DIR / "extracted" / name).exists() for name in EXPECTED_SHA256)
skip_no_source = pytest.mark.skipif(
    not HAS_SOURCE, reason="Braun supplements absent (gitignored); fetched in the milestone-5a CI job"
)


def _row(stim, nostim, pvalue):
    data = {"Ttest_pvalue_InVitroStim": pvalue}
    for i, value in enumerate(stim, start=1):
        data[f"InVitro_PeptideStim_Replicate0{i}"] = value
    for i, value in enumerate(nostim, start=1):
        data[f"InVitro_NoStim_Replicate0{i}"] = value
    return pd.Series(data)


def test_immunogenic_call_follows_paper_rule():
    # P<0.05 AND mean stim >= 3x mean no-stim.
    assert immunogenic_call(_row([100, 90, 110], [0, 0, 0], 0.001)) is True
    assert immunogenic_call(_row([30, 30, 30], [10, 10, 10], 0.01)) is True  # exactly 3-fold
    assert immunogenic_call(_row([10, 10, 10], [5, 5, 5], 0.001)) is False  # significant but <3-fold
    assert immunogenic_call(_row([100, 100, 100], [10, 10, 10], 0.20)) is False  # 10-fold but n.s.
    assert immunogenic_call(_row([100, 100, 100], [0, 0, 0], None)) is None  # unscorable


def _synthetic_extracted(week0_baseline: float = 5.0):
    def iv(pid, peptide, stim, nostim, pvalue, mutation, position, gene="GENEX", change="GENEX|p.A1T"):
        return {
            "Patient_ID": "101",
            "Peptide_ID": pid,
            "Vaccine_Peptide": peptide,
            "Gene_and_Protein_Change": change,
            "Hugo_Symbol": gene,
            "Chromosome": "1",
            "Start_position": position,
            "Variant_Type": "SNV",
            "Mutation_type": mutation,
            "Pool": "A",
            "InVitro_PeptideStim_Replicate01": stim[0],
            "InVitro_PeptideStim_Replicate02": stim[1],
            "InVitro_PeptideStim_Replicate03": stim[2],
            "InVitro_NoStim_Replicate01": nostim[0],
            "InVitro_NoStim_Replicate02": nostim[1],
            "InVitro_NoStim_Replicate03": nostim[2],
            "Ttest_pvalue_InVitroStim": pvalue,
        }

    invitro = pd.DataFrame(
        [
            iv("P1", "MKQEVTIKALKEKIREYEQAL", [100, 90, 110], [0, 0, 0], 0.001, "Driver", 1000),
            iv("P2", "STRDPLSEITKQEKDFLWSHRH", [5, 4, 6], [4, 5, 5], 0.90, "Passenger", 2000),
            iv("P3", "NQRNNVVRNSRTSGYNVRNSRT", [10, 10, 10], [5, 5, 5], 0.001, "Passenger", 3000),
            iv("P4", "KVLHSAILRGCICALVLFRTES", [50, 50, 50], [5, 5, 5], 0.01, "", 4000),
            iv("P5", "GILGFVFTLKVLHSAILRGCIC", [100, 100, 100], [0, 0, 0], None, "Passenger", 5000),
        ]
    )
    exvivo = pd.DataFrame(
        [
            {
                "Patient ID": "101",
                "Vaccine Pool": "A",
                "Treatment": "Vaccine + ipilimumab",
                "Week 0 (mean, background-subtracted)": week0_baseline,
                "Week 3 (mean, background-subtracted)": 200.0,
            }
        ]
    )
    meta = pd.DataFrame(
        [{"ID": "16097-101", "Cohort": "Vaccine + ipilimumab", "Stage": "III"}]
    )
    return {"invitro": invitro, "exvivo": exvivo, "meta": meta}


def _normalize(extracted, tmp_path):
    source = tmp_path / "source.txt"
    source.write_text("synthetic")
    manifest = manifest_from_paths("Braun", "test", "braun_rcc_event_b", "1.0.0", [source])
    adapter = BraunRCCAdapter(tmp_path)
    corpus = adapter.normalize(extracted, manifest)
    return adapter, corpus


def test_normalize_assigns_three_state_labels_and_routes_review(tmp_path):
    adapter, corpus = _normalize(_synthetic_extracted(), tmp_path)
    accepted = validate_corpus(corpus).accepted_corpus

    # 5 assayed peptides in; 3 accepted (P1 positive, P2/P3 negative); P4/P5 held for review.
    assert len(corpus.assays) == 5
    assert len(accepted.assays) == 3
    labels = accepted.assays.response_label.astype(str).str.upper().value_counts().to_dict()
    assert labels == {ResponseLabel.POSITIVE.value: 1, ResponseLabel.TESTED_NEGATIVE.value: 2}

    codes = sorted(issue.code for issue in adapter.review_issues)
    assert codes == ["UNCLASSIFIED_MUTATION", "UNSCORABLE_ASSAY"]

    # No validator issues on clean rows, and no negative inferred from omission.
    assert validate_corpus(corpus).review_queue == ()
    assert (accepted.assays.explicit_assay_inclusion == True).all()  # noqa: E712


def test_de_novo_event_typing_requires_clean_baseline(tmp_path):
    _, de_novo = _normalize(_synthetic_extracted(week0_baseline=5.0), tmp_path)
    events = de_novo.assays.event_type.astype(str).str.upper().unique().tolist()
    assert events == [BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value]
    assert (
        de_novo.assays.relative_to_vaccine.astype(str).str.upper().eq("POST_VACCINE").all()
    )

    # A pre-existing (week-0) pool response above the ELISpot floor blocks a clean de-novo claim.
    _, boosted = _normalize(_synthetic_extracted(week0_baseline=80.0), tmp_path)
    assert (
        boosted.assays.event_type.astype(str).str.upper().eq(BiologicalEvent.UNKNOWN_EVENT.value).all()
    )


def test_every_accepted_label_has_provenance(tmp_path):
    _, corpus = _normalize(_synthetic_extracted(), tmp_path)
    accepted = validate_corpus(corpus).accepted_corpus
    provenance_ids = set(accepted.provenance.provenance_id.astype(str))
    for frame in (accepted.assays, accepted.candidates, accepted.patients, accepted.vaccines):
        assert frame.provenance_id.astype(str).isin(provenance_ids).all()
    assert_quality_gates(accepted)  # gates pass on a clean slice


def test_quality_gates_catch_violations(tmp_path):
    _, corpus = _normalize(_synthetic_extracted(), tmp_path)
    accepted = validate_corpus(corpus).accepted_corpus

    tampered = accepted.assays.copy()
    negatives = tampered.response_label.astype(str).str.upper().eq(ResponseLabel.TESTED_NEGATIVE.value)
    tampered.loc[negatives, "explicit_assay_inclusion"] = False
    broken = validate_corpus(corpus).accepted_corpus
    broken.assays = tampered
    with pytest.raises(AssertionError, match="explicit assay inclusion"):
        assert_quality_gates(broken)


def test_no_event_a_relabelled_as_event_b(tmp_path):
    _, corpus = _normalize(_synthetic_extracted(), tmp_path)
    assert not corpus.assays.event_type.astype(str).str.upper().eq(
        BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value
    ).any()


@skip_no_source
def test_real_source_reconciles_to_paper(tmp_path):
    build = build_braun_corpus(RAW_DIR)
    recon = reconcile_braun(RAW_DIR, build)
    assert recon["summary_reconciles"] is True
    assert recon["accepted"]["positives"] == 61
    assert recon["accepted"]["tested_negatives"] == 68
    assert recon["accepted"]["patients"] == 9
    assert recon["review_queue"]["by_code"].get("UNCLASSIFIED_MUTATION") == 1
    assert_quality_gates(build.result.accepted_corpus)

    audit = corpus_audit(build.result.accepted_corpus, build.review_queue)
    assert audit["event_counts"] == {BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value: 129}
    assert audit["model_readiness"]["event_b_patient_n"] == 9
    assert (
        audit["model_readiness"]["decision"]
        == "EVENT_B_VERTICAL_SLICE_VALIDATED_NOT_YET_SUFFICIENT_FOR_GENERAL_MODEL"
    )


@skip_no_source
def test_real_source_manifest_and_build_are_deterministic():
    first = braun_source_paths(RAW_DIR)
    second = braun_source_paths(RAW_DIR)
    assert [p.name for p in first] == [p.name for p in second]

    build_a = build_braun_corpus(RAW_DIR)
    build_b = build_braun_corpus(RAW_DIR)
    key = ["assay_id", "response_label", "event_type", "provenance_id"]
    left = build_a.result.accepted_corpus.assays.sort_values("assay_id")[key].reset_index(drop=True)
    right = build_b.result.accepted_corpus.assays.sort_values("assay_id")[key].reset_index(drop=True)
    pd.testing.assert_frame_equal(left, right)


@skip_no_source
def test_combined_corpus_separates_events_and_moves_off_zero():
    from event_b.braun_pipeline import load_corpus_from_parquet

    improve_dir = Path("outputs/event_b_corpus")
    if not (improve_dir / "assays.parquet").exists():
        pytest.skip("frozen IMPROVE corpus not materialised (outputs/ is gitignored)")
    build = build_braun_corpus(RAW_DIR)
    combined = combine_corpora(load_corpus_from_parquet(improve_dir), build.result.accepted_corpus)
    audit = corpus_audit(combined, build.review_queue)
    assert audit["event_counts"][BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value] > 0
    assert audit["event_counts"][BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value] == 129
    assert audit["model_readiness"]["event_b_patient_n"] == 9
