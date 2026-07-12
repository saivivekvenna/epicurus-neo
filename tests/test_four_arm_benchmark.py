"""Tests for the four-arm generation x scorer patient-level top-20 benchmark harness.

The harness compares, per patient, four arms that factor the neoantigen pipeline into a
generation stage and a scoring/selection stage:

  1. pvac_prime         standard pVAC candidates       + genuine PRIME  + plain top-k
  2. lossless_prime     lossless-generation union      + genuine PRIME  + plain top-k
  3. lossless_epicurus  lossless-generation union      + Epicurus       + plain top-k
  4. full_epicurus      lossless-generation union      + Epicurus       + route-aware select

Strict label isolation: measured positives are passed separately and only ever used to
score coverage AFTER ranking. Every arm reports an explicit NOT_EVALUABLE status when a
required input is missing, so a cohort without lossless generation or without labels never
silently produces a number.

These tests are pure logic on synthetic frames whose answer is known by construction.
"""

from __future__ import annotations

import pandas as pd

from benchmark.four_arm import (
    FOUR_ARMS,
    REQ_EPICURUS,
    REQ_LABELS,
    REQ_LOSSLESS,
    REQ_PRIME,
    arm_requirements,
    attach_epicurus_score,
    evaluate_eligibility,
    run_arm,
    run_patient,
    stage_attribution,
)


# ---------------------------------------------------------------------------
# Synthetic universe helper
# ---------------------------------------------------------------------------
def _cand(mutation_id, peptide, hla, source, prime, epicurus, **extra):
    """One candidate row with the router-relevant defaults (missense/expressed/valid)."""
    row = {
        "patient_id": "p1",
        "mutation_id": mutation_id,
        "gene_symbol": mutation_id.split("-")[0],
        "mutant_peptide": peptide,
        "hla_allele": hla,
        "mhc_class": "I",
        "source_variant_type": "SNV",
        "expression_tpm": 50.0,
        "binding_percentile_rank": 0.5,
        "n_callers": 2,
        "n_timepoints": 1,
        "candidate_source": source,
        "genuine_prime": prime,
        "epicurus": epicurus,
    }
    row.update(extra)
    return row


def _universe(rows):
    return pd.DataFrame(rows)


ARM_IDS = ["pvac_prime", "lossless_prime", "lossless_epicurus", "full_epicurus"]


# ---------------------------------------------------------------------------
# Eligibility
# ---------------------------------------------------------------------------
def test_four_arms_are_the_canonical_2x2():
    assert [a.arm_id for a in FOUR_ARMS] == ARM_IDS
    by_id = {a.arm_id: a for a in FOUR_ARMS}
    assert by_id["pvac_prime"].generation == "pvac"
    assert by_id["lossless_prime"].generation == "lossless_union"
    assert by_id["lossless_epicurus"].scorer == "epicurus"
    assert by_id["full_epicurus"].selection == "route_aware"


def test_arm_requirements_track_generation_and_scorer():
    by_id = {a.arm_id: a for a in FOUR_ARMS}
    assert REQ_PRIME in arm_requirements(by_id["pvac_prime"])
    assert REQ_LOSSLESS not in arm_requirements(by_id["pvac_prime"])
    assert REQ_LOSSLESS in arm_requirements(by_id["lossless_prime"])
    assert REQ_EPICURUS in arm_requirements(by_id["lossless_epicurus"])
    # labels are required by every arm (no denominator, no metric)
    assert all(REQ_LABELS in arm_requirements(a) for a in FOUR_ARMS)


def test_eligibility_all_inputs_present():
    available = {REQ_LABELS, "pvac_candidates", REQ_LOSSLESS, REQ_PRIME, REQ_EPICURUS, "router_features"}
    elig = evaluate_eligibility(available)
    assert all(elig[a].evaluable for a in ARM_IDS)
    assert all(elig[a].missing == [] for a in ARM_IDS)


def test_eligibility_without_lossless_only_pvac_arm_evaluable():
    available = {REQ_LABELS, "pvac_candidates", REQ_PRIME, REQ_EPICURUS, "router_features"}
    elig = evaluate_eligibility(available)
    assert elig["pvac_prime"].evaluable
    for arm in ("lossless_prime", "lossless_epicurus", "full_epicurus"):
        assert not elig[arm].evaluable
        assert REQ_LOSSLESS in elig[arm].missing


def test_eligibility_without_labels_nothing_is_evaluable():
    available = {"pvac_candidates", REQ_LOSSLESS, REQ_PRIME, REQ_EPICURUS, "router_features"}
    elig = evaluate_eligibility(available)
    assert not any(elig[a].evaluable for a in ARM_IDS)
    assert all(REQ_LABELS in elig[a].missing for a in ARM_IDS)


def test_eligibility_without_epicurus_only_prime_arms_evaluable():
    available = {REQ_LABELS, "pvac_candidates", REQ_LOSSLESS, REQ_PRIME, "router_features"}
    elig = evaluate_eligibility(available)
    assert elig["pvac_prime"].evaluable and elig["lossless_prime"].evaluable
    assert not elig["lossless_epicurus"].evaluable
    assert not elig["full_epicurus"].evaluable


# ---------------------------------------------------------------------------
# Generation recall (the reachability stage)
# ---------------------------------------------------------------------------
def test_generation_recall_excludes_positive_only_in_lossless_from_pvac_arm():
    # HIT lives only in a lossless-recovered candidate; pVAC never generated it.
    uni = _universe([
        _cand("SEEN-1-1", "AAAAAAAAA", "HLA-A*01:01", "pvactools_2025_01", prime=-0.1, epicurus=2.0),
        _cand("MISS-2-2", "CCCCCCCCC", "HLA-A*01:01", "lossless_recovery", prime=-0.2, epicurus=1.0),
    ])
    positives = {"SEEN-1-1", "MISS-2-2"}
    by_id = {a.arm_id: a for a in FOUR_ARMS}

    pvac = run_arm(uni, positives, by_id["pvac_prime"])
    lossless = run_arm(uni, positives, by_id["lossless_prime"])

    assert pvac.generation_recall.n == 1  # only SEEN-1-1 reachable
    assert "MISS-2-2" not in pvac.generation_recall.ids
    assert lossless.generation_recall.n == 2  # union reaches both
    assert set(lossless.generation_recall.ids) == positives


# ---------------------------------------------------------------------------
# Top-20 mutation coverage (ranking stage), granularity = mutation
# ---------------------------------------------------------------------------
def test_hits_at_k_counts_positive_mutations_in_top_k_only():
    # k=2: the positive's only candidate sits at rank 3 -> not covered.
    rows = [
        _cand("A-1-1", "AAAAAAAAA", "HLA-A*01:01", "pvactools_2025_01", prime=-0.01, epicurus=9.0),
        _cand("B-2-2", "CCCCCCCCC", "HLA-A*01:01", "pvactools_2025_01", prime=-0.02, epicurus=8.0),
        _cand("HIT-3-3", "DDDDDDDDD", "HLA-A*01:01", "pvactools_2025_01", prime=-0.03, epicurus=1.0),
    ]
    uni = _universe(rows)
    positives = {"HIT-3-3"}
    by_id = {a.arm_id: a for a in FOUR_ARMS}
    res = run_arm(uni, positives, by_id["pvac_prime"], k=2)
    assert res.hits_at_k == 0
    assert res.recall_at_k == 0.0

    res_k3 = run_arm(uni, positives, by_id["pvac_prime"], k=3)
    assert res_k3.hits_at_k == 1
    assert res_k3.recall_at_k == 1.0


def test_scorer_swap_can_change_coverage():
    # PRIME ranks the decoy top; Epicurus ranks the true positive top. k=1 separates them.
    rows = [
        _cand("DECOY-1-1", "AAAAAAAAA", "HLA-A*01:01", "pvactools_2025_01", prime=-0.01, epicurus=1.0),
        _cand("HIT-2-2", "CCCCCCCCC", "HLA-A*01:01", "pvactools_2025_01", prime=-0.9, epicurus=9.0),
    ]
    uni = _universe(rows)
    positives = {"HIT-2-2"}
    by_id = {a.arm_id: a for a in FOUR_ARMS}
    prime = run_arm(uni, positives, by_id["lossless_prime"], k=1)
    epi = run_arm(uni, positives, by_id["lossless_epicurus"], k=1)
    assert prime.hits_at_k == 0
    assert epi.hits_at_k == 1


def test_plain_topk_tiebreak_is_deterministic():
    # Identical scores -> selection must be stable across row permutations.
    rows = [
        _cand("A-1-1", "AAAAAAAAA", "HLA-A*01:01", "pvactools_2025_01", prime=-0.5, epicurus=1.0),
        _cand("B-2-2", "CCCCCCCCC", "HLA-A*01:01", "pvactools_2025_01", prime=-0.5, epicurus=1.0),
        _cand("C-3-3", "DDDDDDDDD", "HLA-A*01:01", "pvactools_2025_01", prime=-0.5, epicurus=1.0),
    ]
    by_id = {a.arm_id: a for a in FOUR_ARMS}
    uni = _universe(rows)
    rev = _universe(list(reversed(rows)))
    a = run_arm(uni, {"A-1-1", "B-2-2", "C-3-3"}, by_id["pvac_prime"], k=2)
    b = run_arm(rev, {"A-1-1", "B-2-2", "C-3-3"}, by_id["pvac_prime"], k=2)
    assert a.top_k.ids == b.top_k.ids


# ---------------------------------------------------------------------------
# Route-aware selection (arm 4) respects rankability
# ---------------------------------------------------------------------------
def test_route_aware_arm_drops_needs_peptide_generation_rows():
    # A row with no peptide is valid-but-not-rankable; route-aware must not select it.
    rows = [
        _cand("HIT-1-1", "AAAAAAAAA", "HLA-A*01:01", "lossless_recovery", prime=-0.1, epicurus=9.0),
        _cand("NOPEP-2-2", "", "HLA-A*01:01", "lossless_recovery", prime=-0.01, epicurus=99.0),
    ]
    uni = _universe(rows)
    positives = {"HIT-1-1", "NOPEP-2-2"}
    by_id = {a.arm_id: a for a in FOUR_ARMS}
    res = run_arm(uni, positives, by_id["full_epicurus"], k=20)
    assert "HIT-1-1" in res.top_k.ids
    assert "NOPEP-2-2" not in res.top_k.ids
    assert res.rankable_recall.n == 1


# ---------------------------------------------------------------------------
# Stage attribution is additive
# ---------------------------------------------------------------------------
def test_stage_attribution_decomposes_total_additively():
    # Positive reachable only via lossless AND only rankable when Epicurus ranks it up.
    rows = [
        _cand("SEEN-1-1", "AAAAAAAAA", "HLA-A*01:01", "pvactools_2025_01", prime=-0.01, epicurus=1.0),
        _cand("GEN-2-2", "CCCCCCCCC", "HLA-A*01:01", "lossless_recovery", prime=-0.02, epicurus=9.0),
    ]
    uni = _universe(rows)
    positives = {"SEEN-1-1", "GEN-2-2"}
    out = run_patient(uni, positives)
    attr = stage_attribution(out["arms"])
    assert attr["evaluable"] is True
    total = out["arms"]["full_epicurus"].hits_at_k - out["arms"]["pvac_prime"].hits_at_k
    assert attr["total"] == total
    assert attr["generation"] + attr["scorer"] + attr["selection"] == attr["total"]


def test_run_patient_marks_arms_not_evaluable_when_no_lossless_source():
    rows = [
        _cand("A-1-1", "AAAAAAAAA", "HLA-A*01:01", "pvactools_2025_01", prime=-0.1, epicurus=2.0),
    ]
    uni = _universe(rows)
    out = run_patient(uni, {"A-1-1"})
    assert out["arms"]["pvac_prime"].evaluable
    assert not out["arms"]["lossless_prime"].evaluable
    assert REQ_LOSSLESS in out["arms"]["lossless_prime"].missing
    # attribution cannot be computed when an arm is NOT_EVALUABLE
    assert stage_attribution(out["arms"])["evaluable"] is False


def test_attach_epicurus_score_is_monotone_in_prime_within_patient():
    # Holding el/expr fixed, a better (lower) PRIME %rank must yield a higher Epicurus score.
    uni = pd.DataFrame({
        "patient_id": ["p1", "p1"],
        "prime": [0.05, 0.90],
        "el": [0.3, 0.3],
        "expr": [50.0, 50.0],
    })
    scored = attach_epicurus_score(uni)
    assert "epicurus" in scored.columns
    assert scored["epicurus"].iloc[0] > scored["epicurus"].iloc[1]


def test_attach_epicurus_score_tolerates_missing_el():
    uni = pd.DataFrame({
        "patient_id": ["p1", "p1"],
        "prime": [0.05, 0.90],
        "el": [float("nan"), float("nan")],
        "expr": [50.0, 50.0],
    })
    scored = attach_epicurus_score(uni)
    assert scored["epicurus"].notna().all()


def test_run_patient_zero_positives_is_not_evaluable():
    rows = [_cand("A-1-1", "AAAAAAAAA", "HLA-A*01:01", "pvactools_2025_01", prime=-0.1, epicurus=2.0)]
    out = run_patient(_universe(rows), set())
    assert all(not out["arms"][a].evaluable for a in ARM_IDS)
    assert all(REQ_LABELS in out["arms"][a].missing for a in ARM_IDS)
