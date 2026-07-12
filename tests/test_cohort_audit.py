"""Tests for the three-level benchmark cohort eligibility audit.

The benchmark is a hierarchy of three DISTINCT tasks, never pooled:
  1. reachability (raw -> generation, stage-loss attribution)
  2. conditional ranking (only among generated/rankable candidates, within one denominator)
  3. end-to-end patient utility (recognized in final top-20 from common raw inputs) = north star

Each cohort has a fixed role and is eligible for a subset of levels. The four-arm harness only
produces an end-to-end headline where a cohort is Level-3 eligible.
"""

from __future__ import annotations

from benchmark.cohort_audit import (
    COHORTS,
    LEVEL_CONDITIONAL_RANKING,
    LEVEL_END_TO_END,
    LEVEL_REACHABILITY,
    audit_cohort,
    classify_levels,
    level_requirements,
    run_cohort_audit,
)
from benchmark.four_arm import REQ_LABELS, REQ_LOSSLESS, REQ_PRIME


def _cohort(cid):
    return next(c for c in COHORTS if c.cohort_id == cid)


def _by_id(audit):
    return {c["cohort_id"]: c for c in audit["cohorts"]}


# ---------------------------------------------------------------------------
# Level requirements
# ---------------------------------------------------------------------------
def test_reachability_needs_raw_generation_and_labels():
    reqs = level_requirements(LEVEL_REACHABILITY)
    assert REQ_LOSSLESS in reqs and REQ_LABELS in reqs


def test_conditional_ranking_needs_scorer_and_labels_not_raw_generation():
    reqs = level_requirements(LEVEL_CONDITIONAL_RANKING)
    assert REQ_PRIME in reqs and REQ_LABELS in reqs
    assert REQ_LOSSLESS not in reqs  # ranks among ALREADY-generated candidates


def test_end_to_end_needs_raw_generation_scorer_and_labels():
    reqs = level_requirements(LEVEL_END_TO_END)
    assert {REQ_LOSSLESS, REQ_PRIME, REQ_LABELS} <= reqs


# ---------------------------------------------------------------------------
# Per-cohort level eligibility
# ---------------------------------------------------------------------------
def test_sid_is_eligible_for_all_three_levels():
    lv = classify_levels(_cohort("osteosarc_sid"))
    assert lv[LEVEL_REACHABILITY]["eligible"]
    assert lv[LEVEL_CONDITIONAL_RANKING]["eligible"]
    assert lv[LEVEL_END_TO_END]["eligible"]


def test_gartner_is_conditional_ranking_only():
    lv = classify_levels(_cohort("gartner_nci"))
    assert lv[LEVEL_CONDITIONAL_RANKING]["eligible"]
    assert not lv[LEVEL_REACHABILITY]["eligible"]
    assert not lv[LEVEL_END_TO_END]["eligible"]
    assert REQ_LOSSLESS in lv[LEVEL_END_TO_END]["missing"]


def test_zhao_is_eligible_for_no_levels_prime_blocked():
    lv = classify_levels(_cohort("zhao_dc_2026"))
    assert not any(lv[level]["eligible"] for level in
                   (LEVEL_REACHABILITY, LEVEL_CONDITIONAL_RANKING, LEVEL_END_TO_END))
    assert REQ_PRIME in lv[LEVEL_CONDITIONAL_RANKING]["missing"]


def test_rttp_is_eligible_for_no_levels_no_labels():
    lv = classify_levels(_cohort("rttp_sr24_58221"))
    for level in (LEVEL_REACHABILITY, LEVEL_CONDITIONAL_RANKING, LEVEL_END_TO_END):
        assert not lv[level]["eligible"]
        assert REQ_LABELS in lv[level]["missing"]


def test_only_one_cohort_is_end_to_end_eligible():
    audit = run_cohort_audit()
    assert audit["n_end_to_end_eligible"] == 1
    end = [c["cohort_id"] for c in audit["cohorts"]
           if c["levels"]["end_to_end_patient_utility"]["eligible"]]
    assert end == ["osteosarc_sid"]


# ---------------------------------------------------------------------------
# Cohort roles are fixed and explicit
# ---------------------------------------------------------------------------
def test_every_cohort_has_a_role():
    assert all(c.role for c in COHORTS)


def test_cohort_roles_match_the_product_decision():
    role = {c.cohort_id: c.role for c in COHORTS}
    assert "training" in role["cedar_tcell"].lower() or "prior" in role["cedar_tcell"].lower()
    assert "training" in role["zhao_dc_2026"].lower() or "prior" in role["zhao_dc_2026"].lower()
    assert "presentation" in role["cd8_multimer"].lower()
    assert "broad" in role["gartner_nci"].lower()
    assert "prefiltered" in role["improve_srhgroup"].lower() or "subset" in role["improve_srhgroup"].lower()
    assert "end-to-end" in role["osteosarc_sid"].lower()
    assert "deployment" in role["rttp_sr24_58221"].lower()


# ---------------------------------------------------------------------------
# No pooling invariant
# ---------------------------------------------------------------------------
def test_audit_declares_the_no_pooling_invariant():
    audit = run_cohort_audit()
    assert audit["no_pooling"]  # a stated invariant string
    # there must be no cross-cohort aggregate performance metric at the top level
    assert not any("pooled" in k.lower() or "aggregate" in k.lower() for k in audit)


def test_audit_carries_a_per_cohort_denominator_so_cohorts_stay_separate():
    audit = run_cohort_audit()
    for c in audit["cohorts"]:
        assert c["denominator"]


# ---------------------------------------------------------------------------
# Leakage + scorer-only invariants (retained from the arm-level audit)
# ---------------------------------------------------------------------------
def test_cd8_multimer_epicurus_arms_carry_a_leakage_reason():
    c = _by_id(run_cohort_audit())["cd8_multimer"]
    for aid in ("lossless_epicurus", "full_epicurus"):
        assert any(m.startswith("LEAKAGE:") for m in c["arms"][aid]["missing"])
    assert not any(m.startswith("LEAKAGE:") for m in c["arms"]["pvac_prime"]["missing"])


def test_cedar_is_present_as_a_prior_asset_not_a_ranking_cohort():
    lv = classify_levels(_cohort("cedar_tcell"))
    assert not any(lv[level]["eligible"] for level in
                   (LEVEL_REACHABILITY, LEVEL_CONDITIONAL_RANKING, LEVEL_END_TO_END))


def test_audit_cohort_shape_includes_role_and_levels():
    a = audit_cohort(COHORTS[0])
    assert set(a) >= {"cohort_id", "role", "levels", "arms", "denominator", "note"}
