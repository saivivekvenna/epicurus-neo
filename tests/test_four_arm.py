"""Invariant tests for the four-arm harness: arm-universe separation, UNIQUE-mutation hits@20 (a mutation
counts once regardless of how many of its peptides are selected), ≤k slots with n_selected/saturation,
additive stage attribution, eligibility/NOT_EVALUABLE, and label isolation (positives never change ranking).
"""

from __future__ import annotations

import pandas as pd

from benchmark.four_arm import (
    FOUR_ARMS,
    REQ_EPICURUS,
    REQ_LABELS,
    REQ_LOSSLESS,
    REQ_PRIME,
    REQ_PVAC,
    REQ_ROUTER,
    ArmResult,
    Coverage,
    _coverage,
    _generation_rows,
    arm_requirements,
    evaluate_eligibility,
    is_pvac_source,
    run_patient,
    stage_attribution,
)


# ---- UNIQUE-mutation hits@20: 3 peptides from ONE mutation = 1 hit, not 3 --------------------------
def test_coverage_counts_unique_mutations_not_peptides():
    # a mutation appearing in many selected slots is counted ONCE (set intersection)
    cov = _coverage(["M1", "M1", "M1", "M2"], {"M1", "M2", "M3"})
    assert cov.n == 2 and cov.of == 3 and cov.ids == ["M1", "M2"]      # 3/M1-peptides -> 1; M3 not recovered
    assert cov.frac == round(2 / 3, 4)
    # full recovery = 3 distinct mutations present
    assert _coverage(["M1", "M2", "M3", "M3"], {"M1", "M2", "M3"}).n == 3


# ---- arm-universe separation: pvac arm != lossless arms ---------------------------------------------
def test_generation_rows_separate_pvac_from_lossless_universe():
    uni = pd.DataFrame({
        "mutation_id": ["A", "B", "C"],
        "candidate_source": ["pvac", "lossless_recovered", "pvactools_v4"],
    })
    pvac = _generation_rows(uni, "pvac")
    loss = _generation_rows(uni, "lossless_union")
    assert set(pvac["mutation_id"]) == {"A", "C"}          # only pVAC-source rows
    assert set(loss["mutation_id"]) == {"A", "B", "C"}     # lossless union = the whole universe
    assert is_pvac_source("pvac") and is_pvac_source("pvactools_v4") and not is_pvac_source("lossless_recovered")


def test_lossless_prime_and_lossless_epicurus_share_one_universe():
    # both lossless arms draw generation rows from the identical (full) universe; only the scorer differs
    specs = {a.arm_id: a for a in FOUR_ARMS}
    assert specs["lossless_prime"].generation == specs["lossless_epicurus"].generation == "lossless_union"
    assert specs["pvac_prime"].generation == "pvac"
    assert specs["lossless_prime"].scorer == "genuine_prime" and specs["lossless_epicurus"].scorer == "epicurus"


# ---- eligibility / NOT_EVALUABLE --------------------------------------------------------------------
def test_eligibility_requirements_and_not_evaluable():
    # lossless arms need REQ_LOSSLESS; without it they are NOT_EVALUABLE while pvac_prime can still run
    avail = {REQ_LABELS, REQ_PVAC, REQ_PRIME, REQ_ROUTER}
    elig = evaluate_eligibility(avail)
    assert elig["pvac_prime"].evaluable is True
    assert elig["lossless_prime"].evaluable is False and REQ_LOSSLESS in elig["lossless_prime"].missing
    assert elig["lossless_epicurus"].evaluable is False
    # epicurus arms need REQ_EPICURUS
    assert REQ_EPICURUS in arm_requirements(FOUR_ARMS[2])


def test_no_pvac_run_only_disables_pvac_arm():
    # harness bug fix: absence of a genuine pVAC run must NOT disable the lossless arms
    assert REQ_PVAC not in arm_requirements(FOUR_ARMS[1])         # lossless_prime
    assert REQ_PVAC not in arm_requirements(FOUR_ARMS[2])         # lossless_epicurus
    assert REQ_PVAC not in arm_requirements(FOUR_ARMS[3])         # full_epicurus
    assert REQ_PVAC in arm_requirements(FOUR_ARMS[0])             # only pvac_prime needs it
    avail = {REQ_LABELS, REQ_LOSSLESS, REQ_PRIME, REQ_EPICURUS, REQ_ROUTER}   # no REQ_PVAC (no pVAC run)
    elig = evaluate_eligibility(avail)
    assert elig["pvac_prime"].evaluable is False and REQ_PVAC in elig["pvac_prime"].missing
    assert elig["lossless_prime"].evaluable and elig["lossless_epicurus"].evaluable and elig["full_epicurus"].evaluable


# ---- additive stage attribution ---------------------------------------------------------------------
def test_stage_attribution_is_additive():
    def r(arm, hits):
        return ArmResult(arm_id=arm, evaluable=True, top_k=Coverage(hits, 3, []))
    results = {"pvac_prime": r("pvac_prime", 1), "lossless_prime": r("lossless_prime", 2),
               "lossless_epicurus": r("lossless_epicurus", 2), "full_epicurus": r("full_epicurus", 3)}
    a = stage_attribution(results)
    assert a["evaluable"] and a["generation"] == 1 and a["scorer"] == 0 and a["selection"] == 1
    assert a["total"] == a["generation"] + a["scorer"] + a["selection"] == 2


# ---- end-to-end: ≤k slots (saturation) + unique-mutation hits + label isolation ---------------------
def _mini_universe():
    # one recognized mutation M1 with THREE peptides; two decoy mutations. All rankable (peptide+HLA present).
    rows = []
    for i, pep in enumerate(["AAAAAAAAA", "CCCCCCCCC", "DDDDDDDDD"]):        # 3 peptides of M1
        rows.append({"patient_id": "p", "mutation_id": "M1", "candidate_source": "lossless_recovered",
                     "mutant_peptide": pep, "hla_allele": "HLA-A*02:01", "genuine_prime": 0.9 - i * 0.01,
                     "epicurus": 0.9 - i * 0.01})
    for j, mid in enumerate(["D1", "D2"]):
        rows.append({"patient_id": "p", "mutation_id": mid, "candidate_source": "pvac",
                     "mutant_peptide": f"PEP{j}KLMNP", "hla_allele": "HLA-A*02:01",
                     "genuine_prime": 0.5 - j * 0.01, "epicurus": 0.5 - j * 0.01})
    return pd.DataFrame(rows)


def test_run_patient_saturation_and_unique_mutation_and_label_isolation():
    uni = _mini_universe()
    out = run_patient(uni, {"M1"}, k=20)
    lp = out["arms"]["lossless_prime"]
    if lp.evaluable and lp.top_k is not None:
        # M1 has 3 selected peptides but counts as ONE unique-mutation hit
        assert lp.hits_at_k == 1
        # fewer than 20 rankable candidates -> n_selected < 20 (list did NOT saturate)
        assert lp.n_selected is not None and lp.n_selected < 20
        # label isolation: changing the positive set cannot change what/how many rows are selected
        out2 = run_patient(uni, {"D1"}, k=20)
        assert out2["arms"]["lossless_prime"].n_selected == lp.n_selected
