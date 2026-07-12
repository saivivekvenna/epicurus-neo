"""Targeted tests for the pool-size sensitivity diagnostic.

Pure logic only — no PRIME binary, no caches. Verifies the load-bearing invariants:
nestedness (SMALL subset MEDIUM subset LARGE), oracle all-positive retention (variant A keeps every
positive by construction), determinism of the seeded negative sampling, size-matching + potential
positive loss for the label-blind gate (variant B), and the pool<=20 saturation flag.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from event_b.pool_size_sensitivity import (
    GATE_COL, POOL_FRACS, gate_pools, oracle_pools, patient_eligibility,
    patient_metrics, run_cohort, score_arms,
)


def _synth(n_pos=4, n_neg=40, patient="p1", seed=0) -> pd.DataFrame:
    """One-patient synthetic frame with monotone-ish features so scorers are well-defined."""
    rng = np.random.default_rng(seed)
    n = n_pos + n_neg
    labels = ["POSITIVE"] * n_pos + ["TESTED_NEGATIVE"] * n_neg
    return pd.DataFrame({
        "patient_id": [patient] * n,
        "mutant_peptide": [f"PEP{i:04d}" for i in range(n)],
        "hla_allele": ["A0201"] * n,
        "label": labels,
        "prime": rng.uniform(0, 100, n),
        "el": rng.uniform(0, 100, n),
        "expr": rng.uniform(0, 50, n),
    })


def _ids(df):
    return set(df["mutant_peptide"])


# ---- nestedness ----------------------------------------------------------------------------------
def test_oracle_pools_are_nested():
    df = _synth()
    pools = oracle_pools(df, "c", "p1", seed=3)
    assert _ids(pools["SMALL"]) <= _ids(pools["MEDIUM"]) <= _ids(pools["LARGE"])
    assert _ids(pools["LARGE"]) == _ids(df)


def test_gate_pools_are_nested():
    df = _synth()
    pools, _ = gate_pools(df)
    assert _ids(pools["SMALL"]) <= _ids(pools["MEDIUM"]) <= _ids(pools["LARGE"])


# ---- oracle keeps every positive; identical positives across pools -------------------------------
def test_oracle_retains_all_positives_in_every_pool():
    df = _synth(n_pos=5, n_neg=60)
    pos_ids = _ids(df[df["label"] == "POSITIVE"])
    for pool in oracle_pools(df, "c", "p1", seed=7).values():
        assert _ids(pool[pool["label"] == "POSITIVE"]) == pos_ids


def test_oracle_pool_sizes_match_fractions():
    df = _synth(n_pos=4, n_neg=40)
    pools = oracle_pools(df, "c", "p1", seed=1)
    assert (pools["LARGE"]["label"] == "TESTED_NEGATIVE").sum() == 40
    assert (pools["MEDIUM"]["label"] == "TESTED_NEGATIVE").sum() == int(np.ceil(0.5 * 40))
    assert (pools["SMALL"]["label"] == "TESTED_NEGATIVE").sum() == int(np.ceil(0.25 * 40))


# ---- determinism ---------------------------------------------------------------------------------
def test_oracle_sampling_is_deterministic_per_seed():
    df = _synth()
    a = oracle_pools(df, "c", "p1", seed=11)["SMALL"]
    b = oracle_pools(df, "c", "p1", seed=11)["SMALL"]
    assert _ids(a) == _ids(b)


def test_different_seeds_differ_but_stay_nested():
    df = _synth(n_neg=80)
    s5 = oracle_pools(df, "c", "p1", seed=5)
    s6 = oracle_pools(df, "c", "p1", seed=6)
    assert _ids(s5["SMALL"]) != _ids(s6["SMALL"])  # seeds actually change the draw
    assert _ids(s5["SMALL"]) <= _ids(s5["MEDIUM"])  # nesting preserved regardless of seed


# ---- variant B: size match + honest retention reporting ------------------------------------------
def test_gate_pools_size_matched_to_oracle():
    df = _synth(n_pos=4, n_neg=40)
    gate, _ = gate_pools(df)
    orc = oracle_pools(df, "c", "p1", seed=0)
    for pool in POOL_FRACS:
        assert len(gate[pool]) == len(orc[pool])  # same budget, different selection policy


def test_gate_can_drop_positives_when_positives_present_low():
    # positives forced to the WORST EL (highest raw rank) -> gate (keeps low EL) should drop them.
    df = _synth(n_pos=3, n_neg=40, seed=2)
    df.loc[df["label"] == "POSITIVE", "el"] = 999.0  # worst presentation
    pools, _ = gate_pools(df)
    kept = (pools["SMALL"]["label"] == "POSITIVE").sum()
    assert kept < 3  # deployable gate is NOT an oracle


def test_gate_is_label_blind():
    # gate ordering must be identical whether or not labels are flipped (never reads labels).
    df = _synth()
    _, order_a = gate_pools(df)
    df2 = df.copy()
    df2["label"] = "TESTED_NEGATIVE"  # wipe labels
    _, order_b = gate_pools(df2)
    assert np.allclose(order_a, order_b)


# ---- saturation flag -----------------------------------------------------------------------------
def test_saturation_flag_le20():
    small = _synth(n_pos=2, n_neg=10)  # pool size 12 <= 20
    big = _synth(n_pos=2, n_neg=60)    # pool size 62 > 20
    assert patient_metrics(score_arms(small), "genuine_prime")["saturated_le20"] is True
    assert patient_metrics(score_arms(big), "genuine_prime")["saturated_le20"] is False


# ---- metric sanity -------------------------------------------------------------------------------
def test_perfect_ranker_metrics():
    # make genuine_prime perfectly rank positives first: positives get the lowest raw prime %rank.
    df = _synth(n_pos=3, n_neg=30, seed=4)
    df.loc[df["label"] == "POSITIVE", "prime"] = -1.0  # best possible (lower=better) -> top of -prime
    m = patient_metrics(score_arms(df), "genuine_prime")
    assert m["hits@5"] == 3 and m["recall@20"] == 1.0 and m["mrr"] == 1.0
    assert m["median_pos_rank"] == 2.0  # positions 1,2,3


# ---- eligibility ---------------------------------------------------------------------------------
def test_eligibility_excludes_no_positive_and_too_few_negatives():
    frames = pd.concat([
        _synth(n_pos=2, n_neg=40, patient="ok"),
        _synth(n_pos=0, n_neg=40, patient="nopos"),
        _synth(n_pos=2, n_neg=2, patient="fewneg"),
    ], ignore_index=True)
    elig = patient_eligibility(frames)
    assert elig.eligible == ["ok"]
    assert elig.excluded["nopos"] == "no_positive"
    assert "too_few_negatives" in elig.excluded["fewneg"]


# ---- end-to-end smoke on synthetic multi-patient cohort ------------------------------------------
def test_run_cohort_smoke():
    frames = pd.concat([_synth(n_pos=3, n_neg=50, patient=f"p{i}", seed=i) for i in range(4)], ignore_index=True)
    res = run_cohort(frames, "synth", seeds=[0, 1, 2])
    assert res["n_patients_eligible"] == 4
    for pool in POOL_FRACS:
        # oracle retention is 100% -> recall@20 must be achievable and deltas finite
        assert res["variant_A_oracle"][pool]["arms"]["genuine_prime"]["recall@20"] is not None
        # gate retention reported in [0, 1]
        r = res["variant_B_gate"][pool]["positive_retention_mean"]
        assert 0.0 <= r <= 1.0
