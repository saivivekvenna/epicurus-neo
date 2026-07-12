"""Safety-invariant tests for the dynamic upstream gate (Milestone 7).

These lock the properties that make the gate a SAFE rejector: missing evidence keeps, the AND structure
gives high recall, removal is monotone in the threshold, the CP lower bound is correct, cohort identity is
not an input, and the per-patient rails hold.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from event_b.dynamic_gate import (
    GateConfig,
    apply_gate,
    attach_percentiles,
    calibrate_threshold,
    clopper_pearson_lower,
    gate_retention_stats,
    predictor_disagreement,
    within_patient_percentile,
)


def _pool(n_neg=40, n_pos=4, seed=0, patient="p1"):
    """One patient: negatives spread low-ish, positives spread, all three features present."""
    rng = np.random.default_rng(seed)
    rows = []
    for _ in range(n_neg):
        rows.append({"patient_id": patient, "label": "TESTED_NEGATIVE",
                     "el": rng.uniform(0, 5), "prime": rng.uniform(0, 5), "expr": rng.uniform(0, 10)})
    for _ in range(n_pos):
        rows.append({"patient_id": patient, "label": "POSITIVE",
                     "el": rng.uniform(0, 5), "prime": rng.uniform(0, 5), "expr": rng.uniform(0, 10)})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------------------------------
def test_missing_core_axis_forces_keep():
    """A candidate missing a CORE axis (el or prime) must never be vetoed (fail-open on missing evidence)."""
    df = _pool()
    # bottom on prime+expr, but CORE axis el missing -> AND fails -> KEEP
    df.loc[0, ["prime", "expr"]] = [1000.0, -1000.0]
    df.loc[0, "el"] = np.nan
    g = apply_gate(df, GateConfig(t=0.9, el_floor_frac=0.0, pool_floor=0, cov_floor=0.0))
    assert bool(g.loc[0, "dyn_gate_keep"]) is True
    assert g.loc[0, "dyn_veto_eligible"] == False  # noqa: E712


def test_missing_rescue_axis_does_not_block_veto():
    """A missing RESCUE axis (expr) never causes a veto but also cannot rescue: a candidate bad on both
    core axes with expr absent is still vetoed (the gate operates on the core el+prime)."""
    df = _pool()
    df.loc[0, ["el", "prime"]] = [1000.0, 1000.0]  # worst on both core axes
    df.loc[0, "expr"] = np.nan                     # rescue axis absent
    g = apply_gate(df, GateConfig(t=0.9, el_floor_frac=0.0, per_patient_cap=1.0, pool_floor=0, cov_floor=0.0))
    assert bool(g.loc[0, "dyn_veto_eligible"]) is True
    assert bool(g.loc[0, "dyn_gate_keep"]) is False


def test_strong_on_core_axis_survives():
    """A candidate strong (top) on even one CORE axis is never vetoed."""
    df = _pool()
    # bottom on el+expr but BEST on prime (smallest %rank raw) -> survives
    df.loc[1, ["el", "expr"]] = [1000.0, -1000.0]
    df.loc[1, "prime"] = -999.0  # lower_raw_better -> best -> top percentile
    g = apply_gate(df, GateConfig(t=0.9, el_floor_frac=0.0, pool_floor=0, cov_floor=0.0))
    assert bool(g.loc[1, "dyn_gate_keep"]) is True
    assert bool(g.loc[1, "dyn_veto_eligible"]) is False


def test_high_expression_rescues():
    """Rescue-only expr: bad on both core axes but TOP expression -> rescued (kept)."""
    df = _pool()
    df.loc[2, ["el", "prime"]] = [1000.0, 1000.0]  # worst on both core axes
    df.loc[2, "expr"] = 1e9                         # best expression -> top percentile >= t
    g = apply_gate(df, GateConfig(t=0.9, el_floor_frac=0.0, per_patient_cap=1.0, pool_floor=0, cov_floor=0.0))
    assert bool(g.loc[2, "dyn_veto_eligible"]) is True
    assert bool(g.loc[2, "dyn_gate_keep"]) is True
    assert g.loc[2, "dyn_gate_reason"] == "RESCUE_EXPR"


def test_all_core_bad_and_low_expr_is_vetoed():
    """Bad on both core axes AND low expression -> vetoed when rails are off (the 3-way AND)."""
    df = _pool()
    df.loc[3, ["el", "prime"]] = [1000.0, 1000.0]
    df.loc[3, "expr"] = -1000.0
    g = apply_gate(df, GateConfig(t=0.9, el_floor_frac=0.0, per_patient_cap=1.0, pool_floor=0, cov_floor=0.0))
    assert bool(g.loc[3, "dyn_veto_eligible"]) is True
    assert bool(g.loc[3, "dyn_gate_keep"]) is False


def test_positives_never_used_by_gate():
    """Label-blind: shuffling the label column must not change gate keep decisions."""
    df = _pool()
    cfg = GateConfig(t=0.7, el_floor_frac=0.05, pool_floor=0, cov_floor=0.0)
    keep_a = apply_gate(df, cfg)["dyn_gate_keep"].to_numpy()
    df2 = df.copy()
    df2["label"] = np.random.default_rng(1).permutation(df2["label"].to_numpy())
    keep_b = apply_gate(df2, cfg)["dyn_gate_keep"].to_numpy()
    assert np.array_equal(keep_a, keep_b)


def test_cohort_identity_not_an_input():
    """Renaming the patient/cohort must not change decisions (percentiles are within-patient)."""
    df = _pool(patient="gartner:001")
    cfg = GateConfig(t=0.7, el_floor_frac=0.0, pool_floor=0, cov_floor=0.0)
    keep_a = apply_gate(df, cfg)["dyn_gate_keep"].to_numpy()
    df2 = df.copy()
    df2["patient_id"] = "improve:XYZ"
    keep_b = apply_gate(df2, cfg)["dyn_gate_keep"].to_numpy()
    assert np.array_equal(keep_a, keep_b)


def test_removal_monotone_in_t():
    """Higher t => weakly MORE negatives removed (and weakly less retention). Rails off."""
    df = _pool(n_neg=200, n_pos=20, seed=3)
    cfg = lambda t: GateConfig(t=t, el_floor_frac=0.0, per_patient_cap=1.0, pool_floor=0, cov_floor=0.0)  # noqa: E731
    removals = [gate_retention_stats(apply_gate(df, cfg(t)))["negative_removal"] for t in (0.2, 0.4, 0.6, 0.8)]
    assert all(removals[i] <= removals[i + 1] + 1e-9 for i in range(len(removals) - 1))


def test_el_floor_always_kept():
    """The top-M-by-EL candidates are always kept, even at an aggressive t."""
    df = _pool(n_neg=100, n_pos=10, seed=4)
    cfg = GateConfig(t=0.95, el_floor_frac=0.10, per_patient_cap=1.0, pool_floor=0, cov_floor=0.0)
    g = apply_gate(df, cfg)
    s_el = within_patient_percentile(df, "el", higher_better=False)
    n_floor = int(np.ceil(0.10 * len(df)))
    top_el = np.argsort(-s_el, kind="mergesort")[:n_floor]
    assert g["dyn_gate_keep"].to_numpy()[top_el].all()


def test_pool_floor_keeps_all():
    """A pool smaller than pool_floor is never gated."""
    df = _pool(n_neg=5, n_pos=1, seed=5)
    g = apply_gate(df, GateConfig(t=0.95, pool_floor=8))
    assert g["dyn_gate_keep"].all()


def test_per_patient_cap_falls_back():
    """If a veto would remove more than the cap, the patient falls back to keep-all."""
    # all-junk pool: every negative bad on all axes -> t=0.95 would veto nearly all -> cap triggers
    df = _pool(n_neg=50, n_pos=1, seed=6)
    df.loc[df["label"] == "TESTED_NEGATIVE", ["el", "prime"]] = 1000.0
    df.loc[df["label"] == "TESTED_NEGATIVE", "expr"] = -1000.0
    g = apply_gate(df, GateConfig(t=0.95, el_floor_frac=0.0, per_patient_cap=0.5, pool_floor=0, cov_floor=0.0))
    assert g["dyn_gate_keep"].all()  # cap exceeded -> keep-all
    assert (g["dyn_gate_reason"] == "KEEP_CAP_EXCEEDED").any()


def test_clopper_pearson_bounds():
    assert clopper_pearson_lower(0, 10) == 0.0
    assert clopper_pearson_lower(10, 0) == 0.0
    # full retention: lower bound = alpha**(1/n) < 1
    lb_full = clopper_pearson_lower(20, 20, conf=0.95)
    assert 0.0 < lb_full < 1.0
    assert lb_full == pytest.approx(0.05 ** (1 / 20), rel=1e-6)
    # monotone in k
    assert clopper_pearson_lower(18, 20) < clopper_pearson_lower(19, 20) < clopper_pearson_lower(20, 20)
    # matches known reference: 8/10 -> ~0.4931 (one-sided 95%)
    assert clopper_pearson_lower(8, 10, conf=0.95) == pytest.approx(0.4931, abs=1e-3)


def test_calibrate_picks_largest_safe_t():
    """Calibration returns the most aggressive t whose CP lower bound >= target."""
    df = _pool(n_neg=300, n_pos=40, seed=7)
    out = calibrate_threshold(df, target=0.95, base_config=GateConfig(el_floor_frac=0.0, pool_floor=0, cov_floor=0.0))
    # applying chosen t must satisfy the bound on the calibration set
    cfg = GateConfig(t=out["chosen_t"], el_floor_frac=0.0, pool_floor=0, cov_floor=0.0)
    stats = gate_retention_stats(apply_gate(df, cfg))
    assert stats["positive_retention_cp_lb"] >= 0.95 or out["chosen_t"] == 0.0
    # and any strictly larger t on the grid would violate it (largest-safe property)
    larger = [r for r in out["sweep"] if r["t"] > out["chosen_t"]]
    assert all(r["cp_lb"] < 0.95 for r in larger)


def test_disagreement_rescue():
    """With disagreement rescue on, a veto-eligible candidate whose predictors disagree is kept."""
    df = _pool(n_neg=30, n_pos=3, seed=8)
    df.loc[2, ["el", "prime"]] = [1000.0, 1000.0]
    df.loc[2, "expr"] = -1000.0
    dis = np.zeros(len(df))
    dis[2] = 0.9  # high disagreement on that row
    cfg = GateConfig(t=0.9, disagreement_rescue=0.5, el_floor_frac=0.0, per_patient_cap=1.0,
                     pool_floor=0, cov_floor=0.0)
    g = apply_gate(df, cfg, disagreement=dis)
    assert bool(g.loc[2, "dyn_gate_keep"]) is True
    assert g.loc[2, "dyn_gate_reason"] == "RESCUE_PREDICTOR_DISAGREEMENT"


def test_predictor_disagreement_needs_two():
    df = _pool(n_neg=10, n_pos=2)
    df["pred_a"] = df["el"]
    out = predictor_disagreement(df, ["pred_a"])  # only one predictor -> all NaN
    assert np.isnan(out).all()
    df["pred_b"] = df["prime"]
    out2 = predictor_disagreement(df, ["pred_a", "pred_b"])
    assert (~np.isnan(out2)).all()
    assert np.all(out2 >= 0)


def test_attach_percentiles_preserves_missing():
    df = _pool(n_neg=10, n_pos=2)
    df.loc[0, "expr"] = np.nan
    p = attach_percentiles(df)
    assert np.isnan(p.loc[0, "s_expr"])
    assert p["s_el"].notna().all()
