"""Tests for the v2 reselection policy (src/event_b/dynamic_gate_v2.py)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from event_b.dynamic_gate_v2 import (
    K,
    add_v2_features,
    counterfactual_reselect,
    fit_negative_risk,
    reselect,
    utility,
)


def _pool(n=60, seed=0):
    rng = np.random.default_rng(seed)
    prime = np.arange(n).astype(float)  # index 0 = best (lowest %rank) => frozen-Epicurus rank == index
    return pd.DataFrame({
        "patient_id": "p1",
        "label": ["POSITIVE" if i in (25, 26, 27) else "TESTED_NEGATIVE" for i in range(n)],
        "prime": prime, "el": prime + rng.normal(0, 0.1, n), "expr": rng.uniform(0, 10, n),
        "mutant_peptide": [f"PEPTIDE{i:02d}A" for i in range(n)],
    })


def test_risk_in_unit_interval_and_missing_is_low():
    f = _pool()
    m = fit_negative_risk(f)
    r = m.risk(f)
    assert np.all((r >= 0) & (r <= 1))
    f2 = f.copy()
    f2.loc[0, ["expr"]] = np.nan  # expr missing; discord/interaction still present -> still scored
    g = f2.copy()
    g[["el", "prime", "expr"]] = np.nan  # ALL orthogonal inputs unavailable -> risk forced to 0
    gr = m.risk(g)
    assert np.allclose(gr, 0.0)


def test_add_v2_features_discordance():
    f = add_v2_features(_pool())
    assert {"p_el", "p_prime", "p_expr", "discord", "expr_x_pres"} <= set(f.columns)
    assert np.all((f["discord"] >= 0) & (f["discord"] <= 1))


def test_reselect_budget_zero_keeps_all():
    f = _pool()
    m = fit_negative_risk(f)
    keep = reselect(f, m.risk(f), budget_frac=0.0)
    assert keep.all()


def test_reselect_removes_from_threat_zone_only():
    f = _pool()
    m = fit_negative_risk(f)
    keep = reselect(f, m.risk(f), budget_frac=0.2, threat_k=2 * K)
    removed = np.where(~keep)[0]
    # everything removed must be within the top-2K by frozen Epicurus (here = lowest indices, best prime)
    assert removed.max() < 2 * K
    assert (~keep).sum() > 0


class _StubEns:
    """Deterministic q with controllable bounds for testing the counterfactual decision logic."""
    def __init__(self, mean, lcb, ucb):
        self._m, self._l, self._u = mean, lcb, ucb

    def q_with_uncertainty(self, frame, z=1.28):
        return self._m, self._l, self._u


def _q_by_frozen_rank(f, top_q, backfill_q):
    """Build a q array aligned to the ACTUAL frozen-Epicurus ranking (not index order)."""
    from event_b.pool_size_sensitivity import score_arms
    s = score_arms(f).sort_values("frozen_epicurus", ascending=False, kind="mergesort")
    rank_of = {idx: r for r, idx in enumerate(s.index)}
    return np.array([top_q if rank_of[i] < K else backfill_q for i in f.index], float)


def test_counterfactual_abstains_when_no_gain():
    f = _pool()
    mean = _q_by_frozen_rank(f, top_q=0.9, backfill_q=0.1)  # replacements worse than removed -> no gain
    ens = _StubEns(mean, mean - 0.01, mean + 0.01)
    keep = counterfactual_reselect(f, ens, max_budget=8, conservative=True)
    assert keep.all()  # nothing removed


def test_counterfactual_acts_when_replacement_dominates():
    f = _pool()
    mean = _q_by_frozen_rank(f, top_q=0.1, backfill_q=0.9)  # top-20 look negative, backfill positive
    ens = _StubEns(mean, mean - 0.01, mean + 0.01)          # tight bands -> LCB(repl) > UCB(removed)
    keep = counterfactual_reselect(f, ens, max_budget=8, conservative=True)
    assert (~keep).sum() > 0  # it should swap


def test_counterfactual_wide_bands_force_abstention():
    f = _pool()
    mean = _q_by_frozen_rank(f, top_q=0.1, backfill_q=0.9)
    ens = _StubEns(mean, np.zeros(len(f)), np.ones(len(f)))  # wide bands -> conservative gate abstains
    keep = counterfactual_reselect(f, ens, max_budget=8, conservative=True)
    assert keep.all()


def test_utility_penalizes_harm():
    good = np.array([1.0, 0.0, 0.0, 0.0])
    harm = np.array([1.0, -1.0, 0.0, 0.0])
    assert utility(good, lam=1.0, mu=0.5) > utility(harm, lam=1.0, mu=0.5)
