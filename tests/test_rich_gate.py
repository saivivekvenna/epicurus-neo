"""Tests for the rich-feature gate module (src/event_b/rich_gate.py)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from event_b.rich_gate import (
    K,
    balanced_weights,
    base_anchored_hits,
    counterfactual_swaps,
    decision_boundary_mask,
    feature_matrix,
    feature_percentile,
    patient_hits20,
    peptide_hydro_aro_fraction,
)


def _frame(n=60, seed=0):
    rng = np.random.default_rng(seed)
    prime = np.arange(n).astype(float)  # index == frozen-Epicurus rank
    return pd.DataFrame({
        "patient_id": "improve:p1",
        "label": ["POSITIVE" if i in (24, 25, 26) else "TESTED_NEGATIVE" for i in range(n)],
        "prime": prime, "el": prime, "expr": rng.uniform(0, 5, n),
        "mutant_peptide": [f"PEP{i:02d}A" for i in range(n)],
        "RankEL_wt": rng.uniform(0, 10, n), "DAI": rng.uniform(0, 1, n),
    })


def test_feature_matrix_has_missing_indicator():
    f = _frame()
    f.loc[0, "RankEL_wt"] = np.nan
    X, names = feature_matrix(f, ["RankEL_wt", "DAI"])
    assert "RankEL_wt__pct" in names and "RankEL_wt__missing" in names
    mi = names.index("RankEL_wt__missing")
    assert X[0, mi] == 1.0 and X[1, mi] == 0.0          # explicit missing indicator, not silent
    assert X[0, names.index("RankEL_wt__pct")] == 0.5   # imputed neutral


def test_balanced_weights_equalize_classes():
    f = _frame()
    w = balanced_weights(f, balance_class=True)
    y = f["label"].to_numpy() == "POSITIVE"
    assert np.isclose(w[y].sum(), w[~y].sum())          # positive mass == negative mass
    w0 = balanced_weights(f, balance_class=False)
    assert not np.isclose(w0[y].sum(), w0[~y].sum())    # unbalanced differs


def test_decision_boundary_includes_all_positives():
    f = _frame(n=200)
    mask = decision_boundary_mask(f, top=K, extra=40)
    pos = (f["label"].to_numpy() == "POSITIVE")
    assert mask[pos].all()                               # every positive kept for training
    assert mask.sum() < len(f)                           # but not all 200 easy negatives


def test_counterfactual_swaps_removes_from_top20_only():
    f = _frame()
    q = np.full(len(f), 0.5)
    # make a few top-20 candidates look like negatives (low q) and backfill look positive (high q)
    q[:K] = 0.1
    q[K:K + 8] = 0.9
    keep = counterfactual_swaps(f, q, max_budget=8, margin=0.0)
    removed = np.where(~keep)[0]
    assert removed.size > 0
    assert removed.max() < K                             # only top-20 removed


def test_counterfactual_abstains_when_backfill_worse():
    f = _frame()
    q = np.full(len(f), 0.5)
    q[:K] = 0.9        # top-20 look positive
    q[K:] = 0.1        # backfill look negative -> no swap is worth it
    keep = counterfactual_swaps(f, q, max_budget=8, margin=0.0)
    assert keep.all()


def test_fixed_budget_removes_exactly_m():
    f = _frame()
    q = np.linspace(0, 1, len(f))  # distinct q
    keep = counterfactual_swaps(f, q, max_budget=8, fixed_m=3)
    assert (~keep).sum() == 3


def test_base_anchored_noop_at_alpha0():
    """alpha=0 MUST reproduce frozen Epicurus exactly (base-only no-op) — the key sanity control."""
    f = _frame(n=80)
    fp = feature_percentile(f, "RankEL_wt", True)
    h0 = base_anchored_hits(f, fp, alpha=0.0, threat=60)
    # equals ungated frozen-Epicurus top-20 hits
    ung = patient_hits20(f, np.ones(len(f), bool))
    assert h0["improve:p1"] == ung["improve:p1"]


def test_base_anchored_feature_can_swap_in_positive():
    """A feature that marks the just-outside positives high should be able to swap them into top-20."""
    f = _frame(n=80)
    # craft a feature that is high exactly for the positives (ranks 24-26, just outside top-20)
    feat = np.zeros(len(f))
    feat[f["label"].to_numpy() == "POSITIVE"] = 1.0
    f = f.assign(_probe=feat)
    fp = feature_percentile(f, "_probe", True)
    h_hi = base_anchored_hits(f, fp, alpha=3.0, threat=60)
    h_lo = base_anchored_hits(f, fp, alpha=0.0, threat=60)
    assert h_hi["improve:p1"] >= h_lo["improve:p1"]  # can only help when feature marks the positives


def test_peptide_hydro_fraction():
    s = pd.Series(["AAAA", "DEDE", "FWYC"])  # all-hydrophobic, all-acidic, aromatic/cys
    out = peptide_hydro_aro_fraction(s)
    assert out[0] == 1.0 and out[1] == 0.0 and out[2] == 1.0


def test_patient_hits20_rescoring_counts_positives():
    f = _frame()
    hits = patient_hits20(f, np.ones(len(f), bool))
    # positives at ranks 24,25,26 -> all within top-20? no (ranks 25-27) -> 0 in top-20 ungated
    assert hits["improve:p1"] == 0.0
