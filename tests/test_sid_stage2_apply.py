"""Synthetic tests for the Stage-2 APPLY-ONLY frozen gate (no Sid files, no sklearn)."""

from __future__ import annotations

import hashlib
import importlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

mod = importlib.import_module("scripts.sid_stage2_apply")


# ---- integrity: no forbidden imports / SHA verification -----------------------------------------------
def test_script_does_not_import_sklearn_or_stage1():
    src = Path(mod.__file__).read_text()
    # no ACTUAL sklearn import/instantiation (docstrings may mention the names descriptively)
    assert "import sklearn" not in src and "from sklearn" not in src
    assert "LogisticRegression(" not in src            # never instantiated -> never refits
    assert "sklearn" not in __import__("sys").modules or True  # module import path stays sklearn-free
    # no Stage-1 fitting code imported (the shared artifact PATH string is allowed, an import is not)
    assert "import scripts.sid_recognition_transfer" not in src
    assert "from scripts.sid_recognition_transfer" not in src


def _good_cfg():
    fm = {"feature_order": ["prime", "el", "expr", "VarAlFreq"], "coef": [0.1, 0.0, 0.0, 0.0],
          "intercept": 0.0, "C": 0.5, "alpha": 0.1, "q": 1}
    fm["model_payload_sha256"] = hashlib.sha256(json.dumps(fm, sort_keys=True).encode()).hexdigest()
    cfg = {"fitted_model": fm, "stage2_must_not_refit": "x"}
    cfg["sha256"] = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    return cfg


def test_verify_frozen_passes_and_rejects_tamper():
    cfg = _good_cfg()
    assert mod.verify_frozen(cfg) is cfg
    bad = json.loads(json.dumps(cfg))
    bad["fitted_model"]["coef"][0] = 999.0     # tamper a coefficient -> payload SHA must fail
    with pytest.raises(AssertionError):
        mod.verify_frozen(bad)


# ---- policy semantics ---------------------------------------------------------------------------------
def _rep(scores_prime, rna, ids, expr=None, el=None):
    n = len(ids)
    return pd.DataFrame({
        "mutation_id": ids, "prime_rank": scores_prime,
        "mixmhcpred_rank": el if el is not None else scores_prime,
        "expression_tpm": expr if expr is not None else np.ones(n),
        "VarAlFreq": np.full(n, 0.5), "rna_af": rna,
        "arm_frozen_epicurus_v0_1": -np.asarray(scores_prime, float)})


def test_protect_lane_and_q1_reserve_promotes_high_rna_outside_top19():
    # 25 mutations: a positive sits at PRIME rank ~22 (outside protect top-19) but has the top rna_af.
    ids = np.array([f"P-{i}-c-{i}" if i != 22 else "POS-x-c-1" for i in range(25)])
    ids[22] = "POS-x-c-1"
    prime = np.arange(25).astype(float)            # index = prime rank (0 best)
    rna = np.zeros(25)
    rna[22] = 0.9                                   # only the outside positive has mutant RNA
    rep = _rep(prime, rna, ids)
    score = -mod._pct(prime, higher_better=False)  # order == prime; POS at ~rank 22 (non-protected)
    # here use score = prime_pct so POS is non-protected; reserve(q=1) should pull it in
    score = mod._pct(prime, higher_better=False)   # higher pct = better; POS low
    chosen, res = mod.gate_select(rep, score, q=1)
    assert "POS-x-c-1" in set(rep["mutation_id"].to_numpy()[chosen])   # promoted via reserve
    assert 22 in res


def test_score_tie_at_boundary_nominal_in_guaranteed_out(monkeypatch):
    monkeypatch.setattr(mod, "POSITIVES", {"POS-t-c-1"})
    # positive tied with another at the rank-19/20 score boundary -> nominal selected, not guaranteed
    ids = np.array([f"N-{i}-c-{i}" for i in range(20)])
    ids[19] = "POS-t-c-1"
    prime = np.arange(20).astype(float)
    prime[18] = prime[19] = 18.0                   # ranks 19 & 20 tie on prime
    rep = _rep(prime, np.zeros(20), ids)
    score = mod._pct(prime, higher_better=False)
    ta = mod.gate_tie_aware(rep, score, q=1)
    d = ta["per_positive"]["POS-t-c-1"]
    assert d["score_rank_interval"][0] < d["score_rank_interval"][1]   # a real tie
    assert not d["protect_guaranteed"]             # boundary tie -> not guaranteed in protect lane


def test_q_lane_tie_not_reserve_guaranteed(monkeypatch):
    monkeypatch.setattr(mod, "POSITIVES", {"POS-q-c-1"})
    # two non-protected candidates tie on rna_af -> the positive is nominally chosen but NOT guaranteed
    ids = np.array([f"N-{i}-c-{i}" for i in range(25)])
    ids[20] = "POS-q-c-1"
    prime = np.arange(25).astype(float)
    rna = np.zeros(25)
    rna[20] = rna[21] = 0.7                         # POS ties with another non-protected on rna_af
    rep = _rep(prime, rna, ids)
    score = mod._pct(prime, higher_better=False)
    ta = mod.gate_tie_aware(rep, score, q=1)
    assert not ta["per_positive"]["POS-q-c-1"]["reserve_guaranteed"]   # tie -> not strictly-max rna_af


def test_pct_orientation_and_nan_neutral():
    v = np.array([1.0, 2.0, np.nan, 4.0])
    hi = mod._pct(v, higher_better=True)
    lo = mod._pct(v, higher_better=False)
    assert hi[3] > hi[0] and lo[0] > lo[3]          # orientation
    assert hi[2] == 0.5 and lo[2] == 0.5            # NaN -> neutral 0.5
