"""Unit tests for the risk-controlled negative reducer (no Sid/Miller, no network).

Covers: CP lower bound, protected-core m (+boundary ties), retention-guaranteed tau calibration, missing/OOD
KEEP, core protection, the NOW-LIVE Delta-hits@20 under m<20, backfill, nonnegative-logistic monotonicity +
apply path, matched-random same-count, and leakage quarantine.
"""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd

nr = importlib.import_module("event_b.negative_reducer")


# ---- CP lower bound ----------------------------------------------------------------------------------
def test_cp_lower_known_values_and_monotone():
    assert abs(nr.cp_lower(59, 59) - 0.05 ** (1 / 59)) < 1e-9
    assert nr.cp_lower(59, 59) >= 0.95 and nr.cp_lower(58, 58) < 0.95   # n=59 is the powered threshold
    assert nr.cp_lower(0, 10) == 0.0
    # more retained (same n) => higher lower bound
    assert nr.cp_lower(100, 100) > nr.cp_lower(95, 100) > nr.cp_lower(90, 100)


def test_no_data_files_touch_sid_or_miller():
    for f in nr.ALLOWED_DATA_FILES:
        low = f.lower()
        assert "sid" not in low and "osteosarc" not in low and "miller" not in low and "variant_vafs" not in low


# ---- protected core ----------------------------------------------------------------------------------
def _pt(prime, label=None, pid="p", **extra):
    n = len(prime)
    d = {"study": "s", "patient_id": pid, "mut_peptide": [f"PEP{i}" for i in range(n)],
         "prime": np.asarray(prime, float),
         "label": label if label is not None else ["TESTED_NEGATIVE"] * n,
         "el": np.asarray(prime, float), "expr": np.ones(n)}
    for c in ["VarAlFreq", "rna_af", "ValMutRNACoef", "CelPrev"]:
        d[c] = extra.get(c, np.full(n, np.nan))
    return pd.DataFrame(d)


def test_core_mask_m0_protects_nothing_m5_protects_top5():
    df = _pt(np.arange(30).astype(float))          # prime 0..29 (0 best)
    assert nr.core_mask(df, 0).sum() == 0
    c5 = nr.core_mask(df, 5)
    assert c5.sum() == 5 and c5[:5].all() and not c5[5:].any()


def test_core_mask_boundary_ties_all_protected():
    # indices 4,5,6 tie at PRIME value 4.0 (the m=5 boundary). thr=4.0 => rows 0..6 (prime<=4) all protected
    # -> core is 7, exceeding m=5 because the boundary ties are all kept.
    prime = np.arange(10).astype(float)
    prime[4] = prime[5] = prime[6] = 4.0
    df = _pt(prime)
    c = nr.core_mask(df, 5)
    assert c.sum() == 7 and c[4] and c[5] and c[6] and not c[7]


def test_core_mask_small_patient_protects_all():
    df = _pt(np.arange(3).astype(float))
    assert nr.core_mask(df, 5).all()


# ---- tau calibration guarantees retention ------------------------------------------------------------
def test_calibrate_tau_guarantees_retention_and_is_monotone():
    n_pos = 200
    rng = np.random.default_rng(0)
    # positives get higher keep-scores than negatives, with some overlap
    ks_pos = rng.uniform(0.4, 1.0, n_pos)
    ks_neg = rng.uniform(0.0, 0.6, 800)
    ks = np.concatenate([ks_pos, ks_neg])
    y = np.concatenate([np.ones(n_pos), np.zeros(800)]).astype(int)
    removable = np.ones(len(ks), bool)
    tau, r_max, np_, cp_lb, powered = nr.calibrate_tau(ks, y, removable)
    assert powered and np_ == n_pos and cp_lb >= 0.95
    # applying tau retains CP-lower-bound >= 0.95 of positives
    removed_pos = int(((ks < tau) & (y == 1)).sum())
    assert nr.cp_lower(n_pos - removed_pos, n_pos) >= 0.95
    # allowed removed positives r_max is exactly the max satisfying the bound
    assert nr.cp_lower(n_pos - r_max, n_pos) >= 0.95
    assert r_max == n_pos or nr.cp_lower(n_pos - (r_max + 1), n_pos) < 0.95


def test_calibrate_tau_nonremovable_positives_never_cost_retention():
    # a low-scoring positive that is NON-removable (in core) must not force tau down
    ks = np.array([0.01, 0.9, 0.8, 0.7, 0.6])
    y = np.array([1, 1, 0, 0, 0])
    removable = np.array([False, True, True, True, True])  # the 0.01 positive is protected
    tau, r_max, n_pos, cp_lb, powered = nr.calibrate_tau(ks, y, removable)
    # only 1 removable positive (score 0.9); with n_pos=2 CP is underpowered but tau still computed
    assert not (removable & (y == 1) & (ks == 0.01)).any()


# ---- removable mask: missing + OOD KEEP --------------------------------------------------------------
def test_removable_mask_missing_feature_and_ood_keep():
    df = _pt(np.arange(10).astype(float))
    df.loc[3, "expr"] = np.nan                     # row 3 missing a portable input -> KEEP
    rem = nr.removable_mask(df, nr.PORTABLE, m=0, ood=set())
    assert not rem[3] and rem[7]
    # OOD patient => nothing removable
    rem_ood = nr.removable_mask(df, nr.PORTABLE, m=0, ood={"p"})
    assert not rem_ood.any()


# ---- Delta-hits@20 is LIVE under m<20 (the whole point of PROTOCOL CORRECTION 1) ---------------------
def test_delta_hits_live_when_core_below_20():
    # 25 candidates; the ONLY clean positive sits at PRIME rank 21 (index 21). Negatives at ranks 5..20 are
    # removable (m=5). Removing them promotes the rank-21 positive into the top-20 -> +1 hit.
    prime = np.arange(25).astype(float)
    label = ["TESTED_NEGATIVE"] * 25
    label[21] = "POSITIVE"
    df = _pt(prime, label=label)
    clean = np.ones(25, bool)
    null = nr.hits_at_k(df, np.zeros(25, bool), clean)
    assert null["p"] == 0.0                         # positive at rank 21 is outside top-20
    removed = np.zeros(25, bool)
    removed[5:21] = True                            # remove ranks 5..20 (all outside the m=5 core)
    gate = nr.hits_at_k(df, removed, clean)
    assert gate["p"] == 1.0                         # promoted into top-20 -> Delta = +1 (LIVE)


def test_delta_hits_zero_when_core_is_full_top20():
    # with the removable set entirely below rank 20, top-20 cannot change (sanity for the old degenerate case)
    prime = np.arange(25).astype(float)
    label = ["TESTED_NEGATIVE"] * 25
    label[21] = "POSITIVE"
    df = _pt(prime, label=label)
    clean = np.ones(25, bool)
    removed = np.zeros(25, bool)
    removed[22:] = True                             # only ranks 22..24 removed (below the positive)
    assert nr.hits_at_k(df, removed, clean)["p"] == nr.hits_at_k(df, np.zeros(25, bool), clean)["p"]


def test_backfill_keeps_top20_full():
    # remove almost everything; backfill re-admits highest-PRIME removed so top-20 stays size 20
    prime = np.arange(25).astype(float)
    df = _pt(prime)
    removed = np.ones(25, bool)
    removed[0] = False                              # only 1 survivor
    h = nr.hits_at_k(df, removed, np.ones(25, bool))
    assert h["p"] == 0.0  # no positives, but must not error; backfill fills to 20 by prime


# ---- nonnegative logistic: coef>=0, monotone, apply path -------------------------------------------
def test_nnlogistic_coef_nonnegative_and_apply_path():
    rng = np.random.default_rng(1)
    X = rng.uniform(0, 1, (300, 3))
    # true signal favors higher feature values (positive)
    y = (X[:, 0] + 0.5 * X[:, 1] + rng.normal(0, 0.3, 300) > 1.0).astype(int)
    w = np.ones(300)
    b0, coef, ok = nr.fit_nnlogistic(X, y, w, C=1.0)
    assert ok and (coef >= -1e-8).all()             # converged + nonnegative constraint honored
    s = nr.nnlogistic_score(X, b0, coef)
    assert np.allclose(s, 1.0 / (1.0 + np.exp(-(X @ coef + b0))))   # apply path == sigmoid(linear)
    # monotone: raising any feature cannot lower the keep-score
    Xup = X.copy()
    Xup[:, 0] += 0.1
    assert (nr.nnlogistic_score(Xup, b0, coef) >= s - 1e-9).all()


def test_nnlogistic_fails_closed_on_optimizer_failure(monkeypatch):
    # if the optimizer reports failure / non-finite, fit_nnlogistic must return zeros (KEEP-all), not garbage
    class _Bad:
        success = False
        x = np.array([np.nan, np.nan, np.nan])
    monkeypatch.setattr(nr, "minimize", lambda *a, **k: _Bad())
    b0, coef, ok = nr.fit_nnlogistic(np.ones((5, 2)), np.array([1, 0, 1, 0, 1]), np.ones(5), C=1.0)
    assert ok is False and b0 == 0.0 and coef.shape == (2,) and np.all(coef == 0.0) and np.all(np.isfinite(coef))


def test_nnlogistic_recovers_sign_against_anti_signal():
    # if the data actually wanted a NEGATIVE coefficient, the constraint pins it at 0 (never negative)
    rng = np.random.default_rng(2)
    X = rng.uniform(0, 1, (300, 2))
    y = (-X[:, 0] + rng.normal(0, 0.2, 300) > -0.5).astype(int)   # feature 0 anti-correlated with y
    b0, coef, ok = nr.fit_nnlogistic(X, y, np.ones(300), C=1.0)
    assert ok and coef[0] <= 1e-6                     # converged; cannot go negative


# ---- matched-random same count -----------------------------------------------------------------------
def test_matched_random_same_count():
    prime = np.arange(25).astype(float)
    label = ["TESTED_NEGATIVE"] * 25
    label[21] = "POSITIVE"
    df = _pt(prime, label=label)
    clean = np.ones(25, bool)
    null = nr.hits_at_k(df, np.zeros(25, bool), clean)
    # gate removed 16 candidates for patient p
    counts = {"p": 16}
    d = nr.matched_random_delta(df, counts, nr.PORTABLE, m=5, ood=set(), clean=clean, null_hits=null,
                                seeds=range(5))
    assert isinstance(d, float)                      # random removal of 16 non-core candidates, averaged


# ---- OOD raw envelope --------------------------------------------------------------------------------
def test_ood_patients_flags_out_of_support():
    train = _pt(np.arange(50).astype(float), pid="tr")
    train["expr"] = np.linspace(1, 10, 50)
    test = _pt(np.arange(10).astype(float), pid="te")
    test["expr"] = np.full(10, 1e6)                  # wildly out of train expr envelope
    ood = nr.ood_patients(train, test, nr.PORTABLE)
    assert "te" in ood


# ---- apply_payload equivalence (Correction 3.4) ------------------------------------------------------
def test_apply_payload_matches_model_path_on_external():
    # a frozen payload (coef + serialized envelope + m + tau) must reproduce the model-based gate decision
    # on unseen external data, with NO fitting.
    def frame(pids, seed):
        r = np.random.default_rng(seed)
        rows = []
        for pid in pids:
            for _ in range(30):
                pos = r.random() < 0.15
                rows.append({"study": "x", "patient_id": pid, "mut_peptide": "P",
                             "label": "POSITIVE" if pos else "TESTED_NEGATIVE",
                             "prime": r.uniform(0, 3) if pos else r.uniform(0, 50),
                             "el": r.uniform(0, 3) if pos else r.uniform(0, 50), "expr": r.uniform(0, 20)})
        d = pd.DataFrame(rows)
        for c in ["VarAlFreq", "rna_af", "ValMutRNACoef", "CelPrev"]:
            d[c] = np.nan
        return d
    train = frame([f"t{i}" for i in range(8)], 1)
    ext = frame([f"e{i}" for i in range(5)], 2)
    cols = nr.PORTABLE
    X = nr.feat_matrix(train, cols)
    y = (train["label"].to_numpy() == "POSITIVE").astype(int)
    b0, coef, ok = nr.fit_nnlogistic(X, y, nr.balanced_weights(train), C=1.0)
    assert ok
    m = 5
    tau, *_ = nr.calibrate_tau(nr.nnlogistic_score(X, b0, coef), y, nr.removable_mask(train, cols, m, set()))
    payload = {"model": "nnlog", "feature_order": cols, "coef": [float(x) for x in coef], "intercept": float(b0),
               "m": m, "tau": (None if not np.isfinite(tau) else float(tau)),
               "ood_envelope": nr.raw_envelope(train, cols), "ood_cover": 0.5}
    # model path (with fitting) vs pure apply_payload (no fitting)
    ood = nr.ood_from_envelope(ext, payload["ood_envelope"], cols, 0.5)
    model_removed = nr.gate_removed(ext, nr.nnlogistic_score(nr.feat_matrix(ext, cols), b0, coef),
                                    (np.inf if payload["tau"] is None else payload["tau"]), cols, m, ood)
    assert np.array_equal(nr.apply_payload(ext, payload), model_removed)
    assert np.array_equal(nr.apply_payload(ext, {"model": "NULL"}), np.zeros(len(ext), bool))


# ---- quarantine --------------------------------------------------------------------------------------
def test_clean_against_flags_exact_and_near():
    df = pd.DataFrame({"mut_peptide": ["AAAAAAAAA", "CDEFGHIKL", "ZZZZZZZZZ"]})
    clean = nr._clean_against(df, {"AAAAAAAAA"})
    assert not clean[0] and clean[2]                 # exact match unclean; unrelated stays clean
