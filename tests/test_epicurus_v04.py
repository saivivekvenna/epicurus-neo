"""Tests for the Epicurus v0.4 DEVELOPMENT experiment (source-aware tower).

Covers the review-required safeguards on synthetic frames (no PRIME executable needed): the λ=∞ pooled branch is
numerically identical to frozen v0.3; source-label-order invariance; equal source contribution under unequal
counts; unseen-source fallback = shared backbone; multi-init stability; recovery of an opposite per-source
feature weight that pooling cannot fit; preprocessing fit inside train only; per-source intercepts are
rank-inert; label-blind eligibility. Real-data assembly tests are guarded on the caches/split being present.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from event_b.epicurus_v03 import FEATURES, MILRanker, per_patient_metrics, run_oof
from event_b.epicurus_v04 import (
    TowerMILRanker, attrition_report, ext_metrics, rankable_patients,
)

FROZEN = Path("configs/frozen/mil_dev_split_v1.json")
G_PRIME = Path("data/raw/gartner_nci/_cache_gartner_muller_prime.tsv")
_READY = FROZEN.exists() and G_PRIME.exists()


def _synth(n_patients=9, per=30, seed=0, singleton=True):
    """v0.4-schema frame with a clean presentation signal; 3 sources round-robin."""
    rng = np.random.default_rng(seed)
    rows = []
    for pt in range(n_patients):
        src = ["gartner", "improve", "multimer"][pt % 3]
        for i in range(per):
            pos = i < 3
            strong = rng.uniform(2.5, 4.0) if pos else rng.uniform(0.0, 1.5)
            bag = f"{src}:{pt}#b{i}" if singleton else f"{src}:{pt}#B{i // 5}"
            rows.append(dict(
                source=src, patient_id=f"{src}:{pt}", bag_id=bag,
                bag_label="POSITIVE" if pos else "NEGATIVE", eval_positive=int(pos),
                prime_rank=rng.uniform(0, 2), mix_rank=rng.uniform(0, 2), el_strength=strong,  # raw core scores
                f_prime_pct=rng.uniform(0, 1), f_mix_pct=rng.uniform(0, 1),
                f_el_pct=(0.5 + strong / 8) if pos else rng.uniform(0, 0.5), f_pres_abs=strong,
                f_expr=0.0, f_agreto=0.0, f_foreign=0.0, f_bindstab=0.0, f_proc=0.0,
                fold=pt % 5, quarantined=False,
            ))
    return pd.DataFrame(rows)


def _theta_size(member, S=3, n=9, finite_lambda=True):
    has_c = member in ("C", "F")
    has_v = member == "F" and finite_lambda
    return n + 1 + (S if has_c else 0) + (S * n if has_v else 0)


# 1) λ=∞ pooled branch is numerically identical to frozen v0.3 MILRanker
def test_pooled_branch_equals_v03_mil():
    f = _synth(seed=1)
    for C, tau in [(0.1, 1.0), (0.3, 0.5), (1.0, 1.0)]:
        tp = TowerMILRanker(member="P", C=C, tau=tau).fit(f)
        mv = MILRanker(C=C, tau=tau).fit(f)
        assert np.allclose(tp.w0_, mv.coef_, atol=1e-8)
        assert abs(tp.b0_ - mv.intercept_) < 1e-8
        assert np.max(np.abs(tp.raw_score(f) - mv.raw_score(f))) < 1e-8


def test_feature_tower_lambda_inf_reduces_to_calibration():
    f = _synth(seed=2)
    fi = TowerMILRanker(member="F", C=0.3, tau=1.0, lam=np.inf).fit(f)
    c = TowerMILRanker(member="C", C=0.3, tau=1.0).fit(f)
    assert np.max(np.abs(fi.raw_score(f) - c.raw_score(f))) < 1e-8


# 2) source-label-ordering invariance: shuffling the input row order cannot change fitted predictions
def test_source_label_order_invariance():
    f = _synth(seed=3)
    a = TowerMILRanker(member="F", C=0.3, tau=1.0, lam=1.0).fit(f)
    shuffled = f.sample(frac=1.0, random_state=7).reset_index(drop=True)
    b = TowerMILRanker(member="F", C=0.3, tau=1.0, lam=1.0).fit(shuffled)
    assert a.sources_ == b.sources_                       # deterministic (sorted) source encoding
    assert np.max(np.abs(a.raw_score(f) - b.raw_score(f))) < 1e-6


# 3) equal source contribution despite unequal patient/candidate counts
def test_equal_source_contribution_in_objective():
    # gartner: 1 patient, 1 big negative bag (10 children) + 2 small; improve: 3 patients; multimer: 2 patients
    rows = []
    rows += [dict(source="gartner", patient_id="gartner:1", bag_id="gartner:1#big", bag_label="NEGATIVE")] * 10
    rows += [dict(source="gartner", patient_id="gartner:1", bag_id="gartner:1#s1", bag_label="POSITIVE")]
    rows += [dict(source="gartner", patient_id="gartner:1", bag_id="gartner:1#s2", bag_label="NEGATIVE")]
    for p in range(3):
        rows += [dict(source="improve", patient_id=f"improve:{p}", bag_id=f"improve:{p}#b", bag_label="POSITIVE")]
    for p in range(2):
        rows += [dict(source="multimer", patient_id=f"multimer:{p}", bag_id=f"multimer:{p}#b",
                      bag_label="NEGATIVE")]
    lab = pd.DataFrame(rows).reset_index(drop=True)
    codes = pd.factorize(lab["bag_id"])[0]
    order = np.argsort(codes, kind="stable")
    lab_s = lab.iloc[order].reset_index(drop=True)
    codes_s = codes[order]
    boundaries = np.concatenate([[0], np.flatnonzero(np.diff(codes_s)) + 1])
    _, wb = MILRanker._bag_targets(lab_s, boundaries)
    first = lab_s.iloc[boundaries]
    per_src = pd.Series(wb, index=first["source"].to_numpy()).groupby(level=0).sum()
    for s, tot in per_src.items():
        assert np.isclose(tot, 1.0 / 3)                   # each of 3 sources contributes total weight 1/S


# 4) unseen-source fallback uses only the shared backbone
def test_unseen_source_falls_back_to_backbone():
    f = _synth(seed=4)
    m = TowerMILRanker(member="F", C=0.3, tau=1.0, lam=1.0).fit(f)
    g = f.copy(); g["source"] = "BRAND_NEW_COHORT"
    X = m._std(g[FEATURES].to_numpy(float), fit=False)
    pooled = X @ m.w0_ + m.b0_
    assert np.max(np.abs(m.raw_score(g) - pooled)) < 1e-12


# 5) multi-init stability: fixed-seed perturbed inits converge to materially identical predictions
def test_multi_init_stability():
    f = _synth(seed=5)
    base = TowerMILRanker(member="F", C=0.3, tau=1.0, lam=1.0).fit(f)
    s0 = base.raw_score(f)
    size = _theta_size("F")
    for seed in (11, 23, 37):
        rng = np.random.default_rng(seed)
        m = TowerMILRanker(member="F", C=0.3, tau=1.0, lam=1.0).fit(f, init=rng.normal(0, 0.05, size))
        s = m.raw_score(f)
        rho = pd.Series(s0).corr(pd.Series(s), method="spearman")
        assert rho >= 0.999
        # compare on standardized scores (scale-free)
        z0 = (s0 - s0.mean()) / (s0.std() + 1e-9); z = (s - s.mean()) / (s.std() + 1e-9)
        assert np.max(np.abs(z0 - z)) <= 1e-3


# 6) recovers an opposite per-source feature weight that a single pooled weight cannot fit
def test_recovers_opposite_per_source_weight():
    rng = np.random.default_rng(6)
    rows = []
    # two sources; within each patient the positive is flagged by f_expr but with OPPOSITE sign per source
    for src, sign in [("gartner", +1.0), ("improve", -1.0)]:
        for pt in range(6):
            for i in range(20):
                pos = i < 4
                expr = sign * (1.0 if pos else -1.0) + rng.normal(0, 0.05)
                rows.append(dict(source=src, patient_id=f"{src}:{pt}", bag_id=f"{src}:{pt}#b{i}",
                    bag_label="POSITIVE" if pos else "NEGATIVE", eval_positive=int(pos),
                    f_prime_pct=0.0, f_mix_pct=0.0, f_el_pct=0.0, f_pres_abs=0.0,
                    f_expr=expr, f_agreto=0.0, f_foreign=0.0, f_bindstab=0.0, f_proc=0.0,
                    fold=pt % 5, quarantined=False))
    f = pd.DataFrame(rows)
    j = FEATURES.index("f_expr")
    P = TowerMILRanker(member="P", C=1.0, tau=1.0).fit(f)
    F = TowerMILRanker(member="F", C=1.0, tau=1.0, lam=0.03).fit(f)
    wg = F.w0_[j] + F.v_["gartner"][j]; wi = F.w0_[j] + F.v_["improve"][j]
    assert np.sign(wg) != np.sign(wi)                     # opposite effective weights recovered

    def min_sep(model):
        s = model.raw_score(f); f2 = f.assign(_s=s)
        seps = []
        for src, g in f2.groupby("source"):
            ep = g["eval_positive"].to_numpy() == 1
            seps.append(g.loc[ep, "_s"].mean() - g.loc[~ep, "_s"].mean())
        return min(seps)
    assert min_sep(F) > 0.5                               # tower separates BOTH sources
    assert min_sep(P) < 0.1                               # pooled cannot (opposite signs cancel)


# 7) all fitted preprocessing (standardization) is a function of the TRAIN rows only
def test_preprocessing_fit_inside_train_only():
    f = _synth(seed=8)
    m = TowerMILRanker(member="F", C=0.3, tau=1.0, lam=1.0).fit(f)
    mean0, std0 = m.mean_.copy(), m.std_.copy()
    m2 = TowerMILRanker(member="F", C=0.3, tau=1.0, lam=1.0)
    # scoring a totally different held-out frame must not alter the fitted scaler
    held = _synth(seed=999)
    m.raw_score(held)
    assert np.array_equal(m.mean_, mean0) and np.array_equal(m.std_, std0)
    # scaler is deterministic from the labeled train rows
    m2.fit(f)
    assert np.allclose(m.mean_, m2.mean_) and np.allclose(m.std_, m2.std_)


# 8) per-source intercepts are rank-inert within patient (cannot masquerade as ranking lift)
def test_intercepts_are_rank_inert():
    f = _synth(seed=9)
    base = TowerMILRanker(member="C", C=0.3, tau=1.0).fit(f)
    bumped = TowerMILRanker(member="C", C=0.3, tau=1.0)
    bumped.__dict__.update({k: getattr(base, k) for k in
                            ("mean_", "std_", "w0_", "b0_", "v_", "sources_")})
    bumped.c_ = {s: base.c_[s] + (5.0 if s == "gartner" else -3.0) for s in base.sources_}
    a = per_patient_metrics(f, base.raw_score(f))
    b = per_patient_metrics(f, bumped.raw_score(f))
    m = a.merge(b, on=["source", "patient_id"], suffixes=("_a", "_b"))
    assert (m["hits_a"] == m["hits_b"]).all()             # intercept shift changes no within-patient ranking


# 9) OOF integrity: no patient is scored by a model that trained on it
def test_oof_never_scores_a_training_patient():
    f = _synth(n_patients=10, seed=3)
    oof = run_oof(f, lambda: TowerMILRanker(member="F", C=0.3, tau=1.0, lam=1.0))
    assert oof.metrics["patient_id"].is_unique


# 10) label-blind eligibility: flipping every label leaves the rankable set unchanged
def test_rankable_is_label_blind():
    f = _synth(seed=10)
    before = rankable_patients(f)
    flipped = f.copy()
    flipped["eval_positive"] = 1 - flipped["eval_positive"]
    flipped["bag_label"] = np.where(flipped["bag_label"] == "POSITIVE", "NEGATIVE", "POSITIVE")
    after = rankable_patients(flipped)
    assert before == after and len(before) == 9


def test_ext_metrics_best_rank_and_ndcg():
    f = _synth(seed=12)
    m = TowerMILRanker(member="F", C=0.3, tau=1.0, lam=1.0).fit(f)
    mt = ext_metrics(f, m.raw_score(f))
    assert {"hits", "recall", "best_pos_rank", "ndcg"} <= set(mt.columns)
    assert (mt["best_pos_rank"] >= 1).all() and (mt["ndcg"] <= 1.0 + 1e-9).all()


# ---- real-data guards ----
@pytest.mark.skipif(not _READY, reason="v0.4 caches/split not available")
def test_attrition_scored_matches_gate():
    from event_b.epicurus_v04 import assemble_frame
    a = attrition_report(assemble_frame())
    assert a["TOTAL"]["scored_has_positive"] == 118
    assert a["TOTAL"]["feature_bearing"] == 152


@pytest.mark.skipif(not _READY, reason="v0.4 caches/split not available")
def test_test_io_guard_blocks_gartner_test_path():
    from event_b.epicurus_v04 import guard_no_test_io
    with guard_no_test_io():
        with pytest.raises(RuntimeError):
            pd.read_csv("data/raw/gartner_nci/NmersTestingSet.txt", sep="\t")
