"""Tests for the Epicurus v0.5 DEVELOPMENT experiment (context-conditioned pairwise challenger).

Exercises the FROZEN protocol's MATH on synthetic bags — no PRIME executable, no full benchmark, no Q/R fit on
the real corpus. Independent oracles where possible: a central-difference check of the analytic gradient (Q and
R), hand-built `_FitData` for the log-MEAN-exp bag properties, and a from-scratch pairwise-logistic recompute
for the singleton reduction. Also covers: no intercept + strict init-invariance, exact source/patient/positive/
negative-bag weight normalization, the robust HLA-locus parser and rowwise-SD disagreement context, φ-safety +
R(β=0)≡Q, train-only scaling, the deterministic selection tie-break, comparator hyperparameter extraction, and
the fail-closed provenance / Gartner-TEST guards. Real-frame checks are guarded on the caches/split and never fit.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from event_b.epicurus_v03 import FEATURES, AdditiveRanker, MILRanker
from event_b.epicurus_v04 import TowerMILRanker
import event_b.epicurus_v05 as v5
from event_b.epicurus_v05 import ContextPairwiseRanker, _FitData

FROZEN = Path("configs/frozen/mil_dev_split_v1.json")
G_PRIME = Path("data/raw/gartner_nci/_cache_gartner_muller_prime.tsv")
_READY = FROZEN.exists() and G_PRIME.exists()

_AA = list("ACDEFGHIKLMNPQRSTVWY")
_HLA = ["A*02:01", "B*07:02", "C*07:01", "A0201", "HLA-B08:01", "HLA-C*07:02"]


# --------------------------------------------------------------------------------------------------
# synthetic frames + hand-built fit scaffolding
# --------------------------------------------------------------------------------------------------
def _synth(seed=0, n_patients=6, per=8):
    """v0.5-schema frame: 3 sources round-robin; 2 exact positives/patient; Gartner negatives share 2-child
    bags (multi-child log-mean-exp), instance sources are singleton negatives. Carries peptide/HLA for contexts."""
    rng = np.random.default_rng(seed)
    rows = []
    for pt in range(n_patients):
        src = ["gartner", "improve", "multimer"][pt % 3]
        pid = f"{src}:{pt}"
        for i in range(per):
            pos = i < 2
            if src == "gartner" and not pos:
                bag = f"{pid}#B{i // 2}"          # negatives pair up into 2-child bags
            else:
                bag = f"{pid}#b{i}"               # singleton bag
            feats = {f: float(rng.normal()) for f in FEATURES}
            for cf in ("f_prime_pct", "f_mix_pct", "f_el_pct"):
                feats[cf] = float(rng.uniform(0, 1))
            pep = "".join(rng.choice(_AA, size=int(rng.integers(8, 12))))
            rows.append(dict(
                source=src, patient_id=pid, bag_id=bag,
                bag_label="POSITIVE" if pos else "NEGATIVE", eval_positive=int(pos),
                peptide=pep, hla_allele=str(rng.choice(_HLA)),
                fold=pt % 5, quarantined=False, **feats))
    return pd.DataFrame(rows)


def _make_fd(f, has_beta, lam_w=1.0, lam_ctx=10.0):
    """Build a `_FitData` from a frame WITHOUT running the optimizer (train-fit scalers, real design)."""
    lab = f[(f["eval_positive"] == 1) | (f["bag_label"] == "NEGATIVE")].copy()
    lab = v5.add_approved_contexts(lab).reset_index(drop=True)
    X = lab[FEATURES].to_numpy(float)
    fmean, fstd = X.mean(0), X.std(0) + 1e-9
    if has_beta:
        Z = lab[v5.CTX_COLS].to_numpy(float)
        cmean, cstd = Z.mean(0), Z.std(0) + 1e-9
    else:
        cmean = cstd = None
    phi = v5.build_design(lab, has_beta, fmean, fstd, cmean, cstd)
    pen = np.full(len(FEATURES), float(lam_w))
    if has_beta:
        pen = np.concatenate([pen, np.full(len(v5.interaction_names()), float(lam_ctx))])
    return v5._prepare_fit(lab, phi, pen), phi.shape[1]


def _fd_scores(pos_scores, bag_children, w, penalty_len=1):
    """Hand-built 1-feature `_FitData` (θ=[1] ⇒ s_j = φ_j): full control over child scores and pairing."""
    phi_pos = np.array(pos_scores, float).reshape(-1, 1)
    counts = np.array([len(b) for b in bag_children], int)
    boundaries = np.concatenate([[0], np.cumsum(counts)[:-1]]) if len(counts) else np.array([], int)
    phi_child = np.array([c for b in bag_children for c in b], float).reshape(-1, 1)
    pp, pb, pw = [], [], []
    for a in range(len(pos_scores)):
        for b in range(len(bag_children)):
            pp.append(a); pb.append(b); pw.append(w)
    n_pos, n_bag = len(pos_scores), len(bag_children)
    return _FitData(
        phi_pos=phi_pos, phi_child=phi_child, seg_boundaries=boundaries, seg_counts=counts,
        pair_pos=np.array(pp, int), pair_bag=np.array(pb, int), pair_w=np.array(pw, float),
        penalty=np.zeros(penalty_len), pos_patient=np.array(["p"] * n_pos),
        pos_source=np.array(["s"] * n_pos), bag_patient=np.array(["p"] * n_bag),
        bag_source=np.array(["s"] * n_bag))


# --------------------------------------------------------------------------------------------------
# 1) finite-difference check of the analytic gradient (Q and R)
# --------------------------------------------------------------------------------------------------
def test_finite_difference_gradient_Q_and_R():
    for has_beta in (False, True):
        fd, D = _make_fd(_synth(seed=3), has_beta)
        rng = np.random.default_rng(0)
        theta = rng.normal(0, 0.5, D)
        _, grad = v5._loss_grad(theta, fd)
        eps = 1e-6
        num = np.empty(D)
        for i in range(D):
            tp, tm = theta.copy(), theta.copy()
            tp[i] += eps; tm[i] -= eps
            num[i] = (v5._loss_grad(tp, fd)[0] - v5._loss_grad(tm, fd)[0]) / (2 * eps)
        assert np.max(np.abs(num - grad)) < 1e-5, (has_beta, np.max(np.abs(num - grad)))


# --------------------------------------------------------------------------------------------------
# 2) no intercept parameter + strict init-invariant solution (multi-init coefficient/rank equality)
# --------------------------------------------------------------------------------------------------
def test_no_intercept_and_strict_init_invariance():
    f = _synth(seed=4)
    base = ContextPairwiseRanker(member="R", lam_w=1.0, lam_ctx=10.0).fit(f)
    assert base.coef_.shape == (9,) and base.beta_.shape == (16,)
    assert not hasattr(base, "intercept_")
    # θ has exactly D = 9 + 16 slots — no bias term
    assert len(v5.design_column_names(True)) == 25
    base_theta = np.concatenate([base.coef_, base.beta_])
    s0 = base.raw_score(f)
    for seed in (1, 2, 3, 11, 23):
        rng = np.random.default_rng(seed)
        m = ContextPairwiseRanker(member="R", lam_w=1.0, lam_ctx=10.0).fit(f, init=rng.normal(0, 0.5, 25))
        # §9.6 hard gate: perturbed fixed-seed inits ⇒ Spearman = 1.0 AND max|Δcoef| ≤ 1e-6 (unique global min)
        assert np.max(np.abs(base_theta - np.concatenate([m.coef_, m.beta_]))) < 1e-6
        assert pd.Series(s0).corr(pd.Series(m.raw_score(f)), method="spearman") >= 1 - 1e-12
    q = ContextPairwiseRanker(member="Q", lam_w=1.0).fit(f)
    assert q.coef_.shape == (9,) and q.beta_ is None


# --------------------------------------------------------------------------------------------------
# 3) log-MEAN-exp bag semantics: sensitive to any child; invariant to duplicating children + bag weight
# --------------------------------------------------------------------------------------------------
def test_changing_any_child_changes_lme_loss_and_gradient():
    from dataclasses import replace
    fd = _fd_scores(pos_scores=[0.2], bag_children=[[0.1, 0.4, 0.7]], w=1.0)
    theta = np.array([0.5])
    l0, g0 = v5._loss_grad(theta, fd)
    bumped = fd.phi_child.copy(); bumped[1, 0] += 0.9                # change one child's score
    l1, g1 = v5._loss_grad(theta, replace(fd, phi_child=bumped))
    assert abs(l1 - l0) > 1e-8 and abs(g1[0] - g0[0]) > 1e-8


def test_duplicating_identical_children_leaves_logmeanexp_and_weight_unchanged():
    theta = np.array([0.5])
    single = _fd_scores([0.2], [[0.1, 0.6]], w=1.0)
    dup = _fd_scores([0.2], [[0.1, 0.6, 0.1, 0.6]], w=1.0)          # duplicate the two children
    assert abs(v5._loss_grad(theta, single)[0] - v5._loss_grad(theta, dup)[0]) < 1e-12
    # bag weight is per-BAG (no |b| term): duplicating children keeps B_p=1 and pair weight identical
    lab = pd.DataFrame({
        "source": ["gartner"] * 4, "patient_id": ["gartner:1"] * 4,
        "bag_id": ["gartner:1#pos", "gartner:1#B", "gartner:1#B", "gartner:1#s"],
        "bag_label": ["POSITIVE", "NEGATIVE", "NEGATIVE", "NEGATIVE"], "eval_positive": [1, 0, 0, 0]})
    lab_dup = pd.concat([lab, lab[lab["bag_id"] == "gartner:1#B"]], ignore_index=True)  # 2 -> 4 children in #B
    fd1 = v5._prepare_fit(lab, np.zeros((len(lab), 1)), np.zeros(1))
    fd2 = v5._prepare_fit(lab_dup, np.zeros((len(lab_dup), 1)), np.zeros(1))
    assert len(fd1.seg_counts) == len(fd2.seg_counts)              # same number of negative BAGS (2)
    assert np.allclose(np.sort(fd1.pair_w), np.sort(fd2.pair_w))   # pair weights unchanged by child count


# --------------------------------------------------------------------------------------------------
# 4) singleton bags reduce EXACTLY to ordinary pairwise logistic
# --------------------------------------------------------------------------------------------------
def test_singleton_bag_reduces_to_pairwise_logistic():
    pos = [0.3, -0.2]
    negs = [-0.1, 0.5, 0.9]                                        # each a singleton bag
    w = 0.25
    fd = _fd_scores(pos, [[n] for n in negs], w=w)
    theta = np.array([1.3])
    loss, grad = v5._loss_grad(theta, fd)
    # from-scratch pairwise logistic: L = Σ w·softplus(θ·n − θ·a); dL/dθ = Σ w·σ(·)·(n − a)
    exp_loss = exp_grad = 0.0
    for a in pos:
        for n in negs:
            u = theta[0] * n - theta[0] * a
            exp_loss += w * np.logaddexp(0.0, u)
            exp_grad += w * (1 / (1 + np.exp(-u))) * (n - a)
    assert abs(loss - exp_loss) < 1e-12 and abs(grad[0] - exp_grad) < 1e-12


# --------------------------------------------------------------------------------------------------
# 5) exact source / patient / positive / negative-bag weight totals (1/S per source; grand total 1)
# --------------------------------------------------------------------------------------------------
def test_exact_weight_normalization_totals():
    f = _synth(seed=5, n_patients=9, per=8)                        # 3 patients per source
    lab = f[(f["eval_positive"] == 1) | (f["bag_label"] == "NEGATIVE")].reset_index(drop=True)
    fd = v5._prepare_fit(lab, np.zeros((len(lab), 1)), np.zeros(1))
    src_of_pair = fd.pos_source[fd.pair_pos]
    sources = np.unique(src_of_pair)
    S = len(sources)
    total = 0.0
    for s in sources:
        tot_s = fd.pair_w[src_of_pair == s].sum()
        assert np.isclose(tot_s, 1.0 / S), (s, tot_s)             # each source contributes exactly 1/S
        total += tot_s
    assert np.isclose(total, 1.0)                                  # grand total is 1
    # each patient equal within source: patient totals within a source are identical
    pat_of_pair = fd.pos_patient[fd.pair_pos]
    for s in sources:
        pats = np.unique(pat_of_pair[src_of_pair == s])
        tots = [fd.pair_w[pat_of_pair == p].sum() for p in pats]
        assert np.allclose(tots, tots[0])
        assert np.isclose(tots[0], 1.0 / (S * len(pats)))


# --------------------------------------------------------------------------------------------------
# 6) robust HLA-locus parser + exact rowwise-SD disagreement (after mask semantics)
# --------------------------------------------------------------------------------------------------
def test_locus_parser_and_disagreement_context():
    assert v5.locus("A0201") == "A"
    assert v5.locus("HLA-A01:01") == "A"
    assert v5.locus("HLA-A*02:01") == "A"
    assert v5.locus("B*07:02") == "B"
    assert v5.locus("HLA-C07:01") == "C"
    assert v5.locus("DRB1*01:01") == "OTHER"                       # class II -> reference (A)
    # disagreement = population SD (ddof=0) of the leakage-safe masked pcts [f_prime_pct,f_mix_pct,f_el_pct]
    base = pd.DataFrame({"peptide": ["ACDEFGHIK"], "hla_allele": ["A*02:01"],
                         "f_prime_pct": [0.5], "f_mix_pct": [0.5], "f_el_pct": [0.5]})
    out = v5.add_approved_contexts(base)
    assert out["ctx_pred_disagree"].iloc[0] == 0.0                 # masked prime=0.5, all equal -> SD 0
    assert out["ctx_pep_len"].iloc[0] == 9.0
    assert out["ctx_locus_B"].iloc[0] == 0.0 and out["ctx_locus_C"].iloc[0] == 0.0   # locus A -> both dummies 0
    d2 = base.assign(f_prime_pct=[0.2], f_mix_pct=[0.5], f_el_pct=[0.8])
    got = v5.add_approved_contexts(d2)["ctx_pred_disagree"].iloc[0]
    assert abs(got - float(np.std([0.2, 0.5, 0.8]))) < 1e-12       # ddof=0
    b = v5.add_approved_contexts(base.assign(hla_allele=["B*07:02"]))
    assert b["ctx_locus_B"].iloc[0] == 1.0 and b["ctx_locus_C"].iloc[0] == 0.0


# --------------------------------------------------------------------------------------------------
# 7) φ carries no forbidden/source field, and R(β=0) ≡ Q (structurally and at λ_ctx=∞)
# --------------------------------------------------------------------------------------------------
def test_design_has_no_forbidden_fields():
    v5.assert_design_is_safe(v5.design_column_names(True))         # no raise on the real R design
    for bad in ("source", "study_id", "assay", "fold", "outcome", "pool_size", "hla_multiplicity",
                "routes_per_bag", "eval_positive"):
        with pytest.raises(AssertionError):
            v5.assert_design_is_safe(list(FEATURES) + [bad])


def test_R_beta_zero_equals_Q():
    f = _synth(seed=7)
    # structural: with shared scalers, φ_R·[w0;0] == φ_Q·w0 (interactions × 0 vanish)
    lab = v5.add_approved_contexts(
        f[(f["eval_positive"] == 1) | (f["bag_label"] == "NEGATIVE")].copy()).reset_index(drop=True)
    X = lab[FEATURES].to_numpy(float)
    fmean, fstd = X.mean(0), X.std(0) + 1e-9
    Z = lab[v5.CTX_COLS].to_numpy(float)
    cmean, cstd = Z.mean(0), Z.std(0) + 1e-9
    phiQ = v5.build_design(lab, False, fmean, fstd, None, None)
    phiR = v5.build_design(lab, True, fmean, fstd, cmean, cstd)
    w0 = np.random.default_rng(1).normal(size=9)
    assert np.allclose(phiQ @ w0, phiR @ np.concatenate([w0, np.zeros(16)]))
    # registered: λ_ctx=∞ ⇒ β=0 ⇒ R fits the identical 9-dim Q problem
    q = ContextPairwiseRanker(member="Q", lam_w=1.0).fit(f)
    rinf = ContextPairwiseRanker(member="R", lam_w=1.0, lam_ctx=np.inf).fit(f)
    assert rinf.has_beta is False and rinf.beta_ is None
    assert np.max(np.abs(q.coef_ - rinf.coef_)) < 1e-10
    assert np.max(np.abs(q.raw_score(f) - rinf.raw_score(f))) < 1e-10


# --------------------------------------------------------------------------------------------------
# 8) all fitted scaling is a function of the TRAIN fit-rows only (no outer-test leakage)
# --------------------------------------------------------------------------------------------------
def test_train_only_scaling_no_leakage():
    f, held = _synth(seed=8), _synth(seed=999)
    m = ContextPairwiseRanker(member="R", lam_w=1.0, lam_ctx=10.0).fit(f)
    fm, fs, cm, cs = m.feat_mean_.copy(), m.feat_std_.copy(), m.ctx_mean_.copy(), m.ctx_std_.copy()
    m.raw_score(held)                                             # scoring a held-out frame must not mutate
    assert np.array_equal(m.feat_mean_, fm) and np.array_equal(m.feat_std_, fs)
    assert np.array_equal(m.ctx_mean_, cm) and np.array_equal(m.ctx_std_, cs)
    m2 = ContextPairwiseRanker(member="R", lam_w=1.0, lam_ctx=10.0).fit(f)
    assert np.allclose(m.feat_mean_, m2.feat_mean_) and np.allclose(m.ctx_std_, m2.ctx_std_)


# --------------------------------------------------------------------------------------------------
# 9) deterministic selection tie-break (more shrinkage first: larger λ_ctx, then larger λ_w)
# --------------------------------------------------------------------------------------------------
def test_deterministic_selection_tie_break():
    scored = [(1.000, {"lam_w": 0.1, "lam_ctx": 3.0}),
              (1.005, {"lam_w": 0.3, "lam_ctx": 30.0}),
              (1.000, {"lam_w": 1.0, "lam_ctx": 100.0}),
              (0.500, {"lam_w": 10.0, "lam_ctx": np.inf})]
    assert v5._pick_from_scored(scored, eps=0.01) == {"lam_w": 1.0, "lam_ctx": 100.0}   # largest λ_ctx in tie
    tie_inf = [(1.0, {"lam_w": 1.0, "lam_ctx": 100.0}), (1.0, {"lam_w": 1.0, "lam_ctx": np.inf})]
    assert not np.isfinite(v5._pick_from_scored(tie_inf)["lam_ctx"])                    # ∞ wins (most shrunk)
    q = [(1.0, {"lam_w": 0.1}), (1.0, {"lam_w": 10.0}), (0.9, {"lam_w": 1.0})]
    assert v5._pick_from_scored(q) == {"lam_w": 10.0}                                   # Q: larger λ_w in tie


# --------------------------------------------------------------------------------------------------
# 10) comparator hyperparameter extraction uses the correct JSON keys and frozen classes
# --------------------------------------------------------------------------------------------------
def test_comparator_hparam_extraction_keys_and_classes():
    v03 = {"ladder": {
        "rung3_MIL": {"folds": [{"fold": 0, "C": 0.1, "tau": 1.0}, {"fold": 1, "C": 0.3, "tau": 0.5}]},
        "rung2_additive": {"folds": [{"fold": 0, "C": 0.1}, {"fold": 1, "C": 1.0}]}}}
    v04 = {"members": {"F_feature_tower": {"folds": [
        {"fold": 0, "C": 1.0, "tau": 0.5, "lam": 10.0}, {"fold": 1, "C": 0.3, "tau": 1.0, "lam": None}]}}}
    P, A, F = v5.extract_P_hparams(v03), v5.extract_A_hparams(v03), v5.extract_F_hparams(v04)
    assert P[0] == {"C": 0.1, "tau": 1.0} and P[1]["tau"] == 0.5
    assert A[0] == {"C": 0.1} and set(A[1]) == {"C"}               # A carries C only (no tau/lam)
    assert F[0]["lam"] == 10.0 and np.isinf(F[1]["lam"])           # lam None -> ∞
    # reconstruction builds the exact frozen classes
    f = _synth(seed=2)
    pf = {i: {"C": 1.0, "tau": 1.0} for i in range(5)}
    af = {i: {"C": 1.0} for i in range(5)}
    ff = {i: {"C": 1.0, "tau": 1.0, "lam": 1.0} for i in range(5)}
    assert isinstance(v5._reconstruct_oof(f, pf, lambda h: MILRanker(**h), "P").models[0][1], MILRanker)
    assert isinstance(v5._reconstruct_oof(f, af, lambda h: AdditiveRanker(**h), "A").models[0][1], AdditiveRanker)
    assert isinstance(v5._reconstruct_oof(
        f, ff, lambda h: TowerMILRanker(member="F", **h), "F").models[0][1], TowerMILRanker)


# --------------------------------------------------------------------------------------------------
# 10b) reproduction verification is REAL and fail-fast (convex P/A tight; nonconvex F honest tolerance)
# --------------------------------------------------------------------------------------------------
def test_reconstruction_verification_is_real_and_fail_fast():
    import copy
    f = _synth(seed=2)
    v03 = {"ladder": {"rung2_additive": {"folds": [{"fold": i, "C": 1.0} for i in range(5)]}}}
    oofA = v5.reconstruct_A(f, v03)
    stored = [{"fold": fold, "coefficients": {ftr: round(float(c), 4) for ftr, c in zip(FEATURES, m.coef_)}}
              for fold, m in oofA.models]                          # stored = the frozen 4-dp coefficients
    hits = v5.overall_hits(oofA)
    assert v5.verify_convex_reconstruction(oofA, stored, hits)["reproduced"]     # a faithful refit passes
    bad_coef = copy.deepcopy(stored); bad_coef[0]["coefficients"][FEATURES[0]] += 1.0
    with pytest.raises(RuntimeError):
        v5.verify_convex_reconstruction(oofA, bad_coef, hits)      # perturbed coefficient -> fail fast
    with pytest.raises(RuntimeError):
        v5.verify_convex_reconstruction(oofA, stored, hits + 0.5)  # perturbed aggregate hits -> fail fast
    # F: nonconvex, verified to a tolerance on aggregate hits; a gross residual still fails
    v04 = {"members": {"F_feature_tower": {"folds": [
        {"fold": i, "C": 1.0, "tau": 1.0, "lam": 1.0} for i in range(5)]}}}
    oofF = v5.reconstruct_F(f, v04)
    hF = v5.overall_hits(oofF)
    assert v5.verify_f_reconstruction(oofF, hF)["reproduced"] and v5.verify_f_reconstruction(oofF, hF)["nonconvex"]
    with pytest.raises(RuntimeError):
        v5.verify_f_reconstruction(oofF, hF + 1.0)


# --------------------------------------------------------------------------------------------------
# 11) provenance mismatch fails closed; Gartner TEST path access raises
# --------------------------------------------------------------------------------------------------
def test_provenance_fails_closed(tmp_path):
    import hashlib
    target = "src/event_b/epicurus_v05.py"
    good = hashlib.sha256(open(target, "rb").read()).hexdigest()
    bad = tmp_path / "prov_bad.json"
    bad.write_text(json.dumps({"git_head": "x", "input_paths": {"m": target},
                               "inputs_sha256": {"m": "deadbeef"}}))
    with pytest.raises(RuntimeError):
        v5.verify_provenance(bad)
    ok = tmp_path / "prov_ok.json"
    ok.write_text(json.dumps({"git_head": "abc123", "input_paths": {"m": target},
                              "inputs_sha256": {"m": good}}))
    out = v5.verify_provenance(ok)
    assert out["n_inputs_verified"] == 1 and out["git_head"] == "abc123"


def test_gartner_test_path_access_raises():
    from event_b.epicurus_v04 import guard_no_test_io
    with guard_no_test_io():
        with pytest.raises(RuntimeError):
            open("data/raw/gartner_nci/NmersTestingSet.txt")
        with pytest.raises(RuntimeError):
            pd.read_csv("data/raw/gartner_nci/MmpsTestingSet_extract.tsv", sep="\t")


# --------------------------------------------------------------------------------------------------
# 12) selection/OOF integrity smoke on synthetic data (no full benchmark, no real fit)
# --------------------------------------------------------------------------------------------------
def test_oof_integrity_and_selection_smoke():
    f = _synth(seed=11, n_patients=10)
    oofQ = v5.oof_qr(f, "Q")
    assert oofQ.metrics["patient_id"].is_unique                   # no patient scored by a model that saw it
    cfg = v5.select_lambda(f[f["fold"] != 0], "R")
    assert set(cfg) == {"lam_w", "lam_ctx"} and cfg["lam_ctx"] >= cfg["lam_w"]  # registered grid constraint


# ---- real-data guards (never fit Q/R) ----
@pytest.mark.skipif(not _READY, reason="v0.5 caches/split not available")
def test_real_frame_contexts_and_attrition_no_fit():
    from event_b.epicurus_v04 import assemble_frame, attrition_report
    frame = assemble_frame()
    a = attrition_report(frame)
    assert a["TOTAL"]["feature_bearing"] == 152 and a["TOTAL"]["scored_has_positive"] == 118
    ev = v5.add_approved_contexts(frame[~frame["quarantined"]])
    for c in v5.CTX_COLS:
        assert c in ev.columns and ev[c].notna().all()
    assert set(np.unique(ev["ctx_locus_B"])) <= {0.0, 1.0}


@pytest.mark.skipif(not _READY, reason="v0.5 caches/split not available")
def test_context_construction_matches_frozen_audit():
    import sys
    sys.path.insert(0, "scripts")
    import epicurus_v05_context_audit as audit
    from event_b.epicurus_v04 import assemble_frame
    ev = assemble_frame()
    ev = ev[~ev["quarantined"]].reset_index(drop=True)
    a = audit.build_context_columns(ev)                           # frozen pre-fit audit construction
    b = v5.add_approved_contexts(ev)                              # module construction used at fit time
    for c in v5.CTX_COLS:
        assert np.allclose(a[c].to_numpy(float), b[c].to_numpy(float), equal_nan=True), c
