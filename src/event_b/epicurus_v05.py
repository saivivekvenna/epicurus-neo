"""Epicurus v0.5 DEVELOPMENT experiment — deployable context-conditioned pairwise challenger.

DEVELOPMENT ONLY. Implements the FROZEN protocol
`artifacts/milestone_7_decision/epicurus_v05/PREREGISTERED_PROTOCOL.md` verbatim and reuses v0.3's frozen
feature pipeline / evaluation and v0.4's TEST-I/O guard (imported, never reimplemented). The Gartner TEST
holdout is never loaded/scored; the frozen split is used verbatim; no external claim is produced here.

Two new models on top of the frozen ladder (P/A/F, reconstructed by exact frozen code, §2.1):

    Q — shared pairwise   s(x) = w0·x                                 (no intercept; strictly convex)
    R — context ranker     s(x) = w0·x + Σ_{f∈C4}Σ_c β_{f,c}·x_f·z_c   (no intercept; R⊃Q, β=0 ⇒ R≡Q)

The objective (§5) is a weighted within-patient PAIRWISE logistic loss: every exact positive peptide is paired
with every tested-negative BAG in its patient; a negative bag enters through a convex log-MEAN-exp (τ=1) over
ALL its children; instance-source negatives are singleton bags (⇒ ordinary pairwise logistic). Every fitted
coefficient carries a positive L2 (λ_w on w0, λ_ctx on β) ⇒ strictly convex ⇒ unique, init-invariant solution.
Solved deterministically with an analytic gradient (L-BFGS-B, fixed zero init).

This module holds the MATH and reconstruction only (separated from report serialization so unit tests can
exercise it without running the full benchmark — the runner is `scripts/epicurus_v05_dev.py`).
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit

from event_b.epicurus_v03 import (  # frozen v0.3 pipeline / evaluation, reused verbatim
    CORE, FEATURES, K_TOP, AdditiveRanker, MILRanker, OOFResult, _source_balanced_weight, per_patient_metrics,
)
from event_b.epicurus_v04 import TowerMILRanker, _sha256  # frozen v0.4 source-name tower (F) + hashing

# --------------------------------------------------------------------------------------------------
# frozen protocol constants (§4, §5, §7). Changing any of these changes the registered experiment.
# --------------------------------------------------------------------------------------------------
TAU = 1.0                                                              # §5: fixed, NOT on the grid
CTX_COLS = ["ctx_pep_len", "ctx_pred_disagree", "ctx_locus_B", "ctx_locus_C"]   # §4 z (encoding cols)
C4 = ["f_el_pct", "f_pres_abs", "f_prime_pct", "f_mix_pct"]           # §4 core presentation feats for R
DISAGREE_COLS = ["f_prime_pct", "f_mix_pct", "f_el_pct"]              # leakage-safe masked pcts (== audit)
LAM_W_GRID = [0.1, 0.3, 1.0, 3.0, 10.0]                              # §7 λ_w
LAM_CTX_GRID = [3.0, 10.0, 30.0, 100.0, np.inf]                      # §7 λ_ctx (∞ ⇒ β=0 ⇒ R≡Q); λ_ctx ≥ λ_w
SELECT_EPS = 0.01                                                     # §7 selection tolerance

# CORE is exactly the C4 set (guards against a silent feature-list drift renaming a core feature)
assert set(C4) == set(CORE), (C4, CORE)

# forbidden φ tokens (§4): no source/study/assay/fold/outcome/pool-size/multiplicity/route may enter φ
_FORBIDDEN_PHI_TOKENS = ("source", "study", "assay", "fold", "outcome", "pool",
                         "multipl", "route", "label", "eval_positive", "bag")

V05_DIR = Path("artifacts/milestone_7_decision/epicurus_v05")
V05_PROVENANCE = V05_DIR / "PROVENANCE.json"
V03_RESULT = Path("artifacts/milestone_7_decision/epicurus_v03/DEV_RESULT.json")
V04_RESULT = Path("artifacts/milestone_7_decision/epicurus_v04/DEV_RESULT.json")


# ==================================================================================================
# portable, candidate-level contexts (§4) — mirrors the FROZEN pre-fit audit
# (`scripts/epicurus_v05_context_audit.py`) exactly; a consistency test cross-checks the two.
# ==================================================================================================
def locus(allele: str) -> str:
    """HLA class-I locus in {A,B,C} from the normalized allele string; anything else -> 'OTHER' (== A ref).
    Robust across every source format: 'A0201', 'HLA-A01:01', 'HLA-A*02:01'. Identical to the audit `_locus`."""
    a = str(allele)
    m = re.search(r"HLA[-_ ]?([ABC])", a, re.I) or re.match(r"\s*([ABC])[\*0-9]", a, re.I)
    return m.group(1).upper() if m else "OTHER"


def add_approved_contexts(df: pd.DataFrame) -> pd.DataFrame:
    """Attach the four APPROVED candidate-level context columns (§4). Idempotent, OUTCOME-label-blind, and a
    pure function of each candidate's own sequence / allele / leakage-safe masked predictor percentiles — no
    pool enumeration, no labels. Centering/scaling happens later, TRAIN-only, inside the model."""
    if all(c in df.columns for c in CTX_COLS):
        return df
    df = df.copy()
    df["ctx_pep_len"] = df["peptide"].astype(str).str.len().astype(float)          # from the sequence alone
    df["ctx_pred_disagree"] = df[DISAGREE_COLS].std(axis=1, ddof=0).astype(float)   # rowwise SD of masked pcts
    loc = df["hla_allele"].map(locus)
    df["ctx_locus_B"] = (loc == "B").astype(float)   # A(/OTHER) is the reference level; B,C are the encoding
    df["ctx_locus_C"] = (loc == "C").astype(float)
    return df


# ==================================================================================================
# design matrix (§5): φ = x for Q; φ = [x ; vec(x_{C4} ⊗ z)] for R. No intercept, no forbidden field.
# ==================================================================================================
def interaction_names() -> list[str]:
    """The 16 R interaction column names, f (outer) × context (inner) — deterministic order."""
    return [f"{f}__x__{c}" for f in C4 for c in CTX_COLS]


def design_column_names(has_beta: bool) -> list[str]:
    """Ordered φ column names. Q (or R at λ_ctx=∞): the 9 shared features. R: features + 16 interactions."""
    return list(FEATURES) + (interaction_names() if has_beta else [])


def assert_design_is_safe(cols: list[str]) -> None:
    """§4 guard: assert φ contains ONLY the shared features and the registered C4×context interactions — no
    source/study/assay/fold/outcome/pool-size/multiplicity/route field, and no context other than CTX_COLS."""
    assert cols[:len(FEATURES)] == list(FEATURES), ("shared features drifted", cols[:len(FEATURES)])
    for c in cols:
        low = c.lower()
        # the only names that may appear are FEATURES and "<C4feat>__x__<ctx>"
        if "__x__" in c:
            f, ctx = c.split("__x__")
            assert f in C4 and ctx in CTX_COLS, ("illegal interaction column", c)
        else:
            assert c in FEATURES, ("illegal φ column", c)
        # defense in depth: reject any forbidden substring in a non-context, non-feature position
        if c not in FEATURES and "__x__" not in c:
            for tok in _FORBIDDEN_PHI_TOKENS:
                assert tok not in low, ("forbidden token in φ", c, tok)


def _standardize(X: np.ndarray, mean: np.ndarray, std: np.ndarray) -> np.ndarray:
    return (X - mean) / std


def build_design(df: pd.DataFrame, has_beta: bool, feat_mean: np.ndarray, feat_std: np.ndarray,
                 ctx_mean: np.ndarray | None, ctx_std: np.ndarray | None) -> np.ndarray:
    """Standardized φ for the given frame. x = standardized FEATURES; for R the 16 interactions are the
    standardized C4 columns × the standardized context columns (so β is scale-comparable under the shared λ_ctx
    and centered contexts make β=0 ⇒ R≡Q)."""
    x = _standardize(df[FEATURES].to_numpy(float), feat_mean, feat_std)
    if not has_beta:
        return x
    df = add_approved_contexts(df)
    z = _standardize(df[CTX_COLS].to_numpy(float), ctx_mean, ctx_std)
    c4_idx = [FEATURES.index(f) for f in C4]
    inter = np.concatenate([x[:, [j]] * z for j in c4_idx], axis=1)   # f outer, context inner (16 cols)
    return np.concatenate([x, inter], axis=1)


# ==================================================================================================
# θ-independent fit scaffolding (§5): positives, negative bags, within-patient pairs and their weights
# ==================================================================================================
@dataclass
class _FitData:
    phi_pos: np.ndarray            # (n_pos, D)  design rows of the exact positives
    phi_child: np.ndarray          # (n_child, D) design rows of every negative-bag child (bag-contiguous)
    seg_boundaries: np.ndarray     # (n_bag,) first-child index of each negative bag in phi_child
    seg_counts: np.ndarray         # (n_bag,) |b|, children per negative bag
    pair_pos: np.ndarray           # (n_pairs,) positive index of each (positive, negative-bag) pair
    pair_bag: np.ndarray           # (n_pairs,) negative-bag index of each pair
    pair_w: np.ndarray             # (n_pairs,) 1/(S·N_s·P_p·B_p)
    penalty: np.ndarray            # (D,) per-coefficient L2 weight (λ_w on w0, λ_ctx on β)
    pos_patient: np.ndarray        # (n_pos,) patient of each positive (diagnostics / weight-total checks)
    pos_source: np.ndarray         # (n_pos,) source of each positive
    bag_patient: np.ndarray        # (n_bag,) patient of each negative bag
    bag_source: np.ndarray         # (n_bag,) source of each negative bag


def _prepare_fit(lab: pd.DataFrame, phi: np.ndarray, penalty: np.ndarray) -> _FitData:
    """Enumerate the exact positives, the tested-negative bags (Gartner NEGATIVE bags + singleton instance
    negatives), and every within-patient (positive × negative-bag) pair with its exactly-normalized weight.
    All θ-independent; computed once per fit. `phi` is row-aligned to `lab` (both reset to 0..n-1)."""
    lab = lab.reset_index(drop=True)
    pos_mask = (lab["eval_positive"].to_numpy() == 1)
    neg_mask = (lab["bag_label"].to_numpy() == "NEGATIVE")

    pos_rows = np.flatnonzero(pos_mask)
    phi_pos = phi[pos_rows]
    pos_patient = lab["patient_id"].to_numpy()[pos_rows]
    pos_source = lab["source"].to_numpy()[pos_rows]

    # negative bags: children grouped by bag_id, made contiguous so segmented log-mean-exp is a reduceat
    neg = lab.loc[neg_mask]
    neg_codes = pd.factorize(neg["bag_id"])[0]
    order = np.argsort(neg_codes, kind="stable")
    neg_rows = neg.index.to_numpy()[order]
    codes_s = neg_codes[order]
    phi_child = phi[neg_rows]
    boundaries = np.concatenate([[0], np.flatnonzero(np.diff(codes_s)) + 1]) if len(codes_s) else np.array([], int)
    counts = np.diff(np.append(boundaries, len(codes_s))) if len(boundaries) else np.array([], int)
    bag_patient = lab["patient_id"].to_numpy()[neg_rows][boundaries] if len(boundaries) else np.array([])
    bag_source = lab["source"].to_numpy()[neg_rows][boundaries] if len(boundaries) else np.array([])

    # index positives and bags by patient
    pos_by_patient: dict[str, list[int]] = {}
    for a, p in enumerate(pos_patient):
        pos_by_patient.setdefault(p, []).append(a)
    bag_by_patient: dict[str, list[int]] = {}
    bag_src_of: dict[str, str] = {}
    for b, p in enumerate(bag_patient):
        bag_by_patient.setdefault(p, []).append(b)
        bag_src_of[p] = bag_source[b]

    # contributing patients = those with ≥1 positive AND ≥1 negative bag; N_s counts THEM (exact source balance)
    contrib = [p for p in pos_by_patient if p in bag_by_patient]
    n_s: dict[str, int] = {}
    for p in contrib:
        n_s[bag_src_of[p]] = n_s.get(bag_src_of[p], 0) + 1
    n_src = len(n_s)

    pair_pos, pair_bag, pair_w = [], [], []
    for p in contrib:
        aa = pos_by_patient[p]
        bb = bag_by_patient[p]
        w = 1.0 / (n_src * n_s[bag_src_of[p]] * len(aa) * len(bb))   # 1/(S·N_s·P_p·B_p), constant within p
        for a in aa:
            for b in bb:
                pair_pos.append(a); pair_bag.append(b); pair_w.append(w)

    return _FitData(
        phi_pos=phi_pos, phi_child=phi_child, seg_boundaries=boundaries, seg_counts=counts,
        pair_pos=np.asarray(pair_pos, int), pair_bag=np.asarray(pair_bag, int),
        pair_w=np.asarray(pair_w, float), penalty=penalty,
        pos_patient=pos_patient, pos_source=pos_source,
        bag_patient=np.asarray(bag_patient), bag_source=np.asarray(bag_source),
    )


def _loss_grad(theta: np.ndarray, fd: _FitData) -> tuple[float, np.ndarray]:
    """Weighted within-patient pairwise logistic loss with bag-aware log-MEAN-exp negatives (τ=1) and its
    ANALYTIC gradient. Strictly convex once `fd.penalty` is positive on every coefficient."""
    s_pos = fd.phi_pos @ theta
    s_child = fd.phi_child @ theta

    # segmented log-mean-exp over each bag's children (stable): LME_b = seg_max + log(mean exp(s - seg_max))
    rep = lambda a: np.repeat(a, fd.seg_counts)
    seg_max = np.maximum.reduceat(s_child, fd.seg_boundaries)
    e = np.exp(s_child - rep(seg_max))                      # τ = 1
    seg_sum = np.add.reduceat(e, fd.seg_boundaries)
    lme = seg_max + np.log(seg_sum / fd.seg_counts)         # log-MEAN-exp (÷|b| ⇒ duplicate-child invariant)
    p_child = e / rep(seg_sum)                              # softmax weight of each child within its bag

    u = lme[fd.pair_bag] - s_pos[fd.pair_pos]               # margin: negative-bag soft-score minus positive
    loss = float(np.sum(fd.pair_w * np.logaddexp(0.0, u)))  # softplus, stable
    sig = expit(u)                                          # d softplus / du
    wsig = fd.pair_w * sig

    ga = np.bincount(fd.pair_pos, weights=wsig, minlength=fd.phi_pos.shape[0])   # per-positive accumulation
    db = np.bincount(fd.pair_bag, weights=wsig, minlength=len(fd.seg_counts))    # per-bag accumulation
    grad = -(fd.phi_pos.T @ ga) + fd.phi_child.T @ (rep(db) * p_child)

    pen = 0.5 * float(theta @ (fd.penalty * theta))
    grad = grad + fd.penalty * theta
    return loss + pen, grad


# ==================================================================================================
# the model: shared pairwise (Q) and context-conditioned pairwise (R)
# ==================================================================================================
@dataclass
class ContextPairwiseRanker:
    """member ∈ {"Q","R"}. No intercept. λ_w shrinks w0; λ_ctx shrinks β (R only). λ_ctx=∞ ⇒ β=0 ⇒ R≡Q, and
    is fit as the 9-dim Q problem (nested check + runtime reuse). Deterministic analytic-gradient L-BFGS-B."""
    member: str = "R"
    lam_w: float = 1.0
    lam_ctx: float = 10.0            # ignored for Q; np.inf ⇒ β=0 (identical to Q at the same λ_w)
    # fitted state
    feat_mean_: np.ndarray = None
    feat_std_: np.ndarray = None
    ctx_mean_: np.ndarray = None
    ctx_std_: np.ndarray = None
    coef_: np.ndarray = None         # w0 (len 9)
    beta_: np.ndarray = None         # β (len 16) or None when β≡0
    n_iter_: int = 0

    @property
    def has_beta(self) -> bool:
        return bool(self.member == "R" and np.isfinite(self.lam_ctx))

    def _penalty(self) -> np.ndarray:
        pen = np.full(len(FEATURES), float(self.lam_w))
        if self.has_beta:
            pen = np.concatenate([pen, np.full(len(interaction_names()), float(self.lam_ctx))])
        return pen

    def fit(self, train: pd.DataFrame, init: np.ndarray | None = None) -> "ContextPairwiseRanker":
        assert self.member in ("Q", "R"), self.member
        lab = train.loc[(train["eval_positive"] == 1) | (train["bag_label"] == "NEGATIVE")].copy()
        lab = add_approved_contexts(lab).reset_index(drop=True)

        # scalers fit STRICTLY on the training fit-rows (feature std matches v0.3; contexts centered/scaled)
        Xf = lab[FEATURES].to_numpy(float)
        self.feat_mean_, self.feat_std_ = Xf.mean(0), Xf.std(0) + 1e-9
        Zf = lab[CTX_COLS].to_numpy(float)
        self.ctx_mean_, self.ctx_std_ = Zf.mean(0), Zf.std(0) + 1e-9

        cols = design_column_names(self.has_beta)
        assert_design_is_safe(cols)
        phi = build_design(lab, self.has_beta, self.feat_mean_, self.feat_std_, self.ctx_mean_, self.ctx_std_)
        assert phi.shape[1] == len(cols)

        fd = _prepare_fit(lab, phi, self._penalty())
        D = phi.shape[1]
        theta0 = np.zeros(D) if init is None else np.asarray(init, float)
        assert theta0.shape == (D,), (theta0.shape, D)
        # SOLVER TOLERANCE — REVIEW FLAG for Codex. The protocol registers §5 `ftol=1e-9`, but on this
        # O(1)-scale objective (total pair-weight normalizes to 1) `ftol=1e-9` halts on the function plateau
        # ~2e-5 from the unique minimum — which FAILS the registered §9.6 HARD GATE (max|Δcoef| ≤ 1e-6 across
        # perturbed inits) even though the problem is strictly convex (rank/Spearman is 1.0 regardless, so
        # hits@20 is unaffected). Strict convexity guarantees a unique global minimum, so we converge on the
        # gradient: ftol=1e-14 (do not stop on the function plateau early) + gtol=1e-10. Empirically this
        # reaches ≤3e-7 across the registered grid, seeds, and perturbation scales — satisfying the hard gate
        # with margin. maxiter=500 per §5. This is a pre-fit solver-precision decision (not a change to pairs,
        # weights, folds, or the registered grid); flagged for Codex to ratify vs. relaxing the §9.6 threshold.
        res = minimize(_loss_grad, theta0, args=(fd,), jac=True, method="L-BFGS-B",
                       options={"maxiter": 500, "ftol": 1e-14, "gtol": 1e-10})
        self.n_iter_ = int(res.nit)
        self.coef_ = res.x[:len(FEATURES)]
        self.beta_ = res.x[len(FEATURES):].copy() if self.has_beta else None
        return self

    def raw_score(self, df: pd.DataFrame) -> np.ndarray:
        x = _standardize(df[FEATURES].to_numpy(float), self.feat_mean_, self.feat_std_)
        if not self.has_beta:
            return x @ self.coef_
        theta = np.concatenate([self.coef_, self.beta_])
        phi = build_design(df, True, self.feat_mean_, self.feat_std_, self.ctx_mean_, self.ctx_std_)
        return phi @ theta

    def to_dict(self) -> dict:
        d = {"kind": "context_pairwise", "member": self.member, "lam_w": self.lam_w,
             "lam_ctx": (None if not np.isfinite(self.lam_ctx) else self.lam_ctx),
             "w0": {f: round(float(w), 4) for f, w in zip(FEATURES, self.coef_)}}
        if self.beta_ is not None:
            d["beta"] = {n: round(float(b), 4) for n, b in zip(interaction_names(), self.beta_)}
        return d


# ==================================================================================================
# nested selection (§7): inner grouped-OOF by source-balanced hits@20 + the registered tie-break.
# ==================================================================================================
def _grid(member: str) -> list[dict]:
    """Registered grids (§7). Q: λ_w only. R: (λ_w, λ_ctx) with λ_ctx ≥ λ_w; λ_ctx=∞ ⇒ R≡Q (nested check)."""
    if member == "Q":
        return [{"lam_w": lw} for lw in LAM_W_GRID]
    return [{"lam_w": lw, "lam_ctx": lc} for lw in LAM_W_GRID for lc in LAM_CTX_GRID if lc >= lw]


def _effective(member: str, cfg: dict) -> tuple[float, float, bool]:
    """(λ_w, effective λ_ctx, has_beta). member='R' with λ_ctx=∞ collapses to the 9-dim Q problem."""
    lam_w = cfg["lam_w"]
    lam_ctx = cfg.get("lam_ctx", np.inf)
    has_beta = (member == "R") and np.isfinite(lam_ctx)
    return lam_w, (lam_ctx if has_beta else np.inf), has_beta


def _fit_cached(member: str, cfg: dict, sub: pd.DataFrame, cache: dict | None) -> ContextPairwiseRanker:
    """Fit (λ_w, λ_ctx) on `sub`, reusing an identical fit from `cache`. Because λ_ctx=∞ collapses to the Q
    (9-dim) problem, its cache key ignores λ_ctx — so Q and R's λ_ctx=∞ column REUSE one fit (§5/§7)."""
    lam_w, lam_ctx, has_beta = _effective(member, cfg)
    key = (tuple(sorted(int(f) for f in sub["fold"].unique())), round(lam_w, 6), has_beta,
           None if not has_beta else round(lam_ctx, 6))
    if cache is not None and key in cache:
        return cache[key]
    m = ContextPairwiseRanker(member=("R" if has_beta else "Q"), lam_w=lam_w,
                              lam_ctx=(lam_ctx if has_beta else np.inf)).fit(sub)
    if cache is not None:
        cache[key] = m
    return m


def _sb_hits(mt: pd.DataFrame) -> float:
    return float(np.average(mt["hits"], weights=_source_balanced_weight(mt))) if len(mt) else -1e9


def _shrink_key(cfg: dict):
    """Registered preference among near-ties: MORE shrinkage first — larger λ_ctx (∞ wins), then larger λ_w;
    lexical on (λ_ctx, λ_w) is subsumed because a grid point is uniquely (λ_ctx, λ_w)."""
    lam_ctx = cfg.get("lam_ctx", np.inf)
    lc = np.inf if not np.isfinite(lam_ctx) else lam_ctx
    return (-lc, -cfg["lam_w"])


def _pick_from_scored(scored: list[tuple[float, dict]], eps: float = SELECT_EPS) -> dict:
    """(1) maximize inner score; (2) within eps prefer more shrinkage; (3) lexical (subsumed). Deterministic."""
    best = max(s for s, _ in scored)
    cands = [cfg for s, cfg in scored if s >= best - eps]
    return min(cands, key=_shrink_key)


def select_lambda(tr: pd.DataFrame, member: str, k: int = K_TOP, cache: dict | None = None) -> dict:
    """Nested inner grouped-OOF selection INSIDE outer-train (§7): source-balanced mean patient hits@20, then
    the registered tie-break. No outer-test row is ever touched; the registered grid is not altered."""
    inner = sorted(tr["fold"].unique())
    scored = []
    for cfg in _grid(member):
        vals = []
        for vf in inner:
            itr, ite = tr[tr["fold"] != vf], tr[tr["fold"] == vf]
            if itr.empty or ite.empty:
                continue
            m = _fit_cached(member, cfg, itr, cache)
            vals.append(_sb_hits(per_patient_metrics(ite, m.raw_score(ite), k)))
        scored.append((float(np.mean(vals)) if vals else -1e9, cfg))
    return _pick_from_scored(scored)


def oof_qr(frame: pd.DataFrame, member: str, k: int = K_TOP, cache: dict | None = None) -> OOFResult:
    """Outer nested OOF on the 5 frozen folds for Q or R. Selection runs strictly inside each outer-train; the
    outer model is fit at the selected (λ_w, λ_ctx) and scores only its held-out fold."""
    dev = frame[~frame["quarantined"]].copy()
    cache = {} if cache is None else cache
    parts, specs, models = [], [], []
    for f in sorted(dev["fold"].unique()):
        tr, te = dev[dev["fold"] != f], dev[dev["fold"] == f]
        cfg = select_lambda(tr, member, k, cache)
        model = _fit_cached(member, cfg, tr, cache)
        parts.append(per_patient_metrics(te, model.raw_score(te), k))
        lam_ctx = cfg.get("lam_ctx", np.inf)
        specs.append({"fold": int(f), "member": member, "lam_w": cfg["lam_w"],
                      "lam_ctx": (None if not np.isfinite(lam_ctx) else lam_ctx)})
        models.append((int(f), model))
    return OOFResult(pd.concat(parts, ignore_index=True), {"folds": specs}, models)


def overall_hits(oof: OOFResult) -> float:
    """Source-balanced mean patient hits@20 over the OOF metrics — the identical aggregation as v0.3/v0.4."""
    mt = oof.metrics
    return round(float(np.average(mt["hits"], weights=_source_balanced_weight(mt))), 4)


# ==================================================================================================
# frozen comparators P/A/F are REFIT from exact frozen code + stored per-fold hyperparameters (§2.1).
# The result JSONs are aggregate/coefficient records, NOT scoreable state — we reconstruct, then verify.
# ==================================================================================================
def extract_P_hparams(v03: dict) -> dict:
    """P = frozen corrected-v0.3 rung-3 MIL: per-fold (C, τ) from `ladder.rung3_MIL.folds`."""
    return {int(f["fold"]): {"C": f["C"], "tau": f["tau"]} for f in v03["ladder"]["rung3_MIL"]["folds"]}


def extract_A_hparams(v03: dict) -> dict:
    """A = frozen v0.3 rung-2 AdditiveRanker (additive_logistic): per-fold C from `ladder.rung2_additive.folds`."""
    return {int(f["fold"]): {"C": f["C"]} for f in v03["ladder"]["rung2_additive"]["folds"]}


def extract_F_hparams(v04: dict) -> dict:
    """F = frozen v0.4 TowerMILRanker: per-fold (C, τ, λ) from `members.F_feature_tower.folds` (λ None ⇒ ∞)."""
    out = {}
    for f in v04["members"]["F_feature_tower"]["folds"]:
        lam = f["lam"]
        out[int(f["fold"])] = {"C": f["C"], "tau": f["tau"], "lam": (np.inf if lam is None else lam)}
    return out


def _reconstruct_oof(frame: pd.DataFrame, per_fold: dict, make_model, member: str) -> OOFResult:
    """Re-run the exact frozen code at the stored per-fold hyperparameters on each ORIGINAL outer-train fold,
    ZERO retuning. Deterministic ⇒ a refit reproduces the frozen model. Scores only the held-out fold."""
    dev = frame[~frame["quarantined"]].copy()
    parts, specs, models = [], [], []
    for f in sorted(dev["fold"].unique()):
        tr, te = dev[dev["fold"] != f], dev[dev["fold"] == f]
        m = make_model(per_fold[int(f)]).fit(tr)
        parts.append(per_patient_metrics(te, m.raw_score(te)))
        specs.append({"fold": int(f), "member": member, **per_fold[int(f)]})
        models.append((int(f), m))
    return OOFResult(pd.concat(parts, ignore_index=True), {"folds": specs}, models)


def reconstruct_P(frame: pd.DataFrame, v03: dict) -> OOFResult:
    return _reconstruct_oof(frame, extract_P_hparams(v03),
                            lambda h: MILRanker(C=h["C"], tau=h["tau"]), "P")


def reconstruct_A(frame: pd.DataFrame, v03: dict) -> OOFResult:
    return _reconstruct_oof(frame, extract_A_hparams(v03), lambda h: AdditiveRanker(C=h["C"]), "A")


def reconstruct_F(frame: pd.DataFrame, v04: dict) -> OOFResult:
    return _reconstruct_oof(frame, extract_F_hparams(v04),
                            lambda h: TowerMILRanker(member="F", C=h["C"], tau=h["tau"], lam=h["lam"]), "F")


def _coef_max_abs_diff(refit_models: list, stored_folds: list) -> float:
    """Max |refit coef − stored coef| across folds/features. Stored coefs are 4-dp rounded in the frozen JSON."""
    by_fold = {int(f): m for f, m in refit_models}
    worst = 0.0
    for rec in stored_folds:
        m = by_fold[int(rec["fold"])]
        stored = np.array([rec["coefficients"][f] for f in FEATURES])
        worst = max(worst, float(np.max(np.abs(m.coef_ - stored))))
    return worst


def verify_convex_reconstruction(oof: OOFResult, stored_folds: list, stored_overall_hits: float,
                                 tol_coef: float = 2e-3, tol_hits: float = 1e-6) -> dict:
    """FAIL-FAST reproduction check for the CONVEX comparators P and A: per-fold refit coefficients reproduce
    the frozen 4-dp coefficients (tol_coef) AND source-balanced overall hits reproduce the frozen aggregate
    (tol_hits). Any breach raises — no silent substitution of a differently-fit model (§2.1)."""
    coef_diff = _coef_max_abs_diff(oof.models, stored_folds)
    got_hits = overall_hits(oof)
    hits_diff = abs(got_hits - stored_overall_hits)
    ok = (coef_diff <= tol_coef) and (hits_diff <= tol_hits)
    report = {"max_abs_coef_diff": round(coef_diff, 6), "overall_hits_refit": got_hits,
              "overall_hits_frozen": stored_overall_hits, "abs_hits_diff": round(hits_diff, 6),
              "tol_coef": tol_coef, "tol_hits": tol_hits, "reproduced": ok}
    if not ok:
        raise RuntimeError(f"CONVEX comparator reconstruction FAILED (frozen code did not reproduce): {report}")
    return report


def verify_f_reconstruction(oof: OOFResult, stored_overall_hits: float, tol_hits: float = 5e-3) -> dict:
    """HONEST reproduction check for the NONCONVEX comparator F: F rides v0.4's non-convex MIL log-sum-exp with
    a documented multi-init score wobble, so we verify only that source-balanced overall hits reproduce to
    v0.4's own tolerance and REPORT the residual (we do not claim tight coefficient reproduction). A gross
    breach still fails fast."""
    got_hits = overall_hits(oof)
    hits_diff = abs(got_hits - stored_overall_hits)
    ok = hits_diff <= tol_hits
    report = {"overall_hits_refit": got_hits, "overall_hits_frozen": stored_overall_hits,
              "abs_hits_diff": round(hits_diff, 6), "tol_hits": tol_hits, "nonconvex": True, "reproduced": ok,
              "note": "F is nonconvex (v0.4 multi-init wobble); hits verified to v0.4 tolerance, residual reported."}
    if not ok:
        raise RuntimeError(f"NONCONVEX comparator F reconstruction residual exceeds tolerance: {report}")
    return report


# ==================================================================================================
# provenance guard (§10) — recompute every pinned SHA256 and FAIL FAST on any mismatch.
# ==================================================================================================
def verify_provenance(prov_path: str | Path = V05_PROVENANCE) -> dict:
    """Recompute SHA256 of every pinned input and fail fast on any mismatch. Mirrors v0.4's guard, reading the
    v0.5 PROVENANCE.json. Parameterized on the path so it is unit-testable without the full corpus present."""
    prov = json.loads(Path(prov_path).read_text())
    mism = []
    for key, path in prov["input_paths"].items():
        got = _sha256(path)
        want = prov["inputs_sha256"][key]
        if got != want:
            mism.append({"key": key, "path": path, "expected": want, "got": got})
    if mism:
        raise RuntimeError(f"PROVENANCE MISMATCH (frozen inputs changed since preregistration): {mism}")
    return {"git_head": prov["git_head"], "n_inputs_verified": len(prov["input_paths"])}
