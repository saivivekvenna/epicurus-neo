"""Risk-controlled NEGATIVE REDUCER — non-Sid nested engine (CONTRACT.md + PROTOCOL CORRECTION 1).

A label-blind gate that CONFIDENTLY REMOVES tested-negatives at guaranteed high positive retention, then
hands survivors UNCHANGED to genuine PRIME. It never reranks toward positives. Primary target = removable
negatives subject to Clopper-Pearson 95%% lower-bound retention >= 0.95 (powered only at the DEV aggregate /
IMPROVE; small studies abstain from the CP CLAIM but the gate still applies). A bounded protected PRIME core
m in {0,5,10} (NOT a full top-20 lane) keeps Delta-hits@20 live. Models: nonnegative-coefficient logistic
(coef>=0 => monotone keep-score), monotonic shallow HGB, and NULL. Nested: outer leave-one-STUDY-out; inner
patient-grouped CV selects (model, m, tau). Peptide exact/near quarantine recomputed inside each split, for
HIT-COUNTING only (full pool always competes). Accesses NO Sid / Miller file.

This module is import-clean and unit-tested; the runner (`scripts/negative_reducer_run.py`) orchestrates it.
"""

from __future__ import annotations

import hashlib
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import beta as _beta

from event_b.leakage_registry import _kmer_index, canonical_peptide, near_duplicate

# ---- data locations (portable base CSVs for all 3 studies; rich cols joined for IMPROVE from the zip) -----
POOL = Path("artifacts/milestone_7_decision/pool_size_sensitivity")
IMPROVE_ZIP = Path("data/raw/improve/data.zip")
IMPROVE_MEMBER = "data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt"
PRIME_CACHE = Path("data/raw/gartner_nci/_cache_improve_prime.tsv")
ALLOWED_DATA_FILES = {str(POOL / "base_gartner.csv"), str(POOL / "base_multimer.csv"), str(IMPROVE_ZIP),
                      str(PRIME_CACHE)}

STUDIES = ["improve", "gartner", "multimer"]
PORTABLE = ["prime", "el", "expr"]
RICH = PORTABLE + ["VarAlFreq", "rna_af", "ValMutRNACoef", "CelPrev"]
# orientation: True => higher raw value is more recognition-favorable (=> nonneg-coef keep-score monotone up)
HIGHER_BETTER = {"prime": False, "el": False, "expr": True,
                 "VarAlFreq": True, "rna_af": True, "ValMutRNACoef": True, "CelPrev": True}
M_GRID = [0, 5, 10]
C_GRID = [0.5, 1.0, 2.0]
K = 20
NEAR = 0.8
CATASTROPHIC = -0.02
CONF = 0.95
CP_MIN_POS = 59  # min positives for CP-95%-LB>=0.95 to be attainable at 100% retention (0.05^(1/59)=0.9505)


def stable_seed(key) -> int:
    """Process-stable seed from sha256 (Python hash() is per-process randomized)."""
    return int.from_bytes(hashlib.sha256(str(key).encode()).digest()[:4], "big")


# ---- loading -----------------------------------------------------------------------------------------
def _load_improve() -> pd.DataFrame:
    with zipfile.ZipFile(IMPROVE_ZIP) as z, z.open(IMPROVE_MEMBER) as fh:
        d = pd.read_csv(fh, sep="\t")
    d = d[d["response"].isin([0, 1])].copy()
    pr = pd.read_csv(PRIME_CACHE, sep="\t").rename(
        columns={"mutant_peptide": "Mut_peptide", "hla_allele": "HLA_allele"})
    d = d.merge(pr, on=["Mut_peptide", "HLA_allele"], how="left")
    out = pd.DataFrame({
        "study": "improve", "patient_id": "improve:" + d["Patient"].astype(str),
        "mut_peptide": d["Mut_peptide"].astype(str),
        "label": np.where(d["response"] == 1, "POSITIVE", "TESTED_NEGATIVE"),
        "prime": pd.to_numeric(d["prime_rank"], errors="coerce"),
        "el": pd.to_numeric(d["RankEL"], errors="coerce"),
        "expr": pd.to_numeric(d["Expression"], errors="coerce")})
    for c in ["VarAlFreq", "rna_af", "ValMutRNACoef", "CelPrev"]:
        out[c] = pd.to_numeric(d[c], errors="coerce")
    return out


def _load_portable(study: str) -> pd.DataFrame:
    d = pd.read_csv(POOL / f"base_{study}.csv")
    d = d[d["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    out = pd.DataFrame({
        "study": study, "patient_id": study + ":" + d["patient_id"].astype(str),
        "mut_peptide": d["mutant_peptide"].astype(str), "label": d["label"].astype(str),
        "prime": pd.to_numeric(d["prime"], errors="coerce"),
        "el": pd.to_numeric(d["el"], errors="coerce"),
        "expr": pd.to_numeric(d["expr"], errors="coerce")})
    for c in ["VarAlFreq", "rna_af", "ValMutRNACoef", "CelPrev"]:
        out[c] = np.nan
    return out


def load_dev() -> pd.DataFrame:
    """Unified 3-study DEV frame (IMPROVE rich features present; Gartner/multimer portable only)."""
    df = pd.concat([_load_improve(), _load_portable("gartner"), _load_portable("multimer")],
                   ignore_index=True)
    return df.reset_index(drop=True)


# ---- features ----------------------------------------------------------------------------------------
def pct(df: pd.DataFrame, col: str) -> np.ndarray:
    """Within-patient oriented percentile of a raw feature; NaN -> 0.5 (neutral)."""
    v = pd.to_numeric(df[col], errors="coerce")
    v = v if HIGHER_BETTER[col] else -v
    return v.groupby(df["patient_id"]).rank(pct=True).fillna(0.5).to_numpy()


def feat_matrix(df: pd.DataFrame, cols) -> np.ndarray:
    return np.column_stack([pct(df, c) for c in cols])


def feat_present(df: pd.DataFrame, cols) -> np.ndarray:
    """True where EVERY model input is non-NaN (missing evidence => KEEP, never removed)."""
    ok = np.ones(len(df), bool)
    for c in cols:
        ok &= pd.to_numeric(df[c], errors="coerce").notna().to_numpy()
    return ok


def balanced_weights(df: pd.DataFrame) -> np.ndarray:
    """Per-patient inverse-size weights, class-balanced within the fit."""
    w = np.ones(len(df))
    idx = df.index.get_indexer
    for _, g in df.groupby("patient_id"):
        w[idx(g.index)] = 1.0 / len(g)
    y = (df["label"].to_numpy() == "POSITIVE")
    wp, wn = w[y].sum(), w[~y].sum()
    return np.where(y, w * (0.5 / wp if wp else 1.0), w * (0.5 / wn if wn else 1.0))


# ---- models ------------------------------------------------------------------------------------------
def fit_nnlogistic(X: np.ndarray, y: np.ndarray, w: np.ndarray, C: float):
    """Nonnegative-coefficient logistic (coef>=0 via L-BFGS-B; intercept free) on the balanced, L2-penalized
    weighted log-loss. coef>=0 => keep-score monotone nondecreasing in each recognition-favoring percentile.
    Returns (intercept, coef)."""
    n, d = X.shape

    def obj(theta):
        b0, beta = theta[0], theta[1:]
        z = X @ beta + b0
        softplus = np.log1p(np.exp(-np.abs(z))) + np.maximum(z, 0.0)
        loss = float((w * (softplus - y * z)).sum()) + float(beta @ beta) / (2.0 * C)
        p = 1.0 / (1.0 + np.exp(-z))
        g = w * (p - y)
        grad = np.empty(d + 1)
        grad[0] = float(g.sum())
        grad[1:] = X.T @ g + beta / C
        return loss, grad

    bounds = [(None, None)] + [(0.0, None)] * d
    res = minimize(obj, np.zeros(d + 1), jac=True, method="L-BFGS-B", bounds=bounds,
                   options={"maxiter": 500, "ftol": 1e-12})
    if not res.success or not np.all(np.isfinite(res.x)):
        # FAIL CLOSED: a constant keep-score (zero coef) => the gate removes nothing (KEEP-all), never garbage
        return 0.0, np.zeros(d)
    return float(res.x[0]), res.x[1:].copy()


def nnlogistic_score(X: np.ndarray, intercept: float, coef: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-(X @ coef + intercept)))


def hgb_available() -> bool:
    try:
        from sklearn.ensemble import HistGradientBoostingClassifier  # noqa: F401
        import inspect
        from sklearn.ensemble import HistGradientBoostingClassifier as H
        return "monotonic_cst" in inspect.signature(H.__init__).parameters
    except Exception:
        return False


def fit_hgb(X: np.ndarray, y: np.ndarray, w: np.ndarray, ncol: int):
    """Monotonic shallow gradient-boosted trees (depth<=3, monotone-up in every feature)."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    clf = HistGradientBoostingClassifier(max_depth=3, monotonic_cst=[1] * ncol, max_iter=100,
                                         learning_rate=0.1, l2_regularization=1.0, random_state=0)
    clf.fit(X, y, sample_weight=w)
    return clf


# ---- Clopper-Pearson one-sided lower bound -----------------------------------------------------------
def cp_lower(k: int, n: int, conf: float = CONF) -> float:
    """One-sided Clopper-Pearson lower confidence bound on a proportion k/n."""
    if n == 0 or k == 0:
        return 0.0
    if k >= n:
        return float((1.0 - conf) ** (1.0 / n))
    return float(_beta.ppf(1.0 - conf, k, n - k + 1))


# ---- protected PRIME core (m) ------------------------------------------------------------------------
def core_mask(df: pd.DataFrame, m: int) -> np.ndarray:
    """True where a row is in the protected top-m PRIME core (unremovable). Exact-PRIME-score ties at the
    m-boundary are ALL protected (core may exceed m). m=0 => protect nothing."""
    prot = np.zeros(len(df), bool)
    if m <= 0:
        return prot
    idx = df.index.get_indexer
    for _, g in df.groupby("patient_id"):
        loc = idx(g.index)
        p = pd.to_numeric(g["prime"], errors="coerce").to_numpy()
        finite = p[np.isfinite(p)]
        if len(finite) == 0:
            prot[loc] = True  # no PRIME => protect all (cannot rank => cannot safely prune)
            continue
        if len(finite) <= m:
            prot[loc] = True
            continue
        thr = np.sort(finite)[m - 1]     # m-th best PRIME (lower=better); ties at thr all protected
        prot[loc] = (p <= thr) | ~np.isfinite(p)
    return prot


# ---- OOD patient detection (raw-feature envelope; percentiles are always [0,1] so use RAW) ------------
def ood_patients(train: pd.DataFrame, test: pd.DataFrame, cols, cover: float = 0.5) -> set:
    """Patients in `test` whose RAW features fall out of the train [p1,p99] envelope for > `cover` of their
    candidates => OOD => KEEP all (no removal)."""
    env = {}
    for c in cols:
        v = pd.to_numeric(train[c], errors="coerce").to_numpy()
        v = v[np.isfinite(v)]
        env[c] = (np.percentile(v, 1), np.percentile(v, 99)) if len(v) else (-np.inf, np.inf)
    out = np.zeros(len(test), bool)
    for c in cols:
        lo, hi = env[c]
        v = pd.to_numeric(test[c], errors="coerce").to_numpy()
        out |= np.isfinite(v) & ((v < lo) | (v > hi))
    ood = set()
    for pid, g in test.groupby("patient_id"):
        loc = test.index.get_indexer(g.index)
        if out[loc].mean() > cover:
            ood.add(pid)
    return ood


# ---- gate removal + hits@20 --------------------------------------------------------------------------
def removable_mask(df: pd.DataFrame, cols, m: int, ood: set) -> np.ndarray:
    """Candidates the gate is ALLOWED to remove: outside the protected core, all features present, patient
    not OOD. (Threshold then decides which of these are actually removed.)"""
    ok = (~core_mask(df, m)) & feat_present(df, cols)
    ok &= ~df["patient_id"].isin(ood).to_numpy()
    return ok


def calibrate_tau(keepscore: np.ndarray, y: np.ndarray, removable: np.ndarray, conf: float = CONF):
    """Most aggressive keep-score cut tau (remove removable candidates with score < tau) whose CP-lower-bound
    retention over ALL positives >= 0.95. Retention loss can only come from REMOVABLE positives. Returns
    (tau, r_max, n_pos, cp_lb_at_tau, powered)."""
    n_pos = int(y.sum())
    powered = n_pos >= CP_MIN_POS
    pos_rem_scores = np.sort(keepscore[removable & (y == 1)])
    r_max = 0
    for r in range(len(pos_rem_scores) + 1):
        if cp_lower(n_pos - r, n_pos, conf) >= 0.95:
            r_max = r
        else:
            break
    if r_max >= len(pos_rem_scores):
        tau = np.inf  # can remove every removable candidate and still guarantee retention
    else:
        tau = float(pos_rem_scores[r_max])  # remove strictly below this positive's score
    cp_lb = cp_lower(n_pos - min(r_max, len(pos_rem_scores)), n_pos, conf)
    return tau, r_max, n_pos, cp_lb, powered


def gate_removed(df: pd.DataFrame, keepscore: np.ndarray, tau: float, cols, m: int, ood: set) -> np.ndarray:
    rem = removable_mask(df, cols, m, ood)
    return rem & (keepscore < tau)


def _clean_against(df: pd.DataFrame, train_peptides: set) -> np.ndarray:
    """True where a row's mut peptide does NOT exact/near-match the train peptide set (=> countable as hit).
    Quarantine is for HIT-COUNTING only; the row still competes for slots."""
    idx = _kmer_index({p for p in train_peptides if p})
    clean = np.ones(len(df), bool)
    canon = df["mut_peptide"].map(canonical_peptide).to_numpy()
    tr = {p for p in train_peptides if p}
    for i, c in enumerate(canon):
        if c and (c in tr or near_duplicate(c, idx, threshold=NEAR) is not None):
            clean[i] = False
    return clean


def hits_at_k(df: pd.DataFrame, removed: np.ndarray, clean: np.ndarray, k: int = K) -> dict:
    """Per-patient hits@k: PRIME ranks survivors (backfill highest-PRIME removed if <k survivors); count
    positives among the top-k that are leakage-CLEAN."""
    prime = pd.to_numeric(df["prime"], errors="coerce").to_numpy()
    ispos = (df["label"].to_numpy() == "POSITIVE")
    out = {}
    idx = df.index.get_indexer
    for pid, g in df.groupby("patient_id"):
        loc = idx(g.index)
        pr = prime[loc]
        order = np.argsort(np.where(np.isfinite(pr), pr, np.inf), kind="mergesort")
        rem = removed[loc]
        surv = [i for i in order if not rem[i]]
        if len(surv) >= k:
            top = surv[:k]
        else:
            top = surv + [i for i in order if rem[i]][: k - len(surv)]
        top = np.array(top[:k], dtype=int)
        gloc = loc[top]
        out[str(pid)] = float((ispos[gloc] & clean[gloc]).sum())
    return out


def matched_random_delta(df: pd.DataFrame, removed_counts: dict, cols, m: int, ood: set,
                         clean: np.ndarray, null_hits: dict, seeds=range(20)) -> float:
    """Aggregate patient-macro hits@20 delta when removing the SAME per-patient count as the gate, uniformly
    at random among that patient's removable candidates (same m/core). Averaged over seeds."""
    rem_allowed = removable_mask(df, cols, m, ood)
    idx = df.index.get_indexer
    deltas = []
    for s in seeds:
        removed = np.zeros(len(df), bool)
        for pid, g in df.groupby("patient_id"):
            c = int(removed_counts.get(str(pid), 0))
            if c <= 0:
                continue
            loc = idx(g.index)
            elig = loc[rem_allowed[loc]]
            if len(elig) == 0:
                continue
            rng = np.random.default_rng([s, stable_seed(pid)])
            pick = rng.permutation(elig)[: min(c, len(elig))]
            removed[pick] = True
        h = hits_at_k(df, removed, clean)
        deltas.append(np.mean([h[p] - null_hits[p] for p in null_hits]))
    return float(np.mean(deltas))
