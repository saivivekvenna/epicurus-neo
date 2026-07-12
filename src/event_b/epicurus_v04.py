"""Epicurus v0.4 DEVELOPMENT experiment — source-aware tower on the frozen mil_dev_split_v1.

DEVELOPMENT ONLY. Follows `artifacts/milestone_7_decision/epicurus_v04/PREREGISTERED_PROTOCOL.md` exactly and
reuses v0.3's frozen feature pipeline / evaluation verbatim (imported, not reimplemented). The Gartner TEST
holdout is never loaded/scored; an I/O guard makes any TEST path un-openable. No external claim is produced.

Model family (§3 of the protocol), one shared linear instance scorer with MIL log-sum-exp bag aggregation:

    P — pooled          f(x) = w0·x + b0                       (== corrected v0.3; loaded, not retuned)
    C — calibration     f(x) = w0·x + b0 + c_s                 (per-source intercept; rank-inert within patient)
    F — feature tower   f(x) = (w0 + v_s)·x + b0 + c_s         (per-source head v_s; the GATED candidate)

`λ` shrinks the feature heads v_s (relative pooling); `λ=∞` drops them (F→C). Per-source intercepts share the
w0 ridge (fixed, not tuned). Optimization is deterministic L-BFGS-B from a fixed init — NOT asserted convex.
"""

from __future__ import annotations

import builtins
import hashlib
import json
import zipfile
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from event_b.epicurus_v03 import (  # reuse the frozen v0.3 pipeline verbatim
    FEATURES, MILRanker, K_TOP, assemble_frame as _v03_assemble_frame, baseline_score, hits_at_k,
    load_frozen, per_patient_metrics,
)

PROVENANCE = Path("artifacts/milestone_7_decision/epicurus_v04/PROVENANCE.json")
V03_RESULT = Path("artifacts/milestone_7_decision/epicurus_v03/DEV_RESULT.json")
FORBIDDEN_TEST_TOKEN = "TestingSet"   # any path containing this (Nmers/Mmps TestingSet) is a Gartner TEST file


# ==================================================================================================
# guardrails: TEST is never opened (not merely filtered), and pinned provenance is re-verified
# ==================================================================================================
@contextmanager
def guard_no_test_io():
    """Patch every file-opening entrypoint our loaders use so that opening any Gartner TEST path raises.
    Covers builtins.open, pandas.read_csv/read_excel, and zipfile.ZipFile (the C parser bypasses builtins.open,
    so patching read_csv is necessary, not sufficient-by-open alone)."""
    def _check(path):
        try:
            s = str(path)
        except Exception:
            return
        if FORBIDDEN_TEST_TOKEN in s:
            raise RuntimeError(f"BLOCKED Gartner TEST I/O attempt on {s!r} during a v0.4 development run.")

    orig_open, orig_csv, orig_xl, orig_zip = builtins.open, pd.read_csv, pd.read_excel, zipfile.ZipFile

    def open_g(file, *a, **k): _check(file); return orig_open(file, *a, **k)
    def csv_g(fp, *a, **k): _check(fp); return orig_csv(fp, *a, **k)
    def xl_g(io, *a, **k): _check(io); return orig_xl(io, *a, **k)
    def zip_g(file, *a, **k): _check(file); return orig_zip(file, *a, **k)

    builtins.open, pd.read_csv, pd.read_excel, zipfile.ZipFile = open_g, csv_g, xl_g, zip_g
    try:
        yield
    finally:
        builtins.open, pd.read_csv, pd.read_excel, zipfile.ZipFile = orig_open, orig_csv, orig_xl, orig_zip


def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_provenance() -> dict:
    """Recompute SHA256 of every pinned input and fail fast on any mismatch (§11)."""
    prov = json.loads(PROVENANCE.read_text())
    mism = []
    for key, path in prov["input_paths"].items():
        got = _sha256(path)
        want = prov["inputs_sha256"][key]
        if got != want:
            mism.append({"key": key, "path": path, "expected": want, "got": got})
    if mism:
        raise RuntimeError(f"PROVENANCE MISMATCH (frozen inputs changed since preregistration): {mism}")
    return {"git_head": prov["git_head"], "n_inputs_verified": len(prov["input_paths"])}


# ==================================================================================================
# frame assembly (delegates to the frozen v0.3 pipeline, inside the TEST-I/O guard)
# ==================================================================================================
def assemble_frame() -> pd.DataFrame:
    with guard_no_test_io():
        return _v03_assemble_frame()


# ==================================================================================================
# eligibility / attrition — LABEL-BLIND (§8). Uses only core-feature presence, never outcome labels.
# ==================================================================================================
CORE_SCORES = ["prime_rank", "mix_rank", "el_strength"]   # the fail-closed core (see v0.3 assemble_frame)


def rankable_patients(frame: pd.DataFrame) -> set:
    """Patients with ≥1 core-complete candidate. Reads ONLY core-score columns — no eval_positive/bag_label."""
    ok = frame[CORE_SCORES].notna().all(axis=1)
    return set(frame.loc[ok, "patient_id"].unique())


def attrition_report(frame: pd.DataFrame) -> dict:
    """Explain feature-bearing → rankable (label-blind) → scored (has in-pool positive), by source. Rankable and
    scored are computed on the NON-QUARANTINED evaluation pool (the pool the gate uses), so `scored` matches the
    reported n_scored_patients; the drop from recurrent-peptide quarantine is reported separately."""
    patient_fold, _, _ = load_frozen()
    feature_bearing = set(patient_fold)                          # 152 patients in the frozen split
    ev = frame[~frame["quarantined"]]                            # gate evaluates here
    rankable = rankable_patients(ev)                             # label-blind
    scored = {pid for pid, g in ev.groupby("patient_id") if int(g["eval_positive"].sum()) > 0}
    scored_incl_quar = {pid for pid, g in frame.groupby("patient_id") if int(g["eval_positive"].sum()) > 0}

    def src_of(pid: str) -> str:
        return pid.split(":", 1)[0]

    rows = {}
    for src in sorted({src_of(p) for p in feature_bearing}):
        fb = {p for p in feature_bearing if src_of(p) == src}
        rk = {p for p in rankable if src_of(p) == src}
        sc = {p for p in scored if src_of(p) == src}
        lost_q = {p for p in scored_incl_quar if src_of(p) == src} - sc
        fs = ev[ev["source"] == src]
        rows[src] = {
            "feature_bearing": len(fb), "rankable_label_blind": len(rk), "scored_has_positive": len(sc),
            "lost_to_quarantine_only": len(lost_q),
            "core_availability_rate": round(float(fs[CORE_SCORES].notna().all(axis=1).mean()), 4) if len(fs) else None,
            "prime_mask_rate": round(float(fs["prime_masked"].mean()), 4) if len(fs) else None,
            "prime_rank_available_rate": round(float(fs["prime_rank"].notna().mean()), 4) if len(fs) else None,
        }
    rows["TOTAL"] = {"feature_bearing": len(feature_bearing), "rankable_label_blind": len(rankable),
                     "scored_has_positive": len(scored),
                     "lost_to_quarantine_only": len(scored_incl_quar - scored)}
    return rows


# ==================================================================================================
# the source-aware tower
# ==================================================================================================
def _head_source(df: pd.DataFrame) -> pd.Series:
    """Which source indexes the heads/intercepts. Normally the true source; the negative-control diagnostic
    supplies a shuffled `head_source` column to test whether lift is capacity rather than real structure."""
    return df["head_source"] if "head_source" in df.columns else df["source"]


@dataclass
class TowerMILRanker:
    """Partial-pooling MIL ranker. member ∈ {P, C, F}. λ shrinks the feature heads (F only). Balancing/eval
    always use the true `source`; heads use `head_source` (== source unless a negative control overrides it)."""
    member: str = "F"
    C: float = 0.3
    tau: float = 1.0
    lam: float = 1.0
    mean_: np.ndarray = None
    std_: np.ndarray = None
    w0_: np.ndarray = None
    b0_: float = 0.0
    c_: dict = field(default_factory=dict)     # source -> intercept
    v_: dict = field(default_factory=dict)     # source -> head deviation vector
    sources_: list = field(default_factory=list)

    # -- feature scaling fit strictly on the training rows --
    def _std(self, X, fit):
        if fit:
            self.mean_ = X.mean(0); self.std_ = X.std(0) + 1e-9
        return (X - self.mean_) / self.std_

    def _flags(self):
        has_c = self.member in ("C", "F")
        has_v = (self.member == "F") and np.isfinite(self.lam)
        return has_c, has_v

    def fit(self, train, init=None):
        has_c, has_v = self._flags()
        lab = train[train["bag_label"].isin(["POSITIVE", "NEGATIVE"])].copy().reset_index(drop=True)
        self.sources_ = sorted(lab["source"].unique())
        s_index = {s: i for i, s in enumerate(self.sources_)}
        S = len(self.sources_)

        # sort instances so each bag is a contiguous segment (vectorized segmented log-sum-exp)
        codes = pd.factorize(lab["bag_id"])[0]
        order = np.argsort(codes, kind="stable")
        lab = lab.iloc[order].reset_index(drop=True)
        Xs = self._std(lab[FEATURES].to_numpy(float), fit=True)
        codes_s = codes[order]
        boundaries = np.concatenate([[0], np.flatnonzero(np.diff(codes_s)) + 1])
        counts = np.diff(np.append(boundaries, len(codes_s)))
        yb, wb = MILRanker._bag_targets(lab, boundaries)            # reuse v0.3 source×patient×bag balancing
        src_i = _head_source(lab).map(s_index).to_numpy()
        masks = [src_i == s for s in range(S)]                      # S=3 -> cheap per-source vector ops

        n, nbag, tau = Xs.shape[1], len(boundaries), self.tau
        reg = 1.0 / (self.C * nbag)                                 # identical ridge scale to v0.3 (0.5·reg·‖·‖²)
        rep = lambda a: np.repeat(a, counts)

        # theta layout: [w0(n), b0(1), c(S)?, V(S*n)?]
        nc = S if has_c else 0
        nv = S * n if has_v else 0

        def unpack(theta):
            w0 = theta[:n]; b0 = theta[n]
            c = theta[n + 1:n + 1 + nc] if has_c else np.zeros(S)
            V = theta[n + 1 + nc:].reshape(S, n) if has_v else np.zeros((S, n))
            return w0, b0, c, V

        def negll(theta):
            w0, b0, c, V = unpack(theta)
            s = Xs @ w0 + b0
            if has_c:
                cvec = np.empty(len(s))
                for si, m in enumerate(masks): cvec[m] = c[si]
                s = s + cvec
            if has_v:
                for si, m in enumerate(masks):
                    if m.any(): s[m] = s[m] + Xs[m] @ V[si]
            seg_max = np.maximum.reduceat(s, boundaries)
            e = np.exp((s - rep(seg_max)) / tau)
            seg_sum = np.add.reduceat(e, boundaries)
            lse = seg_max + tau * np.log(seg_sum / counts)
            p = 1.0 / (1.0 + np.exp(-lse))
            loss = float(np.sum(wb * -(yb * np.log(p + 1e-12) + (1 - yb) * np.log(1 - p + 1e-12))))
            pen = 0.5 * reg * (w0 @ w0)
            if has_c: pen += 0.5 * reg * (c @ c)
            if has_v: pen += 0.5 * reg * self.lam * float(np.sum(V * V))

            dl_bag = wb * (p - yb)
            sw = e / rep(seg_sum)
            ds = rep(dl_bag) * sw
            grad = np.empty(n + 1 + nc + nv)
            grad[:n] = Xs.T @ ds + reg * w0
            grad[n] = ds.sum()
            if has_c:
                grad[n + 1:n + 1 + nc] = np.array([ds[m].sum() for m in masks]) + reg * c
            if has_v:
                gV = np.stack([(ds[m, None] * Xs[m]).sum(0) if m.any() else np.zeros(n) for m in masks])
                grad[n + 1 + nc:] = (gV + reg * self.lam * V).ravel()
            return loss + pen, grad

        theta0 = np.zeros(n + 1 + nc + nv) if init is None else np.asarray(init, float)
        res = minimize(negll, theta0, jac=True, method="L-BFGS-B", options={"maxiter": 500, "ftol": 1e-9})
        w0, b0, c, V = unpack(res.x)
        self.w0_, self.b0_ = w0, float(b0)
        self.c_ = {s: float(c[i]) for i, s in enumerate(self.sources_)} if has_c else {}
        self.v_ = {s: V[i].copy() for i, s in enumerate(self.sources_)} if has_v else {}
        self._n = n
        return self

    def raw_score(self, df):
        X = self._std(df[FEATURES].to_numpy(float), fit=False)
        s = X @ self.w0_ + self.b0_
        hs = _head_source(df).to_numpy()
        if self.c_:
            s = s + np.array([self.c_.get(x, 0.0) for x in hs])          # unseen source -> 0 intercept
        if self.v_:
            zero = np.zeros(self.w0_.shape[0])
            Vrows = np.stack([self.v_.get(x, zero) for x in hs])         # unseen source -> shared backbone only
            s = s + np.einsum("ij,ij->i", X, Vrows)
        return s

    def effective_weights(self) -> dict:
        """Per-source effective weights w_s = w0 + v_s (identifiable) and deviation norms ‖v_s‖."""
        out = {"shared_w0": {f: round(float(w), 4) for f, w in zip(FEATURES, self.w0_)}, "per_source": {}}
        for s in self.sources_:
            v = self.v_.get(s, np.zeros_like(self.w0_))
            out["per_source"][s] = {
                "w_s": {f: round(float(w), 4) for f, w in zip(FEATURES, self.w0_ + v)},
                "dev_norm": round(float(np.linalg.norm(v)), 4),
                "intercept": round(float(self.c_.get(s, 0.0)), 4),
            }
        return out

    def to_dict(self):
        return {"kind": "tower_MIL", "member": self.member, "C": self.C, "tau": self.tau,
                "lam": (None if np.isinf(self.lam) else self.lam), "sources": list(self.sources_),
                "effective_weights": self.effective_weights()}


# ==================================================================================================
# extra per-patient diagnostics (NON-selected; hits@20 remains the sole gate) — §10
# ==================================================================================================
def _ndcg_at_k(order_pos: np.ndarray, k: int) -> float:
    gains = order_pos[:k].astype(float)
    disc = 1.0 / np.log2(np.arange(2, 2 + len(gains)))
    dcg = float((gains * disc).sum())
    ideal = float((np.sort(order_pos)[::-1][:k] * disc).sum())
    return dcg / ideal if ideal > 0 else 0.0


def ext_metrics(df: pd.DataFrame, score: np.ndarray, k: int = K_TOP) -> pd.DataFrame:
    """Per-patient hits@k, recall@k, best-positive rank (1-indexed), nDCG@k. Diagnostics only."""
    df = df.assign(_s=score)
    rows = []
    for (src, pid), g in df.groupby(["source", "patient_id"]):
        ep = g["eval_positive"].to_numpy()
        npos = int(ep.sum())
        if npos == 0:
            continue
        order = np.argsort(-g["_s"].to_numpy(), kind="stable")
        ranked_pos = ep[order]
        h = int(ranked_pos[:k].sum())
        best_rank = int(np.argmax(ranked_pos) + 1)               # position of the first positive
        rows.append({"source": src, "patient_id": pid, "n_pos": npos, "pool": len(g),
                     "hits": h, "recall": h / npos, "best_pos_rank": best_rank,
                     "ndcg": round(_ndcg_at_k(ranked_pos, k), 4)})
    return pd.DataFrame(rows)
