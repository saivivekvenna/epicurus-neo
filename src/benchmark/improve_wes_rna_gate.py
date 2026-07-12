"""Biology-first WES/RNA gate — pure analysis core (isolated milestone-7 asset).

The intervention under study is a GATE, not a reranker. Starting from the frozen
Epicurus v0.1 base order, a label-blind demotion predicate removes top-k
candidates that fail a biologically necessary presentation/clonality
prerequisite; the freed slots are backfilled strictly by the unchanged base
order. Net Δ(recognized hits@k) per patient = (positives pulled in) − (positives
removed). This module holds the accounting, within-patient partial effects, the
matched-random removal control, and the patient-paired bootstrap. It fits no
model and never consults the outcome to define a rule.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def apply_demotion_gate(
    frame: pd.DataFrame,
    *,
    k: int,
    score_col: str = "score",
    label_col: str = "pos",
    demote_col: str = "demote",
    by: str = "patient_id",
) -> dict:
    """Apply a demotion gate on top of the frozen base order and return per-patient
    Δ(hits@k). A demoted candidate is removed from the eligible pool; the new top-k
    is the first k of the score-sorted survivors (so demoted top-k slots are
    backfilled by base order, and demoted challengers are never pulled in)."""
    deltas: dict = {}
    for pid, g in frame.groupby(by):
        g = g.sort_values(score_col, ascending=False, kind="mergesort")
        lab = g[label_col].to_numpy().astype(int)
        dem = g[demote_col].to_numpy().astype(bool)
        before = int(lab[:k].sum())
        kept = lab[~dem]
        after = int(kept[:k].sum())
        deltas[pid] = after - before
    vals = np.array(list(deltas.values()), dtype=float)
    return {
        "deltas": deltas,
        "mean_delta": float(vals.mean()) if len(vals) else 0.0,
        "total_delta": int(vals.sum()),
        "n_patients": len(vals),
        "n_improved": int((vals > 0).sum()),
        "n_harmed": int((vals < 0).sum()),
    }


def within_patient_bin(
    frame: pd.DataFrame, col: str, *, n_bins: int = 4, by: str = "patient_id"
) -> pd.Series:
    """Within-patient quantile bin (0=low .. n_bins-1=high). Missing values get bin
    -1 (their own stratum) so missingness behaviour is never silently imputed."""
    def _bin(s: pd.Series) -> pd.Series:
        out = pd.Series(-1, index=s.index, dtype=int)
        ok = s.notna()
        if ok.sum() == 0:
            return out
        try:
            q = pd.qcut(s[ok].rank(method="first"), q=min(n_bins, ok.sum()), labels=False)
        except ValueError:
            q = pd.Series(0, index=s[ok].index)
        out.loc[ok.index[ok]] = q.astype(int).to_numpy()
        return out

    return frame.groupby(by)[col].transform(_bin).astype(int)


def partial_effect(
    frame: pd.DataFrame,
    col: str,
    *,
    label_col: str = "pos",
    n_bins: int = 4,
    by: str = "patient_id",
) -> list[dict]:
    """Within-patient partial effect: positive-rate + support per within-patient bin
    (bin -1 = missing). Patient-relative, so it is not confounded by per-patient
    expression scale. No smoothing, no causal claim."""
    b = within_patient_bin(frame, col, n_bins=n_bins, by=by)
    tmp = pd.DataFrame({"bin": b.to_numpy(), "pos": frame[label_col].to_numpy(), "pid": frame[by].to_numpy()})
    rows = []
    for bin_id, g in tmp.groupby("bin"):
        rows.append({
            "bin": int(bin_id),
            "pos_rate": round(float(g["pos"].mean()), 4),
            "n": int(len(g)),
            "n_pos": int(g["pos"].sum()),
            "n_patients": int(g["pid"].nunique()),
        })
    return sorted(rows, key=lambda r: r["bin"])


def matched_random_removal(
    frame: pd.DataFrame,
    *,
    k: int,
    score_col: str = "score",
    label_col: str = "pos",
    demote_col: str = "demote",
    by: str = "patient_id",
    seed: int = 0,
    reps: int = 200,
) -> dict:
    """Control: remove the SAME number of top-k candidates per patient as the real
    gate, but chosen at random, and report the mean Δ(hits@k). A real gate must beat
    this (backfilling from challengers is unfavourable when challenger pos-rate <
    top-k pos-rate, so random removal is not zero)."""
    rng = np.random.default_rng(seed)
    n_removed: dict = {}
    per_patient_means: dict = {}
    for pid, g in frame.groupby(by):
        g = g.sort_values(score_col, ascending=False, kind="mergesort").reset_index(drop=True)
        lab = g[label_col].to_numpy().astype(int)
        topk_idx = np.arange(min(k, len(g)))
        n_dem = int(g[demote_col].to_numpy().astype(bool)[:k].sum())
        n_removed[pid] = n_dem
        before = int(lab[:k].sum())
        if n_dem == 0:
            per_patient_means[pid] = 0.0
            continue
        ds = []
        for _ in range(reps):
            rm = rng.choice(topk_idx, size=n_dem, replace=False)
            mask = np.ones(len(g), dtype=bool)
            mask[rm] = False
            after = int(lab[mask][:k].sum())
            ds.append(after - before)
        per_patient_means[pid] = float(np.mean(ds))
    vals = np.array(list(per_patient_means.values()), dtype=float)
    return {
        "mean_delta": float(vals.mean()) if len(vals) else 0.0,
        "per_patient_means": per_patient_means,
        "n_removed_per_patient": n_removed,
    }


def paired_bootstrap(deltas: dict, *, reps: int = 2000, seed: int = 0, alpha: float = 0.05) -> dict:
    """Patient-equal-weight bootstrap over per-patient Δ. Resamples patients."""
    rng = np.random.default_rng(seed)
    vals = np.array(list(deltas.values()), dtype=float)
    n = len(vals)
    if n == 0:
        return {"mean": 0.0, "lo": 0.0, "hi": 0.0, "n": 0}
    means = np.array([vals[rng.integers(0, n, n)].mean() for _ in range(reps)])
    return {
        "mean": float(vals.mean()),
        "lo": float(np.quantile(means, alpha / 2)),
        "hi": float(np.quantile(means, 1 - alpha / 2)),
        "n": n,
        "frac_gt0": round(float((means > 0).mean()), 4),
    }
