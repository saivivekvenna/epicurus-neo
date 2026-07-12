"""RNA-expression ranking-policy forms + a Pareto/portfolio reserve selector.

The question this answers (label-blind, on development cohorts): should RNA expression enter neoantigen
ranking as (a) a rank penalty, (b) a confidence-only annotation, or (c) soft-saturating / route-dependent
evidence? The forms here are the candidates; `scripts/expression_policy_analysis.py` evaluates them per
development cohort (never pooled) under a strict no-regression requirement against the protected
lossless+PRIME incumbent, and freezes the winner.

All scores are oriented HIGHER = better and operate on WITHIN-PATIENT percentiles so cohorts with
different raw expression/PRIME scales are comparable. No constant here is tuned to any evaluation cohort;
the fixed thresholds (bottom-quartile stratum, unit percentile penalty) are documented defaults.
"""

from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd

# Fixed, documented constants (NOT tuned to any eval cohort).
BOTTOM_STRATUM = 0.25  # within-patient percentile boundary for the "low" expression/presentation stratum
PENALTY = 0.25  # bounded demotion (in percentile units) for the soft-saturating form

EXPRESSION_STRATA = ("expr_absent_or_low", "expr_mid", "expr_high")


def within_patient_percentile(
    frame: pd.DataFrame, col: str, *, higher_better: bool, patient_col: str = "patient_id"
) -> pd.Series:
    """Within-patient percentile of ``col`` oriented so higher = better; NaN -> 0.5."""
    v = pd.to_numeric(frame[col], errors="coerce")
    if not higher_better:
        v = -v
    return v.groupby(frame[patient_col]).rank(pct=True).fillna(0.5)


def _prime_pct(frame: pd.DataFrame) -> pd.Series:
    return within_patient_percentile(frame, "prime", higher_better=False)


def _expr_pct(frame: pd.DataFrame) -> pd.Series:
    return within_patient_percentile(frame, "expr", higher_better=True)


# ---------------------------------------------------------------------------
# (b) confidence-only == the protected PRIME incumbent
# ---------------------------------------------------------------------------
def prime_only_score(frame: pd.DataFrame) -> pd.Series:
    """Protected incumbent: rank by genuine PRIME only; expression does not move the score."""
    return _prime_pct(frame)


def expression_confidence_annotation(frame: pd.DataFrame) -> pd.Series:
    """Confidence-only role: label each candidate's expression stratum WITHOUT changing its rank."""
    ep = _expr_pct(frame)
    strata = np.where(ep <= BOTTOM_STRATUM, EXPRESSION_STRATA[0],
                      np.where(ep >= 1 - BOTTOM_STRATUM, EXPRESSION_STRATA[2], EXPRESSION_STRATA[1]))
    return pd.Series(strata, index=frame.index)


# ---------------------------------------------------------------------------
# (a) rank penalty — expression co-drives the score
# ---------------------------------------------------------------------------
def expr_penalty_score(frame: pd.DataFrame, *, weight: float = 0.5) -> pd.Series:
    """Expression as a co-equal rank driver: blend presentation and expression percentiles."""
    return (1 - weight) * _prime_pct(frame) + weight * _expr_pct(frame)


# ---------------------------------------------------------------------------
# (c) soft-saturating / route-dependent — protects strong presenters
# ---------------------------------------------------------------------------
def soft_saturating_score(frame: pd.DataFrame) -> pd.Series:
    """Demote ONLY candidates that are BOTH low-presentation AND low-expression (bounded, saturating).

    Strong presenters (above the bottom presentation stratum) are never demoted -> the PRIME top is
    protected. High expression yields no reward (saturates). Route-dependent: candidates with multi-source
    support (``n_callers``/``n_timepoints`` >= 2) are exempt from the penalty.
    """
    pp = _prime_pct(frame)
    ep = _expr_pct(frame)
    low_pres = pp <= BOTTOM_STRATUM
    low_expr = ep <= BOTTOM_STRATUM
    exempt = pd.Series(False, index=frame.index)
    for col in ("n_callers", "n_timepoints"):
        if col in frame.columns:
            exempt = exempt | (pd.to_numeric(frame[col], errors="coerce").fillna(0) >= 2)
    penalty = PENALTY * (low_pres & low_expr & ~exempt).astype(float)
    return pp - penalty


# ---------------------------------------------------------------------------
# Pareto / portfolio reserve selector
# ---------------------------------------------------------------------------
def _tie_key(peptide: object, hla: object) -> str:
    return hashlib.md5(f"{peptide}|{hla}".encode()).hexdigest()


def select_portfolio_reserved(
    frame: pd.DataFrame, *, k: int = 20, reserve: int = 2, patient_col: str = "patient_id"
) -> pd.DataFrame:
    """Select k candidates by PRIME, reserving slots for reachability across expression strata and
    predictor disagreement. Deterministic (md5 tiebreak). Operates per patient if multiple present.

    Reserves (each at most 1, only if a fresh candidate exists): the best-PRIME candidate in the bottom
    expression stratum (rescues low-expression strong presenters), and the strongest-presentation
    candidate with maximal presentation-vs-expression disagreement. Remaining slots fill by PRIME.
    NOTE: within a fixed top-k budget reservation can cost recall where expression is itself the signal
    (see the analysis); this selector is exploratory, not the frozen default.
    """
    out = []
    for _, group in frame.groupby(patient_col, sort=True):
        out.append(_select_one(group, k, reserve))
    return pd.concat(out) if out else frame.iloc[0:0]


def _select_one(group: pd.DataFrame, k: int, reserve: int) -> pd.DataFrame:
    g = group.copy()
    g["_pp"] = _prime_pct(g)
    g["_ep"] = _expr_pct(g)
    g["_tk"] = [_tie_key(p, h) for p, h in zip(g.get("mutant_peptide", ""), g.get("hla_allele", ""))]
    by_prime = g.sort_values(["_pp", "_tk"], ascending=[False, True], kind="mergesort")

    chosen = list(by_prime.index[: max(k - reserve, 0)])
    chosen_set = set(chosen)

    # reserve 1: best PRIME among the bottom expression stratum (reachability rescue)
    lo = by_prime[(by_prime["_ep"] <= BOTTOM_STRATUM) & (~by_prime.index.isin(chosen_set))]
    if len(lo) and len(chosen_set) < k:
        chosen_set.add(lo.index[0])

    # reserve 1: strongest presenter with maximal presentation-vs-expression disagreement
    rest = g[~g.index.isin(chosen_set)].copy()
    if len(rest) and len(chosen_set) < k:
        rest["_disag"] = (rest["_pp"] - rest["_ep"]).abs() * (rest["_pp"] >= 0.5)
        rest = rest.sort_values(["_disag", "_tk"], ascending=[False, True], kind="mergesort")
        chosen_set.add(rest.index[0])

    # backfill to k by PRIME
    for idx in by_prime.index:
        if len(chosen_set) >= k:
            break
        chosen_set.add(idx)

    return group.loc[sorted(chosen_set)]


# ---------------------------------------------------------------------------
# No-regression verdict
# ---------------------------------------------------------------------------
def no_regression_verdict(delta_per_patient: np.ndarray, *, n_boot: int = 2000, seed: int = 0) -> dict:
    """Paired verdict for a policy's per-patient hits@k delta vs the protected incumbent.

    ``regresses`` is True when the mean per-patient delta is negative (the policy loses recognized
    candidates on this cohort). A bootstrap CI is reported for context.
    """
    delta = np.asarray(delta_per_patient, dtype=float)
    mean = float(delta.mean()) if delta.size else 0.0
    rng = np.random.default_rng(seed)
    if delta.size:
        boots = np.array([rng.choice(delta, size=delta.size, replace=True).mean() for _ in range(n_boot)])
        lo, hi = float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))
    else:
        lo = hi = 0.0
    return {
        "mean_delta": round(mean, 4),
        "ci95": [round(lo, 4), round(hi, 4)],
        "regresses": bool(mean < 0),
        "n_patients": int(delta.size),
    }
