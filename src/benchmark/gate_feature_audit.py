"""Pure computational core for the gate feature audit (isolated milestone-7 asset).

Motivation. The dynamic gate (``configs/frozen/dynamic_gate_v1.json``) was
FALSIFIED because a label-blind *presentation* gate removes 0% of the
high-presentation decoys that outrank positives: on the stratum of candidates
that already look like strong presenters, presentation features carry no
residual signal (that is the recognition wall). The only way a gate can lift
downstream hits@20 is an ORTHOGONAL feature that separates the top-ranked
TESTED_NEGATIVE decoys from POSITIVES *within that high-presentation stratum*.

This module provides the label-aware measurement of exactly that quantity, plus
coverage/cross-fit helpers. It is import-clean and touches no dynamic_gate file.
It NEVER exposes study identity as a deployable feature — study identity is only
audited as a confound elsewhere in the runner.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

POSITIVE = "POSITIVE"
TESTED_NEGATIVE = "TESTED_NEGATIVE"
UNTESTED = "UNTESTED"

# Only these two labels enter any discrimination metric; UNTESTED is never a
# negative (it is an unmeasured candidate, not a measured non-responder).
SCORABLE = (POSITIVE, TESTED_NEGATIVE)


def high_presentation_mask(
    frame: pd.DataFrame,
    score_col: str,
    *,
    higher_better: bool,
    quantile: float = 0.5,
    top_k: int | None = None,
    by: str = "patient_id",
) -> pd.Series:
    """Boolean mask for the best-presenting candidates within each patient — the
    stratum a presentation gate would retain and where top-ranked decoys outrank
    positives.

    Threshold is either the top ``quantile`` fraction (default) or, when
    ``top_k`` is given, the top-``top_k`` per patient (matching the gate's k=20).
    ``higher_better`` orients the score (PRIME score = higher-better; an EL
    percentile *rank* = lower-better). At least the single best is always kept.
    """
    s = frame[score_col].astype(float)
    if not higher_better:
        s = -s
    # Descending rank within patient: rank 1 = best presenter.
    order = s.groupby(frame[by]).rank(ascending=False, method="first")
    if top_k is not None:
        keep_n = float(top_k)
    else:
        size = frame.groupby(by)[score_col].transform("size")
        keep_n = np.ceil(size * quantile).clip(lower=1)
    return (order <= keep_n).astype(bool)


def _scorable(frame: pd.DataFrame, label_col: str) -> pd.DataFrame:
    return frame[frame[label_col].isin(SCORABLE)]


def conditional_auroc(
    frame: pd.DataFrame,
    feature_col: str,
    *,
    higher_better: bool,
    label_col: str = "label",
    mask: pd.Series | None = None,
    min_per_class: int = 3,
) -> dict:
    """AUROC of ``feature_col`` for POSITIVE-vs-TESTED_NEGATIVE, optionally
    restricted to a stratum ``mask`` (e.g. the high-presentation stratum).

    Returns coverage (non-null fraction over scorable rows), class counts, and
    ``auroc`` oriented so higher = better discrimination of POSITIVES. Returns
    ``auroc=None`` when a class is too small or the feature is all-null.
    """
    sub = frame if mask is None else frame[mask.to_numpy()]
    sub = _scorable(sub, label_col)
    n_scorable = int(len(sub))
    coverage = float(sub[feature_col].notna().mean()) if n_scorable and feature_col in sub else 0.0

    if feature_col not in sub:
        return {"auroc": None, "n": n_scorable, "n_pos": 0, "n_neg": 0, "coverage": 0.0}

    ok = sub[feature_col].notna()
    sub = sub[ok]
    y = (sub[label_col] == POSITIVE).to_numpy().astype(int)
    n_pos, n_neg = int(y.sum()), int((1 - y).sum())
    if n_pos < min_per_class or n_neg < min_per_class:
        return {"auroc": None, "n": n_scorable, "n_pos": n_pos, "n_neg": n_neg, "coverage": round(coverage, 4)}

    x = sub[feature_col].to_numpy(dtype=float)
    if not higher_better:
        x = -x
    auroc = _auroc(y, x)
    return {
        "auroc": None if auroc is None else round(auroc, 4),
        "n": n_scorable,
        "n_pos": n_pos,
        "n_neg": n_neg,
        "coverage": round(coverage, 4),
    }


def grouped_oof_auroc(
    frame: pd.DataFrame,
    feature_cols: list[str],
    *,
    group_col: str = "patient_id",
    label_col: str = "label",
    mask: pd.Series | None = None,
    n_splits: int = 5,
    seed: int = 0,
) -> dict:
    """Cross-fitted (patient-grouped OOF) logistic AUROC over ``feature_cols``.

    Measures whether the feature set JOINTLY carries signal that generalises
    across patients (not memorised study/patient identity). Restricted to
    POSITIVE/TESTED_NEGATIVE, optionally within a stratum ``mask``. Median
    imputation + standardisation are fit on train folds only.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import GroupKFold

    sub = frame if mask is None else frame[mask.to_numpy()]
    sub = _scorable(sub, label_col)
    present = [c for c in feature_cols if c in sub.columns]
    if not present:
        return {"oof_auroc": None, "n": int(len(sub)), "features": [], "reason": "no features present"}

    sub = sub.dropna(subset=[group_col])
    y = (sub[label_col] == POSITIVE).to_numpy().astype(int)
    groups = sub[group_col].to_numpy()
    X = sub[present].to_numpy(dtype=float)

    n_groups = len(np.unique(groups))
    if y.sum() < n_splits or (len(y) - y.sum()) < n_splits or n_groups < n_splits:
        return {"oof_auroc": None, "n": int(len(sub)), "features": present, "reason": "insufficient groups/classes"}

    oof = np.full(len(y), np.nan)
    gkf = GroupKFold(n_splits=n_splits)
    for tr, te in gkf.split(X, y, groups):
        med = np.nanmedian(X[tr], axis=0)
        med = np.where(np.isnan(med), 0.0, med)
        Xtr, Xte = _impute(X[tr], med), _impute(X[te], med)
        mu, sd = Xtr.mean(0), Xtr.std(0)
        sd = np.where(sd == 0, 1.0, sd)
        Xtr, Xte = (Xtr - mu) / sd, (Xte - mu) / sd
        if y[tr].sum() == 0 or y[tr].sum() == len(tr):
            continue
        clf = LogisticRegression(max_iter=1000, C=1.0)
        clf.fit(Xtr, y[tr])
        oof[te] = clf.predict_proba(Xte)[:, 1]

    valid = ~np.isnan(oof)
    if valid.sum() < 2 * n_splits or len(np.unique(y[valid])) < 2:
        return {"oof_auroc": None, "n": int(len(sub)), "features": present, "reason": "degenerate OOF"}
    auroc = _auroc(y[valid], oof[valid])
    return {
        "oof_auroc": None if auroc is None else round(auroc, 4),
        "n": int(len(sub)),
        "n_pos": int(y.sum()),
        "features": present,
    }


def within_patient_variation(frame: pd.DataFrame, col: str, *, group_col: str = "patient_id") -> float:
    """Fraction of patients in which ``col`` takes >1 distinct value across the
    patient's candidates. ~1.0 = candidate-varying (rankable within a patient);
    ~0.0 = patient/sample-constant (context-only, cannot re-order candidates).

    This is the deployability axis for a within-patient gate: a feature that is
    constant across a patient's candidates can shift the whole patient's prior
    but can never remove one decoy while sparing a positive.
    """
    if col not in frame.columns or group_col not in frame.columns:
        return 0.0
    varies = frame.groupby(group_col)[col].nunique(dropna=True) > 1
    return round(float(varies.mean()), 4)


def feature_coverage(frame: pd.DataFrame, columns: list[str]) -> dict:
    """Non-null fraction per column; absent columns report 0.0 (never crash)."""
    out = {}
    n = max(len(frame), 1)
    for c in columns:
        out[c] = round(float(frame[c].notna().sum()) / n, 4) if c in frame.columns else 0.0
    return out


def _impute(x: np.ndarray, med: np.ndarray) -> np.ndarray:
    x = x.copy()
    idx = np.where(np.isnan(x))
    x[idx] = np.take(med, idx[1])
    return x


def _auroc(y: np.ndarray, score: np.ndarray) -> float | None:
    from sklearn.metrics import roc_auc_score

    if len(np.unique(y)) < 2:
        return None
    return float(roc_auc_score(y, score))
