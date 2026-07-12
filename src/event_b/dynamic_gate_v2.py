"""Dynamic gate v2 — budgeted candidate RESELECTION for NET top-20 utility (Milestone 7).

Unlike v1 (a safe pruner that removes only low-ranked negatives and so cannot change top-20 — see
CIRCULARITY_AUDIT.md), v2 estimates NEGATIVE RISK among candidates that OCCUPY OR THREATEN the top-20 and
removes the highest-risk ones so that lower-ranked POSITIVES can backfill. It may sacrifice a positive if
doing so backfills more positives (NET utility, not safe retention).

To avoid re-introducing v1's circularity, the negative-risk model is built to add signal BEYOND the
frozen-Epicurus rank: it uses within-patient percentiles of expression and the EL/PRIME DISCORDANCE and
their interaction — i.e. the part of decoy-ness not already captured by the PRIME-dominated ranker.
EL/PRIME percentiles enter only through discordance/interaction, never as a monotone veto. Study identity
is NEVER a feature. Missing features default to LOW risk (=> keep).

The frozen ranker is applied UNCHANGED to survivors; v2 only chooses the candidate set.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from event_b.dynamic_gate import within_patient_percentile
from event_b.pool_size_sensitivity import score_arms

K = 20


# --------------------------------------------------------------------------------------------------
def add_v2_features(frame: pd.DataFrame) -> pd.DataFrame:
    """Within-patient percentile features designed to be ORTHOGONAL to the PRIME-dominated ranker."""
    f = frame.copy()
    f["p_el"] = within_patient_percentile(f, "el", higher_better=False)
    f["p_prime"] = within_patient_percentile(f, "prime", higher_better=False)
    f["p_expr"] = within_patient_percentile(f, "expr", higher_better=True)
    f["discord"] = np.abs(f["p_el"] - f["p_prime"])          # EL/PRIME disagreement (rank-neutral)
    f["expr_x_pres"] = f["p_expr"] * np.maximum(f["p_el"], f["p_prime"])  # low expr among high presentation
    return f


V2_FEATURES = ["p_expr", "discord", "expr_x_pres"]  # NOT p_el/p_prime directly (avoid ranker circularity)


@dataclass
class RiskModel:
    clf: LogisticRegression
    mu: np.ndarray
    sd: np.ndarray
    features: list[str]

    def risk(self, frame: pd.DataFrame) -> np.ndarray:
        f = add_v2_features(frame)
        X = f[self.features].to_numpy(float)
        present = np.isfinite(X)
        Xf = np.where(present, X, self.mu)  # missing -> mean -> ~neutral; then forced low-risk below
        Z = (Xf - self.mu) / self.sd
        r = self.clf.predict_proba(Z)[:, 1]  # P(TESTED_NEGATIVE)
        # any candidate missing ALL orthogonal features -> lowest risk (never removed on absent evidence)
        r[~present.any(axis=1)] = 0.0
        return r


def fit_negative_risk(train: pd.DataFrame, features: list[str] = V2_FEATURES) -> RiskModel:
    f = add_v2_features(train)
    X = f[features].to_numpy(float)
    mu = np.nanmean(X, axis=0)
    X = np.where(np.isfinite(X), X, mu)
    sd = X.std(axis=0)
    sd = np.where(sd == 0, 1.0, sd)
    y = (train["label"].to_numpy() == "TESTED_NEGATIVE").astype(int)  # predict NEGATIVE
    clf = LogisticRegression(max_iter=2000, class_weight="balanced").fit((X - mu) / sd, y)
    return RiskModel(clf=clf, mu=mu, sd=sd, features=features)


# --------------------------------------------------------------------------------------------------
def reselect(frame: pd.DataFrame, risk: np.ndarray, *, budget_frac: float, threat_k: int = 2 * K) -> np.ndarray:
    """Per patient: among the top `threat_k` by frozen Epicurus (the top-20 + its threat zone), remove up
    to `budget_frac * threat_k` candidates with the HIGHEST negative risk. Returns a keep mask.
    Positives are NOT protected — removal is by risk only (net-utility objective)."""
    keep = np.ones(len(frame), dtype=bool)
    scored = score_arms(frame)
    for pid, gp in scored.groupby("patient_id"):
        local = frame.index.get_indexer(gp.index)
        order = np.argsort(-gp["frozen_epicurus"].to_numpy(), kind="mergesort")
        zone = local[order[:threat_k]]
        n_remove = int(np.floor(budget_frac * len(zone)))
        if n_remove <= 0:
            continue
        zrisk = risk[zone]
        drop = zone[np.argsort(-zrisk, kind="mergesort")[:n_remove]]
        keep[drop] = False
    return keep


def patient_hits20(frame: pd.DataFrame, keep: np.ndarray) -> dict[str, float]:
    """Per-patient hits@20 after frozen-Epicurus reranks the SURVIVORS (percentiles recompute)."""
    sub = frame[keep]
    out = {}
    for pid, gp in sub.groupby("patient_id"):
        r = score_arms(gp).sort_values("frozen_epicurus", ascending=False, kind="mergesort")
        out[str(pid)] = float((r["label"].to_numpy()[:K] == "POSITIVE").sum())
    # patients emptied by the gate score 0
    for pid in frame["patient_id"].astype(str).unique():
        out.setdefault(pid, 0.0)
    return out


@dataclass
class RiskEnsemble:
    models: list[RiskModel]

    def q_with_uncertainty(self, frame: pd.DataFrame, z: float = 1.28) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """q = P(POSITIVE) per candidate: ensemble mean and z-sigma lower/upper bounds (z=1.28 ~ 80%).
        Each member returns P(TESTED_NEGATIVE); q_positive = 1 - that. Missing-all-features -> q=0, wide band."""
        preds = np.column_stack([1.0 - m.risk(frame) for m in self.models])
        mean = preds.mean(axis=1)
        std = preds.std(axis=1)
        return mean, np.clip(mean - z * std, 0, 1), np.clip(mean + z * std, 0, 1)


def fit_risk_ensemble(train: pd.DataFrame, n: int = 12, seed: int = 0) -> RiskEnsemble:
    """Bootstrap ensemble of negative-risk models for calibrated uncertainty (used for LCB/UCB)."""
    models = []
    pids = train["patient_id"].unique()
    for b in range(n):
        rng = np.random.default_rng([seed, b])
        boot_pids = rng.choice(pids, size=len(pids), replace=True)  # patient-level bootstrap
        boot = pd.concat([train[train["patient_id"] == p] for p in boot_pids], ignore_index=True)
        if (boot["label"] == "POSITIVE").sum() == 0 or (boot["label"] == "TESTED_NEGATIVE").sum() == 0:
            continue
        models.append(fit_negative_risk(boot))
    if not models:
        models = [fit_negative_risk(train)]
    return RiskEnsemble(models=models)


def counterfactual_reselect(frame: pd.DataFrame, ens: RiskEnsemble, *, max_budget: int = 8,
                            conservative: bool = True, gain_margin: float = 0.0) -> np.ndarray:
    """Counterfactual top-20 replacement policy. For each patient, the current frozen-Epicurus top-20 are
    removal candidates and ranks 21..20+m are their replacements under UNCHANGED reranking. Remove the m
    lowest-q top-20 candidates ONLY while the replacements are expected to add more recognized hits than
    the removed — conservatively, sum(LCB of replacement q) > sum(UCB of removed q). Choose the largest
    such m in 0..B; m=0 = abstain (weak/OOD signal => do nothing). Returns a keep mask."""
    keep = np.ones(len(frame), dtype=bool)
    mean, lcb, ucb = ens.q_with_uncertainty(frame)
    scored = score_arms(frame)
    for pid, gp in scored.groupby("patient_id"):
        local = frame.index.get_indexer(gp.index)
        order = np.argsort(-gp["frozen_epicurus"].to_numpy(), kind="mergesort")
        ranked_local = local[order]
        top = ranked_local[:K]
        backfill = ranked_local[K:K + max_budget]
        if len(top) == 0 or len(backfill) == 0:
            continue
        # removal candidates: lowest-q within the top-20 first
        rem_order = top[np.argsort(mean[top], kind="mergesort")]
        best_m = 0
        for m in range(1, min(max_budget, len(backfill), len(top)) + 1):
            removed = rem_order[:m]
            replaced = backfill[:m]
            if conservative:
                gain = float(lcb[replaced].sum() - ucb[removed].sum())
            else:
                gain = float(mean[replaced].sum() - mean[removed].sum())
            if gain > gain_margin:
                best_m = m  # keep extending while net expected (conservative) gain clears the margin
            else:
                break
        if best_m > 0:
            keep[rem_order[:best_m]] = False
    return keep


def utility(delta: np.ndarray, *, lam: float, mu: float) -> float:
    """Harm-penalized net utility: mean Δhits@20 - lam*mean(negative patient delta) - mu*fraction harmed."""
    mean_delta = float(np.mean(delta))
    neg = delta[delta < 0]
    mean_neg = float(np.mean(neg)) if len(neg) else 0.0  # negative number
    frac_harmed = float(np.mean(delta < 0))
    return mean_delta + lam * mean_neg - mu * frac_harmed  # lam*mean_neg is a subtraction (mean_neg<0)
