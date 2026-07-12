"""Dynamic upstream gate (Milestone 7) — safe-rejection under extreme asymmetric cost.

A LABEL-BLIND gate that removes tested-negative candidates while retaining recognized positives, then
hands survivors UNCHANGED to the frozen rankers (genuine PRIME / frozen Epicurus v0.1). It is a
selective-prediction / safe-rejection layer, NOT a classifier and NOT a reranker. See
artifacts/milestone_7_decision/dynamic_gate/SPEC.md.

Design invariants (locked by tests/test_dynamic_gate.py):
  * MISSING evidence => KEEP. A missing veto axis abstains; it never votes to remove.
  * AND-of-independent-vetoes: a candidate is veto-eligible only if EVERY present axis is confidently bad
    (percentile < t). Strong on any axis => survives => high recall by construction.
  * Uncertainty / discordance / boundary / OOD => KEEP (Layer 2 rescue overrides a veto).
  * Patient-adaptive controller (Layer 3): honors a top-M-by-EL floor, a per-patient removal cap, and
    pool/coverage floors; otherwise falls back to keep-all for that patient.
  * Cohort / study identity is NEVER an input (features are within-patient percentiles only).

Expression is used ONLY as a rescue axis inside the AND (a low-expression candidate is vetoed only if it
is ALSO low on EL and PRIME); the gate never reweights the ranker => the frozen confidence-only expression
policy (configs/frozen/expression_policy_v1.json) is preserved.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
import pandas as pd

# Core veto axes (within-patient percentile, oriented so higher = better candidate). A candidate is
# veto-eligible only if EVERY core axis is present and confidently bad. Expression is a RESCUE-ONLY axis:
# when present and strong it keeps a candidate, but its absence never causes a veto and it never vetoes on
# its own. On full-coverage cohorts (el+prime+expr all present) this is IDENTICAL to the 3-way AND
# el<t ∧ prime<t ∧ expr<t validated in FEASIBILITY.md; when expr is sparse (e.g. CheckMate) the gate still
# operates on the core el+prime axes rather than disabling itself.
CORE_VETO_AXES: tuple[str, ...] = ("el", "prime")
RESCUE_AXES: tuple[str, ...] = ("expr",)
VETO_AXES: tuple[str, ...] = CORE_VETO_AXES  # back-compat alias
ORIENT_HIGHER_BETTER = {"el": False, "prime": False, "expr": True}  # raw orientation of each feature


# --------------------------------------------------------------------------------------------------
# Within-patient percentiles that PRESERVE missingness (unlike prime_transfer._pct which fills 0.5).
# Missing must stay missing so Layer 1 can treat it as an abstain (=> keep).
# --------------------------------------------------------------------------------------------------
def within_patient_percentile(frame: pd.DataFrame, col: str, higher_better: bool) -> np.ndarray:
    """Per-patient rank(pct=True) of the oriented raw feature. Raw NaN -> NaN (kept missing)."""
    v = pd.to_numeric(frame[col], errors="coerce") if col in frame else pd.Series(np.nan, index=frame.index)
    oriented = v if higher_better else -v
    return oriented.groupby(frame["patient_id"]).rank(pct=True).to_numpy()  # NaN stays NaN


def attach_percentiles(frame: pd.DataFrame, axes: tuple[str, ...] = (*CORE_VETO_AXES, *RESCUE_AXES)) -> pd.DataFrame:
    out = frame.copy()
    for ax in axes:
        out[f"s_{ax}"] = within_patient_percentile(out, ax, ORIENT_HIGHER_BETTER[ax])
    return out


def predictor_disagreement(frame: pd.DataFrame, predictor_cols: list[str]) -> np.ndarray:
    """Spread (max-min) across the available presentation-predictor within-patient percentiles.
    Higher spread = models disagree = uncertainty. NaN where <2 predictors present for a row."""
    present = [c for c in predictor_cols if c in frame]
    if len(present) < 2:
        return np.full(len(frame), np.nan)
    pct = np.column_stack([within_patient_percentile(frame, c, higher_better=False) for c in present])
    with np.errstate(invalid="ignore"):
        n_present = np.sum(~np.isnan(pct), axis=1)
        hi = np.nanmax(pct, axis=1)
        lo = np.nanmin(pct, axis=1)
    spread = hi - lo
    spread[n_present < 2] = np.nan
    return spread


# --------------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class GateConfig:
    """Frozen operating point. `t` is the shared veto percentile threshold (Layer 1). Layer-2 rescue and
    Layer-3 safety rails are the remaining knobs. Defaults are the SAFE end (rescues on, rails permissive)."""
    t: float = 0.6
    veto_axes: tuple[str, ...] = CORE_VETO_AXES  # ALL must be present & bad to veto
    rescue_axes: tuple[str, ...] = RESCUE_AXES   # present & >= t => rescue (keep); absence never vetoes
    disagreement_rescue: float | None = None  # keep if predictor spread > this; None = off
    boundary_margin: float = 0.0              # keep if any axis within this of t; 0 = off
    el_floor_frac: float = 0.05               # always keep this top fraction by EL percentile per patient
    per_patient_cap: float = 0.90             # if a veto would remove > this frac of a pool, keep-all
    pool_floor: int = 8                       # pools smaller than this are never gated
    cov_floor: float = 0.5                    # min mean core-axis coverage to gate a patient
    retention_target: float | None = None     # provenance: the target this t was calibrated to
    calibrated_on: tuple[str, ...] = field(default_factory=tuple)
    version: str = "dynamic_gate_v1"

    def to_json(self) -> dict:
        d = asdict(self)
        d["veto_axes"] = list(self.veto_axes)
        d["rescue_axes"] = list(self.rescue_axes)
        d["calibrated_on"] = list(self.calibrated_on)
        return d

    @staticmethod
    def from_json(d: dict) -> "GateConfig":
        d = dict(d)
        d["veto_axes"] = tuple(d.get("veto_axes", CORE_VETO_AXES))
        d["rescue_axes"] = tuple(d.get("rescue_axes", RESCUE_AXES))
        d["calibrated_on"] = tuple(d.get("calibrated_on", ()))
        return GateConfig(**d)


# --------------------------------------------------------------------------------------------------
# Core gate. Operates per patient. Adds columns; never drops rows (removal = ~dyn_gate_keep).
# --------------------------------------------------------------------------------------------------
def apply_gate(
    frame: pd.DataFrame,
    config: GateConfig,
    *,
    disagreement: np.ndarray | None = None,
    percentiles_attached: bool = False,
) -> pd.DataFrame:
    """Return `frame` + columns: dyn_gate_keep(bool), dyn_gate_reason(str), dyn_veto_eligible(bool),
    dyn_rescued(bool). Label-blind: the `label` column, if present, is never read."""
    all_axes = tuple(dict.fromkeys((*config.veto_axes, *config.rescue_axes)))
    df = frame if percentiles_attached else attach_percentiles(frame, all_axes)
    df = df.copy()
    n = len(df)
    t = config.t

    # ---- Layer 1: AND-of-independent-vetoes over the CORE axes ----
    # Veto-eligible iff EVERY core axis is PRESENT and confidently bad (percentile < t). A missing core
    # axis makes the AND fail => KEEP (fail-open). Strong on any single core axis => KEEP.
    all_present = np.ones(n, dtype=bool)
    all_bad = np.ones(n, dtype=bool)
    n_present = np.zeros(n, dtype=float)
    for ax in all_axes:
        s = pd.to_numeric(df[f"s_{ax}"], errors="coerce").to_numpy()
        n_present += (~np.isnan(s)).astype(float)  # coverage counts core + rescue axes
    for ax in config.veto_axes:
        s = pd.to_numeric(df[f"s_{ax}"], errors="coerce").to_numpy()
        present = ~np.isnan(s)
        all_present &= present
        all_bad &= (present & (s < t))
    veto_eligible = all_present & all_bad

    # ---- Layer 2: rescues override a veto -> KEEP ----
    reason = np.where(veto_eligible, "VETO_LOW_ALL_AXES", "KEEP").astype(object)
    rescued = np.zeros(n, dtype=bool)

    # rescue-only axes: a candidate strong (>= t) on a present rescue axis (e.g. expression) survives.
    for ax in config.rescue_axes:
        s = pd.to_numeric(df[f"s_{ax}"], errors="coerce").to_numpy()
        r = veto_eligible & (~rescued) & (~np.isnan(s)) & (s >= t)
        rescued |= r
        reason[r] = f"RESCUE_{ax.upper()}"

    # boundary rescue: any core axis within margin of t
    if config.boundary_margin > 0:
        near = np.zeros(n, dtype=bool)
        for ax in config.veto_axes:
            s = pd.to_numeric(df[f"s_{ax}"], errors="coerce").to_numpy()
            near |= (~np.isnan(s)) & (np.abs(s - t) < config.boundary_margin)
        r = veto_eligible & near
        rescued |= r
        reason[r] = "RESCUE_NEAR_BOUNDARY"

    # disagreement rescue: predictors disagree
    if config.disagreement_rescue is not None and disagreement is not None:
        dis = np.asarray(disagreement, dtype=float)
        r = veto_eligible & (~rescued) & (~np.isnan(dis)) & (dis > config.disagreement_rescue)
        rescued |= r
        reason[r] = "RESCUE_PREDICTOR_DISAGREEMENT"

    keep = ~veto_eligible | rescued

    # ---- Layer 3: patient-adaptive safety rails ----
    df["_s_el_for_floor"] = pd.to_numeric(df.get("s_el", np.nan), errors="coerce")
    core_present_frac = (n_present / max(len(all_axes), 1))

    keep_out = keep.copy()
    reason_out = reason.copy()
    for pid, gidx in df.groupby("patient_id").groups.items():
        rows = df.index.get_indexer(gidx)
        m = len(rows)
        pool_keep = keep[rows]
        pool_cov = float(np.nanmean(core_present_frac[rows])) if m else 1.0
        # pool floor / coverage floor -> keep-all
        if m < config.pool_floor or pool_cov < config.cov_floor:
            keep_out[rows] = True
            reason_out[rows] = np.where(~pool_keep, "KEEP_POOL_FLOOR", reason[rows])
            continue
        # per-patient removal cap -> keep-all if too aggressive
        removed_frac = float(np.mean(~pool_keep)) if m else 0.0
        if removed_frac > config.per_patient_cap:
            keep_out[rows] = True
            reason_out[rows] = np.where(~pool_keep, "KEEP_CAP_EXCEEDED", reason[rows])
            continue
        # top-M-by-EL floor: always keep the presentation-best few
        if config.el_floor_frac > 0:
            s_el = df["_s_el_for_floor"].to_numpy()[rows]
            n_floor = int(np.ceil(config.el_floor_frac * m))
            if n_floor > 0 and np.any(~np.isnan(s_el)):
                order = np.argsort(-np.where(np.isnan(s_el), -np.inf, s_el), kind="mergesort")
                floor_local = order[:n_floor]
                grows = rows[floor_local]
                rescued_floor = grows[~keep_out[grows]]
                keep_out[grows] = True
                reason_out[rescued_floor] = "KEEP_EL_FLOOR"

    df = df.drop(columns=["_s_el_for_floor"])
    df["dyn_veto_eligible"] = veto_eligible
    df["dyn_rescued"] = rescued
    df["dyn_gate_keep"] = keep_out
    df["dyn_gate_reason"] = reason_out
    df["dyn_gate_config_version"] = config.version
    return df


# --------------------------------------------------------------------------------------------------
# Clopper–Pearson one-sided lower confidence bound on a binomial proportion.
# --------------------------------------------------------------------------------------------------
def clopper_pearson_lower(k: int, n: int, conf: float = 0.95) -> float:
    """One-sided lower bound on p given k successes of n. k=0 -> 0.0; k=n -> alpha**(1/n)."""
    if n <= 0:
        return 0.0
    if k <= 0:
        return 0.0
    if k >= n:
        return float((1.0 - conf) ** (1.0 / n))
    from scipy.stats import beta

    return float(beta.ppf(1.0 - conf, k, n - k + 1))


# --------------------------------------------------------------------------------------------------
# Retention / removal accounting (label used ONLY for measurement, never for gating).
# --------------------------------------------------------------------------------------------------
def gate_retention_stats(gated: pd.DataFrame, conf: float = 0.95) -> dict:
    """Positive retention + negative removal on a frame already passed through apply_gate."""
    keep = gated["dyn_gate_keep"].to_numpy(bool)
    is_pos = (gated["label"].to_numpy() == "POSITIVE")
    is_neg = (gated["label"].to_numpy() == "TESTED_NEGATIVE")
    npos, nneg = int(is_pos.sum()), int(is_neg.sum())
    pos_kept = int(np.sum(keep & is_pos))
    neg_removed = int(np.sum(~keep & is_neg))
    # per-patient retention
    per_patient = []
    for pid, g in gated.groupby("patient_id"):
        pk = g["dyn_gate_keep"].to_numpy(bool)
        gp = (g["label"].to_numpy() == "POSITIVE")
        if gp.sum() == 0:
            continue
        per_patient.append(float(np.sum(pk & gp) / gp.sum()))
    per_patient = np.array(per_patient) if per_patient else np.array([1.0])
    return {
        "n_positives": npos,
        "n_tested_negatives": nneg,
        "positive_retention": pos_kept / npos if npos else float("nan"),
        "positive_retention_cp_lb": clopper_pearson_lower(pos_kept, npos, conf) if npos else float("nan"),
        "positive_retention_per_patient_mean": float(np.mean(per_patient)),
        "positive_retention_per_patient_min": float(np.min(per_patient)),
        "n_patients_losing_a_positive": int(np.sum(per_patient < 1.0)),
        "negative_removal": neg_removed / nneg if nneg else float("nan"),
        "n_kept": int(keep.sum()),
        "n_input": int(len(gated)),
        "kept_fraction": float(keep.mean()),
    }


# --------------------------------------------------------------------------------------------------
# Calibration: pick the most aggressive t whose CP lower bound on calibration positives >= target.
# --------------------------------------------------------------------------------------------------
def _retention_at_t(calib: pd.DataFrame, config: GateConfig, disagreement: np.ndarray | None) -> dict:
    g = apply_gate(calib, config, disagreement=disagreement)
    return gate_retention_stats(g)


def calibrate_threshold(
    calib: pd.DataFrame,
    *,
    target: float,
    base_config: GateConfig | None = None,
    disagreement: np.ndarray | None = None,
    grid: np.ndarray | None = None,
    conf: float = 0.95,
) -> dict:
    """Sweep t; keep the LARGEST t whose calibration-positive retention CP lower bound >= target.
    Uses ONLY calibration positives (never eval-cohort positives). Returns the chosen t and the sweep."""
    base = base_config or GateConfig()
    grid = np.round(np.arange(0.05, 0.951, 0.05), 3) if grid is None else grid
    sweep = []
    chosen_t = 0.0  # degenerate: t=0 vetoes nothing -> retention 1.0 -> always safe
    for t in grid:
        cfg = GateConfig(**{**base.to_json(), "t": float(t), "veto_axes": tuple(base.veto_axes),
                            "calibrated_on": tuple(base.calibrated_on)})
        stats = _retention_at_t(calib, cfg, disagreement)
        row = {"t": float(t), "positive_retention": stats["positive_retention"],
               "cp_lb": stats["positive_retention_cp_lb"], "negative_removal": stats["negative_removal"]}
        sweep.append(row)
        if stats["positive_retention_cp_lb"] >= target:
            chosen_t = max(chosen_t, float(t))
    return {"target": target, "conf": conf, "chosen_t": chosen_t, "sweep": sweep}
