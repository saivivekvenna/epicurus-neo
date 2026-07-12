"""Pool-size sensitivity diagnostic (M7).

Question: does frozen Epicurus vs genuine PRIME (and a pVAC-style+PRIME proxy) improve when the
starting candidate pool shrinks but retains exactly all known positives?

Reuses the FROZEN deterministic scorers only — no model fit, no PRIME binary is ever invoked:

  * genuine PRIME     = -prime  (raw GfellerLab %rank; lower raw = better) -> POOL-INVARIANT per candidate.
  * frozen Epicurus   = event_b.prime_transfer.score_with_frozen (configs/frozen/epicurus_v0_1.json):
                        within-patient-percentile logistic residual -> POOL-DEPENDENT, because the
                        percentiles recompute over whatever pool is scored.
  * pVAC-style+PRIME  = equal-weight mean of within-patient percentile(NetMHCpan-EL) and
                        percentile(PRIME) (both lower-raw-better). pVACseq ranks by best binding; the
                        "+PRIME" grafts immunogenicity. This is a TUNING-FREE APPROXIMATION of the
                        pVACtools binding-first ranking augmented with PRIME, NOT the pVACtools binary
                        (no pVAC scorer exists in this repo). Auxiliary arm, flagged as a proxy.

Two variants of nested pools per patient (SMALL subset MEDIUM subset LARGE; identical positives; identical
per-candidate raw features prime/el/expr):

  A (oracle) : all positives + a deterministic random subsample of negatives (>= 20 seeds). Reranker
               stress test ONLY. Positive retention is 100% BY CONSTRUCTION (an oracle keeps every
               positive) -> diagnostic, NOT north-star validation.
  B (gate)   : top-N by a LABEL-BLIND presentation ordering (within-patient NetMHCpan-EL percentile),
               pool-size matched to A's target counts. Positives compete and MAY be dropped; retention
               is REPORTED, never assumed. Labels are never used to design B.

Cohorts are kept strictly separate (heterogeneous roles; never pooled). Multimer is the frozen model's
TRAINING cohort, so frozen-Epicurus scores there are IN-SAMPLE (optimistic) and flagged as such; genuine
PRIME and the pVAC proxy are unaffected.
"""

from __future__ import annotations

import zlib
from dataclasses import dataclass

import numpy as np
import pandas as pd

from event_b.prime_transfer import _pct, score_with_frozen

# ---- configuration -------------------------------------------------------------------------------
POOL_FRACS = {"LARGE": 1.0, "MEDIUM": 0.5, "SMALL": 0.25}  # fraction of NEGATIVES retained (positives: all)
KS = (5, 10, 20)
BUDGET_FRAC = 0.10  # normalized budget: top 10% of the pool
FIXED_K = 5         # fixed small k (hits@5 doubles as the "fixed k below the smallest pool" report)
MIN_NEG = 4         # a patient needs >= this many tested negatives to define a 25% subsample
ARMS = ("genuine_prime", "frozen_epicurus", "pvac_style_prime")
GATE_COL = "gate_presentation_el_pct"  # label-blind ordering for variant B
IN_SAMPLE_COHORTS = frozenset({"multimer"})  # frozen Epicurus was trained here


# ---- deterministic per-(cohort,patient,seed) RNG -------------------------------------------------
def _rng(cohort: str, patient: str, seed: int) -> np.random.Generator:
    """Stable, decorrelated RNG. crc32 keeps determinism across processes (unlike hash())."""
    tag = zlib.crc32(f"{cohort}:{patient}".encode()) & 0x7FFFFFFF
    return np.random.default_rng(np.random.SeedSequence([int(seed), int(tag)]))


# ---- arm scoring (re-scored WITHIN each pool) ----------------------------------------------------
def score_arms(pool: pd.DataFrame) -> pd.DataFrame:
    """Attach all arm score columns to a pool frame, recomputing every pool-dependent quantity over
    exactly the rows present. Higher = better for every arm."""
    pool = pool.copy()
    pool["genuine_prime"] = -pd.to_numeric(pool["prime"], errors="coerce")
    pool["frozen_epicurus"] = score_with_frozen(pool)
    el_pct = _pct(pool, "el", False)      # lower raw EL rank = better -> higher percentile = better
    prime_pct = _pct(pool, "prime", False)
    pool["pvac_style_prime"] = 0.5 * el_pct + 0.5 * prime_pct
    return pool


def gate_score(pat_df: pd.DataFrame) -> np.ndarray:
    """Label-blind presentation ordering used to define variant-B pools (higher = kept first)."""
    return _pct(pat_df, "el", False)


# ---- nested pool builders (per patient) ----------------------------------------------------------
def oracle_pools(pat_df: pd.DataFrame, cohort: str, patient: str, seed: int) -> dict[str, pd.DataFrame]:
    """Variant A. All positives + prefix subsamples of a single deterministic negative shuffle.
    Prefixes of one shuffle => SMALL subset MEDIUM subset LARGE, guaranteed."""
    pos = pat_df[pat_df["label"] == "POSITIVE"]
    neg = pat_df[pat_df["label"] == "TESTED_NEGATIVE"]
    nneg = len(neg)
    order = _rng(cohort, patient, seed).permutation(nneg)
    neg_shuf = neg.iloc[order]
    pools: dict[str, pd.DataFrame] = {}
    for name, frac in POOL_FRACS.items():
        k = nneg if frac == 1.0 else int(np.ceil(frac * nneg))
        pools[name] = pd.concat([pos, neg_shuf.iloc[:k]], axis=0)
    return pools


def gate_pools(pat_df: pd.DataFrame) -> tuple[dict[str, pd.DataFrame], np.ndarray]:
    """Variant B. Top-N by the label-blind gate ordering, pool size matched to variant A's target
    counts. Prefixes of one ordering => nested. Positives may be dropped."""
    npos = int((pat_df["label"] == "POSITIVE").sum())
    nneg = int((pat_df["label"] == "TESTED_NEGATIVE").sum())
    g = pat_df.assign(**{GATE_COL: gate_score(pat_df)})
    ordered = g.sort_values(GATE_COL, ascending=False, kind="mergesort", na_position="last")
    pools: dict[str, pd.DataFrame] = {}
    for name, frac in POOL_FRACS.items():
        target = len(g) if frac == 1.0 else min(npos + int(np.ceil(frac * nneg)), len(g))
        pools[name] = ordered.iloc[:target]
    return pools, ordered[GATE_COL].to_numpy()


# ---- per-patient ranking metrics -----------------------------------------------------------------
def patient_metrics(pool: pd.DataFrame, score_col: str) -> dict:
    """Rank a scored pool by one arm and compute the diagnostic metrics. NaN scores sort last."""
    r = pool.sort_values(score_col, ascending=False, kind="mergesort", na_position="last")
    is_pos = (r["label"].to_numpy() == "POSITIVE")
    n = len(r)
    npos = int(is_pos.sum())
    ranks = np.flatnonzero(is_pos) + 1  # 1-indexed positions of positives
    out: dict = {"pool_size": n, "n_pos": npos, "saturated_le20": bool(n <= 20)}
    for k in KS:
        out[f"hits@{k}"] = int(is_pos[: min(k, n)].sum())
    kb = max(1, int(np.ceil(BUDGET_FRAC * n)))
    out[f"hits@budget{int(BUDGET_FRAC * 100)}pct"] = int(is_pos[:kb].sum())
    out["budget_k"] = kb
    out[f"hits@fixed{FIXED_K}"] = int(is_pos[: min(FIXED_K, n)].sum())
    k20 = min(20, n)
    out["recall@20"] = out["hits@20"] / npos if npos else float("nan")
    out["precision@20"] = out["hits@20"] / k20 if k20 else float("nan")
    out["mrr"] = float(1.0 / ranks[0]) if ranks.size else 0.0
    out["median_pos_rank"] = float(np.median(ranks)) if ranks.size else float("nan")
    return out


# ---- eligibility ---------------------------------------------------------------------------------
@dataclass
class Eligibility:
    eligible: list[str]
    excluded: dict[str, str]


def patient_eligibility(frame: pd.DataFrame) -> Eligibility:
    eligible, excluded = [], {}
    for pid, gp in frame.groupby("patient_id"):
        npos = int((gp["label"] == "POSITIVE").sum())
        nneg = int((gp["label"] == "TESTED_NEGATIVE").sum())
        if npos < 1:
            excluded[str(pid)] = "no_positive"
        elif nneg < MIN_NEG:
            excluded[str(pid)] = f"too_few_negatives(<{MIN_NEG})"
        else:
            eligible.append(str(pid))
    return Eligibility(sorted(eligible), excluded)


# ---- run one cohort ------------------------------------------------------------------------------
_METRIC_KEYS = (
    [f"hits@{k}" for k in KS]
    + [f"hits@budget{int(BUDGET_FRAC * 100)}pct", f"hits@fixed{FIXED_K}", "recall@20", "precision@20", "mrr", "median_pos_rank"]
)


def _mean(vals: list[float]) -> float | None:
    v = [x for x in vals if x is not None and not (isinstance(x, float) and np.isnan(x))]
    return float(np.mean(v)) if v else None


def run_variant_a(frame: pd.DataFrame, cohort: str, elig: list[str], seeds: list[int]) -> dict:
    """Oracle nested pools over seeds. Returns per-(pool,arm) aggregate + paired epi-minus-prime delta,
    plus the tidy long rows (per seed x pool x arm, patient-averaged)."""
    long_rows: list[dict] = []
    # per-seed patient-mean metric, indexed [pool][arm][metric] -> list over seeds
    agg: dict = {p: {a: {m: [] for m in _METRIC_KEYS} for a in ARMS} for p in POOL_FRACS}
    sat: dict = {p: [] for p in POOL_FRACS}          # per (seed) mean saturation flag
    delta20: dict = {p: [] for p in POOL_FRACS}      # per-seed mean paired (epi - prime) hits@20
    delta_recall: dict = {p: [] for p in POOL_FRACS}
    by_pid = {pid: gp for pid, gp in frame.groupby("patient_id")}
    for seed in seeds:
        seed_metric = {p: {a: {m: [] for m in _METRIC_KEYS} for a in ARMS} for p in POOL_FRACS}
        seed_sat = {p: [] for p in POOL_FRACS}
        seed_d20 = {p: [] for p in POOL_FRACS}
        seed_dr = {p: [] for p in POOL_FRACS}
        for pid in elig:
            pools = oracle_pools(by_pid[pid], cohort, pid, seed)
            for pname, pool in pools.items():
                scored = score_arms(pool)
                per_arm = {a: patient_metrics(scored, a) for a in ARMS}
                seed_sat[pname].append(1.0 if per_arm[ARMS[0]]["saturated_le20"] else 0.0)
                for a in ARMS:
                    for m in _METRIC_KEYS:
                        seed_metric[pname][a][m].append(per_arm[a][m])
                seed_d20[pname].append(per_arm["frozen_epicurus"]["hits@20"] - per_arm["genuine_prime"]["hits@20"])
                seed_dr[pname].append(per_arm["frozen_epicurus"]["recall@20"] - per_arm["genuine_prime"]["recall@20"])
        for pname in POOL_FRACS:
            sat[pname].append(_mean(seed_sat[pname]))
            delta20[pname].append(_mean(seed_d20[pname]))
            delta_recall[pname].append(_mean(seed_dr[pname]))
            for a in ARMS:
                for m in _METRIC_KEYS:
                    pm = _mean(seed_metric[pname][a][m])
                    agg[pname][a][m].append(pm)
                    long_rows.append({"cohort": cohort, "variant": "A_oracle", "pool": pname, "arm": a,
                                      "seed": seed, "metric": m, "value": pm})
    # collapse across seeds
    summary: dict = {}
    for pname in POOL_FRACS:
        summary[pname] = {
            "mean_saturated_le20_frac": _mean(sat[pname]),
            "epicurus_minus_prime_hits@20": {
                "mean": _mean(delta20[pname]),
                "std_over_seeds": float(np.std([x for x in delta20[pname] if x is not None])) if delta20[pname] else None,
                "band95": _percentile_band(delta20[pname]),
            },
            "epicurus_minus_prime_recall@20": {
                "mean": _mean(delta_recall[pname]),
                "band95": _percentile_band(delta_recall[pname]),
            },
            "arms": {a: {m: _mean(agg[pname][a][m]) for m in _METRIC_KEYS} for a in ARMS},
        }
    return {"summary": summary, "long_rows": long_rows}


def _percentile_band(vals: list) -> list | None:
    v = [x for x in vals if x is not None and not (isinstance(x, float) and np.isnan(x))]
    if len(v) < 2:
        return None
    return [round(float(np.percentile(v, 2.5)), 4), round(float(np.percentile(v, 97.5)), 4)]


def run_variant_b(frame: pd.DataFrame, cohort: str, elig: list[str]) -> dict:
    """Deterministic label-blind gate nested pools. Reports positive retention + reranked metrics."""
    long_rows: list[dict] = []
    agg: dict = {p: {a: {m: [] for m in _METRIC_KEYS} for a in ARMS} for p in POOL_FRACS}
    retention: dict = {p: [] for p in POOL_FRACS}         # per-patient positives kept / positives total
    all_retained: dict = {p: True for p in POOL_FRACS}
    sat: dict = {p: [] for p in POOL_FRACS}
    delta20: dict = {p: [] for p in POOL_FRACS}
    for pid in elig:
        gp = frame[frame["patient_id"] == pid]
        npos_total = int((gp["label"] == "POSITIVE").sum())
        pools, _ = gate_pools(gp)
        for pname, pool in pools.items():
            kept_pos = int((pool["label"] == "POSITIVE").sum())
            frac_kept = kept_pos / npos_total if npos_total else float("nan")
            retention[pname].append(frac_kept)
            if kept_pos < npos_total:
                all_retained[pname] = False
            scored = score_arms(pool)
            per_arm = {a: patient_metrics(scored, a) for a in ARMS}
            sat[pname].append(1.0 if per_arm[ARMS[0]]["saturated_le20"] else 0.0)
            delta20[pname].append(per_arm["frozen_epicurus"]["hits@20"] - per_arm["genuine_prime"]["hits@20"])
            for a in ARMS:
                for m in _METRIC_KEYS:
                    agg[pname][a][m].append(per_arm[a][m])
                    long_rows.append({"cohort": cohort, "variant": "B_gate", "pool": pname, "arm": a,
                                      "seed": -1, "metric": m, "value": per_arm[a][m]})
    summary: dict = {}
    for pname in POOL_FRACS:
        summary[pname] = {
            "positive_retention_mean": _mean(retention[pname]),
            "positive_retention_min": float(np.min(retention[pname])) if retention[pname] else None,
            "all_positives_retained": bool(all_retained[pname]),
            "n_patients_losing_a_positive": int(sum(1 for x in retention[pname] if x < 1.0)),
            "mean_saturated_le20_frac": _mean(sat[pname]),
            "epicurus_minus_prime_hits@20_mean": _mean(delta20[pname]),
            "arms": {a: {m: _mean(agg[pname][a][m]) for m in _METRIC_KEYS} for a in ARMS},
        }
    return {"summary": summary, "long_rows": long_rows}


def run_cohort(frame: pd.DataFrame, cohort: str, seeds: list[int]) -> dict:
    frame = frame.copy()
    frame["patient_id"] = frame["patient_id"].astype(str)  # robust to CSV round-trip / numeric ids
    elig = patient_eligibility(frame)
    fr = frame[frame["patient_id"].isin(elig.eligible)].copy()
    a = run_variant_a(fr, cohort, elig.eligible, seeds)
    b = run_variant_b(fr, cohort, elig.eligible)
    return {
        "cohort": cohort,
        "epicurus_in_sample": cohort in IN_SAMPLE_COHORTS,
        "n_patients_eligible": len(elig.eligible),
        "excluded_patients": elig.excluded,
        "n_positives": int((fr["label"] == "POSITIVE").sum()),
        "n_tested_negatives": int((fr["label"] == "TESTED_NEGATIVE").sum()),
        "variant_A_oracle": a["summary"],
        "variant_B_gate": b["summary"],
        "_long_rows": a["long_rows"] + b["long_rows"],
    }
