"""Dynamic upstream gate — evaluation runner (Milestone 7).

Loads the frozen deterministic per-candidate features (prime/el/expr) for the three DEV cohorts
(gartner, improve, multimer) from the pool-size base CSVs, plus the LOCKED external cohort CheckMate 153,
and evaluates the layered safe-rejection gate (src/event_b/dynamic_gate.py):

  1. Leave-one-cohort-out (LOCO) calibration of the veto threshold t at retention targets
     {0.90, 0.95, 0.975, 0.99} (Clopper-Pearson lower bound on calibration positives) -> Pareto frontier.
  2. Gate metrics per cohort: positive retention (mean/min/worst + CP lower bound), negative removal,
     patients losing a positive.
  3. Downstream consequence: gate -> survivors -> UNCHANGED rankers (genuine PRIME, frozen Epicurus v0.1);
     paired recognized hits@20 / recall@20, ungated (LARGE) vs gated, bootstrap CI + verdict.
  4. Baselines: pure-EL gate at matched removal, deterministic-only (Layer 0), random removal matched,
     keep-all, full-retention oracle ceiling.
  5. Freeze ONE config (target 0.95, LOCO across all three dev cohorts) -> configs/frozen/dynamic_gate_v1.json,
     then score CheckMate 153 ONCE with it.

    python -m scripts.dynamic_gate

Writes artifacts/milestone_7_decision/dynamic_gate/{dynamic_gate.json, REPORT.md, pareto.csv, per_patient.csv}.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from benchmark.stats import paired_bootstrap
from event_b.dynamic_gate import (
    GateConfig,
    apply_gate,
    attach_percentiles,
    calibrate_threshold,
    gate_retention_stats,
)
from event_b.leakage_registry import canonical_peptide
from event_b.pool_size_sensitivity import (
    patient_eligibility,
    patient_metrics,
    score_arms,
)

POOL = Path("artifacts/milestone_7_decision/pool_size_sensitivity")
ART = Path("artifacts/milestone_7_decision/dynamic_gate")
CHECKMATE = Path("data/processed/checkmate153.normalized.csv")
FROZEN = Path("configs/frozen/dynamic_gate_v1.json")
DEV_COHORTS = ["gartner", "improve", "multimer"]
IN_SAMPLE = {"multimer"}
TARGETS = [0.90, 0.95, 0.975, 0.99]
DOWNSTREAM_ARMS = ["genuine_prime", "frozen_epicurus"]
K = 20
CONF = 0.95
COMMAND = "python -m scripts.dynamic_gate"

_VERDICT = (
    "**B for the gate as a safe recall-preserving pruner; C for the aggressive top-20 premise on "
    "label-blind presentation features.**\n\n"
    "The layered gate is real, safe, and Pareto-dominates the incumbent EL-percentile gate on positive "
    "retention at matched removal (FEASIBILITY.md). It removes presentation-weak negatives at a calibrated "
    "recall floor, so its honest value is **shrinking the candidate universe handed to expensive "
    "downstream steps (genuine-PRIME scoring, wet-lab) without losing recognized positives** — NOT lifting "
    "top-20.\n\n"
    "It does NOT close the oracle gap. The threshold sweep shows that at any SAFE retention (≈1.0) the "
    "paired downstream Δhits@20 is ≈0: the gate only ever removes negatives that already sat below top-20. "
    "The oracle's large lift (Gartner 0.808→1.652) came from deleting *random* negatives — including the "
    "HIGH-presentation decoys that outrank positives — while keeping positives by construction. A "
    "label-blind gate removes **0%** of those high-presentation decoys (they are indistinguishable from "
    "positives: presentation is the ceiling), so the reranking benefit is unreachable this way. Any "
    "apparent lift arrives only once retention has fallen (denominator effect, bought by dropping "
    "positives) — which the safety bar forbids. This re-derives the project's recognition wall from a new "
    "angle: the gate/oracle decomposition.\n\n"
    "CP retention lower bounds are additionally sample-size-capped (Gartner 46 positives ⇒ max CP-LB "
    "0.937 even at 100% retention), so small cohorts cannot certify ≥0.95 regardless of the gate — an "
    "underpowering limit, not an unsafety.\n\n"
    "**What unlocks the next level:** an orthogonal signal that separates true positives from "
    "high-presentation decoys — **mutant-allele RNA VAF + read support, tumor DNA VAF/depth/purity/CCF, "
    "proteasomal processing, agretopicity** — consumed as KEEP-only rescue axes (never imputed to a "
    "vetoing value). These are absent across the current eval cohorts; the open WES+RNA of Miller IPV "
    "(PRJNA980652) and Gartner reconstruction are the concrete path (see SPEC §7)."
)


# --------------------------------------------------------------------------------------------------
def _load_dev(name: str) -> pd.DataFrame:
    f = pd.read_csv(POOL / f"base_{name}.csv")
    f["patient_id"] = f["patient_id"].astype(str)
    f = f[f["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    elig = patient_eligibility(f)
    return f[f["patient_id"].isin(elig.eligible)].reset_index(drop=True)


def _load_checkmate() -> pd.DataFrame:
    f = pd.read_csv(CHECKMATE)
    f["patient_id"] = ("cm153:" + f["patient_id"].astype(str))
    f = f[f["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    keep = [c for c in ["patient_id", "mutant_peptide", "hla_allele", "label", "prime", "el", "expr"] if c in f]
    f = f[keep]
    elig = patient_eligibility(f)
    return f[f["patient_id"].isin(elig.eligible)].reset_index(drop=True)


def _canon_set(frame: pd.DataFrame) -> set[str]:
    return {canonical_peptide(p) for p in frame.get("mutant_peptide", pd.Series([], dtype=str)).astype(str)} - {""}


def _deleak_calibration(calib: pd.DataFrame, eval_frame: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Drop calibration rows whose canonical peptide exactly matches any eval-cohort peptide."""
    if "mutant_peptide" not in calib:
        return calib, {"dropped_exact_overlap": 0, "kept": len(calib)}
    eval_peps = _canon_set(eval_frame)
    canon = calib["mutant_peptide"].astype(str).map(canonical_peptide)
    drop = canon.isin(eval_peps) & (canon != "")
    return calib[~drop].reset_index(drop=True), {"dropped_exact_overlap": int(drop.sum()), "kept": int((~drop).sum())}


# --------------------------------------------------------------------------------------------------
# Downstream: score a pool with the frozen rankers and return per-patient hits@20 / recall@20.
# --------------------------------------------------------------------------------------------------
def _downstream_per_patient(frame: pd.DataFrame, keep_mask: np.ndarray | None) -> dict[str, dict[str, list]]:
    """For each ranking arm, per-patient hits@20 and recall@20 over the (optionally gated) pool.
    Percentile-based arms recompute within whatever survives -> realistic deployment."""
    out = {a: {"hits": [], "recall": [], "pid": []} for a in DOWNSTREAM_ARMS}
    sub = frame if keep_mask is None else frame[keep_mask]
    for pid, gp in sub.groupby("patient_id"):
        if (gp["label"] == "POSITIVE").sum() == 0:
            # patient's positives all removed -> hits/recall = 0 (charged honestly)
            for a in DOWNSTREAM_ARMS:
                out[a]["hits"].append(0.0)
                out[a]["recall"].append(0.0)
                out[a]["pid"].append(pid)
            continue
        scored = score_arms(gp)
        for a in DOWNSTREAM_ARMS:
            m = patient_metrics(scored, a)
            out[a]["hits"].append(float(m["hits@20"]))
            out[a]["recall"].append(float(m["recall@20"]) if not np.isnan(m["recall@20"]) else 0.0)
            out[a]["pid"].append(pid)
    return out


def _paired_downstream(frame: pd.DataFrame, keep_mask: np.ndarray) -> dict:
    """Paired per-patient hits@20 / recall@20: gated vs ungated (LARGE). Positive delta = gating helps."""
    # patients that lose ALL positives to the gate still appear (hits/recall=0) so the harm is counted.
    ungated = _downstream_per_patient(frame, None)
    # align gated onto the ungated patient set (a patient fully removed by the gate contributes 0s)
    gated_hits_by_pid = {a: {} for a in DOWNSTREAM_ARMS}
    gated_recall_by_pid = {a: {} for a in DOWNSTREAM_ARMS}
    gsub = frame[keep_mask]
    gated = _downstream_per_patient(gsub, None) if len(gsub) else {a: {"hits": [], "recall": [], "pid": []} for a in DOWNSTREAM_ARMS}
    for a in DOWNSTREAM_ARMS:
        for pid, h, r in zip(gated[a]["pid"], gated[a]["hits"], gated[a]["recall"]):
            gated_hits_by_pid[a][pid] = h
            gated_recall_by_pid[a][pid] = r
    res = {}
    for a in DOWNSTREAM_ARMS:
        pids = ungated[a]["pid"]
        ung_h = np.array(ungated[a]["hits"])
        ung_r = np.array(ungated[a]["recall"])
        gat_h = np.array([gated_hits_by_pid[a].get(p, 0.0) for p in pids])
        gat_r = np.array([gated_recall_by_pid[a].get(p, 0.0) for p in pids])
        ci_h = paired_bootstrap(gat_h, ung_h)
        ci_r = paired_bootstrap(gat_r, ung_r)
        res[a] = {
            "ungated_hits@20_mean": round(float(ung_h.mean()), 4),
            "gated_hits@20_mean": round(float(gat_h.mean()), 4),
            "hits@20_delta_gated_minus_ungated": round(float(ci_h.delta), 4),
            "hits@20_delta_ci": [round(ci_h.lo, 4), round(ci_h.hi, 4)],
            "ungated_recall@20_mean": round(float(ung_r.mean()), 4),
            "gated_recall@20_mean": round(float(gat_r.mean()), 4),
            "recall@20_delta": round(float(ci_r.delta), 4),
            "recall@20_delta_ci": [round(ci_r.lo, 4), round(ci_r.hi, 4)],
            "no_regression": bool(ci_h.hi >= 0.0),  # gated not significantly worse than ungated
        }
    return res


# --------------------------------------------------------------------------------------------------
# Baseline gates (all label-blind), each returning a keep-mask matched to the dynamic gate's removal.
# --------------------------------------------------------------------------------------------------
def _el_gate_matched(frame: pd.DataFrame, remove_per_patient: dict[str, int]) -> np.ndarray:
    """Remove the SAME count per patient as the dynamic gate, by lowest within-patient EL percentile."""
    p = attach_percentiles(frame)
    keep = np.ones(len(frame), dtype=bool)
    for pid, gp in p.groupby("patient_id"):
        rows = frame.index.get_indexer(gp.index)
        nrem = remove_per_patient.get(pid, 0)
        if nrem <= 0:
            continue
        s_el = gp["s_el"].to_numpy()
        order = np.argsort(np.where(np.isnan(s_el), np.inf, s_el), kind="mergesort")  # lowest EL first
        keep[rows[order[:nrem]]] = False
    return keep


def _random_gate_matched(frame: pd.DataFrame, remove_per_patient: dict[str, int], seed: int = 0) -> np.ndarray:
    keep = np.ones(len(frame), dtype=bool)
    for pid, gp in frame.groupby("patient_id"):
        rows = frame.index.get_indexer(gp.index)
        neg_local = np.where((gp["label"] == "TESTED_NEGATIVE").to_numpy())[0]
        nrem = min(remove_per_patient.get(pid, 0), len(neg_local))
        if nrem <= 0:
            continue
        rng = np.random.default_rng([seed, abs(hash(pid)) % (2**31)])
        pick = rng.permutation(neg_local)[:nrem]
        keep[rows[pick]] = False
    return keep


def _remove_counts(frame: pd.DataFrame, keep: np.ndarray) -> dict[str, int]:
    out = {}
    for pid, gp in frame.groupby("patient_id"):
        rows = frame.index.get_indexer(gp.index)
        out[pid] = int(np.sum(~keep[rows]))
    return out


# --------------------------------------------------------------------------------------------------
SWEEP_TS = [0.3, 0.5, 0.6, 0.7, 0.75, 0.85, 0.95]


def _threshold_sweep(name: str, frame: pd.DataFrame) -> list[dict]:
    """The decisive diagnostic: as the veto threshold rises, does the gate EVER improve downstream top-20
    at a SAFE retention, or does apparent lift only arrive once positives are being dropped?"""
    ungated = _downstream_per_patient(frame, None)
    rows = []
    for t in SWEEP_TS:
        cfg = GateConfig(t=t, el_floor_frac=0.0, per_patient_cap=1.0, pool_floor=0, cov_floor=0.0)
        g = apply_gate(frame, cfg)
        keep = g["dyn_gate_keep"].to_numpy(bool)
        st = gate_retention_stats(g, conf=CONF)
        gated = _downstream_per_patient(frame[keep], None) if keep.any() else {a: {"hits": [], "pid": []} for a in DOWNSTREAM_ARMS}
        gh = {a: {p: h for p, h in zip(gated[a]["pid"], gated[a]["hits"])} for a in DOWNSTREAM_ARMS}
        row = {"t": t, "negative_removal": round(st["negative_removal"], 4),
               "positive_retention": round(st["positive_retention"], 4),
               "positive_retention_cp_lb": round(st["positive_retention_cp_lb"], 4)}
        for a in DOWNSTREAM_ARMS:
            pids = ungated[a]["pid"]
            ung = np.array(ungated[a]["hits"])
            gat = np.array([gh[a].get(p, 0.0) for p in pids])
            ci = paired_bootstrap(gat, ung)
            row[f"d_hits@20_{a}"] = round(float(ci.delta), 4)
            row[f"d_hits@20_{a}_ci"] = [round(ci.lo, 4), round(ci.hi, 4)]
        rows.append(row)
    return rows


def _oracle_lift_location(frame: pd.DataFrame, t: float = 0.6) -> dict:
    """Fraction of HIGH- vs LOW-presentation negatives the gate removes. The oracle's top-20 lift comes
    from deleting high-presentation decoys (top-EL negatives that outrank positives) — which a label-blind
    presentation gate provably cannot touch."""
    p = attach_percentiles(frame)
    cfg = GateConfig(t=t, el_floor_frac=0.0, per_patient_cap=1.0, pool_floor=0, cov_floor=0.0)
    removed = ~apply_gate(frame, cfg)["dyn_gate_keep"].to_numpy(bool)
    isneg = (frame["label"].to_numpy() == "TESTED_NEGATIVE")
    s_el = p["s_el"].to_numpy()
    hi = isneg & (s_el > 0.75)
    lo = isneg & (s_el < 0.25)
    return {"at_t": t,
            "high_presentation_neg_removed": round(float(np.mean(removed[hi])), 4) if hi.any() else None,
            "low_presentation_neg_removed": round(float(np.mean(removed[lo])), 4) if lo.any() else None,
            "n_high_presentation_neg": int(hi.sum()), "n_low_presentation_neg": int(lo.sum())}


def _evaluate_cohort(name: str, eval_frame: pd.DataFrame, config: GateConfig) -> dict:
    """Apply a (already-calibrated) gate config to one cohort; gate metrics + downstream + baselines."""
    gated = apply_gate(eval_frame, config)
    keep = gated["dyn_gate_keep"].to_numpy(bool)
    stats = gate_retention_stats(gated, conf=CONF)
    remove_counts = _remove_counts(eval_frame, keep)

    downstream = _paired_downstream(eval_frame, keep)

    # baselines at matched removal
    el_keep = _el_gate_matched(eval_frame, remove_counts)
    el_stats = gate_retention_stats(
        eval_frame.assign(dyn_gate_keep=el_keep), conf=CONF) if "label" in eval_frame else {}
    rnd_keep = _random_gate_matched(eval_frame, remove_counts)
    rnd_stats = gate_retention_stats(eval_frame.assign(dyn_gate_keep=rnd_keep), conf=CONF)

    reason_counts = gated.loc[~keep, "dyn_gate_reason"].value_counts().to_dict()
    return {
        "cohort": name,
        "in_sample": name in IN_SAMPLE,
        "config_t": config.t,
        "gate": stats,
        "reason_counts_removed": {str(k): int(v) for k, v in reason_counts.items()},
        "downstream_vs_ungated": downstream,
        "threshold_sweep": _threshold_sweep(name, eval_frame),
        "oracle_lift_location": _oracle_lift_location(eval_frame),
        "baseline_el_gate_matched": {"positive_retention": el_stats.get("positive_retention"),
                                     "positive_retention_cp_lb": el_stats.get("positive_retention_cp_lb"),
                                     "negative_removal": el_stats.get("negative_removal"),
                                     "n_patients_losing_a_positive": el_stats.get("n_patients_losing_a_positive")},
        "baseline_random_matched": {"positive_retention": rnd_stats.get("positive_retention"),
                                    "negative_removal": rnd_stats.get("negative_removal")},
        "oracle_ceiling": {"positive_retention": 1.0, "note": "by construction; diagnostic ceiling only"},
    }


def _pareto(dev: dict[str, pd.DataFrame]) -> list[dict]:
    """LOCO Pareto frontier: for each eval cohort and target, calibrate t on the OTHER dev cohorts, apply."""
    rows = []
    for target in TARGETS:
        for eval_name in DEV_COHORTS:
            calib = pd.concat([dev[c] for c in DEV_COHORTS if c != eval_name], ignore_index=True)
            calib, leak = _deleak_calibration(calib, dev[eval_name])
            cal = calibrate_threshold(calib, target=target,
                                      base_config=GateConfig(retention_target=target,
                                                             calibrated_on=tuple(c for c in DEV_COHORTS if c != eval_name)))
            cfg = GateConfig(t=cal["chosen_t"], retention_target=target,
                             calibrated_on=tuple(c for c in DEV_COHORTS if c != eval_name))
            gated = apply_gate(dev[eval_name], cfg)
            s = gate_retention_stats(gated, conf=CONF)
            rows.append({
                "target": target, "eval_cohort": eval_name, "in_sample": eval_name in IN_SAMPLE,
                "loco_calibrated_t": cal["chosen_t"], "calib_deleak_dropped": leak["dropped_exact_overlap"],
                "positive_retention": round(s["positive_retention"], 4),
                "positive_retention_cp_lb": round(s["positive_retention_cp_lb"], 4),
                "positive_retention_worst_patient": round(s["positive_retention_per_patient_min"], 4),
                "negative_removal": round(s["negative_removal"], 4),
                "n_patients_losing_a_positive": s["n_patients_losing_a_positive"],
                "meets_bar_cp_lb>=0.95": bool(s["positive_retention_cp_lb"] >= 0.95),
            })
    return rows


# --------------------------------------------------------------------------------------------------
def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    dev = {name: _load_dev(name) for name in DEV_COHORTS}

    # ---- Pareto frontier (LOCO) ----
    pareto = _pareto(dev)

    # ---- Freeze ONE deployment config: target 0.95, calibrated LOCO-style on all three dev cohorts ----
    # (deployment sees no eval cohort; calibrating on the full dev union is the honest frozen choice.)
    calib_all = pd.concat(dev.values(), ignore_index=True)
    cal95 = calibrate_threshold(calib_all, target=0.95,
                                base_config=GateConfig(retention_target=0.95, calibrated_on=tuple(DEV_COHORTS)))
    frozen_cfg = GateConfig(t=cal95["chosen_t"], retention_target=0.95, calibrated_on=tuple(DEV_COHORTS))
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(json.dumps({
        "name": "dynamic_gate", "version": "1.0.0", "frozen_on": "2026-07-12",
        "description": "Layered safe-rejection upstream gate: Layer0 deterministic impossibility "
                       "(apply_deterministic_gate) + Layer1 AND-of-core-vetoes (el,prime) + Layer2 "
                       "expression rescue + Layer3 patient-adaptive rails. Removes negatives, hands "
                       "survivors UNCHANGED to genuine PRIME / frozen Epicurus v0.1. Label-blind.",
        "config": frozen_cfg.to_json(),
        "calibration": {"target": 0.95, "conf": CONF, "chosen_t": cal95["chosen_t"],
                        "calibrated_on": DEV_COHORTS, "method": "Clopper-Pearson LB on calibration positives"},
        "invariants": ["missing core axis => KEEP", "strong on any core axis => KEEP",
                       "expression is rescue-only, never a standalone veto, never reweights the ranker",
                       "cohort/study identity is never an input"],
    }, indent=2) + "\n")

    # ---- Evaluate the FROZEN config on each dev cohort (multimer flagged in-sample) ----
    dev_eval = {name: _evaluate_cohort(name, dev[name], frozen_cfg) for name in DEV_COHORTS}

    # ---- LOCKED external test: CheckMate 153, scored ONCE with the frozen config ----
    checkmate = _load_checkmate()
    cm_eval = _evaluate_cohort("checkmate153", checkmate, frozen_cfg)
    cm_eval["expr_coverage"] = round(float(pd.to_numeric(checkmate.get("expr", pd.Series(dtype=float)),
                                                          errors="coerce").notna().mean()), 3)

    report = {
        "experiment": "dynamic_upstream_gate",
        "command": COMMAND,
        "frozen_config": frozen_cfg.to_json(),
        "frozen_config_path": str(FROZEN),
        "pareto_frontier_loco": pareto,
        "dev_cohorts_frozen_config": dev_eval,
        "locked_external_test_checkmate153": cm_eval,
        "verdict": _VERDICT,
        "safety_bar": {
            "1_worst_external_cohort_cp_lb": ">= 0.95 (multimer in-sample excluded)",
            "2_no_downstream_regression": "gated hits@20 delta CI upper bound >= 0 for both rankers",
            "3_removal_material": "negative removal > 0 where the bar allows",
            "4_locked_test_reproduces": "frozen config meets (1)-(2) on CheckMate",
        },
        "caveats": [
            "DEV features are peptide-/presentation-only (prime/el/expr). No raw WES depth, mutant-allele "
            "RNA VAF, purity/CCF, or processing available cross-cohort -> this gate does NOT learn genomics.",
            "multimer is frozen Epicurus' training cohort -> in-sample; excluded from the safety headline.",
            "CheckMate expr coverage is sparse; the gate operates on the core el+prime axes there.",
            "Oracle retention (100%) is a ceiling, never a validation.",
            "Downstream hits@20/recall@20 rise partly because gating shrinks the top-20 denominator; the "
            "no-regression check (gated >= ungated) is the honest safety statement, not the raw lift.",
        ],
    }
    (ART / "dynamic_gate.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.DataFrame(pareto).to_csv(ART / "pareto.csv", index=False)
    _write_per_patient_csv(dev, checkmate, frozen_cfg)
    (ART / "REPORT.md").write_text(_report_md(report))
    print(f"Frozen t={frozen_cfg.t}. Wrote {ART}/dynamic_gate.json, REPORT.md, pareto.csv, per_patient.csv")
    print(f"  frozen config -> {FROZEN}")
    for name in DEV_COHORTS:
        g = dev_eval[name]["gate"]
        print(f"  [{name}{'*IN-SAMPLE' if name in IN_SAMPLE else ''}] retain={g['positive_retention']:.3f} "
              f"cp_lb={g['positive_retention_cp_lb']:.3f} neg_removed={g['negative_removal']:.3f} "
              f"worst_pt={g['positive_retention_per_patient_min']:.3f}")
    g = cm_eval["gate"]
    print(f"  [checkmate153 LOCKED] retain={g['positive_retention']:.3f} cp_lb={g['positive_retention_cp_lb']:.3f} "
          f"neg_removed={g['negative_removal']:.3f} expr_cov={cm_eval['expr_coverage']}")
    return 0


def _write_per_patient_csv(dev, checkmate, cfg):
    rows = []
    for name, frame in [*dev.items(), ("checkmate153", checkmate)]:
        gated = apply_gate(frame, cfg)
        for pid, gp in gated.groupby("patient_id"):
            npos = int((gp["label"] == "POSITIVE").sum())
            kept_pos = int(((gp["label"] == "POSITIVE") & gp["dyn_gate_keep"]).sum())
            nneg = int((gp["label"] == "TESTED_NEGATIVE").sum())
            rem_neg = int(((gp["label"] == "TESTED_NEGATIVE") & ~gp["dyn_gate_keep"]).sum())
            rows.append({"cohort": name, "patient_id": pid, "n_pos": npos, "pos_retained": kept_pos,
                         "retention": kept_pos / npos if npos else np.nan, "n_neg": nneg,
                         "neg_removed": rem_neg, "neg_removal_frac": rem_neg / nneg if nneg else np.nan,
                         "pool_kept": int(gp["dyn_gate_keep"].sum()), "pool_size": len(gp)})
    pd.DataFrame(rows).to_csv(ART / "per_patient.csv", index=False)


def _fmt(x, nd=3):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "—"
    return f"{x:.{nd}f}" if isinstance(x, float) else str(x)


def _report_md(r: dict) -> str:
    L = ["# Dynamic upstream gate — evaluation\n",
         f"`{r['command']}` · frozen config `t={r['frozen_config']['t']}` "
         f"(target 0.95, LOCO-calibrated on {', '.join(r['frozen_config']['calibrated_on'])}).\n",
         "Layered safe-rejection gate (Layer0 deterministic impossibility + Layer1 AND-of-core-vetoes "
         "on el+prime + Layer2 expression rescue + Layer3 patient-adaptive rails). Survivors are handed "
         "UNCHANGED to genuine PRIME / frozen Epicurus v0.1. Label-blind; cohort identity never an input.\n"]

    L.append("\n## Gate metrics — frozen config\n")
    L.append("| cohort | pos retention (CP-LB) | worst-patient | neg removal | pts losing a pos | EL-gate matched retention |")
    L.append("|---|--:|--:|--:|--:|--:|")
    for name, e in r["dev_cohorts_frozen_config"].items():
        g = e["gate"]
        b = e["baseline_el_gate_matched"]
        star = " ⚠️in-sample" if e["in_sample"] else ""
        L.append(f"| {name}{star} | {_fmt(g['positive_retention'])} ({_fmt(g['positive_retention_cp_lb'])}) | "
                 f"{_fmt(g['positive_retention_per_patient_min'])} | {_fmt(g['negative_removal'])} | "
                 f"{g['n_patients_losing_a_positive']} | {_fmt(b['positive_retention'])} |")
    cm = r["locked_external_test_checkmate153"]
    g = cm["gate"]
    b = cm["baseline_el_gate_matched"]
    L.append(f"| **checkmate153 (LOCKED)** | {_fmt(g['positive_retention'])} ({_fmt(g['positive_retention_cp_lb'])}) | "
             f"{_fmt(g['positive_retention_per_patient_min'])} | {_fmt(g['negative_removal'])} | "
             f"{g['n_patients_losing_a_positive']} | {_fmt(b['positive_retention'])} |")
    L.append(f"\n_CheckMate expression coverage {cm['expr_coverage']} → gate acts on core el+prime there._\n")

    L.append("\n## Downstream consequence (gate → UNCHANGED ranker; paired per-patient, gated vs ungated)\n")
    L.append("| cohort | ranker | ungated hits@20 | gated hits@20 | Δ [CI] | no-regression |")
    L.append("|---|---|--:|--:|--:|:--:|")
    for name, e in {**r["dev_cohorts_frozen_config"], "checkmate153": r["locked_external_test_checkmate153"]}.items():
        for a, d in e["downstream_vs_ungated"].items():
            L.append(f"| {name} | {a} | {_fmt(d['ungated_hits@20_mean'])} | {_fmt(d['gated_hits@20_mean'])} | "
                     f"{_fmt(d['hits@20_delta_gated_minus_ungated'])} {d['hits@20_delta_ci']} | "
                     f"{'yes' if d['no_regression'] else 'NO'} |")

    L.append("\n## Pareto frontier — LOCO calibration (negative removal achievable at each retention target)\n")
    L.append("| target | cohort | LOCO t | pos retention (CP-LB) | worst-pt | neg removal | meets CP-LB≥0.95 |")
    L.append("|--:|---|--:|--:|--:|--:|:--:|")
    for row in r["pareto_frontier_loco"]:
        star = "*" if row["in_sample"] else ""
        L.append(f"| {row['target']} | {row['eval_cohort']}{star} | {row['loco_calibrated_t']} | "
                 f"{_fmt(row['positive_retention'])} ({_fmt(row['positive_retention_cp_lb'])}) | "
                 f"{_fmt(row['positive_retention_worst_patient'])} | {_fmt(row['negative_removal'])} | "
                 f"{'yes' if row['meets_bar_cp_lb>=0.95'] else 'no'} |")
    L.append("\n_`*` = multimer, frozen-Epicurus in-sample (optimistic)._\n")

    L.append("\n## Decisive diagnostic — does the gate EVER improve top-20 at a SAFE retention?\n")
    L.append("Threshold sweep (rails off). Watch retention and the paired downstream Δhits@20 together: "
             "where retention stays ≈1.0 the downstream Δ is ≈0 (the gate removes only sub-top-20 "
             "negatives); any positive Δ appears only after retention has already fallen (denominator "
             "effect bought by dropping positives — not a safe win).\n")
    for name, e in {**r["dev_cohorts_frozen_config"], "checkmate153": r["locked_external_test_checkmate153"]}.items():
        L.append(f"\n**{name}** (ungated hits@20: epicurus "
                 f"{_fmt(e['downstream_vs_ungated']['frozen_epicurus']['ungated_hits@20_mean'])}, prime "
                 f"{_fmt(e['downstream_vs_ungated']['genuine_prime']['ungated_hits@20_mean'])})\n")
        L.append("| t | neg removal | pos retention (CP-LB) | Δhits@20 epicurus | Δhits@20 prime |")
        L.append("|--:|--:|--:|--:|--:|")
        for s in e["threshold_sweep"]:
            L.append(f"| {s['t']} | {_fmt(s['negative_removal'])} | {_fmt(s['positive_retention'])} "
                     f"({_fmt(s['positive_retention_cp_lb'])}) | {_fmt(s['d_hits@20_frozen_epicurus'])} "
                     f"{s['d_hits@20_frozen_epicurus_ci']} | {_fmt(s['d_hits@20_genuine_prime'])} "
                     f"{s['d_hits@20_genuine_prime_ci']} |")
        ol = e["oracle_lift_location"]
        L.append(f"\n_Where the oracle lift lives (t={ol['at_t']}): the gate removes "
                 f"**{_fmt(ol['high_presentation_neg_removed'])}** of HIGH-presentation negatives "
                 f"(n={ol['n_high_presentation_neg']}) vs **{_fmt(ol['low_presentation_neg_removed'])}** of "
                 f"LOW-presentation ones (n={ol['n_low_presentation_neg']}). The decoys that outrank "
                 f"positives are the top-EL negatives a label-blind presentation gate cannot touch._\n")

    L.append("\n## Verdict\n" + r["verdict"] + "\n")
    L.append("\n" + "\n".join(f"> {c}" for c in r["caveats"]) + "\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
