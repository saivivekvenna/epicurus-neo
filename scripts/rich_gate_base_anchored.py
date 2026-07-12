"""Base-anchored residual gate — nested IMPROVE Partition eval (Milestone 7, rich-gate v2).

The rich-gate v1 (scripts/rich_gate_experiment.py) FAILED because it used an orthogonal-only q as the
replacement utility, discarding the frozen-Epicurus base rank (NULL/falsified — preserved). v2 fixes this:
the selection utility is BASE-ANCHORED,

    U = base_percentile + alpha * (feature_percentile - 0.5),

so the frozen-Epicurus base ranks the candidates and the orthogonal feature only NUDGES; alpha=0 reproduces
Epicurus exactly (base-only no-op). alpha is selected on TRAINING patients (frozen before the outer
held-out Partition). Epicurus itself is UNCHANGED — it is the base anchor, not replaced.

Predeclared single features only (no multivariate feature-shopping on outer outcomes). Same official 5
patient-disjoint Partitions; matched-random identical-count swap control; swap labels (pos/neg in/out).
DEVELOPMENT evidence only (single study, feature choice followed inspection).

    python -m scripts.rich_gate_base_anchored

Writes artifacts/milestone_7_decision/rich_gate/{base_anchored.json, BASE_ANCHORED_REPORT.md}.
"""

from __future__ import annotations

import json
import zlib
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from benchmark.stats import paired_bootstrap
from event_b.pool_size_sensitivity import patient_eligibility
from event_b.rich_gate import (
    K,
    base_anchored_hits,
    feature_percentile,
    frozen_base_score,
    load_improve_rich,
    peptide_hydro_aro_fraction,
)

POOL = Path("artifacts/milestone_7_decision/pool_size_sensitivity")

ART = Path("artifacts/milestone_7_decision/rich_gate")
THREAT = 60
ALPHA_GRID = [0.0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0]
N_RAND = 100
# Predeclared features (higher_is_better). Primary = the two the user's diagnostic flagged; others = attribution.
FEATURES = [("PropHydroAro", True, "primary"), ("HydroCore", True, "primary"),
            ("HydroAll", True, "attribution"), ("VarAlFreq", True, "attribution"),
            ("SelfSim", False, "attribution"), ("DAI", True, "attribution")]


def _delta(frame, feat_pct, alpha, pids) -> np.ndarray:
    h1 = base_anchored_hits(frame, feat_pct, alpha, THREAT)
    h0 = base_anchored_hits(frame, feat_pct, 0.0, THREAT)
    return np.array([h1[p] - h0[p] for p in pids])


def _select_alpha(train: pd.DataFrame, col: str, hb: bool) -> float:
    """Pick alpha maximizing pooled training-patient net hits@20 delta (fit-free -> frozen for outer)."""
    best_a, best_d = 0.0, 0.0  # default abstain: only adopt alpha if it strictly beats no-op on training
    fp = feature_percentile(train, col, hb)
    pids = sorted(train["patient_id"].unique())
    for a in ALPHA_GRID:
        if a == 0.0:
            continue
        d = float(_delta(train, fp, a, pids).mean())
        if d > best_d:
            best_d, best_a = d, a
    return best_a


def _random_matched(frame, feat_pct, alpha, pids):
    """Matched-random: per patient make the SAME number of swaps the policy made, but swap in RANDOM threat
    candidates (replacing the lowest-base of the Epicurus top-20). >=N_RAND crc32 draws. Fair control."""
    _, swaps = base_anchored_hits(frame, feat_pct, alpha, THREAT, return_swaps=True)
    base = frozen_base_score(frame)
    lab = frame["label"].to_numpy()
    h0 = base_anchored_hits(frame, feat_pct, 0.0, THREAT)
    delta_by_pid = {p: 0.0 for p in pids}
    for pid, gidx in frame.groupby("patient_id").groups.items():
        k = swaps[str(pid)]["n_swaps"]
        if k <= 0:
            continue
        rows = frame.index.get_indexer(gidx)
        bp = pd.Series(base[rows]).rank(pct=True).to_numpy()
        thr = np.argsort(-bp, kind="mergesort")[:THREAT]
        base_top_local = thr[np.argsort(-bp[thr], kind="mergesort")[:K]]
        pool_in = np.array([i for i in thr if i not in set(base_top_local)])  # candidates outside top-20
        drop_lowest = base_top_local[np.argsort(bp[base_top_local], kind="mergesort")[:k]]  # lowest-base of top20
        keep_top = [i for i in base_top_local if i not in set(drop_lowest)]
        seed0 = zlib.crc32(str(pid).encode()) & 0x7FFFFFFF
        acc = 0.0
        for s in range(N_RAND):
            rng = np.random.default_rng([s, seed0])
            add = rng.permutation(pool_in)[:k] if len(pool_in) >= k else pool_in
            newtop = list(keep_top) + list(add)
            acc += float(sum(lab[rows[i]] == "POSITIVE" for i in newtop)) - h0[str(pid)]
        delta_by_pid[str(pid)] = acc / N_RAND
    return np.array([delta_by_pid[p] for p in pids])


def run_feature(frame, col, hb) -> dict:
    parts = sorted(frame["partition"].unique())
    d_all, d_rand_all, ung_all = [], [], []
    per_fold, alphas = {}, []
    swaps_pos_in = swaps_neg_in = swaps_pos_out = 0
    for p in parts:
        ev = frame[frame["partition"] == p].reset_index(drop=True)
        train = frame[frame["partition"] != p].reset_index(drop=True)
        alpha = _select_alpha(train, col, hb)
        alphas.append(alpha)
        fp = feature_percentile(ev, col, hb)
        pids = sorted(ev["patient_id"].unique())
        h0 = base_anchored_hits(ev, fp, 0.0, THREAT)
        d = _delta(ev, fp, alpha, pids)
        d_rand = _random_matched(ev, fp, alpha, pids)
        _, swaps = base_anchored_hits(ev, fp, alpha, THREAT, return_swaps=True)
        swaps_pos_in += sum(s["pos_in"] for s in swaps.values())
        swaps_neg_in += sum(s["neg_in"] for s in swaps.values())
        swaps_pos_out += sum(s["pos_out"] for s in swaps.values())
        d_all.append(d)
        d_rand_all.append(d_rand)
        ung_all.append(np.array([h0[p_] for p_ in pids]))
        per_fold[int(p)] = {"alpha": alpha, "mean_delta": round(float(d.mean()), 4),
                            "random_delta": round(float(d_rand.mean()), 4), "n_patients": len(pids)}
    d = np.concatenate(d_all)
    d_rand = np.concatenate(d_rand_all)
    ung = np.concatenate(ung_all)
    ci = paired_bootstrap(ung + d, ung)
    return {"feature": col, "higher_better": hb, "chosen_alphas": alphas,
            "ungated_hits@20": round(float(ung.mean()), 4), "gated_hits@20": round(float((ung + d).mean()), 4),
            "mean_delta": round(float(d.mean()), 4), "delta_ci": [round(ci.lo, 4), round(ci.hi, 4)],
            "p_better": round(ci.p_better, 3),
            "improved": int((d > 0).sum()), "tied": int((d == 0).sum()), "harmed": int((d < 0).sum()),
            "worst_patient": round(float(d.min()), 4),
            "control_random_matched_delta": round(float(d_rand.mean()), 4),
            "beats_random": bool(d.mean() > d_rand.mean() + 1e-9),
            "swap_labels": {"pos_swapped_in": swaps_pos_in, "neg_swapped_in": swaps_neg_in,
                            "pos_swapped_out": swaps_pos_out},
            "per_fold": per_fold}


def run_leave_source_out(frame, col, hb) -> dict:
    """Transport within IMPROVE: alpha selected on the OTHER tissue sources, applied to the held-out
    source. Tests whether the hydrophobic gain is source-invariant or driven by one tissue."""
    sources = sorted(frame["source"].unique())
    out = {}
    d_all = []
    for s in sources:
        ev = frame[frame["source"] == s].reset_index(drop=True)
        train = frame[frame["source"] != s].reset_index(drop=True)
        alpha = _select_alpha(train, col, hb)
        fp = feature_percentile(ev, col, hb)
        pids = sorted(ev["patient_id"].unique())
        d = _delta(ev, fp, alpha, pids)
        out[s] = {"alpha": alpha, "mean_delta": round(float(d.mean()), 4), "n_patients": len(pids)}
        d_all.append(d)
    d = np.concatenate(d_all)
    return {"per_source": out, "pooled_mean_delta": round(float(d.mean()), 4),
            "all_sources_positive": bool(all(v["mean_delta"] > 0 for v in out.values()))}


def _load_external(name):
    f = pd.read_csv(POOL / f"base_{name}.csv")
    f["patient_id"] = f["patient_id"].astype(str)
    f = f[f["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    el = patient_eligibility(f)
    return f[f["patient_id"].isin(el.eligible)].reset_index(drop=True)


def external_transfer(frame_improve, col="PropHydroAro") -> dict:
    """Freeze alpha on ALL IMPROVE, apply the base-anchored hydrophobicity gate label-blind to external
    pools (Gartner, multimer) using a peptide-derived hydrophobicity proxy. Tests universal transfer."""
    alpha = _select_alpha(frame_improve, col, True)
    out = {"frozen_alpha_on_improve": alpha, "feature": col,
           "note": "external feature = peptide-derived hydro/aromatic fraction (proxy for PropHydroAro)"}
    for name in ["gartner", "multimer"]:
        ext = _load_external(name)
        fp_raw = peptide_hydro_aro_fraction(ext["mutant_peptide"])
        ext = ext.assign(_hydro=fp_raw)
        fp = feature_percentile(ext, "_hydro", True)
        pids = sorted(ext["patient_id"].unique())
        d = _delta(ext, fp, alpha, pids)
        out[name] = {"mean_delta": round(float(d.mean()), 4), "n_patients": len(pids),
                     "improved": int((d > 0).sum()), "tied": int((d == 0).sum()), "harmed": int((d < 0).sum()),
                     "in_sample": name == "multimer"}
    return out


def why_improve_differs(frame) -> dict:
    """Quick diagnostic: is the hydrophobic signal a length/source artifact of IMPROVE?"""
    fp = peptide_hydro_aro_fraction(frame["mutant_peptide"])
    lab = frame["label"].to_numpy()
    rows = []
    for src, gidx in frame.groupby("source").groups.items():
        idx = frame.index.get_indexer(gidx)
        pos = lab[idx] == "POSITIVE"
        rows.append({"source": src, "n": len(idx), "pos": int(pos.sum()),
                     "mean_pep_len": round(float(frame["mutant_peptide"].iloc[idx].str.len().mean()), 2),
                     "hydro_pos": round(float(np.nanmean(fp[idx][pos])), 3),
                     "hydro_neg": round(float(np.nanmean(fp[idx][~pos])), 3)})
    return {"by_source": rows,
            "interpretation": "If positives are systematically MORE hydrophobic than negatives only in "
                              "IMPROVE (a TIL/screened set), the gain reflects IMPROVE's candidate-selection/"
                              "assay regime, not a universal recognition rule."}


def base_only_noop_sanity(frame) -> dict:
    """alpha=0 MUST reproduce Epicurus exactly (delta 0 everywhere)."""
    fp = feature_percentile(frame, "PropHydroAro", True)
    pids = sorted(frame["patient_id"].unique())
    d = _delta(frame, fp, 0.0, pids)
    return {"max_abs_delta_at_alpha0": float(np.abs(d).max()), "passes": bool(np.abs(d).max() == 0.0)}


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    frame = load_improve_rich()
    sanity = base_only_noop_sanity(frame)
    results = {col: run_feature(frame, col, hb) for col, hb, _ in FEATURES}
    lso = run_leave_source_out(frame, "PropHydroAro", True)
    ext = external_transfer(frame, "PropHydroAro")
    why = why_improve_differs(frame)
    transport = _transport_verdict(results, lso, ext)
    report = {"experiment": "rich_gate_v2_base_anchored_IMPROVE",
              "command": "python -m scripts.rich_gate_base_anchored",
              "utility": "U = base_percentile + alpha*(feature_percentile-0.5); Epicurus UNCHANGED (base anchor)",
              "design": "nested 5 patient-disjoint Partitions; alpha selected on training patients, frozen for outer",
              "status_note": "DEVELOPMENT/LOCAL only. rich-gate v1 (orthogonal-only utility) remains NULL/falsified. "
                             "The IMPROVE hydrophobic gain does NOT transfer externally (see external_transfer) "
                             "=> regime-aware abstention required; DO NOT freeze hydrophobicity into product.",
              "base_only_noop_sanity": sanity,
              "feature_classes": {col: cls for col, _, cls in FEATURES},
              "results": results,
              "leave_source_out_within_improve": lso,
              "external_transfer_falsification": ext,
              "why_improve_differs": why,
              "regime_aware_transport_verdict": transport,
              "verdict": _verdict(results, sanity)}
    (ART / "base_anchored.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    (ART / "BASE_ANCHORED_REPORT.md").write_text(_md(report))
    print(f"base-only no-op sanity: {sanity['passes']}")
    for col, _, cls in FEATURES:
        o = results[col]
        print(f"  [{col} {cls}] Δ={o['mean_delta']:+.3f} CI{o['delta_ci']} imp/tie/harm "
              f"{o['improved']}/{o['tied']}/{o['harmed']} randΔ={o['control_random_matched_delta']:+.3f} "
              f"beats_random={o['beats_random']}")
    print(f"  leave-source-out (IMPROVE bladder/melanoma/Basket): pooled Δ={lso['pooled_mean_delta']:+.3f} "
          f"all_sources_positive={lso['all_sources_positive']} {[(s, v['mean_delta']) for s, v in lso['per_source'].items()]}")
    print(f"  EXTERNAL transfer (frozen α={ext['frozen_alpha_on_improve']}): "
          f"gartner Δ={ext['gartner']['mean_delta']:+.3f}  multimer Δ={ext['multimer']['mean_delta']:+.3f}")
    print(f"  TRANSPORT VERDICT: {transport['decision']}")
    return 0


def _transport_verdict(results, lso, ext) -> dict:
    """Regime-aware abstention: the policy may activate ONLY where transport evidence supports it."""
    local_ok = results["PropHydroAro"]["mean_delta"] > 0 and results["PropHydroAro"]["delta_ci"][0] > 0
    source_ok = lso["all_sources_positive"]
    external_ok = ext["gartner"]["mean_delta"] >= 0 and ext["multimer"]["mean_delta"] >= 0
    if local_ok and source_ok and external_ok:
        decision = "ACTIVATE (transport supported everywhere)"
    elif local_ok and source_ok:
        decision = ("REGIME-LOCAL: significant within IMPROVE and across its tissue sources, but FALSIFIED "
                    "externally (Gartner/multimer harmed) => a regime-aware gate must ABSTAIN (no-op) off the "
                    "IMPROVE-like regime. NOT frozen into product.")
    elif local_ok:
        decision = "IMPROVE-ONLY: not even source-invariant within IMPROVE => abstain off-regime."
    else:
        decision = "NULL"
    return {"local_significant": local_ok, "source_invariant_within_improve": source_ok,
            "external_nonharmful": external_ok, "decision": decision,
            "rule": "activate only if local_significant AND source_invariant AND external_nonharmful; else abstain"}


def _verdict(results, sanity) -> str:
    if not sanity["passes"]:
        return "INVALID: base-only no-op sanity failed (alpha=0 should reproduce Epicurus)."
    prim = [results[c] for c in ("PropHydroAro", "HydroCore")]
    best = max(prim, key=lambda o: o["mean_delta"])
    robust = best["mean_delta"] > 0 and best["beats_random"] and best["delta_ci"][0] > 0
    if robust:
        return (f"DEVELOPMENT-POSITIVE (single study): base-anchored gate on {best['feature']} lifts NET "
                f"hits@20 by {best['mean_delta']:+.3f} (CI {best['delta_ci']}, beats matched-random "
                f"{best['control_random_matched_delta']:+.3f}) on nested patient-disjoint IMPROVE folds. Feature "
                f"choice followed inspection => development discovery, NOT external proof. Requires an UNTOUCHED "
                f"external rich cohort to confirm before any freeze.")
    if best["mean_delta"] > 0 and best["beats_random"]:
        return (f"PROMISING but CI spans 0: {best['feature']} Δ {best['mean_delta']:+.3f} (CI {best['delta_ci']}), "
                f"beats matched-random. Under-powered on one study; needs external rich cohort.")
    return (f"NULL: base-anchored single-feature gate does not robustly beat matched-random on IMPROVE "
            f"({best['feature']} Δ {best['mean_delta']:+.3f}).")


def _md(r) -> str:
    L = [f"# Rich-gate v2 — base-anchored residual gate (IMPROVE nested)\n\n`{r['command']}`\n",
         f"Utility: `{r['utility']}`.\n_{r['status_note']}_\n",
         f"Base-only no-op sanity (alpha=0 == Epicurus): **{'PASS' if r['base_only_noop_sanity']['passes'] else 'FAIL'}** "
         f"(max|Δ|={r['base_only_noop_sanity']['max_abs_delta_at_alpha0']}).\n"]
    L.append("\n| feature | class | **Δ hits@20** [CI] | p_better | imp/tie/harm | random-matched Δ | beats random | pos/neg swapped-in | chosen alphas |")
    L.append("|---|---|--:|--:|--:|--:|:--:|--:|--:|")
    for col, _, cls in FEATURES:
        o = r["results"][col]
        s = o["swap_labels"]
        L.append(f"| {col} | {cls} | **{o['mean_delta']:+.3f}** {o['delta_ci']} | {o['p_better']} | "
                 f"{o['improved']}/{o['tied']}/{o['harmed']} | {o['control_random_matched_delta']:+.3f} | "
                 f"{'yes' if o['beats_random'] else 'NO'} | {s['pos_swapped_in']}/{s['neg_swapped_in']} | "
                 f"{o['chosen_alphas']} |")
    L.append("\n## Verdict (local IMPROVE)\n" + r["verdict"] + "\n")

    lso = r["leave_source_out_within_improve"]
    L.append("\n## Transport 1 — leave-source-out WITHIN IMPROVE (PropHydroAro)\n")
    L.append(f"Pooled Δ **{lso['pooled_mean_delta']:+.3f}**; all sources positive: "
             f"**{lso['all_sources_positive']}**. Per source: " +
             ", ".join(f"{s} {v['mean_delta']:+.3f} (α={v['alpha']})" for s, v in lso["per_source"].items()) + ".\n")

    ext = r["external_transfer_falsification"]
    L.append("\n## Transport 2 — EXTERNAL transfer (frozen on IMPROVE → Gartner/multimer)\n")
    L.append(f"Frozen α={ext['frozen_alpha_on_improve']} on IMPROVE, applied label-blind with a "
             f"peptide-derived hydro/aromatic proxy:\n")
    L.append("| external cohort | Δ hits@20 | imp/tie/harm |")
    L.append("|---|--:|--:|")
    for name in ["gartner", "multimer"]:
        e = ext[name]
        star = " ⚠️in-sample" if e["in_sample"] else ""
        L.append(f"| {name}{star} | **{e['mean_delta']:+.3f}** | {e['improved']}/{e['tied']}/{e['harmed']} |")
    L.append("\n_External transfer FALSIFIES a universal hard gate: the IMPROVE hydrophobic gain does not "
             "transport (harms Gartner/multimer). Consistent with the user's independent external run "
             "(Gartner −0.154, multimer −0.577).\n")

    why = r["why_improve_differs"]
    L.append("\n## Why IMPROVE differs (hydrophobicity of positives vs negatives by source)\n")
    L.append("| source | n | pos | mean pep len | hydro(pos) | hydro(neg) |")
    L.append("|---|--:|--:|--:|--:|--:|")
    for row in why["by_source"]:
        L.append(f"| {row['source']} | {row['n']} | {row['pos']} | {row['mean_pep_len']} | "
                 f"{row['hydro_pos']} | {row['hydro_neg']} |")
    L.append(f"\n_{why['interpretation']}_\n")

    tv = r["regime_aware_transport_verdict"]
    L.append("\n## Regime-aware transport verdict\n")
    L.append(f"- local significant (IMPROVE nested): **{tv['local_significant']}**\n"
             f"- source-invariant within IMPROVE: **{tv['source_invariant_within_improve']}**\n"
             f"- external non-harmful: **{tv['external_nonharmful']}**\n\n"
             f"**Decision:** {tv['decision']}\n\n_Rule: {tv['rule']}._\n")
    L.append("\n> Epicurus is UNCHANGED (the base anchor); the feature only nudges. alpha=0 reproduces "
             "Epicurus (no-op sanity). rich-gate v1 (orthogonal-only utility) stays NULL/falsified. The "
             "base-anchored architecture is validated LOCALLY on IMPROVE but hydrophobicity is regime-specific "
             "and MUST NOT be frozen into product; a regime-aware gate abstains off-regime. Next lever = an "
             "untouched external RICH cohort + leave-source-out-supported features.\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
