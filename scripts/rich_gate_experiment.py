"""Rich-feature dynamic gate — nested IMPROVE Partition-held-out evaluation (Milestone 7).

Uses the full 88-column IMPROVE table (orthogonal WES/RNA features) to learn a decision-boundary risk model
and a constrained counterfactual swap policy into UNCHANGED frozen Epicurus. Objective = NET recognized
hits@20. 5 patient-disjoint official Partitions; near-peptide guard (>=0.8) + recurrent-peptide quarantine
between train/eval; within-fold preprocessing. Controls: ungated, matched-random-same-removals, fixed
budget. Ablations: per feature family. DEVELOPMENT evidence only (one study) if positive.

Isolated from the user's scripts/counterfactual_gate_experiment.py (not touched).

    python -m scripts.rich_gate_experiment

Writes artifacts/milestone_7_decision/rich_gate/{rich_gate.json, REPORT.md}.
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
from event_b.leakage_registry import _kmer_index, canonical_peptide, near_duplicate
from event_b.rich_gate import (
    K,
    SAFE_FAMILIES,
    counterfactual_swaps,
    decision_boundary_mask,
    fit_histgbt,
    fit_pairwise,
    load_improve_rich,
    patient_hits20,
    score_arms,
)

ART = Path("artifacts/milestone_7_decision/rich_gate")
MAX_BUDGET = 8
NEAR_THRESH = 0.8
N_RAND = 100       # matched-random draws (stable crc32 seeds)
FIXED_M = 3        # fixed-budget control: always remove the 3 lowest-q of the top-20


def _recurrent_peptides(frame: pd.DataFrame) -> set[str]:
    pp = frame.groupby("mutant_peptide")["patient_id"].nunique()
    return set(pp[pp > 1].index)


def _deleak_train(train: pd.DataFrame, eval_frame: pd.DataFrame, recurrent: set[str]) -> pd.DataFrame:
    """Drop train rows whose peptide is recurrent (cross-patient) OR near-duplicates any eval peptide."""
    eval_canon = {canonical_peptide(p) for p in eval_frame["mutant_peptide"]} - {""}
    eidx = _kmer_index(eval_canon)
    train = train[~train["mutant_peptide"].isin(recurrent)].copy()
    canon = train["mutant_peptide"].map(canonical_peptide)
    near = canon.map(lambda c: bool(c) and near_duplicate(c, eidx, threshold=NEAR_THRESH) is not None)
    return train[~near.to_numpy()].reset_index(drop=True)


def _delta(frame, keep, ung, pids) -> np.ndarray:
    gat = patient_hits20(frame, keep)
    return np.array([gat[p] - ung[p] for p in pids])


def _fold_eval(train_all, ev, families, learner_kind, balance_class=True, margin=0.0):
    ev = ev.reset_index(drop=True)
    recurrent = _recurrent_peptides(pd.concat([train_all, ev], ignore_index=True))
    train = _deleak_train(train_all, ev, recurrent)
    boundary = train[decision_boundary_mask(train)].reset_index(drop=True)
    learner = (fit_histgbt(boundary, families, balance_class=balance_class) if learner_kind == "histgbt"
               else fit_pairwise(boundary, families))
    q = learner.q(ev)
    ung = patient_hits20(ev, np.ones(len(ev), bool))
    pids = sorted(ung)
    keep_cf = counterfactual_swaps(ev, q, max_budget=MAX_BUDGET, margin=margin)
    keep_fx = counterfactual_swaps(ev, q, max_budget=MAX_BUDGET, fixed_m=FIXED_M)
    d_cf = _delta(ev, keep_cf, ung, pids)
    d_fx = _delta(ev, keep_fx, ung, pids)
    # matched-random control: same #removed per patient, sampled from the SAME top-20 removal set,
    # survivors rescored. Stable crc32 seeds; N_RAND draws; localized per patient for speed.
    d_rand = _random_matched_delta(ev, keep_cf, ung, pids)
    return {"pids": pids, "ung": np.array([ung[p] for p in pids]), "d_cf": d_cf, "d_fx": d_fx,
            "d_rand": d_rand, "n_removed": int((~keep_cf).sum())}


def _random_matched_delta(frame, keep_ref, ung, pids):
    """Per-patient matched-random: remove the SAME number the policy removed for that patient, drawn from
    the patient's frozen top-20 (the identical removal-eligible set), rescore survivors. Only patients with
    a removal contribute non-zero; averaged over N_RAND crc32-seeded draws."""
    scored = score_arms(frame)
    delta_by_pid = {p: 0.0 for p in pids}
    for pid, gp in scored.groupby("patient_id"):
        local = frame.index.get_indexer(gp.index)
        n_rm = int(np.sum(~keep_ref[local]))
        if n_rm <= 0:
            continue
        order = np.argsort(-gp["frozen_epicurus"].to_numpy(), kind="mergesort")
        top = local[order[:K]]                      # SAME eligible removal set as the policy (top-20)
        base = ung[str(pid)]
        seed0 = zlib.crc32(str(pid).encode()) & 0x7FFFFFFF
        acc = 0.0
        for s in range(N_RAND):
            rng = np.random.default_rng([s, seed0])
            drop = rng.permutation(top)[:n_rm]
            keep = np.ones(len(frame), bool)
            keep[drop] = False
            sub = frame[keep & frame.index.isin(gp.index)]
            r = score_arms(sub).sort_values("frozen_epicurus", ascending=False, kind="mergesort")
            acc += float((r["label"].to_numpy()[:K] == "POSITIVE").sum()) - base
        delta_by_pid[str(pid)] = acc / N_RAND
    return np.array([delta_by_pid[p] for p in pids])


def _summ(d_cf, d_rand, ung, d_fx=None) -> dict:
    ci = paired_bootstrap(ung + d_cf, ung)  # paired gated vs ungated per patient
    out = {"n_patients": len(d_cf), "ungated_hits@20": round(float(ung.mean()), 4),
           "gated_hits@20": round(float((ung + d_cf).mean()), 4),
           "mean_delta": round(float(d_cf.mean()), 4), "delta_ci": [round(ci.lo, 4), round(ci.hi, 4)],
           "p_better": round(ci.p_better, 3),
           "improved": int((d_cf > 0).sum()), "tied": int((d_cf == 0).sum()), "harmed": int((d_cf < 0).sum()),
           "worst_patient": round(float(d_cf.min()), 4),
           "control_random_matched_delta": round(float(d_rand.mean()), 4),
           "beats_random": bool(d_cf.mean() > d_rand.mean() + 1e-9)}
    if d_fx is not None:
        out["fixed_budget_delta"] = round(float(d_fx.mean()), 4)
    return out


def run(families, learner_kind, frame, balance_class=True) -> dict:
    parts = sorted(frame["partition"].unique())
    d_cf_all, d_fx_all, d_rand_all, ung_all, nrem = [], [], [], [], 0
    per_fold = {}
    for p in parts:
        ev = frame[frame["partition"] == p].reset_index(drop=True)
        train_all = frame[frame["partition"] != p].reset_index(drop=True)
        r = _fold_eval(train_all, ev, families, learner_kind, balance_class)
        d_cf_all.append(r["d_cf"])
        d_fx_all.append(r["d_fx"])
        d_rand_all.append(r["d_rand"])
        ung_all.append(r["ung"])
        nrem += r["n_removed"]
        per_fold[int(p)] = {"n_patients": len(r["pids"]), "mean_delta": round(float(r["d_cf"].mean()), 4),
                            "random_delta": round(float(r["d_rand"].mean()), 4)}
    d_cf, d_fx = np.concatenate(d_cf_all), np.concatenate(d_fx_all)
    d_rand, ung = np.concatenate(d_rand_all), np.concatenate(ung_all)
    return {"learner": learner_kind, "families": families, "balance_class": balance_class,
            "candidates_removed_total": nrem, "overall": _summ(d_cf, d_rand, ung, d_fx), "per_fold": per_fold}


def run_single_feature(frame, col, higher_better=True) -> dict:
    """Predeclared, UNtrained one-feature guard policy (same counterfactual swap rule). If the learned
    model can't beat this, it is overcomplication. Applied on all patients (no fitting -> no folds)."""
    from event_b.dynamic_gate import within_patient_percentile
    q = within_patient_percentile(frame, col, higher_better=higher_better)
    q = np.where(np.isnan(q), 0.5, q)
    ung = patient_hits20(frame, np.ones(len(frame), bool))
    pids = sorted(ung)
    keep = counterfactual_swaps(frame, q, max_budget=MAX_BUDGET, margin=0.0)
    d = _delta(frame, keep, ung, pids)
    d_rand = _random_matched_delta(frame, keep, ung, pids)
    out = _summ(d, d_rand, np.array([ung[p] for p in pids]))
    out["feature"] = col
    out["candidates_removed_total"] = int((~keep).sum())
    return out


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    frame = load_improve_rich()
    print(f"IMPROVE rich: {len(frame)} rows, {frame['patient_id'].nunique()} patients, "
          f"{(frame['label'] == 'POSITIVE').sum()} positives, partitions {sorted(frame['partition'].unique())}")

    results = {"primary_histgbt_all_families": run(SAFE_FAMILIES, "histgbt", frame, balance_class=True),
               "pairwise_all_families": run(SAFE_FAMILIES, "pairwise", frame),
               "histgbt_UNBALANCED_ablation": run(SAFE_FAMILIES, "histgbt", frame, balance_class=False)["overall"],
               "guard_single_feature_PropHydroAro": run_single_feature(frame, "PropHydroAro", higher_better=True)}
    # ablations: each family alone, HistGBT (attribution)
    ablations = {}
    for fam in SAFE_FAMILIES:
        ablations[f"only_{fam}"] = run([fam], "histgbt", frame)["overall"]
    results["ablations_single_family_histgbt"] = ablations

    report = {"experiment": "rich_feature_dynamic_gate_IMPROVE",
              "command": "python -m scripts.rich_gate_experiment",
              "design": "nested 5 patient-disjoint Partition folds; decision-boundary risk learner; "
                        "counterfactual swaps into UNCHANGED frozen Epicurus; NET hits@20.",
              "status_note": "DEVELOPMENT evidence only (single study IMPROVE); NOT external proof.",
              "leakage_controls": ["recurrent-peptide quarantine", f"near-peptide guard >={NEAR_THRESH}",
                                   "within-fold preprocessing", "patient-disjoint partitions",
                                   "excluded label/identity/TME/pipeline-score columns"],
              "results": results,
              "verdict": _verdict(results)}
    (ART / "rich_gate.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    (ART / "REPORT.md").write_text(_md(report))
    print("Wrote rich_gate.json + REPORT.md")
    for name in ["primary_histgbt_all_families", "pairwise_all_families"]:
        o = results[name]["overall"]
        print(f"  [{name}] ung={o['ungated_hits@20']} gated={o['gated_hits@20']} "
              f"Δ={o['mean_delta']:+.3f} CI{o['delta_ci']} imp/tie/harm {o['improved']}/{o['tied']}/{o['harmed']} "
              f"randΔ={o['control_random_matched_delta']:+.3f} fixedΔ={o.get('fixed_budget_delta')} beats_random={o['beats_random']}")
    g = results["guard_single_feature_PropHydroAro"]
    print(f"  [guard PropHydroAro] Δ={g['mean_delta']:+.3f} CI{g['delta_ci']} randΔ={g['control_random_matched_delta']:+.3f} beats_random={g['beats_random']}")
    u = results["histgbt_UNBALANCED_ablation"]
    print(f"  [histgbt UNBALANCED] Δ={u['mean_delta']:+.3f} removed={results['histgbt_UNBALANCED_ablation'].get('n_patients')}")
    return 0


def _verdict(results) -> str:
    o = results["primary_histgbt_all_families"]["overall"]
    beats = o["beats_random"] and o["mean_delta"] > 0
    sig = o["delta_ci"][0] > 0
    if sig and beats:
        return (f"DEVELOPMENT-POSITIVE: rich-feature counterfactual gate lifts NET hits@20 by "
                f"{o['mean_delta']:+.3f} (CI {o['delta_ci']}, beats matched-random {o['control_random_matched_delta']:+.3f}) "
                f"on nested patient-disjoint IMPROVE folds. Single study => development evidence, not external proof.")
    if o["mean_delta"] > 0 and beats:
        return (f"PROMISING but not significant: Δ {o['mean_delta']:+.3f} (CI {o['delta_ci']} spans 0), beats "
                f"matched-random. Needs power / external cohort.")
    return (f"NULL/negative on IMPROVE nested folds: Δ {o['mean_delta']:+.3f} (CI {o['delta_ci']}); "
            f"random-matched {o['control_random_matched_delta']:+.3f}. Rich features do not convert to net hits@20 here.")


def _md(r) -> str:
    L = [f"# Rich-feature dynamic gate — IMPROVE nested Partition eval\n\n`{r['command']}`\n",
         f"_{r['status_note']}_\n", f"Leakage controls: {', '.join(r['leakage_controls'])}.\n"]
    L.append("\n## Primary + pairwise (all safe families)\n")
    L.append("| model | ungated | gated | **Δ hits@20** [CI] | p_better | imp/tie/harm | random-matched Δ | beats random |")
    L.append("|---|--:|--:|--:|--:|--:|--:|:--:|")
    for name in ["primary_histgbt_all_families", "pairwise_all_families"]:
        o = r["results"][name]["overall"]
        L.append(f"| {name.replace('_all_families','')} | {o['ungated_hits@20']} | {o['gated_hits@20']} | "
                 f"**{o['mean_delta']:+.3f}** {o['delta_ci']} | {o['p_better']} | "
                 f"{o['improved']}/{o['tied']}/{o['harmed']} | {o['control_random_matched_delta']:+.3f} | "
                 f"{'yes' if o['beats_random'] else 'NO'} |")
    L.append("\n## Single-family ablations (HistGBT)\n")
    L.append("| family | Δ hits@20 [CI] | imp/tie/harm | random-matched Δ | beats random |")
    L.append("|---|--:|--:|--:|:--:|")
    for fam, o in r["results"]["ablations_single_family_histgbt"].items():
        L.append(f"| {fam} | {o['mean_delta']:+.3f} {o['delta_ci']} | {o['improved']}/{o['tied']}/{o['harmed']} | "
                 f"{o['control_random_matched_delta']:+.3f} | {'yes' if o['beats_random'] else 'NO'} |")
    L.append("\n## Verdict\n" + r["verdict"] + "\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
