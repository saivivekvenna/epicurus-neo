"""Dynamic gate v2 — budgeted reselection: outer leave-one-STUDY-out evaluation (Milestone 7).

Tests whether a negative-risk reselection policy (src/event_b/dynamic_gate_v2.py) lifts NET patient
hits@20 after gate->UNCHANGED frozen Epicurus, and whether any lift TRANSFERS across studies. Mirrors the
outer benchmark of two independent experiments (do not modify those):
  * scripts/hard_decoy_gate_experiment.py (sequence-only; collapsed OOD)
  * scripts/dynamic_utility_gate_experiment.py (HistGBT; Gartner +0.077, no transfer, 25mer confound)

Primary metric: cross-fitted mean patient hits@20 (paired delta vs ungated). The DECISIVE control is
random-matched-pool pruning: v2 must beat removing the SAME number of candidates at random, or it is only
shrinking the pool (denominator), not selecting.

Budget selected by INNER patient-group CV on the training studies (harm-penalized utility); applied ONCE to
the held-out study. Study identity is never a model input. No CheckMate use (consumed). Multimer is
frozen-Epicurus IN-SAMPLE (flagged).

    python -m scripts.dynamic_gate_v2

Writes artifacts/milestone_7_decision/dynamic_gate/{dynamic_gate_v2.json, V2_REPORT.md}.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from event_b.dynamic_gate import GateConfig, apply_gate
from event_b.dynamic_gate_v2 import (
    K,
    counterfactual_reselect,
    fit_negative_risk,
    fit_risk_ensemble,
    patient_hits20,
    reselect,
    utility,
)
from event_b.leakage_registry import canonical_peptide
from event_b.pool_size_sensitivity import patient_eligibility, score_arms

POOL = Path("artifacts/milestone_7_decision/pool_size_sensitivity")
ART = Path("artifacts/milestone_7_decision/dynamic_gate")
COHORTS = ["gartner", "improve", "multimer"]
IN_SAMPLE = {"multimer"}
BUDGETS = [0.05, 0.1, 0.15, 0.2, 0.3, 0.4]
THREAT_K = 2 * K
LAM, MU = 1.0, 0.5  # harm-penalized utility weights for budget selection
INDEP_BENCHMARK = {  # scripts/dynamic_utility_gate_experiment.py (user-provided, for side-by-side)
    "gartner": {"ungated": 0.8077, "gated": 0.8846, "delta": +0.0769, "improved_tied_harmed": "2/24/0"},
    "improve": {"ungated": 1.0714, "gated": 1.0429, "delta": -0.0286, "improved_tied_harmed": "4/60/6"},
    "multimer": {"ungated": 0.9231, "gated": 0.8077, "delta": -0.1154, "improved_tied_harmed": "0/24/2"},
}


def load(name):
    f = pd.read_csv(POOL / f"base_{name}.csv")
    f["patient_id"] = (name + ":" + f["patient_id"].astype(str))
    f = f[f["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    el = patient_eligibility(f)
    f = f[f["patient_id"].isin(el.eligible)].reset_index(drop=True)
    f["cohort"] = name
    return f


def _ungated_hits(frame) -> dict[str, float]:
    return patient_hits20(frame, np.ones(len(frame), bool))


def _delta(frame, keep, ung) -> np.ndarray:
    gat = patient_hits20(frame, keep)
    pids = sorted(ung)
    return np.array([gat[p] - ung[p] for p in pids])


def _inner_cv_select_budget(train: pd.DataFrame) -> float:
    """Pick the budget maximizing harm-penalized utility over inner patient-group folds of the train set."""
    pids = train["patient_id"].unique()
    folds = {p: i % 3 for i, p in enumerate(sorted(pids))}
    best_b, best_u = 0.0, -1e9
    for b in BUDGETS:
        deltas = []
        for fold in range(3):
            tr = train[train["patient_id"].map(folds) != fold]
            va = train[train["patient_id"].map(folds) == fold]
            if va.empty or (tr["label"] == "POSITIVE").sum() == 0:
                continue
            rm = fit_negative_risk(tr)
            ung = _ungated_hits(va)
            keep = reselect(va, rm.risk(va), budget_frac=b, threat_k=THREAT_K)
            deltas.append(_delta(va, keep, ung))
        if not deltas:
            continue
        d = np.concatenate(deltas)
        u = utility(d, lam=LAM, mu=MU)
        if u > best_u:
            best_u, best_b = u, b
    return best_b


def _deleak(train: pd.DataFrame, eval_frame: pd.DataFrame) -> pd.DataFrame:
    ep = {canonical_peptide(p) for p in eval_frame["mutant_peptide"].astype(str)} - {""}
    canon = train["mutant_peptide"].astype(str).map(canonical_peptide)
    return train[~(canon.isin(ep) & (canon != ""))].reset_index(drop=True)


def _random_matched(frame, keep_v2, ung, seed=0) -> np.ndarray:
    """Remove the SAME number per patient as v2, but at random from the threat zone. Decisive control."""
    keep = np.ones(len(frame), bool)
    scored = score_arms(frame)
    for pid, gp in scored.groupby("patient_id"):
        local = frame.index.get_indexer(gp.index)
        n_rm = int(np.sum(~keep_v2[local]))
        if n_rm <= 0:
            continue
        order = np.argsort(-gp["frozen_epicurus"].to_numpy(), kind="mergesort")
        zone = local[order[:THREAT_K]]
        rng = np.random.default_rng([seed, abs(hash(pid)) % (2**31)])
        keep[rng.permutation(zone)[:n_rm]] = False
    return keep


def _summ(delta: np.ndarray) -> dict:
    return {"mean_delta_hits@20": round(float(np.mean(delta)), 4),
            "improved": int(np.sum(delta > 0)), "tied": int(np.sum(delta == 0)), "harmed": int(np.sum(delta < 0)),
            "worst_patient_delta": round(float(np.min(delta)), 4),
            "utility_lam1_mu0.5": round(utility(delta, lam=LAM, mu=MU), 4)}


def run_loco() -> dict:
    frames = {c: load(c) for c in COHORTS}
    out = {}
    for held in COHORTS:
        train = pd.concat([frames[c] for c in COHORTS if c != held], ignore_index=True)
        train = _deleak(train, frames[held])
        budget = _inner_cv_select_budget(train)
        rm = fit_negative_risk(train)
        ev = frames[held]
        ung = _ungated_hits(ev)
        keep_v2 = reselect(ev, rm.risk(ev), budget_frac=budget, threat_k=THREAT_K)
        d_v2 = _delta(ev, keep_v2, ung)
        # counterfactual replacement policy (uncertainty-aware; auto-abstains)
        ens = fit_risk_ensemble(train)
        keep_cf = counterfactual_reselect(ev, ens, max_budget=8, conservative=True)
        d_cf = _delta(ev, keep_cf, ung)
        # controls
        d_rand = np.mean([_delta(ev, _random_matched(ev, keep_v2, ung, seed=s), ung) for s in range(5)], axis=0)
        v1keep = apply_gate(ev, GateConfig()).__getitem__("dyn_gate_keep").to_numpy(bool)
        d_v1 = _delta(ev, v1keep, ung)
        # budget Pareto on the held-out (diagnostic only; budget was chosen by inner CV)
        pareto = []
        for b in BUDGETS:
            kb = reselect(ev, rm.risk(ev), budget_frac=b, threat_k=THREAT_K)
            db = _delta(ev, kb, ung)
            pareto.append({"budget_frac": b, **_summ(db),
                           "pos_retention": round(_pos_ret(ev, kb), 4), "neg_removal": round(_neg_rem(ev, kb), 4)})
        out[held] = {
            "in_sample": held in IN_SAMPLE,
            "n_patients": int(ev["patient_id"].nunique()),
            "selected_budget_frac_inner_cv": budget,
            "ungated_hits@20_mean": round(float(np.mean(list(ung.values()))), 4),
            "v2_fixed_budget": {**_summ(d_v2), "pos_retention": round(_pos_ret(ev, keep_v2), 4),
                                "neg_removal": round(_neg_rem(ev, keep_v2), 4)},
            "v2_counterfactual_abstaining": {**_summ(d_cf),
                                             "candidates_removed": int((~keep_cf).sum()),
                                             "pos_retention": round(_pos_ret(ev, keep_cf), 4)},
            "control_random_matched_pool": _summ(d_rand),
            "baseline_v1_AND_gate": _summ(d_v1),
            "independent_histgbt_benchmark": INDEP_BENCHMARK.get(held),
            "budget_pareto_heldout": pareto,
            "mean_peptide_length": round(float(ev["mutant_peptide"].astype(str).str.len().mean()), 1),
        }
    return out


def _pos_ret(frame, keep) -> float:
    isp = frame["label"].to_numpy() == "POSITIVE"
    return float(np.sum(keep & isp) / isp.sum()) if isp.sum() else float("nan")


def _neg_rem(frame, keep) -> float:
    isn = frame["label"].to_numpy() == "TESTED_NEGATIVE"
    return float(np.sum(~keep & isn) / isn.sum()) if isn.sum() else float("nan")


_VERDICT = (
    "**C (data-limited), no freeze — the learned selection adds nothing over pool reduction, and no policy "
    "transfers.** The DECISIVE control settles it: on Gartner, removing the same number of candidates AT "
    "RANDOM from the top-20 threat zone (random-matched Δ ≈ **+0.06-0.07**) does as well as or BETTER than "
    "v2's learned negative-risk selection (fixed-budget Δ +0.038; counterfactual Δ +0.077). So Gartner's "
    "only positive number is a DENOMINATOR effect (shrinking the pool raises hits@20 mechanically) that "
    "random matches — the negative-risk-at-top-ranks mechanism itself buys nothing. On the well-powered "
    "IMPROVE (deployable 9.4-mer regime) both v2 policies are NEGATIVE (fixed −0.016, counterfactual "
    "−0.115) and on multimer negative (−0.053 / −0.210). The apparent Gartner utility in the independent "
    "HistGBT run (+0.077) is the same denominator+25-mer-regime confound.\n\n"
    "**Counterfactual/abstention finding.** The uncertainty-aware counterfactual policy (remove a top-20 "
    "candidate only when replacement-q LCB > removed-q UCB) was SUPPOSED to auto-abstain when the signal is "
    "weak/OOD, but it did NOT — it removed 77/139/43 candidates and harmed IMPROVE/multimer. Reason: a "
    "bootstrap-logistic ensemble is OVER-CONFIDENT (LCB≈UCB), so the conservative gate collapses to the "
    "mean and never abstains. The design is right but only as safe as its uncertainty calibration: genuine "
    "OOD abstention needs conformal / distributional uncertainty, not bootstrap variance. With correct "
    "abstention the best achievable here is to DO NOTHING (= ungated), because there is no transferable "
    "signal to act on.\n\n"
    "Five independent angles agree (this discordance/expression fixed-budget; the counterfactual; the "
    "independent HistGBT direct-utility; the independent sequence-only that collapsed OOD; the "
    "backfill/diversity probe): no source-invariant negative-risk signal exists among top ranks in the "
    "current cohorts — the recognition wall holds at the top-20 too. NO v2 is frozen. Deliverable = best "
    "Pareto policy behind a feature-availability / peptide-length-regime OOD router with CALIBRATED "
    "abstention (keep-all off-regime; never route by study label), and the exact data limitation: a valid "
    "v2 needs minimal-peptide-regime (8-11mer) cohorts carrying orthogonal WES/RNA features (Miller IPV "
    "PRJNA980652; Gartner reconstruction restricted to class-I minimal epitopes)."
)


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    loco = run_loco()
    report = {"experiment": "dynamic_gate_v2_reselection", "command": "python -m scripts.dynamic_gate_v2",
              "objective": "NET patient hits@20 after gate->unchanged frozen Epicurus (may sacrifice positives)",
              "outer": "leave-one-STUDY-out; budget by inner patient-group CV; study identity never a feature",
              "decisive_control": "random-matched-pool: v2 must beat removing the same count at random",
              "loco": loco, "verdict": _VERDICT,
              "caveats": [
                  "multimer is frozen-Epicurus IN-SAMPLE (optimistic).",
                  "Gartner is mostly 25-mers (regime confound) — not a deployable class-I minimal-epitope set.",
                  "No CheckMate use (consumed locked v1 evidence).",
                  "If v2 ~ random-matched-pool, the 'gain' is pool-size reduction (denominator), not selection.",
              ]}
    (ART / "dynamic_gate_v2.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    (ART / "V2_REPORT.md").write_text(_report_md(report))
    print("Wrote dynamic_gate_v2.json + V2_REPORT.md")
    for held, e in loco.items():
        star = "*" if e["in_sample"] else ""
        print(f"  [{held}{star} len={e['mean_peptide_length']}] ung={e['ungated_hits@20_mean']} "
              f"v2fixedΔ={e['v2_fixed_budget']['mean_delta_hits@20']:+.3f} "
              f"CF_abstainΔ={e['v2_counterfactual_abstaining']['mean_delta_hits@20']:+.3f} "
              f"(CF removed {e['v2_counterfactual_abstaining']['candidates_removed']}) "
              f"randmatchΔ={e['control_random_matched_pool']['mean_delta_hits@20']:+.3f}  "
              f"v1Δ={e['baseline_v1_AND_gate']['mean_delta_hits@20']:+.3f}")
    return 0


def _report_md(r: dict) -> str:
    L = ["# Dynamic gate v2 — budgeted reselection (outer leave-one-study-out)\n",
         f"`{r['command']}` · objective: {r['objective']}.\n",
         "v2 removes the highest negative-risk candidates from each patient's top-20 threat zone so lower "
         "positives can backfill; the frozen ranker is applied UNCHANGED to survivors. Risk model uses "
         "expression + EL/PRIME discordance + interaction (NOT the PRIME-dominated rank directly). Budget "
         "chosen by inner patient-group CV; study identity never an input.\n"]
    L.append("\n## Outer leave-one-study-out\n")
    L.append("`random-matched Δ` removes the SAME count at random from the threat zone — the decisive "
             "control. v2 must beat it or the 'gain' is pool reduction, not selection.\n")
    L.append("| held-out | pep len | ungated | v2 fixed Δ | v2 counterfactual Δ (removed) | **random-matched Δ** | v1 AND Δ | indep HistGBT Δ |")
    L.append("|---|--:|--:|--:|--:|--:|--:|--:|")
    for held, e in r["loco"].items():
        star = " ⚠️IS" if e["in_sample"] else ""
        v2 = e["v2_fixed_budget"]
        cf = e["v2_counterfactual_abstaining"]
        ib = e["independent_histgbt_benchmark"]
        ibd = f"{ib['delta']:+.3f}" if ib else "—"
        L.append(f"| {held}{star} | {e['mean_peptide_length']} | {e['ungated_hits@20_mean']} | "
                 f"{v2['mean_delta_hits@20']:+.3f} ({v2['improved']}/{v2['tied']}/{v2['harmed']}) | "
                 f"{cf['mean_delta_hits@20']:+.3f} ({cf['candidates_removed']}) | "
                 f"**{e['control_random_matched_pool']['mean_delta_hits@20']:+.3f}** | "
                 f"{e['baseline_v1_AND_gate']['mean_delta_hits@20']:+.3f} | {ibd} |")
    L.append("\n_⚠️IS = multimer, frozen-Epicurus in-sample. `random-matched Δ` removes the same count at "
             "random from the threat zone — if v2 ≈ this, the 'gain' is pool reduction, not selection._\n")
    L.append("\n## Budget Pareto (held-out; diagnostic — budget was fixed by inner CV)\n")
    for held, e in r["loco"].items():
        L.append(f"\n**{held}** (ungated {e['ungated_hits@20_mean']}):\n")
        L.append("| budget | Δhits@20 | imp/tie/harm | worst-pt | pos retention | neg removal |")
        L.append("|--:|--:|--:|--:|--:|--:|")
        for p in e["budget_pareto_heldout"]:
            L.append(f"| {p['budget_frac']} | {p['mean_delta_hits@20']:+.3f} | "
                     f"{p['improved']}/{p['tied']}/{p['harmed']} | {p['worst_patient_delta']:+.3f} | "
                     f"{p['pos_retention']} | {p['neg_removal']} |")
    L.append("\n## Verdict\n" + r["verdict"] + "\n")
    L.append("\n" + "\n".join(f"> {c}" for c in r["caveats"]) + "\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
