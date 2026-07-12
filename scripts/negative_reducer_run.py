"""Risk-controlled negative reducer — NESTED LEAVE-ONE-STUDY-OUT runner (CONTRACT.md + CORRECTIONS 1-2).

Outer = leave-one-STUDY-out over {IMPROVE, Gartner, multimer}. Inner = patient-grouped CV within outer-train
selecting (model in {NULL, nonnegative-logistic C-grid}, protected core m in {0,5,10}) by inner-OOF
patient-macro Delta-hits@20; tau is the most aggressive cut whose pooled OOF CP-95%-LB retention >= 0.95.
Monotonic HGB is a DIAGNOSTIC only (Correction 2.3): reported, never selected/frozen. tau is always
calibrated on OUT-OF-FOLD scores (Correction 2.2). The selected recipe is refit on full outer-train to score
the untouched held-out study; retention is MEASURED out-of-study (CP claim only where powered). Aggregate +
IMPROVE CP retention both enforced (2.4). The deployment recipe is chosen by the SAME inner_select run once
on ALL DEV (2.5). Freezes a fitted apply-only payload iff §5 eligibility passes, else NULL. Reads NO Sid /
Miller file.

    PYTHONPATH=src python -m scripts.negative_reducer_run          # full run + freeze + repro check
    PYTHONPATH=src python -m scripts.negative_reducer_run --sha     # print frozen SHA only (repro subprocess)
"""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
import scipy
import sklearn
from sklearn.model_selection import GroupKFold

from event_b.leakage_registry import canonical_peptide
from event_b import negative_reducer as nr

ART = Path("artifacts/milestone_7_decision/negative_reducer")
FROZEN = Path("configs/frozen/negative_reducer_v1.json")
MODEL_RANK = {"NULL": 0, "nnlog": 1, "hgb": 2}


def _is_pos(df):
    return (df["label"].to_numpy() == "POSITIVE").astype(int)


# ---- uniform model interface (nnlog serializable / frozen-eligible; hgb diagnostic-only) --------------
def fit_model(df, cols, model, C):
    X, y, w = nr.feat_matrix(df, cols), _is_pos(df), nr.balanced_weights(df)
    if model == "nnlog":
        b0, coef = nr.fit_nnlogistic(X, y, w, C)
        return {"kind": "nnlog", "intercept": b0, "coef": coef}
    return {"kind": "hgb", "clf": nr.fit_hgb(X, y, w, len(cols))}


def score_model(fm, df, cols):
    X = nr.feat_matrix(df, cols)
    if fm["kind"] == "nnlog":
        return nr.nnlogistic_score(X, fm["intercept"], fm["coef"])
    return fm["clf"].predict_proba(X)[:, 1]


def make_oof(train: pd.DataFrame, cols, model: str, C):
    """Inner patient-grouped OOF keep-scores + per-row leakage-clean mask. Fold models use FOLD-LOCAL
    balanced weights (Correction 2.1 — no validation-fold class-total leakage)."""
    X = nr.feat_matrix(train, cols)
    y = _is_pos(train)
    groups = train["patient_id"].to_numpy()
    n_groups = len(np.unique(groups))
    oof = np.full(len(train), np.nan)
    clean = np.ones(len(train), bool)
    gkf = GroupKFold(n_splits=min(4, n_groups))
    for tr, va in gkf.split(X, y, groups):
        w_tr = nr.balanced_weights(train.iloc[tr])                    # fold-local weights
        if model == "nnlog":
            b0, coef = nr.fit_nnlogistic(X[tr], y[tr], w_tr, C)
            oof[va] = nr.nnlogistic_score(X[va], b0, coef)
        else:
            clf = nr.fit_hgb(X[tr], y[tr], w_tr, len(cols))
            oof[va] = clf.predict_proba(X[va])[:, 1]
        trainpep = set(train.iloc[tr]["mut_peptide"].map(canonical_peptide))
        clean[va] = nr._clean_against(train.iloc[va].reset_index(drop=True), trainpep)
    return oof, clean


def objective_for_m(train, cols, oof, oof_clean, m):
    """Inner-OOF patient-macro Delta-hits@20 for protected core m, with OOF CP-calibrated tau."""
    y = _is_pos(train)
    removable = nr.removable_mask(train, cols, m, ood=set())
    tau, r_max, n_pos, cp_lb, powered = nr.calibrate_tau(oof, y, removable)
    removed = removable & (oof < tau)
    null = nr.hits_at_k(train, np.zeros(len(train), bool), oof_clean)
    gate = nr.hits_at_k(train, removed, oof_clean)
    delta = float(np.mean([gate[p] - null[p] for p in null]))
    return {"m": m, "tau": (None if not np.isfinite(tau) else float(tau)), "tau_inf": bool(not np.isfinite(tau)),
            "delta": round(delta, 5), "neg_removed": int((removed & (y == 0)).sum()),
            "cp_lb": round(cp_lb, 4), "powered": bool(powered), "removed_total": int(removed.sum())}


def _tie(cand):
    """max-key: delta desc, then NULL > larger m > simpler model > less aggressive (fewer removed)."""
    return (cand["delta"], 1 if cand["model"] == "NULL" else 0, cand.get("m", 0),
            -MODEL_RANK[cand["model"]], -cand.get("removed_total", 0))


def inner_select(train: pd.DataFrame, cols):
    """Nested selection over the SELECTABLE/freezable space {NULL, nnlog x C} x m (HGB excluded, 2.3).
    Returns (chosen, selectable_table, hgb_diagnostic_table)."""
    cands = [{"model": "NULL", "C": None, "m": 0, "delta": 0.0, "tau": None, "tau_inf": False,
              "neg_removed": 0, "removed_total": 0, "cp_lb": 1.0, "powered": True}]
    for C in nr.C_GRID:
        oof, oof_clean = make_oof(train, cols, "nnlog", C)
        for m in nr.M_GRID:
            cands.append({"model": "nnlog", "C": C, **objective_for_m(train, cols, oof, oof_clean, m)})
    hgb_diag = []
    if nr.hgb_available():
        oof, oof_clean = make_oof(train, cols, "hgb", None)
        for m in nr.M_GRID:
            hgb_diag.append({"model": "hgb", "C": None, **objective_for_m(train, cols, oof, oof_clean, m)})
    chosen = max(cands, key=_tie)
    return chosen, cands, hgb_diag


def paired_bootstrap(gate: dict, null: dict, n=5000, seed=12345):
    pids = sorted(null)
    d = np.array([gate[p] - null[p] for p in pids], float)
    rng = np.random.default_rng(seed)
    means = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(n)]) if len(d) else np.zeros(1)
    return {"mean_delta": round(float(d.mean()) if len(d) else 0.0, 4),
            "ci95": [round(float(np.percentile(means, 2.5)), 4), round(float(np.percentile(means, 97.5)), 4)],
            "p_gt_0": round(float((means > 0).mean()), 4), "n_patients": len(pids)}


def apply_to_test(train, test, cols, choice):
    """Correction 2.2: tau from patient-grouped OOF on outer-train for the selected (model,C,m); full model
    refit on outer-train only to SCORE the untouched held-out study. Returns metrics + the removal mask
    (removal depends only on test FEATURES + train-derived tau, never on test labels)."""
    model, m = choice["model"], choice.get("m", 0)
    if model == "NULL":
        removed = np.zeros(len(test), bool)
        tau = np.inf
        ood = set()
    else:
        oof, _ = make_oof(train, cols, model, choice["C"])
        rem_train = nr.removable_mask(train, cols, m, ood=set())
        tau, *_ = nr.calibrate_tau(oof, _is_pos(train), rem_train)     # OOF-calibrated, not in-sample
        fm = fit_model(train, cols, model, choice["C"])
        ood = nr.ood_patients(train, test, cols)
        s_test = score_model(fm, test, cols)
        removed = nr.gate_removed(test, s_test, tau, cols, m, ood)

    trainpep = set(train["mut_peptide"].map(canonical_peptide))
    clean = nr._clean_against(test, trainpep)
    null = nr.hits_at_k(test, np.zeros(len(test), bool), clean)
    gate = nr.hits_at_k(test, removed, clean)

    y = _is_pos(test)
    n_pos, n_neg = int(y.sum()), int((y == 0).sum())
    pos_removed, neg_removed = int((removed & (y == 1)).sum()), int((removed & (y == 0)).sum())
    retention = (n_pos - pos_removed) / n_pos if n_pos else 1.0
    cp_lb = nr.cp_lower(n_pos - pos_removed, n_pos)
    powered = n_pos >= nr.CP_MIN_POS

    counts, idx = {}, test.index.get_indexer
    for pid, g in test.groupby("patient_id"):
        counts[str(pid)] = int(removed[idx(g.index)].sum())
    rand_delta = (nr.matched_random_delta(test, counts, cols, m, ood, clean, null)
                  if (neg_removed + pos_removed) > 0 else 0.0)

    prime = pd.to_numeric(test["prime"], errors="coerce").to_numpy()
    rank_improved = []
    for pid, g in test.groupby("patient_id"):
        loc = idx(g.index)
        order = np.argsort(np.where(np.isfinite(prime[loc]), prime[loc], np.inf), kind="mergesort")
        base_rank = {int(order[r]): r for r in range(len(order))}
        surv = [i for i in order if not removed[loc][i]]
        surv_rank = {int(surv[r]): r for r in range(len(surv))}
        for i in np.where(y[loc] == 1)[0]:
            if i in surv_rank:
                rank_improved.append(base_rank[i] - surv_rank[i])
    macro_delta = float(np.mean([gate[p] - null[p] for p in null]))

    return {
        "n_pos": n_pos, "n_neg": n_neg, "pos_removed": pos_removed, "neg_removed": neg_removed,
        "raw_retention": round(retention, 4), "cp_lb_retention": round(cp_lb, 4), "cp_powered": bool(powered),
        "cp_claim": ("retention CP-95%%-LB=%.4f" % cp_lb) if powered
                    else "CP-underpowered: claim abstained (gate still applied)",
        "neg_removal_frac": round(neg_removed / n_neg, 4) if n_neg else 0.0,
        "macro_delta_hits20": round(macro_delta, 4), "matched_random_delta": round(float(rand_delta), 4),
        "beats_random": bool(macro_delta > rand_delta + 1e-9),
        "pooled_hits_gate": int(sum(gate.values())), "pooled_hits_null": int(sum(null.values())),
        "mean_pos_rank_improvement": round(float(np.mean(rank_improved)) if rank_improved else 0.0, 3),
        "paired_bootstrap": paired_bootstrap(gate, null),
        "tau": (None if not np.isfinite(tau) else float(tau)),
        "removed_idx": [int(i) for i in np.where(removed)[0]],       # audit / blindness tests
        "_gate_hits": gate, "_null_hits": null,
    }


def serialize_recipe(df, cols, choice):
    """Full-DEV freeze payload for a SELECTABLE recipe (NULL or nnlog): OOF-calibrate tau, refit on all DEV."""
    model, m = choice["model"], choice.get("m", 0)
    if model == "NULL":
        return {"model": "NULL", "m": 0, "tau": None, "note": "no-gate"}
    oof, _ = make_oof(df, cols, model, choice["C"])
    tau, *_ = nr.calibrate_tau(oof, _is_pos(df), nr.removable_mask(df, cols, m, ood=set()))
    fm = fit_model(df, cols, model, choice["C"])
    payload = {"feature_order": list(cols), "model": model, "C": choice["C"], "m": m,
               "tau": (None if not np.isfinite(tau) else float(tau)),
               "coef": [float(x) for x in fm["coef"]], "intercept": float(fm["intercept"]),
               "percentile_policy": "within-patient rank(pct=True), oriented; NaN->0.5",
               "direction_policy": {c: ("higher_raw_better" if nr.HIGHER_BETTER[c] else "lower_raw_better") for c in cols},
               "protected_core": "top-m genuine-PRIME per patient unremovable; boundary ties all protected",
               "removal_rule": "remove non-core, feature-present, in-support candidates with keepscore<tau (tau=None => remove none)",
               "ood_policy": "raw-feature [p1,p99] envelope; >50% out-of-support candidates => KEEP patient"}
    payload["model_payload_sha256"] = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return payload


def _file_sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest() if Path(p).exists() else "MISSING"


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def run_selection(df, cols):
    """Full nested LOSO + eligibility + frozen payload. Pure over df (no file writes). Deterministic."""
    folds, agg_gate, agg_null = {}, {}, {}
    agg_pos = agg_pos_removed = 0
    for held in nr.STUDIES:
        train = df[df["study"] != held].reset_index(drop=True)
        test = df[df["study"] == held].reset_index(drop=True)
        chosen, cands, hgb_diag = inner_select(train, cols)
        res = apply_to_test(train, test, cols, chosen)
        agg_gate.update(res.pop("_gate_hits"))
        agg_null.update(res.pop("_null_hits"))
        agg_pos += res["n_pos"]
        agg_pos_removed += res["pos_removed"]
        folds[held] = {"selected": {k: chosen[k] for k in ("model", "C", "m", "delta", "tau")},
                       "inner_candidates": cands, "hgb_diagnostic": hgb_diag, "outer_test": res}

    agg_delta = float(np.mean([agg_gate[p] - agg_null[p] for p in agg_null]))
    agg_boot = paired_bootstrap(agg_gate, agg_null)
    per_study_rand = {h: folds[h]["outer_test"]["matched_random_delta"] for h in nr.STUDIES}
    npat = {h: sum(1 for p in agg_null if p.startswith(h + ":")) for h in nr.STUDIES}
    agg_rand = float(sum(per_study_rand[h] * npat[h] for h in nr.STUDIES) / sum(npat.values()))
    worst_delta = min(folds[h]["outer_test"]["macro_delta_hits20"] for h in nr.STUDIES)

    agg_cp = nr.cp_lower(agg_pos - agg_pos_removed, agg_pos)
    improve_cp = folds["improve"]["outer_test"]["cp_lb_retention"]
    every_noncat = all(folds[h]["outer_test"]["macro_delta_hits20"] >= nr.CATASTROPHIC for h in nr.STUDIES)
    agg_gain = agg_delta > 0 and agg_delta > agg_rand + 1e-9
    any_removal = any(folds[h]["outer_test"]["neg_removal_frac"] > 0 for h in nr.STUDIES)
    cp_ok = agg_cp >= 0.95 and improve_cp >= 0.95                      # Correction 2.4: aggregate AND IMPROVE
    eligible = bool(every_noncat and agg_gain and any_removal and cp_ok)

    # Correction 2.5: deployment recipe = same inner_select run once on ALL DEV
    deploy_choice, deploy_cands, deploy_hgb = inner_select(df, cols)
    if eligible:
        payload = serialize_recipe(df, cols, deploy_choice)
        frozen = {"name": "negative_reducer", "version": "1.0.0", "frozen": payload,
                  "full_dev_selection": {k: deploy_choice[k] for k in ("model", "C", "m", "delta")},
                  "stage2_must_not_refit": "Stage 2 applies this payload only (score->threshold->remove); never refit."}
    else:
        frozen = {"name": "negative_reducer", "version": "1.0.0", "frozen": "NULL",
                  "reason": "nested LOSO evidence did not pass CONTRACT §5 (see eligibility)",
                  "full_dev_selection": {k: deploy_choice[k] for k in ("model", "C", "m", "delta")},
                  "stage2_must_not_refit": "N/A (NULL frozen; nothing to apply)"}
    frozen["sha256"] = hashlib.sha256(json.dumps(frozen, sort_keys=True, default=str).encode()).hexdigest()

    result = {
        "contract": "artifacts/milestone_7_decision/negative_reducer/CONTRACT.md",
        "corrections": ["PROTOCOL CORRECTION 1", "PROTOCOL CORRECTION 2"],
        "provenance": {"git_commit": _git_commit(), "sklearn": sklearn.__version__, "scipy": scipy.__version__,
                       "pythonhashseed": os.environ.get("PYTHONHASHSEED", "unset")},
        "hgb_available": nr.hgb_available(), "n_dev_rows": int(len(df)),
        "n_dev_patients": int(df["patient_id"].nunique()),
        "outer_folds": {h: {kk: vv for kk, vv in folds[h].items()} for h in folds},
        "deploy_full_dev_candidates": deploy_cands, "deploy_full_dev_hgb_diagnostic": deploy_hgb,
        "aggregate": {"macro_delta_hits20": round(agg_delta, 4), "matched_random_delta": round(agg_rand, 4),
                      "paired_bootstrap": agg_boot, "worst_study_delta": round(worst_delta, 4),
                      "agg_pos": agg_pos, "agg_pos_removed": agg_pos_removed,
                      "agg_cp_lb_retention": round(agg_cp, 4), "improve_cp_lb_retention": round(improve_cp, 4)},
        "eligibility": {"every_study_noncatastrophic": every_noncat, "aggregate_gain_beats_random": agg_gain,
                        "any_negative_removal": any_removal, "aggregate_cp_ge_0_95": bool(agg_cp >= 0.95),
                        "improve_cp_ge_0_95": bool(improve_cp >= 0.95), "ELIGIBLE": eligible},
        "frozen": frozen, "verdict": "NON-NULL GATE FROZEN" if eligible else "NULL FROZEN (honest negative)",
    }
    return result


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    df = nr.load_dev()
    cols = nr.PORTABLE
    result = run_selection(df, cols)
    result["provenance"]["data_file_sha256"] = {p: _file_sha(Path(p)) for p in sorted(nr.ALLOWED_DATA_FILES)}
    if "--sha" in argv:                                               # repro subprocess: print SHA, write nothing
        print(result["frozen"]["sha256"])
        return 0

    # cross-PYTHONHASHSEED reproducibility (Correction 2.6)
    env = {**os.environ, "PYTHONHASHSEED": "1"}
    try:
        repro = subprocess.check_output([sys.executable, "-m", "scripts.negative_reducer_run", "--sha"],
                                        env={**env, "PYTHONPATH": "src"}, text=True).strip()
    except Exception as e:
        repro = f"repro-failed:{e}"
    result["provenance"]["repro_frozen_sha256_hashseed1"] = repro
    result["provenance"]["reproducible"] = bool(repro == result["frozen"]["sha256"])

    ART.mkdir(parents=True, exist_ok=True)
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(json.dumps(result["frozen"], indent=2, default=str) + "\n")
    (ART / "nested_loso.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    (ART / "REPORT.md").write_text(render_md(result))
    for h in nr.STUDIES:
        o = result["outer_folds"][h]["outer_test"]
        s = result["outer_folds"][h]["selected"]
        print(f"[{h}] {s['model']} C={s['C']} m={s['m']} inner Δ={s['delta']:+.4f} -> retention={o['raw_retention']} "
              f"neg_rm={o['neg_removal_frac']} Δhits20={o['macro_delta_hits20']:+.4f} rand={o['matched_random_delta']:+.4f} "
              f"beats={o['beats_random']} {o['cp_claim']}")
    a = result["aggregate"]
    print("AGG Δhits20=%+.4f rand=%+.4f boot=%s aggCP=%.4f improveCP=%.4f worst=%+.4f -> %s"
          % (a["macro_delta_hits20"], a["matched_random_delta"], a["paired_bootstrap"]["ci95"],
             a["agg_cp_lb_retention"], a["improve_cp_lb_retention"], a["worst_study_delta"], result["verdict"]))
    print("frozen sha:", result["frozen"]["sha256"][:16], "reproducible:", result["provenance"]["reproducible"])
    return 0


def render_md(r):
    L = ["# Risk-controlled negative reducer — nested LOSO (non-Sid)\n",
         f"_DEV {r['n_dev_patients']} patients / {r['n_dev_rows']} rows. NO Sid/Miller. HGB available "
         f"{r['hgb_available']} (diagnostic-only). git {r['provenance']['git_commit'][:10]}; "
         f"sklearn {r['provenance']['sklearn']}; scipy {r['provenance']['scipy']}; "
         f"reproducible={r['provenance'].get('reproducible')}._\n",
         "\n## Outer leave-one-study-out\n",
         "| held-out | chosen (model,C,m) | inner Δ | test retention | CP-LB (pow) | neg removed | "
         "Δhits@20 | matched-rand | beats | pos-rank↑ |",
         "|---|---|--:|--:|--:|--:|--:|--:|:--:|--:|"]
    for h in nr.STUDIES:
        s, o = r["outer_folds"][h]["selected"], r["outer_folds"][h]["outer_test"]
        L.append(f"| {h} | {s['model']},{s['C']},m={s['m']} | {s['delta']:+.3f} | {o['raw_retention']} | "
                 f"{o['cp_lb_retention']} ({'y' if o['cp_powered'] else 'n'}) | {o['neg_removal_frac']} | "
                 f"{o['macro_delta_hits20']:+.3f} | {o['matched_random_delta']:+.3f} | "
                 f"{'y' if o['beats_random'] else 'n'} | {o['mean_pos_rank_improvement']:+.2f} |")
    a = r["aggregate"]
    L.append(f"\n## Aggregate\nΔhits@20 = **{a['macro_delta_hits20']:+.4f}** (matched-random "
             f"{a['matched_random_delta']:+.4f}); bootstrap {a['paired_bootstrap']['mean_delta']:+.4f} "
             f"CI{a['paired_bootstrap']['ci95']} p>0={a['paired_bootstrap']['p_gt_0']}; worst study "
             f"{a['worst_study_delta']:+.4f}. Aggregate CP-LB retention {a['agg_cp_lb_retention']} "
             f"({a['agg_pos_removed']}/{a['agg_pos']} pos removed); IMPROVE CP-LB {a['improve_cp_lb_retention']}.\n")
    L.append("\n## §5 eligibility\n" + "\n".join(f"- {k}: **{v}**" for k, v in r["eligibility"].items()))
    fr = r["frozen"]
    if fr.get("frozen") == "NULL":
        L.append(f"\n## FROZEN: **NULL** — {fr['reason']}\nfull-DEV selection {fr['full_dev_selection']}; "
                 f"SHA `{fr['sha256'][:16]}`. Honest negative.\n")
    else:
        p = fr["frozen"]
        L.append(f"\n## FROZEN: non-null gate {p['model']} m={p['m']} tau={p['tau']}\n"
                 f"coef={p.get('coef')} intercept={p.get('intercept')}; SHA `{fr['sha256'][:16]}`.\n")
    L.append(f"\n**Verdict: {r['verdict']}.** Sid/Miller locked; no application performed.\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
