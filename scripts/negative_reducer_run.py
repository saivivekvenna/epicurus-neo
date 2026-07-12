"""Risk-controlled negative reducer — NESTED LEAVE-ONE-STUDY-OUT runner (CONTRACT + CORRECTIONS 1-3).

Outer = leave-one-STUDY-out over {IMPROVE, Gartner, multimer}. Inner = patient-grouped CV within outer-train
selecting (model in {NULL, nonnegative-logistic C-grid}, protected core m in {0,5,10}) by inner-OOF
patient-macro Delta-hits@20; tau = most aggressive cut whose pooled OOF CP-95%-LB retention >= 0.95.
Monotonic HGB is DIAGNOSTIC-only (never selected/frozen). tau always OUT-OF-FOLD, never in-sample. Any
non-converged fit => KEEP-ALL / ineligible (Correction 3.1). Non-null freeze also requires a valid all-DEV
deploy recipe (3.2). CP recomputed exactly (3.3). Payload is fully applyable incl OOD envelope + versions,
with a pure apply_payload equivalence path (3.4). Reproducibility + manifest are fail-closed (3.5/3.6).
Reads NO Sid / Miller file. Freezes a fitted apply-only payload iff eligible, else NULL.

    PYTHONPATH=src python -m scripts.negative_reducer_run          # full run + freeze + repro + manifest
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

# Correction 3.6: exact expected corpus from CONTRACT §1 (patients, positives, negatives).
EXPECTED_MANIFEST = {
    "studies": ["improve", "gartner", "multimer"],
    "per_study": {"improve": (70, 467, 17053), "gartner": (26, 46, 3722), "multimer": (26, 34, 8069)},
    "total": (122, 547, 28844),
    "required_cols": ["study", "patient_id", "mut_peptide", "label", "prime", "el", "expr"],
    "labels": {"POSITIVE", "TESTED_NEGATIVE"},
}


def _is_pos(df):
    return (df["label"].to_numpy() == "POSITIVE").astype(int)


# ---- uniform model interface (nnlog serializable/frozen-eligible; hgb diagnostic-only) ----------------
def fit_model(df, cols, model, C):
    """Returns a model dict with an 'ok' convergence flag. nnlog carries coef/intercept; hgb carries clf."""
    X, y, w = nr.feat_matrix(df, cols), _is_pos(df), nr.balanced_weights(df)
    if model == "nnlog":
        b0, coef, ok = nr.fit_nnlogistic(X, y, w, C)
        return {"kind": "nnlog", "intercept": b0, "coef": coef, "ok": bool(ok)}
    return {"kind": "hgb", "clf": nr.fit_hgb(X, y, w, len(cols)), "ok": True}


def score_model(fm, df, cols):
    X = nr.feat_matrix(df, cols)
    if fm["kind"] == "nnlog":
        return nr.nnlogistic_score(X, fm["intercept"], fm["coef"])
    return fm["clf"].predict_proba(X)[:, 1]


def make_oof(train: pd.DataFrame, cols, model: str, C):
    """Inner patient-grouped OOF keep-scores + per-row leakage-clean mask. Fold models use FOLD-LOCAL
    balanced weights. Returns (oof, clean, ok); ok=False if ANY fold fit failed to converge (=> the whole
    (model,C) candidate is ineligible)."""
    X = nr.feat_matrix(train, cols)
    y = _is_pos(train)
    groups = train["patient_id"].to_numpy()
    n_groups = len(np.unique(groups))
    oof = np.full(len(train), np.nan)
    clean = np.ones(len(train), bool)
    ok = True
    for tr, va in GroupKFold(n_splits=min(4, n_groups)).split(X, y, groups):
        w_tr = nr.balanced_weights(train.iloc[tr])                    # fold-local weights
        if model == "nnlog":
            b0, coef, succ = nr.fit_nnlogistic(X[tr], y[tr], w_tr, C)
            ok = ok and succ
            oof[va] = nr.nnlogistic_score(X[va], b0, coef)
        else:
            clf = nr.fit_hgb(X[tr], y[tr], w_tr, len(cols))
            oof[va] = clf.predict_proba(X[va])[:, 1]
        trainpep = set(train.iloc[tr]["mut_peptide"].map(canonical_peptide))
        clean[va] = nr._clean_against(train.iloc[va].reset_index(drop=True), trainpep)
    return oof, clean, ok


def objective_for_m(train, cols, oof, oof_clean, m):
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
    """Selection over {NULL, nnlog x C} x m (HGB excluded). Non-converged (model,C) candidates are skipped
    and logged. Returns (chosen, selectable_table, hgb_diagnostic_table, failures)."""
    cands = [{"model": "NULL", "C": None, "m": 0, "delta": 0.0, "tau": None, "tau_inf": False,
              "neg_removed": 0, "removed_total": 0, "cp_lb": 1.0, "powered": True}]
    failures = []
    for C in nr.C_GRID:
        oof, oof_clean, ok = make_oof(train, cols, "nnlog", C)
        if not ok:
            failures.append({"model": "nnlog", "C": C, "reason": "optimizer_non_convergence"})
            continue
        for m in nr.M_GRID:
            cands.append({"model": "nnlog", "C": C, **objective_for_m(train, cols, oof, oof_clean, m)})
    hgb_diag = []
    if nr.hgb_available():
        oof, oof_clean, _ = make_oof(train, cols, "hgb", None)
        for m in nr.M_GRID:
            hgb_diag.append({"model": "hgb", "C": None, **objective_for_m(train, cols, oof, oof_clean, m)})
    chosen = max(cands, key=_tie)
    return chosen, cands, hgb_diag, failures


def paired_bootstrap(gate: dict, null: dict, n=5000, seed=12345):
    pids = sorted(null)
    d = np.array([gate[p] - null[p] for p in pids], float)
    rng = np.random.default_rng(seed)
    means = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(n)]) if len(d) else np.zeros(1)
    return {"mean_delta": round(float(d.mean()) if len(d) else 0.0, 4),
            "ci95": [round(float(np.percentile(means, 2.5)), 4), round(float(np.percentile(means, 97.5)), 4)],
            "p_gt_0": round(float((means > 0).mean()), 4), "n_patients": len(pids)}


def apply_to_test(train, test, cols, choice):
    """tau from patient-grouped OOF on outer-train; full model refit on outer-train ONLY to SCORE the
    untouched held-out study. Any failed fit (OOF fold or full refit) => KEEP-ALL (Correction 3.1). Removal
    depends only on test FEATURES + train-derived tau, never on test labels."""
    model, m = choice["model"], choice.get("m", 0)
    fit_failed = False
    if model == "NULL":
        removed, tau, ood = np.zeros(len(test), bool), np.inf, set()
    else:
        oof, _, oof_ok = make_oof(train, cols, model, choice["C"])
        rem_train = nr.removable_mask(train, cols, m, ood=set())
        tau, *_ = nr.calibrate_tau(oof, _is_pos(train), rem_train)     # OOF-calibrated, not in-sample
        fm = fit_model(train, cols, model, choice["C"])
        if not (oof_ok and fm["ok"]):
            removed, ood, fit_failed = np.zeros(len(test), bool), set(), True   # KEEP-ALL, never threshold a bad model
        else:
            ood = nr.ood_patients(train, test, cols)
            removed = nr.gate_removed(test, score_model(fm, test, cols), tau, cols, m, ood)

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
        "fit_failed": bool(fit_failed),
        "raw_retention": round(retention, 4), "cp_lb_retention": round(cp_lb, 4), "cp_powered": bool(powered),
        "cp_claim": (f"retention CP-95%-LB={cp_lb:.4f}") if powered
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
    """Full-DEV freeze payload for a SELECTABLE recipe. OOF-calibrate tau, refit on all DEV. Returns None if
    the full-DEV fit fails (=> caller freezes NULL). Payload is fully applyable (Correction 3.4): coef/
    intercept + exact OOD envelope + cover + versions; model-payload SHA covers all of it."""
    model, m = choice["model"], choice.get("m", 0)
    if model == "NULL":
        return {"model": "NULL", "m": 0, "tau": None, "note": "no-gate"}
    oof, _, oof_ok = make_oof(df, cols, model, choice["C"])
    tau, *_ = nr.calibrate_tau(oof, _is_pos(df), nr.removable_mask(df, cols, m, ood=set()))
    fm = fit_model(df, cols, model, choice["C"])
    if not (oof_ok and fm["ok"]):
        return None                                              # fail closed => NULL freeze
    payload = {"feature_order": list(cols), "model": model, "C": choice["C"], "m": m,
               "tau": (None if not np.isfinite(tau) else float(tau)),
               "coef": [float(x) for x in fm["coef"]], "intercept": float(fm["intercept"]),
               "ood_envelope": {c: nr.raw_envelope(df, cols)[c] for c in cols}, "ood_cover": 0.5,
               "sklearn_version": sklearn.__version__, "scipy_version": scipy.__version__,
               "percentile_policy": "within-patient rank(pct=True), oriented; NaN->0.5",
               "direction_policy": {c: ("higher_raw_better" if nr.HIGHER_BETTER[c] else "lower_raw_better") for c in cols},
               "protected_core": "top-m genuine-PRIME per patient unremovable; boundary ties all protected",
               "removal_rule": "apply_payload: remove non-core, feature-present, in-support candidates with keepscore<tau (tau=None => none)"}
    payload["model_payload_sha256"] = hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return payload


def _file_sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest() if Path(p).exists() else "MISSING"


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def validate_manifest(df: pd.DataFrame, file_hashes: dict, expected: dict = EXPECTED_MANIFEST):
    """Correction 3.6: fail-closed corpus check BEFORE selection. Raises ValueError on any mismatch."""
    if set(df["study"].unique()) != set(expected["studies"]):
        raise ValueError(f"manifest: studies {sorted(df['study'].unique())} != {expected['studies']}")
    if not set(expected["required_cols"]).issubset(df.columns):
        raise ValueError(f"manifest: missing columns {set(expected['required_cols']) - set(df.columns)}")
    if set(df["label"].unique()) != expected["labels"]:
        raise ValueError(f"manifest: labels {sorted(df['label'].unique())} != {sorted(expected['labels'])}")
    if not (df.index.is_unique and list(df.index) == list(range(len(df)))):
        raise ValueError("manifest: index is not a unique 0..n-1 reset index")
    tot_pat = tot_pos = tot_neg = 0
    for s, (ep, epos, eneg) in expected["per_study"].items():
        sub = df[df["study"] == s]
        pat, pos, neg = sub["patient_id"].nunique(), int((sub["label"] == "POSITIVE").sum()), int((sub["label"] == "TESTED_NEGATIVE").sum())
        if (pat, pos, neg) != (ep, epos, eneg):
            raise ValueError(f"manifest: {s} (pat,pos,neg)={(pat, pos, neg)} != {(ep, epos, eneg)}")
        tot_pat, tot_pos, tot_neg = tot_pat + pat, tot_pos + pos, tot_neg + neg
    if (tot_pat, tot_pos, tot_neg) != tuple(expected["total"]):
        raise ValueError(f"manifest: totals {(tot_pat, tot_pos, tot_neg)} != {tuple(expected['total'])}")
    missing = [k for k, v in file_hashes.items() if v == "MISSING"]
    if missing:
        raise ValueError(f"manifest: data files missing hash: {missing}")
    return True


def deploy_validate(df, cols, choice, loso_eligible):
    """Correction 3.2: a non-null freeze also requires a VALID all-DEV deploy recipe: model!=NULL, delta>0,
    finite/applyable tau, successful full-DEV fit (serialize_recipe not None), and >=1 real removal on DEV.
    Returns (deploy_ok, payload_or_None, removed_negatives)."""
    if not loso_eligible or choice["model"] == "NULL" or choice.get("delta", 0.0) <= 0:
        return False, None, 0
    payload = serialize_recipe(df, cols, choice)
    if payload is None or payload.get("model") in (None, "NULL") or payload.get("tau") is None:
        return False, payload, 0
    removed = nr.apply_payload(df, payload)
    removed_neg = int((removed & (_is_pos(df) == 0)).sum())
    return (removed_neg >= 1), payload, removed_neg


def finalize_guard(result: dict) -> dict:
    """Correction 3.5: refuse to publish a NON-NULL config if it is not cross-PYTHONHASHSEED reproducible.
    Downgrade to NULL (reason repro_mismatch) BEFORE any artifact write."""
    fr = result["frozen"]
    if fr.get("frozen") != "NULL" and not result.get("provenance", {}).get("reproducible", False):
        null = {"name": "negative_reducer", "version": "1.0.0", "frozen": "NULL",
                "reason": "repro_mismatch: cross-PYTHONHASHSEED frozen SHA mismatch/failure; refusing to publish non-null",
                "full_dev_selection": fr.get("full_dev_selection"),
                "stage2_must_not_refit": "N/A (NULL frozen; nothing to apply)"}
        null["sha256"] = hashlib.sha256(json.dumps(null, sort_keys=True, default=str).encode()).hexdigest()
        result["frozen"] = null
        result["eligibility"]["ELIGIBLE"] = False
        result["verdict"] = "NULL FROZEN (repro fail-closed)"
    return result


def run_selection(df, cols):
    """Full nested LOSO + eligibility + frozen decision. Pure over df (no file writes). Deterministic."""
    folds, agg_gate, agg_null = {}, {}, {}
    agg_pos = agg_pos_removed = 0
    improve_np = improve_pr = 0
    for held in nr.STUDIES:
        train = df[df["study"] != held].reset_index(drop=True)
        test = df[df["study"] == held].reset_index(drop=True)
        chosen, cands, hgb_diag, failures = inner_select(train, cols)
        res = apply_to_test(train, test, cols, chosen)
        agg_gate.update(res.pop("_gate_hits"))
        agg_null.update(res.pop("_null_hits"))
        agg_pos += res["n_pos"]
        agg_pos_removed += res["pos_removed"]
        if held == "improve":
            improve_np, improve_pr = res["n_pos"], res["pos_removed"]
        folds[held] = {"selected": {k: chosen[k] for k in ("model", "C", "m", "delta", "tau")},
                       "inner_candidates": cands, "inner_failures": failures, "hgb_diagnostic": hgb_diag,
                       "outer_test": res}

    agg_delta = float(np.mean([agg_gate[p] - agg_null[p] for p in agg_null]))
    agg_boot = paired_bootstrap(agg_gate, agg_null)
    npat = {h: sum(1 for p in agg_null if p.startswith(h + ":")) for h in nr.STUDIES}
    agg_rand = float(sum(folds[h]["outer_test"]["matched_random_delta"] * npat[h] for h in nr.STUDIES) / sum(npat.values()))
    worst_delta = min(folds[h]["outer_test"]["macro_delta_hits20"] for h in nr.STUDIES)

    agg_cp = nr.cp_lower(agg_pos - agg_pos_removed, agg_pos)          # exact (Correction 3.3)
    improve_cp = nr.cp_lower(improve_np - improve_pr, improve_np)     # exact from raw ints, NOT rounded field
    every_noncat = all(folds[h]["outer_test"]["macro_delta_hits20"] >= nr.CATASTROPHIC for h in nr.STUDIES)
    agg_gain = agg_delta > 0 and agg_delta > agg_rand + 1e-9
    any_removal = any(folds[h]["outer_test"]["neg_removal_frac"] > 0 for h in nr.STUDIES)
    cp_ok = agg_cp >= 0.95 and improve_cp >= 0.95
    loso_eligible = bool(every_noncat and agg_gain and any_removal and cp_ok)

    # Correction 3.2/2.5: deployment recipe = same inner_select on ALL DEV; validate it before non-null freeze
    deploy_choice, deploy_cands, deploy_hgb, deploy_fail = inner_select(df, cols)
    deploy_ok, payload, deploy_removed_neg = deploy_validate(df, cols, deploy_choice, loso_eligible)
    eligible = bool(loso_eligible and deploy_ok)

    if eligible:
        frozen = {"name": "negative_reducer", "version": "1.0.0", "frozen": payload,
                  "full_dev_selection": {k: deploy_choice[k] for k in ("model", "C", "m", "delta")},
                  "deploy_removed_negatives_on_dev": deploy_removed_neg,
                  "stage2_must_not_refit": "Stage 2 applies this payload only via apply_payload; never refit."}
    else:
        reason = ("nested LOSO evidence did not pass CONTRACT §5 (see eligibility)" if not loso_eligible
                  else "LOSO passed but all-DEV deploy recipe invalid (NULL / no-removal / failed-fit / non-finite tau)")
        frozen = {"name": "negative_reducer", "version": "1.0.0", "frozen": "NULL", "reason": reason,
                  "full_dev_selection": {k: deploy_choice[k] for k in ("model", "C", "m", "delta")},
                  "deploy_ok": deploy_ok, "stage2_must_not_refit": "N/A (NULL frozen; nothing to apply)"}
    frozen["sha256"] = hashlib.sha256(json.dumps(frozen, sort_keys=True, default=str).encode()).hexdigest()

    return {
        "contract": "artifacts/milestone_7_decision/negative_reducer/CONTRACT.md",
        "corrections": ["PROTOCOL CORRECTION 1", "PROTOCOL CORRECTION 2", "PROTOCOL CORRECTION 3"],
        "provenance": {"git_commit": _git_commit(), "sklearn": sklearn.__version__, "scipy": scipy.__version__,
                       "pythonhashseed": os.environ.get("PYTHONHASHSEED", "unset")},
        "hgb_available": nr.hgb_available(), "n_dev_rows": int(len(df)),
        "n_dev_patients": int(df["patient_id"].nunique()),
        "outer_folds": folds, "deploy_full_dev_candidates": deploy_cands,
        "deploy_full_dev_hgb_diagnostic": deploy_hgb, "deploy_full_dev_failures": deploy_fail,
        "aggregate": {"macro_delta_hits20": round(agg_delta, 4), "matched_random_delta": round(agg_rand, 4),
                      "paired_bootstrap": agg_boot, "worst_study_delta": round(worst_delta, 4),
                      "agg_pos": agg_pos, "agg_pos_removed": agg_pos_removed,
                      "agg_cp_lb_retention": round(agg_cp, 4), "improve_cp_lb_retention": round(improve_cp, 4),
                      "improve_n_pos": improve_np, "improve_pos_removed": improve_pr},
        "eligibility": {"every_study_noncatastrophic": every_noncat, "aggregate_gain_beats_random": agg_gain,
                        "any_negative_removal": any_removal, "aggregate_cp_ge_0_95": bool(agg_cp >= 0.95),
                        "improve_cp_ge_0_95": bool(improve_cp >= 0.95), "loso_eligible": loso_eligible,
                        "deploy_recipe_valid": deploy_ok, "ELIGIBLE": eligible},
        "frozen": frozen, "verdict": "NON-NULL GATE FROZEN" if eligible else "NULL FROZEN (honest negative)",
    }


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    df = nr.load_dev()
    cols = nr.PORTABLE
    file_hashes = {p: _file_sha(Path(p)) for p in sorted(nr.ALLOWED_DATA_FILES)}
    validate_manifest(df, file_hashes)                               # Correction 3.6: abort on wrong corpus
    result = run_selection(df, cols)
    result["provenance"]["data_file_sha256"] = file_hashes
    if "--sha" in argv:                                              # repro subprocess: print SHA, write nothing
        print(result["frozen"]["sha256"])
        return 0

    # cross-PYTHONHASHSEED reproducibility (Correction 2.6/3.5)
    try:
        repro = subprocess.check_output([sys.executable, "-m", "scripts.negative_reducer_run", "--sha"],
                                        env={**os.environ, "PYTHONHASHSEED": "1", "PYTHONPATH": "src"},
                                        text=True).strip()
    except Exception as e:
        repro = f"repro-failed:{e}"
    result["provenance"]["repro_frozen_sha256_hashseed1"] = repro
    result["provenance"]["reproducible"] = bool(repro == result["frozen"]["sha256"])
    result = finalize_guard(result)                                 # fail-closed before writes

    ART.mkdir(parents=True, exist_ok=True)
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(json.dumps(result["frozen"], indent=2, default=str) + "\n")
    (ART / "nested_loso.json").write_text(json.dumps(result, indent=2, default=str) + "\n")
    (ART / "REPORT.md").write_text(render_md(result))
    for h in nr.STUDIES:
        o, s = result["outer_folds"][h]["outer_test"], result["outer_folds"][h]["selected"]
        print(f"[{h}] {s['model']} C={s['C']} m={s['m']} innerΔ={s['delta']:+.4f} -> ret={o['raw_retention']} "
              f"neg_rm={o['neg_removal_frac']} Δh20={o['macro_delta_hits20']:+.4f} rand={o['matched_random_delta']:+.4f} "
              f"beats={o['beats_random']} {o['cp_claim']}")
    a = result["aggregate"]
    print("AGG Δh20=%+.4f rand=%+.4f boot=%s aggCP=%.4f impCP=%.4f worst=%+.4f -> %s"
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
             f"{a['worst_study_delta']:+.4f}. Aggregate CP-LB {a['agg_cp_lb_retention']} "
             f"({a['agg_pos_removed']}/{a['agg_pos']} pos removed); IMPROVE CP-LB {a['improve_cp_lb_retention']} "
             f"({a['improve_pos_removed']}/{a['improve_n_pos']}).\n")
    L.append("\n## §5 eligibility\n" + "\n".join(f"- {k}: **{v}**" for k, v in r["eligibility"].items()))
    fr = r["frozen"]
    if fr.get("frozen") == "NULL":
        L.append(f"\n## FROZEN: **NULL** — {fr['reason']}\nfull-DEV selection {fr['full_dev_selection']}; "
                 f"SHA `{fr['sha256'][:16]}`. Honest negative.\n")
    else:
        p = fr["frozen"]
        L.append(f"\n## FROZEN: non-null gate {p['model']} m={p['m']} tau={p['tau']}\n"
                 f"coef={p.get('coef')} intercept={p.get('intercept')}; removed "
                 f"{fr['deploy_removed_negatives_on_dev']} DEV negatives; SHA `{fr['sha256'][:16]}`.\n")
    L.append(f"\n**Verdict: {r['verdict']}.** Sid/Miller locked; no application performed.\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
