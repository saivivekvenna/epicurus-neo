"""Runner tests for the nested-LOSO negative reducer — SYNTHETIC frames only (no cohort data, no network).

Covers Correction-2 audit fixes: model fit/score (nnlog + hgb), NULL behavior, HGB excluded from the
selectable/freezable set, OUT-OF-FOLD tau (not in-sample), outer-test LABEL BLINDNESS (test labels cannot
change the removal mask or tau), aggregate-CP eligibility recording, and deterministic full-DEV recipe
selection (identical frozen SHA across two calls).
"""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd
import pytest

nr = importlib.import_module("event_b.negative_reducer")
run = importlib.import_module("scripts.negative_reducer_run")


def synth_dev(seed=0):
    """Three synthetic 'studies' with the STUDIES names, patient-grouped, portable features, a signal so a
    keep-score model is meaningful. No cohort files touched."""
    rng = np.random.default_rng(seed)
    rows = []
    for study, npat, ncand in [("improve", 12, 60), ("gartner", 10, 40), ("multimer", 10, 40)]:
        for pi in range(npat):
            pid = f"{study}:{pi}"
            for ci in range(ncand):
                pos = rng.random() < 0.12
                # positives tend to have LOWER prime/el (better) + higher expr
                prime = rng.uniform(0, 3) if pos else rng.uniform(0, 60)
                el = rng.uniform(0, 3) if pos else rng.uniform(0, 60)
                expr = rng.uniform(5, 20) if pos else rng.uniform(0, 15)
                rows.append({"study": study, "patient_id": pid, "mut_peptide": f"{study}{pi}_{ci}PEP",
                             "label": "POSITIVE" if pos else "TESTED_NEGATIVE",
                             "prime": prime, "el": el, "expr": expr,
                             "VarAlFreq": np.nan, "rna_af": np.nan, "ValMutRNACoef": np.nan, "CelPrev": np.nan})
    df = pd.DataFrame(rows)
    # guarantee >=1 positive per study for grouped CV stability
    for study in nr.STUDIES:
        if not ((df.study == study) & (df.label == "POSITIVE")).any():
            df.loc[df.index[(df.study == study)][0], "label"] = "POSITIVE"
    return df.reset_index(drop=True)


def test_fit_and_score_nnlog_and_hgb():
    df = synth_dev()
    for model, C in [("nnlog", 1.0), ("hgb", None)]:
        fm = run.fit_model(df, nr.PORTABLE, model, C)
        s = run.score_model(fm, df, nr.PORTABLE)
        assert s.shape == (len(df),) and np.all((s >= 0) & (s <= 1))
    fm = run.fit_model(df, nr.PORTABLE, "nnlog", 1.0)
    assert (np.asarray(fm["coef"]) >= -1e-8).all()   # frozen-eligible model stays nonnegative/monotone


def test_inner_select_excludes_hgb_from_selectable():
    df = synth_dev()
    train = df[df.study != "multimer"].reset_index(drop=True)
    chosen, cands, hgb_diag, failures = run.inner_select(train, nr.PORTABLE)
    assert chosen["model"] in {"NULL", "nnlog"}                 # HGB never selectable/freezable
    assert all(c["model"] in {"NULL", "nnlog"} for c in cands)
    assert nr.hgb_available() == (len(hgb_diag) > 0)            # HGB present only as a diagnostic table
    assert all(h["model"] == "hgb" for h in hgb_diag)
    assert isinstance(failures, list)


def test_null_removes_nothing():
    df = synth_dev()
    train = df[df.study != "multimer"].reset_index(drop=True)
    test = df[df.study == "multimer"].reset_index(drop=True)
    res = run.apply_to_test(train, test, nr.PORTABLE, {"model": "NULL", "C": None, "m": 0})
    assert res["neg_removed"] == 0 and res["pos_removed"] == 0 and res["removed_idx"] == []
    assert res["macro_delta_hits20"] == 0.0


def test_outer_test_label_blindness():
    # permuting the held-out study's LABELS must not change the removal mask or tau (removal uses features +
    # train-derived tau only). Metrics may change; the DECISION must not.
    df = synth_dev()
    train = df[df.study != "multimer"].reset_index(drop=True)
    test = df[df.study == "multimer"].reset_index(drop=True)
    choice = {"model": "nnlog", "C": 1.0, "m": 5}
    base = run.apply_to_test(train, test, nr.PORTABLE, choice)
    test_perm = test.copy()
    rng = np.random.default_rng(7)
    test_perm["label"] = rng.permutation(test_perm["label"].to_numpy())
    perm = run.apply_to_test(train, test_perm, nr.PORTABLE, choice)
    assert base["removed_idx"] == perm["removed_idx"]           # decision is label-blind
    assert base["tau"] == perm["tau"]


def test_tau_is_out_of_fold_not_in_sample():
    # apply_to_test must calibrate tau on patient-grouped OOF scores, NOT on the full-refit in-sample scores.
    df = synth_dev()
    train = df[df.study != "multimer"].reset_index(drop=True)
    test = df[df.study == "multimer"].reset_index(drop=True)
    choice = {"model": "nnlog", "C": 1.0, "m": 5}
    res = run.apply_to_test(train, test, nr.PORTABLE, choice)
    # recompute the OOF-calibrated tau independently -> must match
    oof, _, _ = run.make_oof(train, nr.PORTABLE, "nnlog", 1.0)
    exp_tau, *_ = nr.calibrate_tau(oof, run._is_pos(train), nr.removable_mask(train, nr.PORTABLE, 5, set()))
    exp = None if not np.isfinite(exp_tau) else float(exp_tau)
    assert res["tau"] == exp
    # and it is NOT the in-sample tau (full-refit scores) unless they coincide by chance
    fm = run.fit_model(train, nr.PORTABLE, "nnlog", 1.0)
    in_s = run.score_model(fm, train, nr.PORTABLE)
    ins_tau, *_ = nr.calibrate_tau(in_s, run._is_pos(train), nr.removable_mask(train, nr.PORTABLE, 5, set()))
    ins = None if not np.isfinite(ins_tau) else float(ins_tau)
    # they need not be equal; assert apply used OOF (already checked) — this documents the distinction
    assert res["tau"] == exp and (ins is None or isinstance(ins, float))


def test_make_oof_uses_fold_local_weights():
    # OOF fold models must fit on fold-local weights: with a study whose class balance differs, using global
    # weights vs local weights yields different OOF scores. Assert make_oof matches an explicit fold-local recompute.
    df = synth_dev(3)
    train = df[df.study != "gartner"].reset_index(drop=True)
    oof, _, ok = run.make_oof(train, nr.PORTABLE, "nnlog", 1.0)
    assert ok and np.isfinite(oof).all()                      # all folds converged + scored out-of-fold
    from sklearn.model_selection import GroupKFold
    X = nr.feat_matrix(train, nr.PORTABLE)
    y = run._is_pos(train)
    groups = train["patient_id"].to_numpy()
    exp = np.full(len(train), np.nan)
    for tr, va in GroupKFold(n_splits=min(4, len(np.unique(groups)))).split(X, y, groups):
        w_tr = nr.balanced_weights(train.iloc[tr])            # fold-local
        b0, coef, _ = nr.fit_nnlogistic(X[tr], y[tr], w_tr, 1.0)
        exp[va] = nr.nnlogistic_score(X[va], b0, coef)
    assert np.allclose(oof, exp)


def test_run_selection_records_aggregate_and_improve_cp():
    df = synth_dev()
    r = run.run_selection(df, nr.PORTABLE)
    e = r["eligibility"]
    assert "aggregate_cp_ge_0_95" in e and "improve_cp_ge_0_95" in e
    a = r["aggregate"]
    assert "agg_cp_lb_retention" in a and "improve_cp_lb_retention" in a and "agg_pos" in a
    # loso eligibility is a strict AND of the recorded conditions
    assert e["loso_eligible"] == bool(e["every_study_noncatastrophic"] and e["aggregate_gain_beats_random"]
                                      and e["any_negative_removal"] and e["aggregate_cp_ge_0_95"]
                                      and e["improve_cp_ge_0_95"])
    # final ELIGIBLE also requires a valid deploy recipe (Correction 3.2)
    assert e["ELIGIBLE"] == bool(e["loso_eligible"] and e["deploy_recipe_valid"])


def test_run_selection_improve_cp_is_exact_not_rounded():
    # Correction 3.3: improve_cp_ge_0_95 must be computed from EXACT CP on raw n_pos/pos_removed, not the
    # rounded display field.
    df = synth_dev()
    r = run.run_selection(df, nr.PORTABLE)
    a = r["aggregate"]
    exact = nr.cp_lower(a["improve_n_pos"] - a["improve_pos_removed"], a["improve_n_pos"])
    assert r["eligibility"]["improve_cp_ge_0_95"] == bool(exact >= 0.95)


def test_cp_rounding_boundary_hazard():
    # concrete case the exact check defends against: retain 373/385 positives -> true CP 0.949988 < 0.95,
    # but round(.,4) == 0.9500 would WRONGLY pass. Eligibility must use the exact value.
    v = nr.cp_lower(373, 385)
    assert v < 0.95 and round(v, 4) == 0.95
    assert (v >= 0.95) is False               # exact check correctly FAILS
    assert (round(v, 4) >= 0.95) is True      # rounded check would WRONGLY pass -> why we keep exact


def test_full_refit_failure_keeps_all(monkeypatch):
    # Correction 3.1 end-to-end: OOF fits succeed (aggressive tau possible) but the FULL refit fails ->
    # apply_to_test must KEEP-ALL (remove nothing), never threshold the degenerate constant-0.5 model.
    df = synth_dev()
    train = df[df.study != "multimer"].reset_index(drop=True)
    test = df[df.study == "multimer"].reset_index(drop=True)
    orig = run.fit_model
    monkeypatch.setattr(run, "fit_model",
                        lambda d, c, m, C: {**orig(d, c, m, C), "ok": False} if m == "nnlog" else orig(d, c, m, C))
    res = run.apply_to_test(train, test, nr.PORTABLE, {"model": "nnlog", "C": 1.0, "m": 5})
    assert res["fit_failed"] and res["removed_idx"] == [] and res["neg_removed"] == 0


def test_deploy_validate_rejects_null_and_accepts_removing_recipe():
    # Correction 3.2: NULL / delta<=0 recipes are invalid; a real removing recipe is valid.
    df = synth_dev()
    assert run.deploy_validate(df, nr.PORTABLE, {"model": "NULL", "C": None, "m": 0, "delta": 0.0}, True) == (False, None, 0)
    assert run.deploy_validate(df, nr.PORTABLE, {"model": "nnlog", "C": 1.0, "m": 0, "delta": 0.5}, False)[0] is False
    ok, payload, rn = run.deploy_validate(df, nr.PORTABLE, {"model": "nnlog", "C": 1.0, "m": 0, "delta": 0.5}, True)
    assert ok and payload is not None and rn >= 1 and payload["tau"] is not None


def test_finalize_guard_downgrades_non_null_on_repro_mismatch():
    # Correction 3.5: a non-null frozen config that is not reproducible must be downgraded to NULL.
    res = {"frozen": {"frozen": {"model": "nnlog", "m": 0}, "full_dev_selection": {"model": "nnlog"}},
           "provenance": {"reproducible": False}, "eligibility": {"ELIGIBLE": True}, "verdict": "NON-NULL GATE FROZEN"}
    out = run.finalize_guard(res)
    assert out["frozen"]["frozen"] == "NULL" and "repro_mismatch" in out["frozen"]["reason"]
    assert out["eligibility"]["ELIGIBLE"] is False and "sha256" in out["frozen"]
    # a reproducible non-null config is left intact
    res2 = {"frozen": {"frozen": {"model": "nnlog"}}, "provenance": {"reproducible": True},
            "eligibility": {"ELIGIBLE": True}, "verdict": "NON-NULL GATE FROZEN"}
    assert run.finalize_guard(res2)["frozen"]["frozen"] == {"model": "nnlog"}


def test_validate_manifest_accepts_matching_and_rejects_violations():
    # Correction 3.6: build a tiny synthetic corpus + matching expected override; then break each invariant.
    def corpus():
        rows = []
        spec = {"improve": (2, 1, 3), "gartner": (2, 1, 2), "multimer": (2, 1, 2)}
        for study, (npat, pos, neg) in spec.items():
            labels = ["POSITIVE"] * pos + ["TESTED_NEGATIVE"] * neg   # round-robin over npat patients
            for i, lab in enumerate(labels):
                rows.append({"study": study, "patient_id": f"{study}:{i % npat}", "mut_peptide": f"{study}P{i}",
                             "label": lab, "prime": 1.0, "el": 1.0, "expr": 1.0})
        return pd.DataFrame(rows).reset_index(drop=True)

    df = corpus()
    exp = {"studies": ["improve", "gartner", "multimer"],
           "per_study": {s: (df[df.study == s].patient_id.nunique(),
                             int((df[df.study == s].label == "POSITIVE").sum()),
                             int((df[df.study == s].label == "TESTED_NEGATIVE").sum())) for s in
                         ["improve", "gartner", "multimer"]},
           "total": (df.patient_id.nunique(),
                     int((df.label == "POSITIVE").sum()), int((df.label == "TESTED_NEGATIVE").sum())),
           "required_cols": ["study", "patient_id", "mut_peptide", "label", "prime", "el", "expr"],
           "labels": {"POSITIVE", "TESTED_NEGATIVE"}}
    good_hashes = {"f": "abc"}
    assert run.validate_manifest(df, good_hashes, exp) is True
    # MISSING hash -> raise
    with pytest.raises(ValueError):
        run.validate_manifest(df, {"f": "MISSING"}, exp)
    # wrong study set -> raise
    bad = df.copy()
    bad.loc[bad.study == "multimer", "study"] = "OTHER"
    with pytest.raises(ValueError):
        run.validate_manifest(bad.reset_index(drop=True), good_hashes, exp)
    # third label -> raise
    bad2 = df.copy()
    bad2.loc[0, "label"] = "UNKNOWN"
    with pytest.raises(ValueError):
        run.validate_manifest(bad2, good_hashes, exp)
    # non-reset index -> raise
    with pytest.raises(ValueError):
        run.validate_manifest(df.iloc[1:], good_hashes, exp)


def test_full_dev_recipe_selection_is_deterministic():
    df = synth_dev()
    r1 = run.run_selection(df, nr.PORTABLE)
    r2 = run.run_selection(df, nr.PORTABLE)
    assert r1["frozen"]["sha256"] == r2["frozen"]["sha256"]       # deterministic freeze
    assert r1["frozen"]["full_dev_selection"] == r2["frozen"]["full_dev_selection"]
    assert r1["frozen"].get("frozen") == r2["frozen"].get("frozen")
