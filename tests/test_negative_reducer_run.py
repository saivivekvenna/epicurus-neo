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
    chosen, cands, hgb_diag = run.inner_select(train, nr.PORTABLE)
    assert chosen["model"] in {"NULL", "nnlog"}                 # HGB never selectable/freezable
    assert all(c["model"] in {"NULL", "nnlog"} for c in cands)
    assert nr.hgb_available() == (len(hgb_diag) > 0)            # HGB present only as a diagnostic table
    assert all(h["model"] == "hgb" for h in hgb_diag)


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
    oof, _ = run.make_oof(train, nr.PORTABLE, "nnlog", 1.0)
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
    oof, _ = run.make_oof(train, nr.PORTABLE, "nnlog", 1.0)
    assert np.isfinite(oof).all()                              # all rows scored out-of-fold
    from sklearn.model_selection import GroupKFold
    X = nr.feat_matrix(train, nr.PORTABLE)
    y = run._is_pos(train)
    groups = train["patient_id"].to_numpy()
    exp = np.full(len(train), np.nan)
    for tr, va in GroupKFold(n_splits=min(4, len(np.unique(groups)))).split(X, y, groups):
        w_tr = nr.balanced_weights(train.iloc[tr])            # fold-local
        b0, coef = nr.fit_nnlogistic(X[tr], y[tr], w_tr, 1.0)
        exp[va] = nr.nnlogistic_score(X[va], b0, coef)
    assert np.allclose(oof, exp)


def test_run_selection_records_aggregate_and_improve_cp():
    df = synth_dev()
    r = run.run_selection(df, nr.PORTABLE)
    e = r["eligibility"]
    assert "aggregate_cp_ge_0_95" in e and "improve_cp_ge_0_95" in e
    a = r["aggregate"]
    assert "agg_cp_lb_retention" in a and "improve_cp_lb_retention" in a and "agg_pos" in a
    # eligibility is a strict AND of the recorded conditions
    assert e["ELIGIBLE"] == bool(e["every_study_noncatastrophic"] and e["aggregate_gain_beats_random"]
                                 and e["any_negative_removal"] and e["aggregate_cp_ge_0_95"]
                                 and e["improve_cp_ge_0_95"])


def test_full_dev_recipe_selection_is_deterministic():
    df = synth_dev()
    r1 = run.run_selection(df, nr.PORTABLE)
    r2 = run.run_selection(df, nr.PORTABLE)
    assert r1["frozen"]["sha256"] == r2["frozen"]["sha256"]       # deterministic freeze
    assert r1["frozen"]["full_dev_selection"] == r2["frozen"]["full_dev_selection"]
    assert r1["frozen"].get("frozen") == r2["frozen"].get("frozen")
