"""Epicurus v0.5 DEVELOPMENT runner — deployable context-conditioned pairwise challenger.

DEVELOPMENT ONLY. Executes the FROZEN protocol
`artifacts/milestone_7_decision/epicurus_v05/PREREGISTERED_PROTOCOL.md` verbatim:

  * provenance re-verified (fail-fast); Gartner TEST never opened (I/O guard);
  * frozen comparators P (MIL) / A (AdditiveRanker) / F (source-name tower) are REFIT from exact frozen code at
    the stored per-fold hyperparameters (ZERO retuning) and reproduction-verified (P/A tight fail-fast, F to
    v0.4's nonconvex tolerance with the residual reported);
  * Q (shared pairwise, no context) and R (portable-context pairwise, the GATED candidate) are nested-selected
    inside each outer-train on the registered grid, then scored out-of-fold on the 5 frozen folds;
  * primary gate = R vs GENUINE PRIME (raw unmasked prime_rank), ACCEPT iff source-balanced patient-paired
    bootstrap CI_lower > 0; every registered diagnostic reported (none used for selection).

The A-vs-Q and A-vs-P contrasts are DESCRIPTIVE, never causal (A vs Q changes objective form AND Gartner
negative aggregation/bag discipline; Q vs P mixes objective + exact-witness supervision). Emits DEV_RESULT.json,
DEV_REPORT.md, and a frozen-config JSON. Nothing is committed by this script.

NOTE: the compute (event_b.epicurus_v05) is separated from this serialization layer so the math is unit-tested
without running the full benchmark.
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

from event_b.epicurus_v03 import (  # noqa: E402  frozen v0.3 evaluation, reused verbatim
    FEATURES, baseline_score, evaluate_model, load_frozen, paired_bootstrap, per_patient_metrics,
)
from event_b.epicurus_v04 import (  # noqa: E402
    assemble_frame, attrition_report, ext_metrics, guard_no_test_io,
)
from event_b import epicurus_v05 as v5  # noqa: E402

OUT = Path("artifacts/milestone_7_decision/epicurus_v05")
FROZEN = Path("configs/frozen/epicurus_v0_5_dev.json")
V03_RESULT = Path("artifacts/milestone_7_decision/epicurus_v03/DEV_RESULT.json")
V04_RESULT = Path("artifacts/milestone_7_decision/epicurus_v04/DEV_RESULT.json")


# ==================================================================================================
# leakage assertions (label-free) — mirror v0.4 §10.7
# ==================================================================================================
def assert_leakage_safe(frame: pd.DataFrame) -> dict:
    dev = frame[~frame["quarantined"]]
    pf = dev.groupby("patient_id")["fold"].nunique()
    assert (pf == 1).all(), "a patient spans >1 fold"
    pep_folds = dev.groupby("canonical_peptide")["fold"].nunique()
    assert (pep_folds == 1).all(), "a non-quarantined peptide crosses folds"
    with guard_no_test_io():
        from event_b.nci_crosswalk import build_train_crosswalk
        train_pats = set("gartner:" + build_train_crosswalk().instances["patient_id"].astype(str))
    g_dev = set(dev.loc[dev["source"] == "gartner", "patient_id"])
    assert g_dev <= train_pats, "a Gartner dev patient is not in the TRAIN crosswalk"
    return {"patient_fold_functional": True, "peptide_fold_functional": True,
            "gartner_dev_patients": len(g_dev), "gartner_test_patients_present": 0}


def assert_io_guard_blocks_test() -> dict:
    """Positively confirm the guard raises on both registered Gartner TEST paths (never merely filtered)."""
    prov = json.loads(v5.V05_PROVENANCE.read_text())
    blocked = []
    with guard_no_test_io():
        for p in prov["forbidden_gartner_test_paths"]:
            try:
                open(p)
                raise AssertionError(f"TEST path was NOT blocked: {p}")
            except RuntimeError:
                blocked.append(p)
    return {"forbidden_paths_blocked": blocked}


# ==================================================================================================
# frozen comparator reconstruction (§2.1) + real, fail-fast reproduction verification (§9.7)
# ==================================================================================================
def reconstruct_and_verify(frame: pd.DataFrame) -> dict:
    v03 = json.loads(V03_RESULT.read_text())
    v04 = json.loads(V04_RESULT.read_text())
    oofP = v5.reconstruct_P(frame, v03)
    oofA = v5.reconstruct_A(frame, v03)
    oofF = v5.reconstruct_F(frame, v04)
    repP = v5.verify_convex_reconstruction(oofP, v03["ladder"]["rung3_MIL"]["folds"],
                                           v03["ladder"]["rung3_MIL"]["eval"]["overall_hits_model"])
    repA = v5.verify_convex_reconstruction(oofA, v03["ladder"]["rung2_additive"]["folds"],
                                           v03["ladder"]["rung2_additive"]["eval"]["overall_hits_model"])
    repF = v5.verify_f_reconstruction(oofF, v04["members"]["F_feature_tower"]["eval"]["overall_hits_model"])
    return {"oofP": oofP, "oofA": oofA, "oofF": oofF,
            "reproduction": {"P": repP, "A": repA, "F": repF}}


# ==================================================================================================
# paired mechanism contrasts (identical patients) — descriptive; the ONLY gate is R vs genuine PRIME
# ==================================================================================================
def mechanism_contrast(a, b, an: str, bn: str) -> dict:
    m = (a.metrics[["source", "patient_id", "hits"]].rename(columns={"hits": f"hits_{an}"})
         .merge(b.metrics[["source", "patient_id", "hits"]].rename(columns={"hits": f"hits_{bn}"}),
                on=["source", "patient_id"]))
    return paired_bootstrap(m, an, bn)


def _gate(evR: dict) -> dict:
    """Primary gate (§7): candidate = R vs GENUINE PRIME. ACCEPT iff paired-bootstrap CI_lower > 0."""
    vp = evR["vs_prime"]
    # strongest presentation baseline = the presentation-type baseline with the higher overall hits
    strongest = "presentation" if evR["overall_hits_presentation"] >= evR["overall_hits_mix"] else "mix"
    vpr = evR[f"vs_{strongest}"]
    beats = vp["ci_lo"] > 0
    no_reg = (vpr["delta"] >= 0) and (vpr["ci_hi"] >= 0)          # design §5.4: point ≥ 0 AND CI not entirely < 0
    return {"beats_genuine_prime": beats, "strongest_presentation": strongest,
            "no_regression_vs_strongest_presentation": no_reg,
            "verdict": "ACCEPT" if beats else "REJECT",
            "delta_vs_prime": vp, "delta_vs_strongest_presentation": vpr,
            "note": "Development-only; ACCEPT licenses ONE pre-registered look at the SEMI-CONSUMED Gartner "
                    "TEST. No external claim without an untouched cohort."}


def _modal_cfg(oof) -> dict:
    keys = [(f["lam_w"], f["lam_ctx"]) for f in oof.spec["folds"]]
    lam_w, lam_ctx = max(set(keys), key=keys.count)
    return {"lam_w": lam_w, "lam_ctx": (np.inf if lam_ctx is None else lam_ctx)}


# ==================================================================================================
# registered diagnostics (§9) — all reported, NONE used for selection
# ==================================================================================================
def _oof_hits_frame(oof, name: str) -> pd.DataFrame:
    return oof.metrics[["source", "patient_id", "hits"]].rename(columns={"hits": f"hits_{name}"})


def context_alias_reaudit(frame: pd.DataFrame) -> dict:
    """§9.3 transferability guard: re-run the §4 patient-grouped, source-balanced source-classification alias
    audit at fit time on the approved contexts, reusing the FROZEN pre-fit audit code verbatim."""
    sys.path.insert(0, "scripts")
    import epicurus_v05_context_audit as audit  # frozen pre-fit audit (hash-pinned)
    ev = v5.add_approved_contexts(frame[~frame["quarantined"]].reset_index(drop=True))
    w = audit._sb_pb_weights(ev)
    per_ctx, encoding = {}, {"ctx_pep_len": ["ctx_pep_len"], "ctx_pred_disagree": ["ctx_pred_disagree"],
                            "ctx_hla_locus": ["ctx_locus_B", "ctx_locus_C"]}
    for name, cols in encoding.items():
        per_ctx[name] = {"alias": audit._alias_metrics(ev, cols, w),
                         "within_patient_variation": audit._within_patient_variation(ev, cols[0])}
    block = audit._alias_metrics(ev, list(v5.CTX_COLS), w)
    return {"per_context": per_ctx, "approved_block_alias": block,
            "note": "near-perfect aliasing ⇒ NO transferability claim for that context (descriptive)."}


def context_ablations(frame: pd.DataFrame, base_delta: float) -> dict:
    """§9.3 leave-one-context-out ablation of R: zero each context (⇒ its 4 interaction columns vanish) and
    re-select+refit R, reporting the change in the R vs PRIME delta. Diagnostic only. `base_delta` is the
    already-computed full-R delta vs PRIME (not recomputed here)."""
    out = {"full_delta_vs_prime": base_delta, "leave_one_context_out": {}}
    for ctx in v5.CTX_COLS:
        fr = v5.add_approved_contexts(frame).copy()
        fr[ctx] = 0.0                                            # constant ⇒ centered to 0 ⇒ interactions vanish
        e = evaluate_model(fr, v5.oof_qr(fr, "R"))[0]
        out["leave_one_context_out"][ctx] = {"delta_vs_prime_without": e["vs_prime"]["delta"],
                                              "drop_vs_full": round(base_delta - e["vs_prime"]["delta"], 4)}
    return out


def shuffled_context_control(frame: pd.DataFrame) -> dict:
    """§9.4 negative control: deterministically shuffle each candidate's contexts WITHIN its patient (breaks the
    context↔candidate link, preserves the per-patient marginal), refit R. Comparable lift ⇒ capacity artifact."""
    fr = v5.add_approved_contexts(frame).copy()
    rng = np.random.default_rng(20260711)
    for _, idx in fr.groupby("patient_id").groups.items():
        idx = np.asarray(idx)
        perm = rng.permutation(len(idx))
        for c in v5.CTX_COLS:
            fr.loc[idx, c] = fr.loc[idx, c].to_numpy()[perm]
    e_sh = evaluate_model(fr, v5.oof_qr(fr, "R"))[0]
    e_true = evaluate_model(frame, v5.oof_qr(frame, "R"))[0]
    return {"true_context_delta_vs_prime": e_true["vs_prime"]["delta"],
            "shuffled_context_delta_vs_prime": e_sh["vs_prime"]["delta"],
            "note": "comparable shuffled lift ⇒ gain is capacity, not portable context structure. Diagnostic."}


def loso_transfer(frame: pd.DataFrame, modalR: dict) -> dict:
    """§9.5 leave-one-source-out stress test (DESCRIPTIVE, n=3): train Q and R on the other two sources at the
    modal config, score the held-out source. Not a gate."""
    dev = frame[~frame["quarantined"]].copy()
    modalQ = {"lam_w": modalR["lam_w"]}
    out = {}
    for s in sorted(dev["source"].unique()):
        tr, te = dev[dev["source"] != s], dev[dev["source"] == s]
        R = v5.ContextPairwiseRanker(member="R", **modalR).fit(tr)
        Q = v5.ContextPairwiseRanker(member="Q", **modalQ).fit(tr)
        mR = per_patient_metrics(te, R.raw_score(te))
        mQ = per_patient_metrics(te, Q.raw_score(te))
        mP = per_patient_metrics(te, baseline_score(te, "prime"))
        out[s] = {"n_patients": int(len(mR)),
                  "R_mean_hits": round(float(mR["hits"].mean()), 3) if len(mR) else None,
                  "Q_mean_hits": round(float(mQ["hits"].mean()), 3) if len(mQ) else None,
                  "prime_mean_hits": round(float(mP["hits"].mean()), 3) if len(mP) else None}
    return {"per_held_out_source": out, "note": "descriptive (n=3); no transfer claim."}


def multi_init_stability(frame: pd.DataFrame, modalR: dict) -> dict:
    """§9.6 HARD GATE (convexity/determinism): ≥2 perturbed fixed-seed inits on fold-0 train at the modal R
    config ⇒ Spearman=1.0 and max|Δcoef| ≤ 1e-6 (unique global minimum)."""
    dev = frame[~frame["quarantined"]].copy()
    tr, te = dev[dev["fold"] != 0], dev[dev["fold"] == 0]
    base = v5.ContextPairwiseRanker(member="R", **modalR).fit(tr)
    base_theta = np.concatenate([base.coef_, base.beta_]) if base.beta_ is not None else base.coef_
    dcoef, rho = [], []
    for seed in (11, 23, 37):
        m = v5.ContextPairwiseRanker(member="R", **modalR).fit(
            tr, init=np.random.default_rng(seed).normal(0, 0.05, base_theta.shape[0]))
        theta = np.concatenate([m.coef_, m.beta_]) if m.beta_ is not None else m.coef_
        dcoef.append(float(np.max(np.abs(base_theta - theta))))
        rho.append(float(pd.Series(base.raw_score(te)).corr(pd.Series(m.raw_score(te)), method="spearman")))
    return {"max_abs_coef_delta": round(max(dcoef), 9), "min_spearman": round(min(rho), 9),
            "pass": bool(max(dcoef) <= 1e-6 and min(rho) >= 1 - 1e-9)}


def prime_masking_availability(frame: pd.DataFrame) -> dict:
    """§9.2 PRIME availability/masking rate per source (label-blind), plus the label-blind attrition report."""
    ev = frame[~frame["quarantined"]]
    per_source = {}
    for s, g in ev.groupby("source"):
        per_source[s] = {"n_rows": int(len(g)),
                         "prime_rank_available_rate": round(float(g["prime_rank"].notna().mean()), 4),
                         "prime_mask_rate": round(float(g["prime_masked"].mean()), 4)}
    return {"per_source": per_source, "attrition": attrition_report(frame)}


def diagnostics(frame: pd.DataFrame, oofQ, oofR, oofP, oofA, oofF, evR: dict) -> dict:
    modalR = _modal_cfg(oofR)                                     # deterministic from oofR — computed once
    d = {"selected_per_fold": oofR.spec["folds"]}
    # (1) contrasts (CIs) — R vs PRIME is the gate; the rest are descriptive
    d["contrasts"] = {
        "R_minus_Q": mechanism_contrast(oofR, oofQ, "R", "Q"),                 # isolates context (both new)
        "Q_minus_P": mechanism_contrast(oofQ, oofP, "Q", "P"),                 # objective + exact-witness supv.
        "A_minus_Q": {**mechanism_contrast(oofA, oofQ, "A", "Q"),
                      "interpretation": "DESCRIPTIVE, not causal: A vs Q changes objective form AND Gartner "
                                        "negative aggregation/bag discipline (not a pure objective isolation)."},
        "A_minus_P": {**mechanism_contrast(oofA, oofP, "A", "P"),
                      "interpretation": "DESCRIPTIVE: supervision granularity (exact-witness pointwise vs bag-MIL)."},
        "R_minus_F": mechanism_contrast(oofR, oofF, "R", "F"),                 # portable vs source-name ceiling
    }
    # (2) per-source ext metrics (R vs baselines) + PRIME masking/availability
    d["per_source_ext"] = _per_source_ext(frame, oofR)
    d["prime_masking_availability"] = prime_masking_availability(frame)
    # (3-6) transferability / capacity / transfer / determinism
    d["context_alias_reaudit"] = context_alias_reaudit(frame)
    d["context_ablations"] = context_ablations(frame, evR["vs_prime"]["delta"])
    d["shuffled_context_control"] = shuffled_context_control(frame)
    d["loso_transfer"] = loso_transfer(frame, modalR)
    d["multi_init_stability"] = multi_init_stability(frame, modalR)
    # (9) interpretation flag: is any gain recognition signal or presentation reweighting?
    d["interpretation"] = {
        "note": "context × presentation-feature interactions are presentation REWEIGHTING, not a new recognition "
                "axis. Flag if pred_disagree interactions dominate the effective β.",
        "effective_beta_final": _effective_beta(frame, modalR)}
    return d


def _per_source_ext(frame: pd.DataFrame, oofR) -> dict:
    dev = frame[~frame["quarantined"]].copy()
    by_fold = dict(oofR.models)
    dev["_s"] = np.nan
    for f, g in dev.groupby("fold"):
        dev.loc[g.index, "_s"] = by_fold[int(f)].raw_score(g)
    mm = ext_metrics(dev, dev["_s"].to_numpy())
    for base in ("prime", "presentation"):
        b = ext_metrics(dev, baseline_score(dev, base)).add_suffix(f"_{base}")
        mm = mm.merge(b.rename(columns={f"source_{base}": "source", f"patient_id_{base}": "patient_id"}),
                      on=["source", "patient_id"])
    res = {}
    for src, g in mm.groupby("source"):
        res[src] = {"patients": int(len(g)), "hits": round(float(g["hits"].mean()), 3),
                    "recall": round(float(g["recall"].mean()), 3),
                    "best_pos_rank": round(float(g["best_pos_rank"].mean()), 2),
                    "ndcg": round(float(g["ndcg"].mean()), 3),
                    "hits_prime": round(float(g["hits_prime"].mean()), 3),
                    "hits_presentation": round(float(g["hits_presentation"].mean()), 3)}
    return res


def _effective_beta(frame: pd.DataFrame, modalR: dict) -> dict:
    dev = frame[~frame["quarantined"]]
    m = v5.ContextPairwiseRanker(member="R", **modalR).fit(dev)
    if m.beta_ is None:
        return {"note": "selected R has λ_ctx=∞ (β≡0): no context interactions carry weight."}
    return {n: round(float(b), 4) for n, b in zip(v5.interaction_names(), m.beta_)}


# ==================================================================================================
# main (NOT executed during CHECKPOINT 1 — writes DEV_RESULT/DEV_REPORT and a frozen config)
# ==================================================================================================
def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    prov = v5.verify_provenance()
    print("provenance OK:", prov)
    io_guard = assert_io_guard_blocks_test()
    print("Gartner TEST I/O guard OK:", io_guard)

    frame = v5.add_approved_contexts(assemble_frame())
    _, _, k = load_frozen()
    leak = assert_leakage_safe(frame)
    print("leakage assertions OK:", leak)

    rec = reconstruct_and_verify(frame)
    oofP, oofA, oofF = rec["oofP"], rec["oofA"], rec["oofF"]
    print("frozen comparator reproduction:", rec["reproduction"])

    cache: dict = {}
    oofQ = v5.oof_qr(frame, "Q", cache=cache)
    oofR = v5.oof_qr(frame, "R", cache=cache)           # the GATED candidate
    evR, jR = evaluate_model(frame, oofR)
    evQ, _ = evaluate_model(frame, oofQ)
    evP, _ = evaluate_model(frame, oofP)
    evA, _ = evaluate_model(frame, oofA)
    evF, _ = evaluate_model(frame, oofF)
    gate = _gate(evR)
    print("R gate:", gate["verdict"], evR["vs_prime"])

    diag = diagnostics(frame, oofQ, oofR, oofP, oofA, oofF, evR)

    result = {
        "protocol": str(OUT / "PREREGISTERED_PROTOCOL.md"),
        "provenance": prov, "io_guard": io_guard, "registered_candidate": "R_context_pairwise", "k": k,
        "scored_patients": evR["n_scored_patients"], "leakage_assertions": leak,
        "frozen_comparator_reproduction": rec["reproduction"],
        "members": {
            "P_pooled_frozen_v03": {"eval": evP, "folds": oofP.spec["folds"]},
            "A_additive_frozen_v03": {"eval": evA, "folds": oofA.spec["folds"]},
            "F_feature_tower_frozen_v04": {"eval": evF, "folds": oofF.spec["folds"]},
            "Q_shared_pairwise": {"eval": evQ, "folds": oofQ.spec["folds"]},
            "R_context_pairwise": {"eval": evR, "gate": gate, "folds": oofR.spec["folds"]},
        },
        "registered_gate_R": gate, "diagnostics": diag, "verdict": gate["verdict"],
        "preservation": "v0.1 remains the frozen model of record; v0.5 is "
                        f"{'ACCEPTED_DEVELOPMENT' if gate['verdict'] == 'ACCEPT' else 'REJECTED_DEVELOPMENT'}. "
                        "Gartner TEST NOT opened; no external claim.",
    }
    (OUT / "DEV_RESULT.json").write_text(json.dumps(result, indent=2, default=str))

    final = v5.ContextPairwiseRanker(member="R", **_modal_cfg(oofR)).fit(frame[~frame["quarantined"]])
    FROZEN.write_text(json.dumps({
        "name": "epicurus_v0_5_dev", "kind": "context_conditioned_pairwise",
        "status": "ACCEPTED_DEVELOPMENT" if gate["verdict"] == "ACCEPT" else "REJECTED_DEVELOPMENT",
        "supersedes_frozen": False, "features": FEATURES, "contexts": v5.CTX_COLS,
        "model": final.to_dict(), "registered_gate": gate,
        "protocol": str(OUT / "PREREGISTERED_PROTOCOL.md"),
        "note": "Development-only deployable challenger; not generalizable without an untouched cohort."},
        indent=2, default=str))

    _write_report(result)
    print(json.dumps({"verdict": gate["verdict"], "R_vs_prime": evR["vs_prime"],
                      "contrasts": {k2: v2 for k2, v2 in diag["contrasts"].items()}}, indent=2, default=str))
    print(f"\nwrote {OUT/'DEV_RESULT.json'}, {OUT/'DEV_REPORT.md'}, {FROZEN}")


def _char(d: dict) -> str:
    if d["ci_lo"] > 0:
        return "significantly BETTER (CI>0)"
    if d["ci_hi"] < 0:
        return "significantly WORSE (CI<0)"
    return "statistically TIED (CI spans 0)"


def _write_report(r: dict) -> None:
    g = r["registered_gate_R"]
    evR = r["members"]["R_context_pairwise"]["eval"]
    evQ = r["members"]["Q_shared_pairwise"]["eval"]
    evP = r["members"]["P_pooled_frozen_v03"]["eval"]
    evF = r["members"]["F_feature_tower_frozen_v04"]["eval"]
    c = r["diagnostics"]["contrasts"]
    L = [f"# Epicurus v0.5 DEVELOPMENT — context-conditioned pairwise challenger — verdict: **{r['verdict']}**\n",
         "Preregistered: `PREREGISTERED_PROTOCOL.md`. DEVELOPMENT ONLY — Gartner TEST not opened; no external "
         "claim. The ONLY gate is R vs GENUINE PRIME (raw unmasked prime_rank).\n",
         f"Provenance verified ({r['provenance']['n_inputs_verified']} inputs, git "
         f"{r['provenance']['git_head'][:10]}). {r['scored_patients']} scored patients "
         "(source-balanced, patient-paired bootstrap).\n",
         "## Frozen comparator reproduction (§2.1)\n",
         f"- P (convex): {r['frozen_comparator_reproduction']['P']}",
         f"- A (convex): {r['frozen_comparator_reproduction']['A']}",
         f"- F (nonconvex, honest tolerance): {r['frozen_comparator_reproduction']['F']}\n",
         "## Registered gate (candidate = R)\n",
         f"- vs **genuine PRIME**: Δhits@20 = {g['delta_vs_prime']['delta']} "
         f"CI[{g['delta_vs_prime']['ci_lo']}, {g['delta_vs_prime']['ci_hi']}] → beats PRIME: "
         f"**{g['beats_genuine_prime']}** ({_char(g['delta_vs_prime'])})",
         f"- vs **strongest presentation** ({g['strongest_presentation']}): Δ = "
         f"{g['delta_vs_strongest_presentation']['delta']} "
         f"CI[{g['delta_vs_strongest_presentation']['ci_lo']}, {g['delta_vs_strongest_presentation']['ci_hi']}] "
         f"→ no regression: **{g['no_regression_vs_strongest_presentation']}**\n",
         "## Members (OOF hits@20)\n",
         "| member | overall hits | Δ vs PRIME (CI) |", "|---|--:|---|",
         f"| P — pooled (frozen v0.3) | {evP['overall_hits_model']} | {evP['vs_prime']['delta']} "
         f"[{evP['vs_prime']['ci_lo']}, {evP['vs_prime']['ci_hi']}] |",
         f"| Q — shared pairwise | {evQ['overall_hits_model']} | {evQ['vs_prime']['delta']} "
         f"[{evQ['vs_prime']['ci_lo']}, {evQ['vs_prime']['ci_hi']}] |",
         f"| **R — context pairwise** | {evR['overall_hits_model']} | {evR['vs_prime']['delta']} "
         f"[{evR['vs_prime']['ci_lo']}, {evR['vs_prime']['ci_hi']}] |",
         f"| F — source-name tower (frozen v0.4) | {evF['overall_hits_model']} | {evF['vs_prime']['delta']} "
         f"[{evF['vs_prime']['ci_lo']}, {evF['vs_prime']['ci_hi']}] |\n",
         "## Contrasts (paired; descriptive — the only gate is R vs PRIME)\n",
         f"- **R − Q** (isolates context): {c['R_minus_Q']['delta']} "
         f"CI[{c['R_minus_Q']['ci_lo']}, {c['R_minus_Q']['ci_hi']}] → {_char(c['R_minus_Q'])}",
         f"- **Q − P** (objective + exact-witness supervision, not a pure isolation): {c['Q_minus_P']['delta']} "
         f"CI[{c['Q_minus_P']['ci_lo']}, {c['Q_minus_P']['ci_hi']}] → {_char(c['Q_minus_P'])}",
         f"- **A − Q** (DESCRIPTIVE — objective form AND Gartner bag discipline, not a pure objective isolation): "
         f"{c['A_minus_Q']['delta']} CI[{c['A_minus_Q']['ci_lo']}, {c['A_minus_Q']['ci_hi']}]",
         f"- **A − P** (DESCRIPTIVE — supervision granularity): {c['A_minus_P']['delta']} "
         f"CI[{c['A_minus_P']['ci_lo']}, {c['A_minus_P']['ci_hi']}]",
         f"- **R − F** (portable context vs source-name ceiling): {c['R_minus_F']['delta']} "
         f"CI[{c['R_minus_F']['ci_lo']}, {c['R_minus_F']['ci_hi']}] → {_char(c['R_minus_F'])}\n",
         f"## Verdict\n\n**{r['verdict']}.** {r['preservation']}\n"]
    (OUT / "DEV_REPORT.md").write_text("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
