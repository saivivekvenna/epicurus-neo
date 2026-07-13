"""Within-patient DECISION-PROBLEM benchmark, corrected for label ascertainment.

The north-star question is within-patient top-k on a patient's candidate universe. But a complete
candidate universe is NOT a completely assayed universe, so the harness now separates two claims
that the original M7 code conflated:

    MODE A - full-universe POSITIVE-UNLABELED retrieval (allowed, ascertainment-limited):
        rank the whole candidate list and ask how many experimentally-confirmed positives land in
        the patient's top-k. Non-positives are a mix of tested-negative and UNTESTED, so this is a
        positive-retrieval metric, NOT discrimination against confirmed negatives.

    MODE B - tested-subset discrimination (genuine, non-circular):
        AUROC/AP and supervised training on POSITIVE vs TESTED_NEGATIVE ONLY. FAIL-CLOSED: if a
        corpus has no tested negatives (e.g. the Müller min file), this is BLOCKED, not fabricated.

    MODE C - selection-bias sensitivity diagnostics:
        how tested vs untested candidates differ on presentation/expression, and where positives
        fall in the presentation ranking. Does not pretend to correct unknown selection.

Leakage control throughout: GroupKFold by patient_id; fold-local impute+scale; ln_NumTested is a
per-patient covariate and is never a feature. NOT a PRIME head-to-head (no genuine PRIME scores;
never call a MixMHCpred/NetMHCpan-EL baseline "PRIME").
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import GroupKFold

from benchmark.scorecard import scorecard
from event_b.neoranking_corpus import (
    PRESENTATION_BASELINE,
    PRESENTATION_FEATURES,
    RECOGNITION_FEATURES,
    load_neoranking_nci,
    oriented_feature_matrix,
    shared_peptide_diagnostics,
)
from event_b.zhao_features import _seq_features

SEQ_FEATURE_KEYS = ["gravy", "nonanchor_gravy", "aromatic_frac", "net_charge", "n_positive", "n_negative"]

BENCHMARK_VERSION = "decision-benchmark-2.0.0-ascertainment-corrected"
N_SPLITS = 5
K_POLICY = (5, 10, 20, 50)
PRIMARY_K = 20
SEED = 0
TESTED_LABELS = ("POSITIVE", "TESTED_NEGATIVE")


def _seq_matrix(peptides) -> np.ndarray:
    rows = [[_seq_features(p).get(key, 0.0) for key in SEQ_FEATURE_KEYS] for p in peptides]
    return np.asarray(rows, dtype=float) if rows else np.empty((0, len(SEQ_FEATURE_KEYS)))


def _zscore(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    sd = v.std()
    return (v - v.mean()) / sd if sd > 0 else v * 0.0


def _oof_scores(X: np.ndarray, y: np.ndarray, groups: np.ndarray) -> np.ndarray:
    """Patient-grouped out-of-fold decision scores; impute (train median) + scale (train stats)
    fit strictly inside each training fold, so no eval-fold statistic leaks."""
    oof = np.zeros(len(y), dtype=float)
    n = min(N_SPLITS, len(np.unique(groups)))
    if n < 2:
        return oof
    for tr, ev in GroupKFold(n_splits=n).split(X, y, groups):
        Xtr, Xev = X[tr].copy(), X[ev].copy()
        medians = np.nanmedian(Xtr, axis=0)
        medians = np.where(np.isnan(medians), 0.0, medians)
        Xtr = np.where(np.isnan(Xtr), medians, Xtr)
        Xev = np.where(np.isnan(Xev), medians, Xev)
        mu, sigma = Xtr.mean(axis=0), Xtr.std(axis=0)
        sigma = np.where(sigma == 0, 1.0, sigma)
        Xtr = (Xtr - mu) / sigma
        Xev = (Xev - mu) / sigma
        if len(np.unique(y[tr])) < 2:
            continue
        model = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced")
        model.fit(Xtr, y[tr])
        oof[ev] = model.decision_function(Xev)
    return oof


def _supervised_oof(X, y_full, groups, tested_mask) -> tuple[np.ndarray, bool]:
    """Patient-grouped OOF that NEVER trains on UNTESTED rows.

    For each outer patient fold, the model is fit ONLY on the training patients' tested rows
    (POSITIVE vs TESTED_NEGATIVE) but predicts EVERY evaluation-patient candidate (so full-universe
    top-k still covers untested candidates at eval time). Impute/scale are fit on the tested train
    rows only. Returns (oof_scores, trained_any); trained_any is False when no fold had both tested
    classes, i.e. the learned arm is not computable and must be BLOCKED (never fabricated).
    """
    X = np.asarray(X, dtype=float)
    y_full = np.asarray(y_full, dtype=int)
    tested_mask = np.asarray(tested_mask, dtype=bool)
    oof = np.full(len(y_full), np.nan, dtype=float)
    n = min(N_SPLITS, len(np.unique(groups)))
    if n < 2:
        return oof, False
    trained_any = False
    for tr, ev in GroupKFold(n_splits=n).split(X, y_full, groups):
        tr_tested = tr[tested_mask[tr]]
        ytr = y_full[tr_tested]
        if len(np.unique(ytr)) < 2:
            continue  # cannot fit a genuine pos-vs-tested-neg model in this fold
        Xtr, Xev = X[tr_tested].copy(), X[ev].copy()
        medians = np.nanmedian(Xtr, axis=0)
        medians = np.where(np.isnan(medians), 0.0, medians)
        Xtr = np.where(np.isnan(Xtr), medians, Xtr)
        Xev = np.where(np.isnan(Xev), medians, Xev)
        mu, sigma = Xtr.mean(axis=0), Xtr.std(axis=0)
        sigma = np.where(sigma == 0, 1.0, sigma)
        model = LogisticRegression(max_iter=2000, C=1.0, class_weight="balanced")
        model.fit((Xtr - mu) / sigma, ytr)
        oof[ev] = model.decision_function((Xev - mu) / sigma)
        trained_any = True
    return oof, trained_any


def _card(frame: pd.DataFrame, arm: str, k: int, baseline: str = "presentation_only") -> dict:
    card = scorecard(
        frame, score_col=arm, baseline_col=baseline,
        group_col="patient_id", k=k, label_col="label",
        ascending=False, baseline_ascending=False,
    )
    entry = card[f"hits@{k}"]
    return {
        "hits_delta": entry["delta_vs_baseline"],
        "hits_value": entry["value"],
        "delta_ci": entry["delta_ci"],
        "p_better": entry["p_better"],
        "capture_delta": card["capture_fraction"]["delta_vs_baseline"],
        "verdict": card["verdict"],
    }


def _top_recall(frame: pd.DataFrame, arm: str, k: int) -> tuple[int, int]:
    got = tot = 0
    for _, g in frame.groupby("patient_id"):
        ranked = g.sort_values(arm, ascending=False, kind="mergesort")
        got += int(ranked["y"].to_numpy()[:k].sum())
        tot += int(g["y"].sum())
    return got, tot


def fail_closed_auroc(y: np.ndarray, score: np.ndarray, *, tested_negative_available: bool) -> dict:
    """AUROC/AP that REFUSES to score unless genuine tested negatives are present."""
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    if not tested_negative_available:
        return {"status": "BLOCKED",
                "reason": "no tested negatives; a discrimination metric would be positives-vs-"
                          "unlabeled (ascertainment-biased) and is not fabricated."}
    finite = np.isfinite(score)
    y, score = y[finite], score[finite]
    if len(np.unique(y)) < 2:
        return {"status": "BLOCKED", "reason": "only one label class in the tested subset."}
    return {"status": "executed", "auroc": float(roc_auc_score(y, score)),
            "ap": float(average_precision_score(y, score)), "n": int(len(y)),
            "positives": int(y.sum()), "negatives": int((y == 0).sum()),
            "dropped_nonfinite": int((~finite).sum())}


def _tested_discrimination(frame, arms, *, tested_negative_available) -> dict:
    """MODE B: AUROC/AP on tested-pos vs tested-neg ONLY, over precomputed arm columns (each already
    trained tested-only via _supervised_oof). Fail-closed when there are no tested negatives."""
    tmask = frame["label"].isin(TESTED_LABELS).to_numpy()
    tf = frame[tmask]
    n_pos = int((tf["label"] == "POSITIVE").sum())
    n_neg = int((tf["label"] == "TESTED_NEGATIVE").sum())
    if not tested_negative_available or n_neg == 0 or n_pos == 0:
        return {"status": "BLOCKED",
                "reason": "no genuine tested negatives in this corpus (see label audit); "
                          "tested-pos-vs-tested-neg AUROC/AP and supervised-negative training are "
                          "BLOCKED and not fabricated.",
                "tested_positives": n_pos, "tested_negatives": n_neg}
    yt = tf["y"].to_numpy(int)
    per_arm = {
        arm: fail_closed_auroc(yt, tf[arm].to_numpy(float), tested_negative_available=True)
        for arm in ["presentation_only", *arms]
    }
    return {"status": "executed", "tested_positives": n_pos, "tested_negatives": n_neg,
            "per_arm_auroc_ap": per_arm,
            "note": "AUROC/AP on genuinely tested positives vs tested negatives only; every learned "
                    "arm was trained tested-only (UNTESTED rows never entered any supervised fit), "
                    "patient-held-out."}


def _incremental_over_presentation_ensemble(frame, tested_negative_available) -> dict:
    """The decisive ablation: does adding expression/VAF beat a learned PRESENTATION ENSEMBLE?

    Compares full_model and recognition_only against presentation_ensemble (five presentation
    predictors, learned tested-only) — NOT against a single raw EL score — on both tested-only
    AUROC/AP and paired within-patient top-k. Fail-closed when the ensemble or tested negatives are
    unavailable. This is the only comparison that supports an 'orthogonal recognition adds signal'
    claim; full_model beating one raw predictor could be mere presentation ensembling.
    """
    if not tested_negative_available or "presentation_ensemble" not in frame.columns:
        return {"status": "BLOCKED",
                "reason": "needs a learned presentation_ensemble and genuine tested negatives; "
                          "not computable here (see label audit)."}
    tmask = frame["label"].isin(TESTED_LABELS).to_numpy()
    tf = frame[tmask]
    yt = tf["y"].to_numpy(int)
    base = fail_closed_auroc(yt, tf["presentation_ensemble"].to_numpy(float), tested_negative_available=True)
    out = {"status": "executed",
           "baseline_arm": "presentation_ensemble",
           "presentation_ensemble_auroc": base.get("auroc"),
           "presentation_ensemble_ap": base.get("ap"),
           "challengers": {}}
    for arm in ["full_model", "recognition_only", "recognition_residual"]:
        if arm not in frame.columns:
            continue
        disc = fail_closed_auroc(yt, tf[arm].to_numpy(float), tested_negative_available=True)
        topk = {f"hits@{k}": _card(frame, arm, k, baseline="presentation_ensemble") for k in K_POLICY}
        primary = topk[f"hits@{PRIMARY_K}"]
        out["challengers"][arm] = {
            "auroc_tested": disc.get("auroc"),
            "delta_auroc_vs_ensemble": (disc.get("auroc") - base.get("auroc"))
            if disc.get("auroc") is not None and base.get("auroc") is not None else None,
            "ap_tested": disc.get("ap"),
            f"hits@{PRIMARY_K}_delta_vs_ensemble": primary["hits_delta"],
            f"hits@{PRIMARY_K}_delta_ci": primary["delta_ci"],
            "p_better_vs_ensemble": primary["p_better"],
            "verdict_vs_ensemble": primary["verdict"],
            "topk_vs_ensemble": {kk: topk[kk]["verdict"] for kk in topk},
        }
    # Explicit lift decomposition so the interpretation is not left to prose: how much of
    # full_model's top-k lift over the single raw score is presentation ENSEMBLING (ensemble vs raw
    # EL) versus the expression/VAF ADDITION (full_model vs ensemble).
    raw = _top_recall(frame, "presentation_only", PRIMARY_K)
    ens = _top_recall(frame, "presentation_ensemble", PRIMARY_K)
    full = _top_recall(frame, "full_model", PRIMARY_K)
    out["topk_decomposition"] = {
        f"raw_el_top{PRIMARY_K}": raw,
        f"presentation_ensemble_top{PRIMARY_K}": ens,
        f"full_model_top{PRIMARY_K}": full,
        "ensembling_gain_over_raw_el": ens[0] - raw[0],
        "expression_vaf_gain_over_ensemble": full[0] - ens[0],
    }

    fm = out["challengers"].get("full_model", {})
    verdict = fm.get("verdict_vs_ensemble")
    delta = fm.get(f"hits@{PRIMARY_K}_delta_vs_ensemble") or 0.0
    dauroc = fm.get("delta_auroc_vs_ensemble")
    dauroc_txt = f"{dauroc:+.3f}" if dauroc is not None else "n/a"
    if verdict == "ACCEPT":
        out["conclusion"] = (
            "full_model beats the LEARNED presentation ensemble at within-patient top-k: expression/"
            "VAF add orthogonal signal that survives a fair (ensembled) presentation baseline.")
    elif verdict == "REJECT" or delta < 0:
        out["conclusion"] = (
            f"full_model does not beat the learned presentation ensemble (hits@{PRIMARY_K} delta "
            f"{delta:+.3f}); no orthogonal recognition gain over a fair baseline here.")
    else:
        out["conclusion"] = (
            f"full_model shows a LARGE but statistically UNESTABLISHED incremental over the learned "
            f"presentation ensemble (hits@{PRIMARY_K} delta {delta:+.3f} [{full[0]}/{full[1]} vs "
            f"{ens[0]}/{ens[1]}], ΔAUROC {dauroc_txt}); the paired within-patient CI spans zero. "
            f"This is an underpowered/inconclusive expression+VAF signal (few tested positives), NOT "
            f"evidence of no effect, and NOT attributable to presentation ensembling (the ensemble "
            f"alone moves top-{PRIMARY_K} by only {ens[0] - raw[0]:+d} over raw EL, vs "
            f"{full[0] - ens[0]:+d} added by expression/VAF).")
    return out


def _selection_bias_diagnostics(frame, presentation_col, feature_cols) -> dict:
    """MODE C: how tested vs untested candidates differ, and where positives sit in presentation rank."""
    f = frame.copy()
    f["pres_rank"] = f.groupby("patient_id")[presentation_col].rank(ascending=False, method="first")
    state = {s: f[f["label"] == s] for s in ("POSITIVE", "TESTED_NEGATIVE", "UNTESTED")}
    medians = {
        s: {c: (float(g[c].median()) if c in g and len(g) else None) for c in feature_cols}
        for s, g in state.items() if len(g)
    }
    pos = state["POSITIVE"]
    bands = {"rank_1_5": (1, 5), "rank_6_20": (6, 20), "rank_21_100": (21, 100),
             "rank_over_100": (101, 10**12)}
    return {
        "label_state_counts": frame["label"].value_counts().to_dict(),
        "feature_medians_by_state": medians,
        "positive_presentation_rank_bands": {
            n: int(pos["pres_rank"].between(lo, hi).sum()) for n, (lo, hi) in bands.items()
        },
        "note": "If UNTESTED candidates are systematically weaker binders / lower expression than "
                "tested ones, the positives were presentation-selected for assay, so full-universe "
                "AUROC is ascertainment-inflated. This diagnostic does NOT correct unknown selection.",
    }


def run_benchmark_core(
    frame: pd.DataFrame,
    *,
    presentation_baseline: dict,
    presentation_features: dict,
    recognition_features: dict,
    baseline_label: str,
    corpus_meta: dict,
    reconciliation: dict,
    honesty_note: str,
    tested_negative_available: bool,
) -> tuple[dict, pd.DataFrame]:
    groups = frame["patient_id"].to_numpy()
    y = frame["y"].to_numpy(dtype=int)

    X_pres = oriented_feature_matrix(frame, presentation_features)
    X_reco = oriented_feature_matrix(frame, recognition_features)
    X_seq = _seq_matrix(frame["mutant_peptide"].tolist())
    X_reco_seq = np.hstack([X_reco, X_seq])
    X_all = np.hstack([X_pres, X_reco])

    frame = frame.copy()
    frame["presentation_only"] = oriented_feature_matrix(frame, presentation_baseline)[:, 0]

    # Learned arms are trained ONLY on tested rows (pos vs tested-neg). If a corpus has no tested
    # negatives (e.g. Müller min), an ordinary supervised model would have to treat UNTESTED as
    # negative — forbidden — so those arms are BLOCKED (never fabricated), leaving only the
    # training-free presentation ranker for Mode A retrieval.
    tested_mask = frame["label"].isin(TESTED_LABELS).to_numpy()
    both_tested_classes = bool((y[tested_mask] == 1).any() and (y[tested_mask] == 0).any())
    learned_available = tested_negative_available and both_tested_classes
    blocked_reason = None
    arms: list[str] = []
    if learned_available:
        # presentation_ensemble = learned combination of the presentation predictors ALONE. It is
        # the fair baseline for the recognition claim: full_model must beat THIS, not a single raw
        # EL score, or its gain is presentation ensembling rather than orthogonal recognition.
        specs = {
            "presentation_ensemble": X_pres,
            "recognition_only": X_reco,
            "recognition_seq_only": X_reco_seq,
            "full_model": X_all,
        }
        for name, matrix in specs.items():
            scores, ok = _supervised_oof(matrix, y, groups, tested_mask)
            if ok:
                frame[name] = scores
                arms.append(name)
        if "recognition_only" in frame:
            frame["recognition_residual"] = _zscore(frame["presentation_only"].to_numpy()) + _zscore(
                np.nan_to_num(frame["recognition_only"].to_numpy(), nan=0.0)
            )
            arms.insert(1, "recognition_residual")
    else:
        blocked_reason = (
            "learned recognition/full-model arms are BLOCKED: no genuine tested negatives, so an "
            "ordinary supervised model would have to treat UNTESTED rows as negatives. Only the "
            "training-free presentation ranker is reported for Mode A here (use a dedicated PU "
            "method, or the Gartner/multimer tested corpora, for learned arms)."
        )

    n_untested = int((~tested_mask).sum())
    fully_tested = n_untested == 0
    semantics = (
        "full TESTED candidate universe (zero UNTESTED rows); every non-positive is a confirmed "
        "tested-negative, so this top-k is over a fully-assayed universe."
        if fully_tested
        else "positive-UNLABELED retrieval over the full candidate universe; non-positives mix "
        "tested-negative and UNTESTED, so this is positive retrieval, NOT discrimination vs "
        "confirmed negatives."
    )
    scored_arms = ["presentation_only", *arms]
    # Rank NaN (rows an arm could not score) last so they never occupy a top-k slot.
    for arm in arms:
        frame[arm] = frame[arm].fillna(frame[arm].min() - 1.0 if frame[arm].notna().any() else -1e9)

    mode_a = {
        "semantics": semantics,
        "fully_tested_universe": fully_tested,
        "learned_arms_blocked": blocked_reason,
        "raw_topk_recall": {
            arm: {f"top{k}": _top_recall(frame, arm, k) for k in K_POLICY} for arm in scored_arms
        },
        "scorecards_vs_baseline": {
            k: {arm: _card(frame, arm, k) for arm in arms} for k in K_POLICY
        },
    }

    # MODE B: tested-subset discrimination (fail-closed).
    mode_b = _tested_discrimination(
        frame, arms, tested_negative_available=tested_negative_available
    )

    # DECISIVE ABLATION: recognition vs a learned presentation ENSEMBLE (not one raw score).
    incremental = _incremental_over_presentation_ensemble(frame, tested_negative_available)

    # MODE C: selection-bias sensitivity.
    feature_cols = [c for c in [*presentation_features, *recognition_features] if c in frame.columns]
    mode_c = _selection_bias_diagnostics(frame, "presentation_only", feature_cols)

    report = {
        "benchmark_version": BENCHMARK_VERSION,
        "corpus": corpus_meta,
        "reconciliation": reconciliation,
        "tested_negative_available": tested_negative_available,
        "baseline": baseline_label,
        "primary_k": PRIMARY_K,
        "k_policy": list(K_POLICY),
        "presentation_features": presentation_features,
        "recognition_features": recognition_features,
        "shared_peptide_diagnostics": shared_peptide_diagnostics(frame),
        "mode_a_positive_unlabeled_retrieval": mode_a,
        "mode_b_tested_discrimination": mode_b,
        "incremental_recognition_over_presentation_ensemble": incremental,
        "mode_c_selection_bias": mode_c,
        "honesty_note": honesty_note,
    }
    return report, frame


def run_decision_benchmark(path=None) -> tuple[dict, pd.DataFrame]:
    corpus = load_neoranking_nci() if path is None else load_neoranking_nci(path)
    return run_benchmark_core(
        corpus.frame,
        presentation_baseline={PRESENTATION_BASELINE: +1},
        presentation_features=PRESENTATION_FEATURES,
        recognition_features=RECOGNITION_FEATURES,
        baseline_label="presentation_only (Score_EL, NetMHCpan eluted-ligand — NOT PRIME)",
        corpus_meta=corpus.provenance,
        reconciliation=corpus.reconciliation,
        tested_negative_available=False,
        honesty_note=(
            "Müller NCI min file: VALIDATED=0 is UNTESTED (no per-peptide assay indicator). Mode A "
            "(positive-unlabeled top-k) is valid but ascertainment-limited; Mode B (tested AUROC/AP "
            "+ supervised) is BLOCKED here. Not a PRIME head-to-head; Score_EL is NetMHCpan EL, not "
            "PRIME."
        ),
    )


def run_multimer_benchmark() -> tuple[dict, pd.DataFrame]:
    from event_b.cd8_multimer_corpus import (
        MULTIMER_PRESENTATION_BASELINE,
        MULTIMER_PRESENTATION_FEATURES,
        MULTIMER_RECOGNITION_FEATURES,
        load_cd8_multimer,
    )

    corpus = load_cd8_multimer()
    return run_benchmark_core(
        corpus.frame,
        presentation_baseline=MULTIMER_PRESENTATION_BASELINE,
        presentation_features=MULTIMER_PRESENTATION_FEATURES,
        recognition_features=MULTIMER_RECOGNITION_FEATURES,
        baseline_label="presentation_only (EL %Rank, oriented — NOT PRIME)",
        corpus_meta=corpus.provenance,
        reconciliation=corpus.reconciliation,
        tested_negative_available=corpus.reconciliation.get("tested_negative_available", False),
        honesty_note=(
            "INDEPENDENT pMHC-multimer cohort. All candidates multimer-tested (mmc1 reconciliation) "
            "=> NO = genuine TESTED_NEGATIVE, so Mode B discrimination is valid here. RF classifier "
            "score (the paper's own trained model) is excluded. Not a PRIME head-to-head."
        ),
    )


def run_gartner_benchmark() -> tuple[dict, pd.DataFrame]:
    """The corrected NCI substrate WITH genuine tested negatives (Screening Status three-state)."""
    from event_b.gartner_nci_corpus import (
        PRESENTATION_BASELINE as GB,
        PRESENTATION_FEATURES as GF,
        RECOGNITION_FEATURES as GR,
        load_gartner_nci,
    )

    corpus = load_gartner_nci()
    return run_benchmark_core(
        corpus.frame,
        presentation_baseline=GB,
        presentation_features=GF,
        recognition_features=GR,
        baseline_label="presentation_only (NetMHCpan EL %rank, oriented — NOT PRIME)",
        corpus_meta=corpus.provenance,
        reconciliation=corpus.reconciliation,
        tested_negative_available=True,
        honesty_note=(
            "Gartner NCI (same cohort as Müller) WITH genuine three-state ascertainment. Mode A is "
            "positive-unlabeled over the full universe (incl. UNTESTED); Mode B AUROC/AP + supervised "
            "are on genuinely tested pos-vs-neg only. This is the honest presentation-discrimination "
            "number the Müller min file could not support. Not a PRIME head-to-head."
        ),
    )


def render_decision_markdown(report: dict) -> str:
    rec = report["reconciliation"]
    a = report["mode_a_positive_unlabeled_retrieval"]
    b = report["mode_b_tested_discrimination"]
    c = report["mode_c_selection_bias"]
    k = report["primary_k"]
    lines = [
        "# Within-patient decision-problem benchmark (ascertainment-corrected)",
        "",
        f"**Corpus:** {report['corpus']['citation']}",
        f"- observed {rec['observed']}; label states {rec.get('label_state_counts', {})}; "
        f"tested-negative available: **{report['tested_negative_available']}**.",
        f"- Baseline: {report['baseline']}. Leakage: GroupKFold by patient.",
        "",
        f"> {report['honesty_note']}",
        "",
        "## Mode A — full-universe top-k retrieval",
        "",
        f"_{a['semantics']}_",
        "",
        (f"**Learned arms BLOCKED:** {a['learned_arms_blocked']}" if a.get("learned_arms_blocked")
         else "_Learned arms trained tested-only (UNTESTED never in any fit); predict full universe._"),
        "",
        "| Arm | " + " | ".join(f"top{kk}" for kk in report["k_policy"]) + " |",
        "|---|" + "---|" * len(report["k_policy"]),
    ]
    for arm, byk in a["raw_topk_recall"].items():
        cells = " | ".join(f"{byk[f'top{kk}'][0]}/{byk[f'top{kk}'][1]}" for kk in report["k_policy"])
        lines.append(f"| {arm} | {cells} |")
    lines += ["", f"### hits@{k} vs baseline (paired patient bootstrap)", "",
              "| Arm | Δ | 95% CI | P(better) | verdict |", "|---|---:|---|---:|---|"]
    for arm, e in a["scorecards_vs_baseline"][k].items():
        ci = e["delta_ci"]
        lines.append(f"| {arm} | {e['hits_delta']:+.3f} | [{ci[0]:+.3f}, {ci[1]:+.3f}] | "
                     f"{e['p_better']:.3f} | {e['verdict']} |")
    lines += ["", "## Mode B — tested-subset discrimination (genuine, fail-closed)", ""]
    if b.get("status") == "BLOCKED":
        lines += [f"**BLOCKED.** {b['reason']} (tested pos={b.get('tested_positives')}, "
                  f"tested neg={b.get('tested_negatives')})."]
    else:
        lines += [f"Tested positives={b['tested_positives']}, tested negatives={b['tested_negatives']}.",
                  "", "| Arm | AUROC (tested) | AP (tested) |", "|---|---:|---:|"]
        for arm, d in b["per_arm_auroc_ap"].items():
            if d.get("status") == "executed":
                lines.append(f"| {arm} | {d['auroc']:.4f} | {d['ap']:.4f} |")
    inc = report.get("incremental_recognition_over_presentation_ensemble", {})
    lines += ["", "## Decisive ablation — recognition vs a LEARNED presentation ensemble", ""]
    if inc.get("status") == "BLOCKED":
        lines.append(f"**BLOCKED.** {inc['reason']}")
    elif inc.get("status") == "executed":
        lines.append(f"Baseline = presentation_ensemble (5 presentation predictors, learned "
                     f"tested-only): AUROC {inc['presentation_ensemble_auroc']:.4f}.")
        lines += ["", "| Challenger | AUROC (tested) | ΔAUROC vs ensemble | "
                  f"Δhits@{k} vs ensemble | P(better) | verdict |", "|---|---:|---:|---:|---:|---|"]
        for arm, d in inc["challengers"].items():
            da = d["delta_auroc_vs_ensemble"]
            lines.append(
                f"| {arm} | {d['auroc_tested']:.4f} | {da:+.4f} | "
                f"{d[f'hits@{k}_delta_vs_ensemble']:+.3f} | {d['p_better_vs_ensemble']:.3f} | "
                f"{d['verdict_vs_ensemble']} |")
        lines += ["", f"**{inc['conclusion']}**"]
    lines += ["", "## Mode C — selection-bias sensitivity", "",
              f"- Label states: {c['label_state_counts']}",
              f"- Positive presentation-rank bands: {c['positive_presentation_rank_bands']}",
              f"- Feature medians by state: {c['feature_medians_by_state']}",
              f"- {c['note']}", ""]
    return "\n".join(lines)
