"""Epicurus v0.5 — PRE-FIT context-feasibility + source-ALIAS audit (DESIGN/PROVENANCE, not model code).

Read-only, OUTCOME-label-blind, NO Q/R model is fit and NO ranking metric is computed. Its job is to
prove — before the v0.5 preregistration is finalized — that each candidate scoring-time context variable is
  (a) computable from information available while scoring a brand-new patient's candidate pool, from the
      candidate itself (peptide sequence, its HLA allele, leakage-safe already-computed predictor values) and
      NOT from source-specific pool enumeration / denominator construction, and
  (b) NOT a source alias — quantified by an ACTUAL patient-grouped, source-balanced source-classification
      audit (macro balanced-accuracy / macro one-vs-rest AUROC against source-balanced chance), NOT by range
      overlap alone (range overlap is retained only as a descriptive statistic, it is NOT the gate).

Source labels are used here descriptively (to quantify aliasing), never to select or fit a Q/R model. Outcome
labels (`eval_positive`, `bag_label`) are NEVER read.  Runs inside the v0.4 Gartner-TEST I/O guard
(`guard_no_test_io`); Gartner TEST is never opened.

FROZEN ALIAS DECISION RULE (fixed before any Q/R fit): a candidate context is APPROVED for R's interaction
block iff (i) it is candidate-level and inference-computable without source-specific pool enumeration, and
(ii) its STANDALONE patient-grouped source-balanced macro one-vs-rest AUROC <= ALIAS_AUROC_MAX (0.70). Any
context that requires pool enumeration / a source-specific denominator, or is a deterministic transform of
already-modeled features, is rejected on structural grounds regardless of its alias number.

Emits `artifacts/milestone_7_decision/epicurus_v05/CONTEXT_FEASIBILITY_AUDIT.json`.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold

from event_b.epicurus_v04 import assemble_frame, guard_no_test_io  # reuse the guarded frame verbatim

OUT = Path("artifacts/milestone_7_decision/epicurus_v05/CONTEXT_FEASIBILITY_AUDIT.json")
FORBIDDEN_LABEL_COLS = ("eval_positive", "bag_label")  # never used in this audit
ALIAS_AUROC_MAX = 0.70   # FROZEN: standalone macro OVR AUROC above this => rejected as a source alias
CV_SEED = 0              # FROZEN: deterministic StratifiedGroupKFold shuffle seed
DISAGREE_COLS = ["f_prime_pct", "f_mix_pct", "f_el_pct"]  # leakage-safe masked predictor percentiles


# --------------------------------------------------------------------------------------------------
# descriptive statistics (retained, but NOT the alias gate)
# --------------------------------------------------------------------------------------------------
def _overlap(a: np.ndarray, b: np.ndarray) -> float:
    """Fraction of the pooled [min,max] range shared by both groups (descriptive only, NOT the gate)."""
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    lo = max(a.min(), b.min()); hi = min(a.max(), b.max())
    inter = max(0.0, hi - lo)
    union = max(a.max(), b.max()) - min(a.min(), b.min())
    return float(inter / union) if union > 0 else 1.0


def _per_source_stats(df: pd.DataFrame, col: str) -> dict:
    out = {}
    for src, g in df.groupby("source"):
        v = pd.to_numeric(g[col], errors="coerce").to_numpy(float)
        v = v[np.isfinite(v)]
        if len(v) == 0:
            out[src] = {"n": 0}
            continue
        out[src] = {"n": int(len(v)), "min": round(float(v.min()), 4), "p50": round(float(np.median(v)), 4),
                    "max": round(float(v.max()), 4), "mean": round(float(v.mean()), 4),
                    "std": round(float(v.std()), 4)}
    return out


def _within_patient_variation(df: pd.DataFrame, col: str) -> dict:
    """Mean within-patient std / #distinct — >0 means the variable varies WITHIN a patient (candidate-level;
    can carry a within-patient ranking effect); ~0 means patient-constant (interaction-only / patient tag)."""
    stds, ndist = [], []
    for _, g in df.groupby("patient_id"):
        v = pd.to_numeric(g[col], errors="coerce").to_numpy(float)
        v = v[np.isfinite(v)]
        if len(v) >= 2:
            stds.append(float(v.std())); ndist.append(int(len(np.unique(np.round(v, 6)))))
    return {"mean_within_patient_std": round(float(np.mean(stds)), 4) if stds else 0.0,
            "mean_within_patient_distinct": round(float(np.mean(ndist)), 2) if ndist else 1.0,
            "patient_constant": (float(np.mean(stds)) < 1e-9) if stds else True}


def _pairwise_overlap(df: pd.DataFrame, cols: list[str]) -> dict:
    """Descriptive range overlap on the first numeric encoding column (retained for the record, NOT the gate)."""
    col = cols[0]
    srcs = sorted(df["source"].unique())
    out = {}
    for i, a in enumerate(srcs):
        for b in srcs[i + 1:]:
            va = pd.to_numeric(df.loc[df["source"] == a, col], errors="coerce").to_numpy(float)
            vb = pd.to_numeric(df.loc[df["source"] == b, col], errors="coerce").to_numpy(float)
            out[f"{a}|{b}"] = round(_overlap(va, vb), 4)
    return out


# --------------------------------------------------------------------------------------------------
# the ACTUAL source-alias audit: patient-grouped, source-balanced source classification
# --------------------------------------------------------------------------------------------------
def _sb_pb_weights(df: pd.DataFrame) -> np.ndarray:
    """Source-balanced x patient-balanced row weights: each source totals weight 1, each patient equal within
    source, rows within a patient share the patient weight. Prevents Gartner's ~290k rows from dominating and
    treats patients (not rows) as the effective independent unit. Mean-normalized for solver stability."""
    nsrc = df["source"].nunique()
    npat_src = df.groupby("source")["patient_id"].transform("nunique").astype(float)
    nrow_pat = df.groupby("patient_id")["peptide"].transform("size").astype(float)
    w = 1.0 / (nsrc * npat_src * nrow_pat)
    return (w / w.mean()).to_numpy(float)


def _alias_metrics(df: pd.DataFrame, cols: list[str], w: np.ndarray) -> dict:
    """Predict SOURCE from the context column(s) with GROUP-held-out patients (StratifiedGroupKFold, stratified
    on source, grouped on patient) and source-balanced x patient-balanced weights. Report macro balanced
    accuracy (mean per-source recall; source-balanced chance = 1/3) and macro OVR AUROC (chance = 0.5)."""
    X = df[cols].to_numpy(float)
    y = df["source"].to_numpy()
    groups = df["patient_id"].to_numpy()
    classes = sorted(df["source"].unique())
    proba = np.zeros((len(df), len(classes)))
    pred = np.empty(len(df), dtype=object)
    sgkf = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=CV_SEED)
    for tr, te in sgkf.split(X, y, groups):
        clf = LogisticRegression(max_iter=2000, class_weight="balanced", C=1.0)
        clf.fit(X[tr], y[tr], sample_weight=w[tr])
        p = clf.predict_proba(X[te])
        # align columns to global class order
        idx = [list(clf.classes_).index(c) for c in classes]
        proba[te] = p[:, idx]
        pred[te] = clf.classes_[p.argmax(1)]
    yoh = pd.get_dummies(pd.Series(y))[classes].to_numpy()
    return {"macro_balanced_accuracy": round(float(balanced_accuracy_score(y, pred)), 4),
            "macro_ovr_auroc": round(float(roc_auc_score(yoh, proba, average="macro", sample_weight=w)), 4),
            "source_balanced_chance": {"macro_balanced_accuracy": round(1.0 / len(classes), 4),
                                       "macro_ovr_auroc": 0.5}}


# --------------------------------------------------------------------------------------------------
# candidate contexts
# --------------------------------------------------------------------------------------------------
def _locus(allele: str) -> str:
    """PREREGISTERED identifiable encoding: HLA class-I locus in {A,B,C} from the normalized allele string.
    Handles every source format ('A0201', 'HLA-A01:01', 'HLA-A*02:01'); anything else -> 'OTHER' (reference)."""
    a = str(allele)
    m = re.search(r"HLA[-_ ]?([ABC])", a, re.I) or re.match(r"\s*([ABC])[\*0-9]", a, re.I)
    return m.group(1).upper() if m else "OTHER"


def build_context_columns(ev: pd.DataFrame) -> pd.DataFrame:
    """Attach every CANDIDATE context to the (non-quarantined) eval frame. Approved contexts are candidate-level
    and inference-derivable from the candidate itself; rejected ones depend on source-specific pool enumeration
    or are already-modeled transforms (kept here only so the audit records WHY they were rejected)."""
    ev = ev.copy()
    # --- APPROVED candidate-level contexts (per-candidate, no pool enumeration) ---
    ev["ctx_pep_len"] = ev["peptide"].astype(str).str.len().astype(float)          # from the sequence alone
    ev["ctx_pred_disagree"] = ev[DISAGREE_COLS].std(axis=1, ddof=0).astype(float)  # rowwise SD of masked pcts
    loc = ev["hla_allele"].map(_locus)
    ev["ctx_locus_B"] = (loc == "B").astype(float)   # A is the reference level; loc_B, loc_C are the encoding
    ev["ctx_locus_C"] = (loc == "C").astype(float)
    # --- REJECTED: pool-enumeration / denominator-dependent (structurally rejected; alias number corroborates) ---
    pool = ev.groupby("patient_id")["peptide"].transform("size").astype(float)
    ev["ctx_log_pool"] = np.log1p(pool)
    ev["ctx_hla_multiplicity"] = ev.groupby("patient_id")["hla_allele"].transform("nunique").astype(float)
    comp_raw = ev.groupby(["patient_id", "hla_allele"])["peptide"].transform("size").astype(float)
    ev["ctx_hla_competition_raw"] = comp_raw
    ev["ctx_hla_comp_frac"] = comp_raw / pool
    ev["ctx_routes_per_bag"] = ev.groupby("bag_id")["peptide"].transform("size").astype(float)
    # --- REJECTED: already-modeled transforms ---
    ev["ctx_pred_consensus"] = ev[DISAGREE_COLS].mean(axis=1).astype(float)  # deterministic mean of modeled feats
    ev["ctx_prime_masked"] = ev["prime_masked"].astype(float)
    return ev


# name -> (encoding columns, pool_enumeration_dependent, structural_reason_if_rejected, description)
CANDIDATES = {
    "ctx_pep_len": (["ctx_pep_len"], False, None,
                    "peptide length (per-candidate, from sequence)"),
    "ctx_pred_disagree": (["ctx_pred_disagree"], False, None,
                          "rowwise SD of the leakage-safe masked predictor percentiles "
                          "[f_prime_pct,f_mix_pct,f_el_pct] (per-candidate)"),
    "ctx_hla_locus": (["ctx_locus_B", "ctx_locus_C"], False, None,
                      "HLA class-I locus A/B/C from the normalized allele (per-candidate; A=reference)"),
    "ctx_log_pool": (["ctx_log_pool"], True,
                     "pool-size has no source-invariant denominator (Gartner is combinatorially "
                     "peptide/HLA-expanded; IMPROVE/multimer denominators built differently)",
                     "log1p(patient candidate-pool size)"),
    "ctx_hla_multiplicity": (["ctx_hla_multiplicity"], True,
                             "#distinct HLA in the pool is enumeration-dependent, not documented genotype "
                             "breadth (no completeness flag across sources)",
                             "distinct HLA alleles observed in the patient pool"),
    "ctx_hla_comp_frac": (["ctx_hla_comp_frac"], True,
                          "still depends on source-specific pool enumeration (denominator = pool size)",
                          "HLA competition as within-patient pool fraction"),
    "ctx_hla_competition_raw": (["ctx_hla_competition_raw"], True,
                                "raw count is pool-size-driven; near-disjoint across sources",
                                "raw HLA competition count"),
    "ctx_routes_per_bag": (["ctx_routes_per_bag"], True,
                           "bag cardinality = routes/mutation; a categorical bag-vs-instance ascertainment "
                           "proxy (Gartner 1-74 vs instance sources ~1)",
                           "bag cardinality = routes/mutation"),
    "ctx_pred_consensus": (["ctx_pred_consensus"], False,
                           "redundant: deterministic mean of already-modeled presentation features; carries no "
                           "portable information beyond x",
                           "mean of the masked predictor percentiles (confidence transform)"),
    "ctx_prime_masked": (["ctx_prime_masked"], False,
                         "already modeled (f_prime_pct->0.5 on masked rows) and its rate is source-skewed; "
                         "re-using it as context risks re-introducing leaked-PRIME structure",
                         "PRIME leakage-mask indicator"),
}


def main() -> None:
    with guard_no_test_io():
        frame = assemble_frame()
    ev = frame[~frame["quarantined"]].copy().reset_index(drop=True)   # OUTCOME labels never read below
    ev = build_context_columns(ev)
    w = _sb_pb_weights(ev)

    report = {
        "experiment": "epicurus_v0.5_context_feasibility_and_alias_audit",
        "note": "PRE-FIT, OUTCOME-label-blind, no Q/R model, no ranking metric. Source used descriptively for "
                "aliasing only. The GATE is the patient-grouped source-balanced source-classification alias "
                "audit; range overlap is descriptive only.",
        "n_rows_eval_pool": int(len(ev)),
        "n_patients": int(ev["patient_id"].nunique()),
        "sources": sorted(ev["source"].unique()),
        "per_source_patient_counts": {s: int(g["patient_id"].nunique()) for s, g in ev.groupby("source")},
        "forbidden_label_cols_untouched": list(FORBIDDEN_LABEL_COLS),
        "alias_audit_spec": {
            "task": "predict source from the context encoding column(s)",
            "cv": "StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=%d) — stratify=source, group=patient"
                  % CV_SEED,
            "weighting": "source-balanced x patient-balanced (each source totals 1; patients equal within source; "
                         "rows within a patient share) — patients are the effective independent unit",
            "classifier": "LogisticRegression(lbfgs, class_weight='balanced', C=1.0) — deterministic",
            "metrics": "macro balanced accuracy (chance 1/3) and macro one-vs-rest AUROC (chance 0.5)",
            "frozen_rule": "APPROVE iff candidate-level & not pool-enumeration-dependent AND standalone "
                           "macro_ovr_auroc <= %.2f" % ALIAS_AUROC_MAX,
        },
        "contexts": {},
    }

    approved, rejected = [], {}
    for name, (cols, pool_dep, structural_reason, desc) in CANDIDATES.items():
        alias = _alias_metrics(ev, cols, w)
        wp = _within_patient_variation(ev, cols[0])
        is_alias = alias["macro_ovr_auroc"] > ALIAS_AUROC_MAX
        if structural_reason is not None:
            verdict, reason = "REJECT", structural_reason + (
                "; also a source alias (AUROC %.3f > %.2f)" % (alias["macro_ovr_auroc"], ALIAS_AUROC_MAX)
                if is_alias else "")
        elif is_alias:
            verdict, reason = "REJECT", "source alias: standalone macro AUROC %.3f > %.2f" % (
                alias["macro_ovr_auroc"], ALIAS_AUROC_MAX)
        else:
            verdict, reason = "APPROVE", None
        report["contexts"][name] = {
            "description": desc,
            "encoding_cols": cols,
            "verdict": verdict,
            "reject_reason": reason,
            "inference_computable": True,
            "pool_enumeration_dependent": pool_dep,
            "alias": alias,
            "within_patient_variation": wp,
            "per_source": _per_source_stats(ev, cols[0]),
            "pairwise_range_overlap_DESCRIPTIVE_ONLY": _pairwise_overlap(ev, cols),
        }
        (approved.append(name) if verdict == "APPROVE" else rejected.setdefault(name, reason))

    # residual JOINT aliasing of the approved block (diagnostic; not a per-context gate)
    block_cols = [c for n in approved for c in CANDIDATES[n][0]]
    report["approved_block_alias"] = {"contexts": approved, "encoding_cols": block_cols,
                                      "alias": _alias_metrics(ev, block_cols, w) if block_cols else None}
    report["approved_contexts"] = approved
    report["rejected_contexts"] = rejected

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps({"approved": approved, "rejected": list(rejected),
                      "block_alias": report["approved_block_alias"]["alias"]}, indent=2))


if __name__ == "__main__":
    main()
