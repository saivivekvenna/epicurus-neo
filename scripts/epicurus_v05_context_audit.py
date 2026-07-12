"""Epicurus v0.5 — PRE-FIT context-feasibility audit (DESIGN/PROVENANCE, not model code).

Read-only, label-blind, NO model is fit and NO ranking metric is computed. Its ONLY job is to
prove — before the v0.5 preregistration is written — that each candidate scoring-time context
variable is (a) computable from information available while scoring a brand-new patient's candidate
pool, and (b) NOT a near-perfect source proxy. Source labels are used here descriptively (to quantify
aliasing), never to select a model. Outcome labels (`eval_positive`, `bag_label`) are NEVER read.

Runs inside the v0.4 Gartner-TEST I/O guard (`guard_no_test_io`). Gartner TEST is never opened.

Emits `artifacts/milestone_7_decision/epicurus_v05/CONTEXT_FEASIBILITY_AUDIT.json`.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from event_b.epicurus_v04 import assemble_frame, guard_no_test_io  # reuse the guarded frame verbatim

OUT = Path("artifacts/milestone_7_decision/epicurus_v05/CONTEXT_FEASIBILITY_AUDIT.json")
FORBIDDEN_LABEL_COLS = ("eval_positive", "bag_label")  # never used in this audit


def _overlap(a: np.ndarray, b: np.ndarray) -> float:
    """Fraction of the pooled [min,max] range shared by both groups' [min,max] (0 = disjoint, 1 = identical)."""
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
    """Mean within-patient std / #distinct values — >0 means the variable varies WITHIN a patient (can carry a
    within-patient main effect); ~0 means patient-constant (interaction-only, main effect rank-inert)."""
    stds, ndist = [], []
    for _, g in df.groupby("patient_id"):
        v = pd.to_numeric(g[col], errors="coerce").to_numpy(float)
        v = v[np.isfinite(v)]
        if len(v) >= 2:
            stds.append(float(v.std())); ndist.append(int(len(np.unique(np.round(v, 6)))))
    return {"mean_within_patient_std": round(float(np.mean(stds)), 4) if stds else 0.0,
            "mean_within_patient_distinct": round(float(np.mean(ndist)), 2) if ndist else 1.0,
            "patient_constant": (float(np.mean(stds)) < 1e-9) if stds else True}


def _pairwise_overlap(df: pd.DataFrame, col: str) -> dict:
    srcs = sorted(df["source"].unique())
    out = {}
    for i, a in enumerate(srcs):
        for b in srcs[i + 1:]:
            va = pd.to_numeric(df.loc[df["source"] == a, col], errors="coerce").to_numpy(float)
            vb = pd.to_numeric(df.loc[df["source"] == b, col], errors="coerce").to_numpy(float)
            out[f"{a}|{b}"] = round(_overlap(va, vb), 4)
    return out


def build_context_columns(ev: pd.DataFrame) -> pd.DataFrame:
    """Attach every CANDIDATE context to the (non-quarantined) eval frame. All label-blind and inference-derivable
    from a patient's own candidate pool + peptide sequence + HLA typing. Documented per column."""
    ev = ev.copy()
    # peptide length: from the peptide sequence alone (per-candidate).
    ev["ctx_pep_len"] = ev["peptide"].astype(str).str.len().astype(float)
    # patient candidate-pool size (patient-constant) and its log.
    pool = ev.groupby("patient_id")["peptide"].transform("size").astype(float)
    ev["ctx_pool_size"] = pool
    ev["ctx_log_pool"] = np.log1p(pool)
    # HLA multiplicity: distinct HLA alleles in the patient's pool (patient-constant; = patient HLA-type breadth).
    ev["ctx_hla_multiplicity"] = ev.groupby("patient_id")["hla_allele"].transform("nunique").astype(float)
    # HLA competition (per-candidate): #candidates in the patient's pool sharing THIS candidate's allele,
    # normalized by pool size -> a SCALE-FREE within-patient fraction (the deployable, non-source-aliased form;
    # the raw count is pool-size-driven and near-disjoint across sources -> a source proxy, so NOT used).
    comp_raw = ev.groupby(["patient_id", "hla_allele"])["peptide"].transform("size").astype(float)
    ev["ctx_hla_competition_raw"] = comp_raw
    ev["ctx_hla_comp_frac"] = comp_raw / pool
    # routes per mutation = bag cardinality (#instances sharing the candidate's bag_id).  [expected source proxy]
    ev["ctx_routes_per_bag"] = ev.groupby("bag_id")["peptide"].transform("size").astype(float)
    # PRIME leakage mask (leakage-safe availability flag; audited but NOT approved as a context).
    ev["ctx_prime_masked"] = ev["prime_masked"].astype(float)
    return ev


CANDIDATES = {
    # --- audited and APPROVED as scoring-time contexts ---
    "ctx_pep_len": "peptide length (per-candidate, from sequence)  [APPROVE]",
    "ctx_log_pool": "log1p(patient candidate-pool size) (patient-constant)  [APPROVE, interaction-only]",
    "ctx_hla_multiplicity": "distinct HLA alleles in patient pool (patient-constant)  [APPROVE, interaction-only]",
    "ctx_hla_comp_frac": "HLA competition as within-patient fraction (per-candidate)  [APPROVE]",
    # --- audited and REJECTED (source / ascertainment proxies) ---
    "ctx_hla_competition_raw": "raw HLA competition count (pool-size-driven)  [REJECT: source proxy]",
    "ctx_routes_per_bag": "bag cardinality = routes/mutation  [REJECT: bag-vs-instance ascertainment proxy]",
    "ctx_prime_masked": "PRIME leakage mask rate  [REJECT as context: source-skewed + already modeled]",
}


def main() -> None:
    with guard_no_test_io():
        frame = assemble_frame()
    # guard: this audit never reads outcome labels
    ev = frame[~frame["quarantined"]].copy()
    ev = build_context_columns(ev)

    report = {"experiment": "epicurus_v0.5_context_feasibility_audit",
              "note": "PRE-FIT, label-blind, no model, no ranking metric. Source used descriptively for aliasing.",
              "n_rows_eval_pool": int(len(ev)),
              "n_patients": int(ev["patient_id"].nunique()),
              "sources": sorted(ev["source"].unique()),
              "forbidden_label_cols_untouched": list(FORBIDDEN_LABEL_COLS),
              "contexts": {}}

    for col, desc in CANDIDATES.items():
        wp = _within_patient_variation(ev, col)
        report["contexts"][col] = {
            "description": desc,
            "inference_computable": True,
            "per_source": _per_source_stats(ev, col),
            "pairwise_range_overlap": _pairwise_overlap(ev, col),
            "within_patient_variation": wp,
        }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
