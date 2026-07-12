"""v2 feasibility probe: is there BACKFILL headroom, and can a label-blind lever exploit it?

The reframed objective is NET top-20 utility: a gate may remove positives if removing higher-ranked
NEGATIVES backfills MORE positives into the top-20. That is only possible if:
  (H1) positives sit just OUTSIDE the top-20 (backfill supply), AND negatives sit INSIDE it (removal targets);
  (H2) some LABEL-BLIND signal can pick which in-top-20 candidates to remove so positives (not other
       negatives) backfill.

We measure, per cohort, on frozen-Epicurus ranking:
  A. Headroom: mean #positives at ranks 21-40 (backfillable) and #negatives at ranks 1-20 (removable);
     the oracle top-20 ceiling min(npos,20) vs current hits@20.
  B. Diversity lever (label-blind, deployable): does capping redundancy free top-20 slots for positives?
     - exact-peptide dedup (keep best-scoring row per peptide);
     - per-HLA-allele cap (where alleles exist: improve/multimer);
     rerank survivors with frozen Epicurus, recount hits@20.
  C. Negative-risk lever probe: among the top-40, is predictor DISAGREEMENT (Gartner 5 predictors) or
     low EXPRESSION separable between positives and negatives? (mean feature, pos vs neg.)

    python -m scripts.dynamic_gate_v2_feasibility

Prints only. multimer is frozen-Epicurus IN-SAMPLE (flagged). No files written; this is a go/no-go probe.
"""

from __future__ import annotations

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from event_b.pool_size_sensitivity import patient_eligibility, score_arms

POOL = Path("artifacts/milestone_7_decision/pool_size_sensitivity")
K = 20


def load(name):
    f = pd.read_csv(POOL / f"base_{name}.csv")
    f["patient_id"] = f["patient_id"].astype(str)
    f = f[f["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    el = patient_eligibility(f)
    return f[f["patient_id"].isin(el.eligible)].reset_index(drop=True)


def ranked(gp):
    return score_arms(gp).sort_values("frozen_epicurus", ascending=False, kind="mergesort")


def hits20(r):
    return int((r["label"].to_numpy()[:K] == "POSITIVE").sum())


def headroom(name, f):
    per = []
    for _, gp in f.groupby("patient_id"):
        r = ranked(gp)
        lab = r["label"].to_numpy()
        npos = int((lab == "POSITIVE").sum())
        h = int((lab[:K] == "POSITIVE").sum())
        neg_in_top = int((lab[:K] == "TESTED_NEGATIVE").sum())
        pos_backfill = int((lab[K:2 * K] == "POSITIVE").sum())  # ranks 21-40
        ceil = min(npos, K)
        per.append((h, ceil, neg_in_top, pos_backfill))
    per = np.array(per, float)
    print(f"\n=== {name}: {f.patient_id.nunique()} pts")
    print(f"  ungated hits@20 mean={per[:,0].mean():.3f}  oracle ceiling min(npos,20) mean={per[:,1].mean():.3f}  "
          f"gap={per[:,1].mean()-per[:,0].mean():.3f}")
    print(f"  removable negatives in top20 mean={per[:,2].mean():.2f}  positives in rank21-40 (backfill supply) "
          f"mean={per[:,3].mean():.2f}")
    return per


def dedup_best(gp, by):
    """Keep the single best-frozen-Epicurus row per group key `by` (label-blind)."""
    s = score_arms(gp)
    return s.sort_values("frozen_epicurus", ascending=False, kind="mergesort").drop_duplicates(by, keep="first")


def allele_cap(gp, cap):
    s = score_arms(gp).sort_values("frozen_epicurus", ascending=False, kind="mergesort")
    if "hla_allele" not in s or s["hla_allele"].nunique() <= 1:
        return None
    s["_rk"] = s.groupby("hla_allele").cumcount()
    return s[s["_rk"] < cap]


def diversity(name, f):
    base_h, dedup_h, cap2_h, cap3_h = [], [], [], []
    for _, gp in f.groupby("patient_id"):
        base_h.append(hits20(ranked(gp)))
        d = dedup_best(gp, "mutant_peptide")
        dedup_h.append(hits20(d.sort_values("frozen_epicurus", ascending=False, kind="mergesort")))
        for cap, store in [(2, cap2_h), (3, cap3_h)]:
            c = allele_cap(gp, cap)
            store.append(hits20(c.sort_values("frozen_epicurus", ascending=False, kind="mergesort")) if c is not None else hits20(ranked(gp)))
    b = np.mean(base_h)
    print(f"  diversity Δhits@20: exact-peptide-dedup {np.mean(dedup_h)-b:+.3f}  "
          f"allele-cap2 {np.mean(cap2_h)-b:+.3f}  allele-cap3 {np.mean(cap3_h)-b:+.3f}")


def negrisk(name, f):
    """Among top-40 by frozen Epicurus, is expression separable pos vs neg? (disagreement needs 5-pred loader.)"""
    pos_expr, neg_expr = [], []
    for _, gp in f.groupby("patient_id"):
        r = ranked(gp).head(2 * K)
        lab = r["label"].to_numpy()
        ex = pd.to_numeric(r["expr"], errors="coerce").to_numpy()
        pos_expr.extend(ex[lab == "POSITIVE"])
        neg_expr.extend(ex[lab == "TESTED_NEGATIVE"])
    pe, ne = np.array(pos_expr, float), np.array(neg_expr, float)
    pe, ne = pe[np.isfinite(pe)], ne[np.isfinite(ne)]
    if len(pe) and len(ne):
        print(f"  top-40 expression: positives median={np.median(pe):.2f} (n={len(pe)})  "
              f"negatives median={np.median(ne):.2f} (n={len(ne)})  → separable if pos>neg")


for name in ["gartner", "improve", "multimer"]:
    f = load(name)
    star = " (IN-SAMPLE)" if name == "multimer" else ""
    headroom(name + star, f)
    diversity(name, f)
    negrisk(name, f)
