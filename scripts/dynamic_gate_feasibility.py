"""Feasibility probe: can ANY label-blind signal rescue the positives the EL gate loses?

Uses the base cohort CSVs already written by scripts.pool_size_sensitivity
(columns: patient_id, mutant_peptide, hla_allele, label, prime, el, expr).

Question 1: where do positives sit on within-patient percentiles of el / prime / expr?
Question 2: the EL gate at a given keep-fraction loses the low-EL positives.
            Do those lost positives sit HIGHER on prime or expr (i.e. rescuable)?
Question 3: AND-of-vetoes — veto iff el_pct<t AND prime_pct<t AND expr_pct<t.
            What negative-removal / positive-retention does that give vs a pure-EL top-K gate
            at MATCHED removal? (the honest head-to-head)
"""
import numpy as np
import pandas as pd
from pathlib import Path

ART = Path("artifacts/milestone_7_decision/pool_size_sensitivity")

def pct(frame, col, higher_better):
    v = pd.to_numeric(frame[col], errors="coerce")
    if not higher_better:
        v = -v
    return v.groupby(frame["patient_id"]).rank(pct=True).to_numpy()  # NaN stays NaN here

def load(name):
    f = pd.read_csv(ART / f"base_{name}.csv")
    f["patient_id"] = f["patient_id"].astype(str)
    f = f[f["label"].isin(["POSITIVE","TESTED_NEGATIVE"])].copy()
    # keep only patients with >=1 pos and >=4 neg (matches MIN_NEG eligibility)
    keep=[]
    for pid,g in f.groupby("patient_id"):
        if (g.label=="POSITIVE").sum()>=1 and (g.label=="TESTED_NEGATIVE").sum()>=4:
            keep.append(pid)
    f=f[f.patient_id.isin(keep)].copy()
    f["el_pct"]=pct(f,"el",False)
    f["prime_pct"]=pct(f,"prime",False)
    f["expr_pct"]=pct(f,"expr",True)
    # coverage
    f["el_miss"]=pd.to_numeric(f["el"],errors="coerce").isna()
    f["prime_miss"]=pd.to_numeric(f["prime"],errors="coerce").isna()
    f["expr_miss"]=pd.to_numeric(f["expr"],errors="coerce").isna()
    return f

def q1_positions(f,name):
    pos=f[f.label=="POSITIVE"]
    print(f"\n=== {name}: {f.patient_id.nunique()} pts, {len(pos)} pos, {(f.label=='TESTED_NEGATIVE').sum()} neg")
    for c in ["el_pct","prime_pct","expr_pct"]:
        v=pos[c].dropna()
        print(f"  positive {c}: median={v.median():.3f} q25={v.quantile(.25):.3f} "
              f"frac_in_bottom_half={np.mean(v<0.5):.3f} frac_in_bottom_quartile={np.mean(v<0.25):.3f} nan={pos[c].isna().mean():.2f}")
    print(f"  feature missingness: el={f.el_miss.mean():.2f} prime={f.prime_miss.mean():.2f} expr={f.expr_miss.mean():.2f}")

def q2_rescue(f,name,keep_frac=0.5):
    """EL gate keeps top keep_frac by el_pct per patient. Among positives it DROPS, how do they
    look on prime_pct/expr_pct? If they're high there, they're rescuable."""
    dropped=[]
    for pid,g in f.groupby("patient_id"):
        g=g.sort_values("el_pct",ascending=False,na_position="last")
        k=int(np.ceil(keep_frac*len(g)))
        kept=set(g.index[:k])
        dp=g[(g.label=="POSITIVE") & (~g.index.isin(kept))]
        dropped.append(dp)
    dp=pd.concat(dropped) if dropped else f.iloc[:0]
    npos=(f.label=="POSITIVE").sum()
    print(f"\n  [{name}] EL-gate keep_frac={keep_frac}: drops {len(dp)}/{npos} positives ({len(dp)/npos:.2%})")
    if len(dp):
        print(f"     dropped-positive prime_pct: median={dp.prime_pct.median():.3f} frac>0.5={np.mean(dp.prime_pct>0.5):.3f} frac>0.75={np.mean(dp.prime_pct>0.75):.3f}")
        print(f"     dropped-positive expr_pct:  median={dp.expr_pct.median():.3f} frac>0.5={np.mean(dp.expr_pct>0.5):.3f} frac>0.75={np.mean(dp.expr_pct>0.75):.3f}")
        rescuable=np.mean((dp.prime_pct>0.5)|(dp.expr_pct>0.5))
        print(f"     rescuable by (prime_pct>0.5 OR expr_pct>0.5): {rescuable:.3f}")

def q3_andveto(f,name):
    """AND-of-vetoes at a single global percentile threshold t (fillna->keep, i.e. missing never vetoed).
    veto iff el_pct<t AND prime_pct<t AND expr_pct<t. Report neg removal & pos retention.
    Compare to pure-EL top-K gate at the SAME removal."""
    print(f"\n  [{name}] AND-of-vetoes (missing=KEEP) vs matched pure-EL gate:")
    el=f.el_pct.fillna(1.0)
    pr=f.prime_pct.fillna(1.0)
    ex=f.expr_pct.fillna(1.0)
    ispos=(f.label=="POSITIVE").to_numpy()
    print(f"     {'t':>5} {'AND:negRemoved':>15} {'AND:posRetained':>16} | {'EL@matched:posRetained':>22}")
    for t in [0.25,0.4,0.5,0.6,0.7,0.75,0.85]:
        veto=(el<t)&(pr<t)&(ex<t)
        keep=~veto.to_numpy()
        neg_removed=np.mean(~keep[~ispos])
        pos_retained=np.mean(keep[ispos])
        # matched pure-EL gate: remove same TOTAL count globally by lowest el_pct
        n_remove=int((~keep).sum())
        order=np.argsort(el.to_numpy(),kind="mergesort")  # ascending: lowest el first removed
        elremove=np.zeros(len(f),bool)
        elremove[order[:n_remove]]=True
        elkeep=~elremove
        el_pos_ret=np.mean(elkeep[ispos])
        print(f"     {t:>5} {neg_removed:>15.3f} {pos_retained:>16.3f} | {el_pos_ret:>22.3f}")

for name in ["gartner","improve","multimer"]:
    f=load(name)
    q1_positions(f,name)
    for kf in [0.5,0.25]:
        q2_rescue(f,name,kf)
    q3_andveto(f,name)
