"""Recognition-transfer STAGE 1 — NON-SID selection + freeze (CONTRACT.md; §5 leakage-clean).

Selects ONE config/null using ONLY non-Sid IMPROVE, evaluated OUT-OF-FOLD over the 5 official
patient-disjoint Partitions with PEPTIDE-LEAKAGE QUARANTINE (held-out rows whose mutant peptide exact/near-
matches outer-train are dropped from the leakage-clean metric; selection is on leakage-clean only). Declared
config = (arm, C, α, q): anchored logistic `score = prime_pct + α·OOF_pred_pct`, then a promote-side q-slot
mutant-RNA reserve on that score (q=0 pure anchored; α=0 pure PRIME/reserve; α=0,q=0 null; α>0,q>0 combined).
Controls: matched-random reserve, per-cancer-cohort transport, external Gartner/multimer transport.
Freezes the §6 winner (config + SHA-256). **Accesses NO Sid file/label.** Stage 2 does not exist yet.

    python -m scripts.sid_recognition_transfer

Writes artifacts/milestone_7_decision/sid_recognition_transfer/{stage1_nonsid.json, STAGE1_REPORT.md}.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from event_b.leakage_registry import _kmer_index, canonical_peptide, near_duplicate

IMPROVE_ZIP = Path("data/raw/improve/data.zip")
IMPROVE_MEMBER = "data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt"
PRIME_CACHE = Path("data/raw/gartner_nci/_cache_improve_prime.tsv")
POOL = Path("artifacts/milestone_7_decision/pool_size_sensitivity")
ALLOWED_DATA_FILES = {str(IMPROVE_ZIP), str(PRIME_CACHE)}
ART = Path("artifacts/milestone_7_decision/sid_recognition_transfer")
FROZEN = Path("configs/frozen/sid_recognition_gate_v1.json")
K = 20
# Preregistered CONTRACT §4 grid (exact).
ALPHA_GRID = [0.0, 0.10, 0.20, 0.25, 0.30, 0.50]
Q_GRID = [0, 1, 2, 3, 4]
C_GRID = [0.5, 1.0, 2.0]
CATASTROPHIC = -0.10
NEAR = 0.8
ARMS = {"core_deployable": ["prime", "el", "expr", "VarAlFreq"],
        "improve_rich_partial_bridge": ["prime", "el", "expr", "rna_var", "rna_af", "ValMutRNACoef",
                                        "VarAlFreq", "CelPrev"]}
ARM_RANK = {"core_deployable": 0, "improve_rich_partial_bridge": 1}  # simpler core preferred in ties


def stable_seed(pid) -> int:
    """Process-stable per-patient seed (Python hash() is randomized across processes)."""
    return int.from_bytes(hashlib.sha256(str(pid).encode()).digest()[:4], "big")


def tiebreak_key(cfg, hits):
    """Conservative deterministic tie-break: prefer NULL, then lower q, lower alpha, simpler (core) arm,
    lower C. Used with min() over (config, hits)."""
    arm, C, a, q = cfg
    is_null = (a == 0 and q == 0)
    return (-hits, 0 if is_null else 1, q, a, ARM_RANK[arm], C)


def load_improve_full() -> pd.DataFrame:
    with zipfile.ZipFile(IMPROVE_ZIP) as z, z.open(IMPROVE_MEMBER) as fh:
        d = pd.read_csv(fh, sep="\t")
    d = d[d["response"].isin([0, 1])].copy()
    pr = pd.read_csv(PRIME_CACHE, sep="\t").rename(
        columns={"mutant_peptide": "Mut_peptide", "hla_allele": "HLA_allele"})
    d = d.merge(pr, on=["Mut_peptide", "HLA_allele"], how="left")
    out = pd.DataFrame({
        "patient_id": "improve:" + d["Patient"].astype(str), "source": d["cohort"].astype(str),
        "partition": d["Partition"].astype(int), "mut_peptide": d["Mut_peptide"].astype(str),
        "label": np.where(d["response"] == 1, "POSITIVE", "TESTED_NEGATIVE"),
        "prime": pd.to_numeric(d["prime_rank"], errors="coerce"),
        "el": pd.to_numeric(d["RankEL"], errors="coerce"),
        "expr": pd.to_numeric(d["Expression"], errors="coerce")})
    for c in ["rna_var", "rna_af", "ValMutRNACoef", "VarAlFreq", "CelPrev"]:
        out[c] = pd.to_numeric(d[c], errors="coerce")
    return out.reset_index(drop=True)


def leaked_mask(df: pd.DataFrame) -> np.ndarray:
    """Per row: True if (held-out) mutant peptide exact/near(>=0.8)-matches its outer-train partition."""
    canon = df["mut_peptide"].map(canonical_peptide)
    leaked = np.zeros(len(df), bool)
    for p in sorted(df["partition"].unique()):
        te = (df["partition"] == p).to_numpy()
        train = {canon.iloc[i] for i in np.where(~te)[0]} - {""}
        idx = _kmer_index(train)
        for i in np.where(te)[0]:
            c = canon.iloc[i]
            if c and (c in train or near_duplicate(c, idx, threshold=NEAR) is not None):
                leaked[i] = True
    return leaked


# ---- percentiles / OOF logistic ----------------------------------------------------------------------
def _pct(df, col, higher_better):
    v = pd.to_numeric(df[col], errors="coerce")
    v = v if higher_better else -v
    return v.groupby(df["patient_id"]).rank(pct=True).fillna(0.5).to_numpy()


def _feat(df, cols):
    return np.column_stack([_pct(df, c, c not in {"prime", "el"}) for c in cols])


def _bal(df):
    w = np.ones(len(df))
    for _, idx in df.groupby("patient_id").groups.items():
        w[df.index.get_indexer(idx)] = 1.0 / len(idx)
    y = (df["label"].to_numpy() == "POSITIVE")
    wp, wn = w[y].sum(), w[~y].sum()
    return np.where(y, w * (0.5 / wp if wp else 1), w * (0.5 / wn if wn else 1))


def oof_pred(df, cols, C):
    pred = np.full(len(df), np.nan)
    y = (df["label"].to_numpy() == "POSITIVE").astype(int)
    X = _feat(df, cols)
    for p in sorted(df["partition"].unique()):
        te = (df["partition"] == p).to_numpy()
        tr = ~te
        if y[tr].sum() and (y[tr] == 0).sum():
            clf = LogisticRegression(max_iter=2000, C=C).fit(X[tr], y[tr], sample_weight=_bal(df.iloc[tr]))
            pred[te] = clf.predict_proba(X[te])[:, 1]
    return pred


# ---- config -> per-patient hits@20 over the FULL pool; clean_mask ONLY gates hit-counting -------------
def config_hits(df, prime_pct, oof, alpha, q, clean_mask, *, random_reserve=False, seed=0):
    """Rank the FULL patient pool (baseline & challenger face identical pools); a selected positive counts
    only if it is leakage-CLEAN (clean_mask True). Leaked candidates still compete for the 20 slots."""
    score = prime_pct + alpha * pd.Series(oof).groupby(df["patient_id"].to_numpy()).rank(pct=True).fillna(0.5).to_numpy()
    rna = pd.to_numeric(df["rna_af"], errors="coerce").to_numpy()
    ispos = (df["label"].to_numpy() == "POSITIVE")
    out = {}
    for pid, g in df.groupby("patient_id"):
        loc = df.index.get_indexer(g.index)          # FULL pool for this patient
        s, r = score[loc], rna[loc]
        n = len(loc)
        k = min(K, n)
        order_s = np.argsort(-s, kind="mergesort")
        prot = list(order_s[: max(k - q, 0)])
        ps = set(prot)
        if q and not random_reserve:
            res = [i for i in np.argsort(-np.where(np.isfinite(r), r, -1), kind="mergesort")
                   if i not in ps and np.isfinite(r[i]) and r[i] > 0][: max(k - len(prot), 0)]
        elif q:
            rng = np.random.default_rng([seed, stable_seed(pid)])
            elig = [i for i in range(n) if i not in ps and np.isfinite(r[i]) and r[i] > 0]
            res = list(rng.permutation(elig))[: max(k - len(prot), 0)]
        else:
            res = []
        chosen = prot + res
        chosen += [i for i in order_s if i not in set(chosen)][: k - len(chosen)]
        chosen = chosen[:k]
        gloc = loc[chosen]
        out[str(pid)] = float((ispos[gloc] & clean_mask[gloc]).sum())  # leaked positives do not count
    return out


def _cohort_delta(df, hits, base):
    src = {p: df[df.patient_id == p].source.iloc[0] for p in hits}
    coh = sorted(set(src.values()))
    return {c: round(sum(hits.get(p, 0) - base.get(p, 0) for p in set(hits) | set(base) if src.get(p) == c), 4)
            for c in coh}


def fit_predict(train, test, cols, C):
    y = (train["label"].to_numpy() == "POSITIVE").astype(int)
    clf = LogisticRegression(max_iter=2000, C=C).fit(_feat(train, cols), y, sample_weight=_bal(train))
    return clf.predict_proba(_feat(test, cols))[:, 1]


def serialize_fitted_model(df, cols, C, alpha, q):
    """Fit the selected family logistic on ALL non-Sid IMPROVE and serialize everything Stage 2 needs to
    APPLY it WITHOUT refitting: feature order, coefficients, intercept, C, sklearn version, the percentile/
    direction policy, the anchor (alpha) and reserve (q), and a payload hash. Stage 2 must only apply these
    coefficients (linear predictor -> sigmoid), never refit after Sid is available."""
    import sklearn
    y = (df["label"].to_numpy() == "POSITIVE").astype(int)
    clf = LogisticRegression(max_iter=2000, C=C).fit(_feat(df, cols), y, sample_weight=_bal(df))
    payload = {
        "feature_order": list(cols),
        "coef": [float(x) for x in clf.coef_[0]],          # FULL precision (rounding could flip top-20 ties)
        "intercept": float(clf.intercept_[0]),
        "C": C, "alpha": alpha, "q": q,
        "sklearn_version": sklearn.__version__,
        "percentile_policy": "within-patient rank(pct=True) of the oriented raw feature; NaN->0.5",
        "direction_policy": {c: ("lower_raw_better" if c in {"prime", "el"} else "higher_raw_better") for c in cols},
        "reserve": "promote-side: protect top (20-q) by score, fill q by mutant-RNA rna_af (higher=better), "
                   "backfill by score; q disabled where rna_af absent",
        "apply": "score = within_patient_pct(-prime) + alpha * within_patient_pct(sigmoid(coef.feat_pct + intercept))",
    }
    payload["model_payload_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()
    return payload


def configs_for(arm):  # per-family config space INCLUDING the null (alpha=0, q=0)
    return [(arm, C, a, q) for C in C_GRID for a in ALPHA_GRID for q in Q_GRID]


def _clean_total(sub, oof, alpha, q, clean_sub):
    oof_use = oof if alpha > 0 else np.zeros(len(sub))
    return sum(config_hits(sub, _pct(sub, "prime", False), oof_use, alpha, q, clean_sub).values())


def nested_cv_family(df, keep_clean_full, arm):
    """Nested CV for ONE family/arm. Inner selection (incl. null, conservative tie-break) uses an inner
    leakage mask RECOMPUTED within the outer-train subset (so outer-test peptide membership cannot leak
    into inner selection). Outer-test uses the full-data outer-fold mask."""
    parts = sorted(df["partition"].unique())
    part_arr = df["partition"].to_numpy()
    per_fold, outer_hits, outer_null, outer_rand = [], {}, {}, {}
    for f in parts:
        inner = df[part_arr != f].reset_index(drop=True)
        clean_in = ~leaked_mask(inner)                              # recomputed inside outer-train
        oof_in = {C: oof_pred(inner, ARMS[arm], C) for C in C_GRID}
        scored = [((arm, C, a, q), _clean_total(inner, oof_in[C], a, q, clean_in))
                  for C in C_GRID for a in ALPHA_GRID for q in Q_GRID]
        best = min(scored, key=lambda kv: tiebreak_key(kv[0], kv[1]))[0]
        _, C, a, q = best
        test = df[part_arr == f].reset_index(drop=True)
        clean_te = keep_clean_full[part_arr == f]
        pp_te = _pct(test, "prime", False)
        pred = fit_predict(inner, test, ARMS[arm], C) if a > 0 else np.zeros(len(test))
        h = config_hits(test, pp_te, pred, a, q, clean_te)
        h0 = config_hits(test, pp_te, np.zeros(len(test)), 0.0, 0, clean_te)
        # nested matched-random comparator for the q reserve (same #slots, random RNA order)
        if q > 0:
            hr = {p: float(np.mean([config_hits(test, pp_te, pred, a, q, clean_te, random_reserve=True,
                                                seed=s)[p] for s in range(20)])) for p in h}
        else:
            hr = dict(h)
        outer_hits.update(h)
        outer_null.update(h0)
        outer_rand.update(hr)
        per_fold.append({"fold": int(f), "chosen": {"arm": arm, "C": C, "alpha": a, "q": q},
                         "outer_test_hits": int(sum(h.values())), "null_hits": int(sum(h0.values())),
                         "matched_random_hits": round(float(sum(hr.values())), 2)})
    return per_fold, outer_hits, outer_null, outer_rand


def _paired_bootstrap(a: dict, b: dict, n: int = 5000, seed: int = 12345) -> dict:
    """Patient-level paired bootstrap CI of mean(a-b). Fixed seed; returns delta + 95% CI + P(>0)."""
    pids = sorted(a)
    d = np.array([a[p] - b.get(p, 0) for p in pids], float)
    rng = np.random.default_rng(seed)
    means = np.array([d[rng.integers(0, len(d), len(d))].mean() for _ in range(n)])
    return {"mean_delta": round(float(d.mean()), 4),
            "ci95": [round(float(np.percentile(means, 2.5)), 4), round(float(np.percentile(means, 97.5)), 4)],
            "p_gt_0": round(float((means > 0).mean()), 4), "n_patients": len(pids), "seed": seed}


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    df = load_improve_full()
    keep_all = np.ones(len(df), bool)
    keep_clean = ~leaked_mask(df)
    n_leaked = int((~keep_clean).sum())

    prime_pct = _pct(df, "prime", False)
    null_hits = config_hits(df, prime_pct, np.zeros(len(df)), 0.0, 0, keep_clean)
    null_tot = {"clean": sum(null_hits.values()),
                "all": sum(config_hits(df, prime_pct, np.zeros(len(df)), 0.0, 0, keep_all).values())}

    # (B) NESTED CV PER FAMILY -> choose an eligible family by nested outer evidence.
    families = {}
    for arm in ARMS:
        pf, oh, on, orand = nested_cv_family(df, keep_clean, arm)
        nt, nn = int(sum(oh.values())), int(sum(on.values()))
        coh = _cohort_delta(df, oh, on)
        families[arm] = {"nested_total": nt, "nested_null": nn, "nested_delta": nt - nn,
                         "per_cohort_delta": coh, "chosen_config_per_fold": pf,
                         "nested_matched_random_total": round(float(sum(orand.values())), 2),
                         "beats_nested_matched_random": bool(nt > sum(orand.values()) + 1e-9),
                         "paired_bootstrap_delta_vs_null": _paired_bootstrap(oh, on),
                         "eligible": bool(nt > nn and min(coh.values()) >= CATASTROPHIC
                                          and all(v >= 0 for v in coh.values()))}
    elig_fams = [a for a in ARMS if families[a]["eligible"]]
    sel_family = min(elig_fams, key=lambda a: (-families[a]["nested_total"], ARM_RANK[a])) if elig_fams else None

    # (A) Deployment params by FULL non-Sid CV ONLY within the selected family.
    grid = []
    if sel_family is not None:
        for C in C_GRID:
            oof = oof_pred(df, ARMS[sel_family], C)
            for a in ALPHA_GRID:
                for q in Q_GRID:
                    if a == 0 and q == 0:
                        continue
                    clean = config_hits(df, prime_pct, oof, a, q, keep_clean)
                    rand = np.mean([sum(config_hits(df, prime_pct, oof, a, q, keep_clean, random_reserve=True,
                                                    seed=s).values()) for s in range(20)]) if q else None
                    coh = _cohort_delta(df, clean, null_hits)
                    grid.append({"arm": sel_family, "C": C, "alpha": a, "q": q,
                                 "hits_clean": int(sum(clean.values())),
                                 "delta_clean": int(sum(clean.values()) - null_tot["clean"]),
                                 "per_cohort_delta_clean": coh, "worst_cohort": round(min(coh.values()), 4),
                                 "matched_random_reserve_hits_clean": round(float(rand), 2) if rand is not None else None,
                                 "beats_random": bool(rand is None or sum(clean.values()) > rand + 1e-9)})
    eligible = [g for g in grid if g["delta_clean"] > 0 and g["worst_cohort"] >= CATASTROPHIC
                and all(v >= 0 for v in g["per_cohort_delta_clean"].values()) and g["beats_random"]]
    winner = min(eligible, key=lambda g: tiebreak_key((g["arm"], g["C"], g["alpha"], g["q"]),
                                                      g["hits_clean"])) if eligible else None
    ext = external_transport(df, winner) if winner else None
    fitted = (serialize_fitted_model(df, ARMS[winner["arm"]], winner["C"], winner["alpha"], winner["q"])
              if winner else None)
    cfg = {"name": "sid_recognition_gate", "version": "1.0.0", "stage": "1_non_sid_frozen",
           "frozen": winner or {"arm": "null", "alpha": 0.0, "q": 0},
           "fitted_model": fitted,
           "stage2_must_not_refit": "Stage 2 applies fitted_model.coef/intercept only (linear predictor -> "
                                    "sigmoid); it MUST NOT refit any model after Sid is available.",
           "null_hits_clean": null_tot["clean"], "n_leaked_rows_quarantined": n_leaked,
           "selected_on": "IMPROVE official partitions OOF, LEAKAGE-CLEAN only; matched-random; transport",
           "sid_never_consulted": True,
           "declared_multiplicity_configs": len(ARMS) * len(C_GRID) * len(ALPHA_GRID) * len(Q_GRID),
           "nested_cv_per_family": families,
           "selected_family_by_nested_evidence": sel_family,
           "deployment_params_by": "full non-Sid CV within the selected family only",
           "external_transport_gartner_multimer": ext,
           "improve_to_sid_feature_map": {"prime": "prime_rank", "el": "mixmhcpred_rank(EL proxy)",
                                          "expr": "gene TPM", "VarAlFreq": "WES tumor VAF",
                                          "rna_af": "variant_vafs_long tumor-RNA VAF", "CelPrev": "MISSING->neutral"},
           "limitations": ("patient-level paired bootstrap CI for the nested per-family delta is COMPUTED "
                           "and reported (see nested_cv_per_family[*].paired_bootstrap_delta_vs_null; fixed "
                           "seed 12345), NOT gated. Stage 1 was NARROWED to the exact runner vs the original "
                           "contract (see CONTRACT.md 'STAGE-1 v2 PROTOCOL CORRECTION'): PRIME base only "
                           "(MixMHCpred base deferred), 2 anchored full-feature arms (not the physchem/"
                           "expression/absence families), 5 IMPROVE official partitions (Gartner/multimer are "
                           "an ANCHORED-COMPONENT-ONLY external transport check, q disabled there). Sid "
                           "previously inspected -> Stage 2 is exploratory confirmation, not pristine validation.")}
    cfg["sha256"] = hashlib.sha256(json.dumps(cfg, sort_keys=True).encode()).hexdigest()
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(json.dumps(cfg, indent=2) + "\n")

    report = {"stage": 1, "n_rows": int(len(df)), "n_leaked_quarantined": n_leaked,
              "null_hits": null_tot, "grid": grid, "frozen": cfg, "sha256": cfg["sha256"],
              "data_files_read": sorted(ALLOWED_DATA_FILES),
              "checkpoint": "PRE-SID: no Sid file/label accessed; STOP for audit before Stage 2."}
    (ART / "stage1_nonsid.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    (ART / "STAGE1_REPORT.md").write_text(_md(report))

    print(f"leakage: quarantined {n_leaked}/{len(df)} held-out rows (exact/near>=0.8). null hits clean="
          f"{null_tot['clean']} all={null_tot['all']}")
    top = sorted(grid, key=lambda g: -g["hits_clean"])[:6]
    for g in top:
        print(f"  {g['arm'][:5]} C={g['C']} a={g['alpha']} q={g['q']}: clean {g['hits_clean']} "
              f"(Δ{g['delta_clean']:+d}) cohorts{g['per_cohort_delta_clean']} rand={g['matched_random_reserve_hits_clean']} "
              f"beats_rand={g['beats_random']}")
    for arm, fam in families.items():
        bs = fam["paired_bootstrap_delta_vs_null"]
        print(f"  NESTED[{arm[:5]}]: total {fam['nested_total']} vs null {fam['nested_null']} "
              f"(Δ{fam['nested_delta']:+d}) cohorts{fam['per_cohort_delta']} eligible={fam['eligible']} "
              f"bootΔ={bs['mean_delta']:+.3f} CI{bs['ci95']} p>0={bs['p_gt_0']} "
              f"matched_random={fam['nested_matched_random_total']} beats_rand={fam['beats_nested_matched_random']}")
    print(f"selected family (nested): {sel_family}")
    print(f"FROZEN: {cfg['frozen'].get('arm')} alpha={cfg['frozen'].get('alpha')} q={cfg['frozen'].get('q')} "
          f"C={cfg['frozen'].get('C')} clean_hits={cfg['frozen'].get('hits_clean')} sha={cfg['sha256'][:16]}")
    print(f"external transport (winner): {ext}")
    print("CHECKPOINT:", report["checkpoint"])
    return 0


def _load_external(name):
    f = pd.read_csv(POOL / f"base_{name}.csv")
    f = f[f["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    f["patient_id"] = f["patient_id"].astype(str)
    return f.reset_index(drop=True)


def external_transport(df, winner):
    if winner is None:
        return None
    cols = ARMS[winner["arm"]]
    y = (df["label"].to_numpy() == "POSITIVE").astype(int)
    clf = LogisticRegression(max_iter=2000, C=winner["C"]).fit(_feat(df, cols), y, sample_weight=_bal(df))
    out = {}
    for name in ["gartner", "multimer"]:
        ext = _load_external(name)
        feats = [( _pct(ext, c, c not in {"prime", "el"}) if c in ext.columns else np.full(len(ext), 0.5))
                 for c in cols]
        pred = clf.predict_proba(np.column_stack(feats))[:, 1]
        pp = _pct(ext, "prime", False)
        km = np.ones(len(ext), bool)
        base = config_hits(ext.assign(rna_af=np.nan), pp, np.zeros(len(ext)), 0.0, 0, km)
        alt = config_hits(ext.assign(rna_af=np.nan), pp, pred, winner["alpha"], 0, km)  # q=0 externally (no RNA/VAF signal)
        out[name] = {"base": int(sum(base.values())), "gated": int(sum(alt.values())),
                     "delta": int(sum(alt.values()) - sum(base.values())),
                     "note": "ANCHORED-COMPONENT-ONLY: q reserve DISABLED externally (no comparable mutant-RNA "
                             "signal); VAF/RNA features neutral where absent. Not the full frozen policy."}
    return out


def _md(r):
    L = [f"# Recognition transfer — STAGE 1 (non-Sid freeze, leakage-clean)\n\n_frozen SHA-256 "
         f"`{r['sha256'][:16]}`; data: {', '.join(r['data_files_read'])} (NO Sid). "
         f"Quarantined {r['n_leaked_quarantined']}/{r['n_rows']} peptide-leaked held-out rows; "
         f"selection on leakage-clean only. null hits clean={r['null_hits']['clean']}._\n"]
    L.append("\n## Nested CV per family (inner mask recomputed within outer-train; null included, "
             "conservative tie-break)\n")
    L.append("| family | nested total | null | Δ | per-cohort Δ | paired-boot Δ [CI95] p>0 | "
             "matched-random (beats) | eligible | chosen (fold:α,q) |")
    L.append("|---|--:|--:|--:|---|---|---|:--:|---|")
    for arm, fam in r["frozen"]["nested_cv_per_family"].items():
        ch = [(x["fold"], x["chosen"]["alpha"], x["chosen"]["q"]) for x in fam["chosen_config_per_fold"]]
        bs = fam["paired_bootstrap_delta_vs_null"]
        L.append(f"| {arm} | {fam['nested_total']} | {fam['nested_null']} | {fam['nested_delta']:+d} | "
                 f"{fam['per_cohort_delta']} | {bs['mean_delta']:+.3f} {bs['ci95']} {bs['p_gt_0']} | "
                 f"{fam['nested_matched_random_total']} ({'y' if fam['beats_nested_matched_random'] else 'n'}) | "
                 f"{'yes' if fam['eligible'] else 'no'} | {ch} |")
    L.append(f"\n**Selected family (by nested evidence): {r['frozen']['selected_family_by_nested_evidence']}**. "
             "Deployment params chosen by full non-Sid CV within that family only:\n")
    L.append("\n| arm | C | α | q | clean hits (Δ) | per-cohort Δ | matched-random | beats rand |")
    L.append("|---|--:|--:|--:|--:|---|--:|:--:|")
    for g in sorted(r["grid"], key=lambda x: -x["hits_clean"])[:14]:
        L.append(f"| {g['arm']} | {g['C']} | {g['alpha']} | {g['q']} | {g['hits_clean']} ({g['delta_clean']:+d}) | "
                 f"{g['per_cohort_delta_clean']} | {g['matched_random_reserve_hits_clean']} | "
                 f"{'y' if g['beats_random'] else 'n'} |")
    f = r["frozen"]["frozen"]
    L.append(f"\n## FROZEN (§6, leakage-clean)\narm **{f.get('arm')}** α={f.get('alpha')} q={f.get('q')} "
             f"C={f.get('C')} clean_hits={f.get('hits_clean')}; multiplicity {r['frozen']['declared_multiplicity_configs']} "
             f"configs; external transport {r['frozen']['external_transport_gartner_multimer']}. "
             f"SHA-256 `{r['sha256']}`.\n\n_{r['frozen']['limitations']}_\n\n> {r['checkpoint']}\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
