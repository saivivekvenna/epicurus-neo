"""Pre-registered incremental TRANSFER test + external-validation lane (north-star protocol, steps 2/4-6).

ONE frozen residual — presentation (NetMHCpan-EL) + tumor-context (expression) ON TOP OF genuine PRIME,
learned only from the independent CD8 multimer cohort's measured POSITIVE/TESTED_NEGATIVE labels — is
applied, frozen, to TWO untouched external cohorts: Gartner NCI (the step-2 transfer test) and IMPROVE
(the step-4/5 external-validation lane). No tuning against either eval cohort.

All three cohorts use genuine GfellerLab PRIME %rank from the installed tool (lower = better), so there is
no reliance on any cohort's ambiguously-oriented precomputed PRIME column. Features are within-patient
PERCENTILES (higher = better) of genuine PRIME, NetMHCpan-EL, expression — a scale-free, decision-relevant
representation that harmonizes TPM vs decile encodings. Design fixed in PROTOCOL.md.

Leakage fail-closed: multimer training peptides that exact/near-match any eval-cohort peptide, or are in
PRIME's training set (PRIME is a feature), are dropped; patient IDs are namespaced-disjoint; PRIME's
synthetic random-proteome negatives are never training data.
"""

from __future__ import annotations

import json
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from benchmark.metrics import hits_at_k
from event_b.decision_benchmark import _card, fail_closed_auroc
from event_b.leakage_registry import _kmer_index, canonical_peptide, near_duplicate
from event_b.minimal_epitope_expansion import DEFAULT_MUT_POS0, DEFAULT_RULES, score_mutations_with_prime
from event_b.prime_adapter import score_prime
from event_b.prime_training import prime_leakage_mask

CACHE = Path("data/raw/gartner_nci")
HLA_FILE = CACHE / "HLA_allotypes.txt"
IMPROVE_ZIP = Path("data/raw/improve/data.zip")
IMPROVE_MEMBER = "data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt"
PRIMARY_RULE = DEFAULT_RULES[0]  # classI_8_11_mut
K = 20
# (raw column, higher_is_better) — oriented so higher percentile = better candidate.
FEATS = [("prime", False), ("el", False), ("expr", True)]


# --------------------------------------------------------------------------------------------------
# Genuine PRIME scoring helpers (cached; PRIME %rank is deterministic per (peptide, allele)).
# --------------------------------------------------------------------------------------------------
def _score_prime_pairs(pairs: pd.DataFrame, cache_name: str) -> pd.DataFrame:
    """Return unique (mutant_peptide, hla_allele) -> prime_rank, cached to a TSV."""
    cache = CACHE / cache_name
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    uniq = pairs.rename(columns={pairs.columns[0]: "mutant_peptide", pairs.columns[1]: "hla_allele"})
    uniq = uniq[["mutant_peptide", "hla_allele"]].astype(str).drop_duplicates()
    res = score_prime(uniq, peptide_col="mutant_peptide", hla_col="hla_allele").scored
    out = res.drop_duplicates(["mutant_peptide", "hla_allele"])[["mutant_peptide", "hla_allele", "prime_rank"]]
    out.to_csv(cache, sep="\t", index=False)
    return out


# --------------------------------------------------------------------------------------------------
# Cohort prep -> unified columns: patient_id, mutant_peptide, label, prime, el, expr (+ comparators).
# --------------------------------------------------------------------------------------------------
def _multimer() -> pd.DataFrame:
    from event_b.cd8_multimer_corpus import load_cd8_multimer

    m = load_cd8_multimer().frame.copy().reset_index(drop=True)
    pr = _score_prime_pairs(m[["mutant_peptide", "hla_allele"]], "_cache_multimer_prime2.tsv")
    m = m.merge(pr, on=["mutant_peptide", "hla_allele"], how="left")
    m = m[m["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    return pd.DataFrame({
        "patient_id": m["patient_id"].astype(str),
        "mutant_peptide": m["mutant_peptide"].astype(str),
        "hla_allele": m["hla_allele"].astype(str),
        "label": m["label"].to_numpy(),
        "prime": pd.to_numeric(m["prime_rank"], errors="coerce"),
        "el": pd.to_numeric(m["EL (%Rank score)"], errors="coerce"),
        "expr": pd.to_numeric(m["RNA expression (TPM)"], errors="coerce"),
    })


def _gartner() -> pd.DataFrame:
    from event_b.gartner_nci_corpus import load_gartner_nci

    g = load_gartner_nci().frame.copy().reset_index(drop=True)
    hla = pd.read_csv(HLA_FILE, sep="\t", dtype=str)
    hla["patient_num"] = hla["Patient"].astype(str)
    g["patient_num"] = g["patient_id"].astype(str).str.extract(r"(\d+)\s*$")[0]
    g = g.merge(hla[["patient_num", "Alleles"]], on="patient_num", how="left")
    cache = CACHE / "_cache_gartner_prime.tsv"
    if cache.exists():
        g = g.merge(pd.read_csv(cache, sep="\t"), on="candidate_id", how="left")
    else:
        sc = score_mutations_with_prime(g.assign(hla_allele=g["Alleles"]), peptide_col="mutant_peptide",
                                        allele_col="hla_allele", default_mut_pos0=DEFAULT_MUT_POS0,
                                        rule=PRIMARY_RULE)
        g["prime_rank"] = pd.to_numeric(sc["prime_rank"], errors="coerce").to_numpy()
        g[["candidate_id", "prime_rank"]].to_csv(cache, sep="\t", index=False)
    g["prime"] = pd.to_numeric(g["prime_rank"], errors="coerce")
    g["el"] = pd.to_numeric(g["netmhcpan_el_rank"], errors="coerce")
    g["expr"] = pd.to_numeric(g["expr_decile"], errors="coerce")
    g["cmp_mhcflurry"] = -pd.to_numeric(g["mhcflurry_rank"], errors="coerce")
    g["cmp_pres_ens"] = _pct_mean(g, {"netmhcpan_el_rank": False, "netmhcpan_ba_rank": False,
                                      "mixmhcpred_rank": False, "mhcflurry_rank": False,
                                      "hlathena_rank": False})
    return g


def _improve() -> pd.DataFrame:
    with ZipFile(IMPROVE_ZIP) as z, z.open(IMPROVE_MEMBER) as fh:
        imp = pd.read_csv(fh, sep="\t")
    imp = imp[imp["response"].isin([0, 1])].copy().reset_index(drop=True)
    imp["patient_id"] = "improve:" + imp["Patient"].astype(str)
    imp["mutant_peptide"] = imp["Mut_peptide"].astype(str)
    imp["hla_allele"] = imp["HLA_allele"].astype(str)
    pr = _score_prime_pairs(imp[["mutant_peptide", "hla_allele"]], "_cache_improve_prime.tsv")
    imp = imp.merge(pr, on=["mutant_peptide", "hla_allele"], how="left")
    out = pd.DataFrame({
        "patient_id": imp["patient_id"],
        "mutant_peptide": imp["mutant_peptide"],
        "hla_allele": imp["hla_allele"],
        "label": np.where(imp["response"] == 1, "POSITIVE", "TESTED_NEGATIVE"),
        "prime": pd.to_numeric(imp["prime_rank"], errors="coerce"),
        "el": pd.to_numeric(imp["RankEL"], errors="coerce"),
        "expr": pd.to_numeric(imp["Expression"], errors="coerce"),
    })
    out["cmp_pres_ens"] = _pct_mean(imp.assign(_el=imp["RankEL"], _ba=imp["RankBA"], patient_id=imp["patient_id"]),
                                    {"_el": False, "_ba": False})
    return out


# --------------------------------------------------------------------------------------------------
def _pct(frame: pd.DataFrame, col: str, higher_better: bool) -> np.ndarray:
    v = pd.to_numeric(frame[col], errors="coerce")
    if not higher_better:
        v = -v
    return v.groupby(frame["patient_id"]).rank(pct=True).fillna(0.5).to_numpy()


def _pct_mean(frame: pd.DataFrame, spec: dict) -> np.ndarray:
    return np.mean(np.vstack([_pct(frame, c, hi) for c, hi in spec.items()]), axis=0)


def _feature_matrix(frame: pd.DataFrame) -> np.ndarray:
    return np.column_stack([_pct(frame, col, hi) for col, hi in FEATS])


def _leakage_exclude(train: pd.DataFrame, eval_peptides: set[str]) -> tuple[pd.DataFrame, dict]:
    canon = [canonical_peptide(p) for p in train["mutant_peptide"].astype(str)]
    eindex = _kmer_index({canonical_peptide(p) for p in eval_peptides} - {""})
    near_e = np.array([bool(c) and (near_duplicate(c, eindex, threshold=0.8) is not None) for c in canon])
    prime_leak = np.asarray(prime_leakage_mask(train["mutant_peptide"].astype(str).tolist(), near=True), dtype=bool)
    drop = near_e | prime_leak
    stats = {"near_or_exact_eval": int(near_e.sum()), "prime_training_leak": int(prime_leak.sum()),
             "dropped_total": int(drop.sum()), "kept": int((~drop).sum())}
    return train[~drop].copy(), stats


def _arm_metrics(frame: pd.DataFrame, arm: str) -> dict:
    hits = hits_at_k(frame, group_col="patient_id", score_col=arm, label_col="label", k=K, ascending=False)
    got = tot = 0
    for _, gp in frame.groupby("patient_id"):
        r = gp.sort_values(arm, ascending=False, kind="mergesort")
        got += int((r["label"].to_numpy()[:K] == "POSITIVE").sum())
        tot += int((gp["label"].to_numpy() == "POSITIVE").sum())
    tested = frame[frame["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])]
    yt = (tested["label"] == "POSITIVE").astype(int).to_numpy()
    au = fail_closed_auroc(yt, tested[arm].to_numpy(float), tested_negative_available=True)
    return {"mean_hits@20": round(float(np.mean(hits)), 4), "recall_top20": f"{got}/{tot}",
            "recall_frac": round(got / tot, 4) if tot else None,
            "tested_auroc": au.get("auroc"), "tested_ap": au.get("ap")}


def _evaluate(frame: pd.DataFrame, clf, mu, sd, extra_arms: list[str], _prescored: bool = False) -> dict:
    """Score epicurus_residual + comparators on one frozen eval cohort; metrics + paired bootstrap."""
    frame = frame.copy()
    if not _prescored:
        frame["epicurus_residual"] = clf.decision_function((_feature_matrix(frame) - mu) / sd)
    frame["genuine_prime"] = -frame["prime"]
    frame["netmhcpan_el"] = -frame["el"]
    arms = ["epicurus_residual", "genuine_prime", "netmhcpan_el", *extra_arms]
    for a in arms:
        frame[a] = pd.to_numeric(frame[a], errors="coerce")
        frame[a] = frame[a].fillna(frame[a].min() - 1.0 if frame[a].notna().any() else -1e9)
    per_arm = {a: _arm_metrics(frame, a) for a in arms}
    paired = {}
    for a in arms:
        if a == "epicurus_residual":
            continue
        c = _card(frame, "epicurus_residual", K, baseline=a)
        paired[a] = {"residual_minus_arm_hits@20": c["hits_delta"], "delta_ci": c["delta_ci"],
                     "p_residual_better": c["p_better"], "verdict": c["verdict"]}
    return {"n_candidates": int(len(frame)), "n_patients": int(frame["patient_id"].nunique()),
            "label_counts": frame["label"].value_counts().to_dict(),
            "per_arm_metrics": per_arm, "residual_vs_comparators_top20_paired_bootstrap": paired}


FROZEN_SPEC = Path("configs/frozen/epicurus_v0_1.json")


def load_frozen_spec() -> dict:
    return json.loads(FROZEN_SPEC.read_text())


def score_with_frozen(frame: pd.DataFrame, spec: dict | None = None) -> np.ndarray:
    """Apply the IMMUTABLE frozen Epicurus formula (configs/frozen/epicurus_v0_1.json) to a frame with
    columns patient_id, prime, el, expr. Deterministic — no retraining. Higher = better. Use this to
    validate the frozen formula on any newly acquired cohort mapped to the unified schema."""
    f = (spec or load_frozen_spec())["formula"]
    orient = f["orientation"]
    X = np.column_stack([_pct(frame, c, orient[c] == "higher_raw_better") for c in f["features"]])
    Z = (X - np.asarray(f["standardizer_mean"])) / np.asarray(f["standardizer_std"])
    coef = np.asarray([f["coefficients"][c] for c in f["features"]])
    return Z @ coef + float(f["intercept"])


def external_validate(frame: pd.DataFrame, *, extra_arms: list[str] | None = None) -> dict:
    """Run the frozen Epicurus formula + genuine PRIME + NetMHCpan-EL on a unified cohort frame
    (patient_id, hla_allele, label, prime, el, expr [+ any cmp_* comparator columns]) and return
    patient top-20 / recall / tested AUROC-AP + paired bootstrap. One call per newly acquired cohort."""
    frame = frame.copy()
    frame["epicurus_residual"] = score_with_frozen(frame)
    return _evaluate(frame, None, None, None, extra_arms=extra_arms or [], _prescored=True)


def fit_frozen_residual() -> dict:
    """Train the ONE frozen residual on the multimer cohort (leakage-clean vs BOTH eval cohorts)."""
    g = _gartner()
    imp = _improve()
    mm_raw = _multimer()
    eval_peptides = set(g["mutant_peptide"].astype(str).str.upper()) | set(imp["mutant_peptide"].astype(str).str.upper())
    mm, leak = _leakage_exclude(mm_raw, eval_peptides)
    Xtr = _feature_matrix(mm)
    ytr = (mm["label"] == "POSITIVE").astype(int).to_numpy()
    mu, sd = Xtr.mean(axis=0), Xtr.std(axis=0)
    sd = np.where(sd == 0, 1.0, sd)
    clf = LogisticRegression(max_iter=2000, class_weight="balanced").fit((Xtr - mu) / sd, ytr)
    meta = {
        "training_cohort": "cd8_multimer (independent; measured POSITIVE vs TESTED_NEGATIVE)",
        "features": [f[0] for f in FEATS],
        "orientation": {f[0]: ("higher_raw_better" if f[1] else "lower_raw_better") for f in FEATS},
        "coefficients": {f[0]: round(float(c), 5) for f, c in zip(FEATS, clf.coef_[0])},
        "intercept": round(float(clf.intercept_[0]), 5),
        "standardizer_mean": [round(float(x), 5) for x in mu],
        "standardizer_std": [round(float(x), 5) for x in sd],
        "train_positives": int(ytr.sum()), "train_tested_negatives": int((ytr == 0).sum()),
        "leakage_excluded_vs_eval": leak,
        "synthetic_negatives_used_as_data": False,
    }
    return {"clf": clf, "mu": mu, "sd": sd, "meta": meta, "gartner": g, "improve": imp}


def run_transfer_and_validation() -> dict:
    fit = fit_frozen_residual()
    clf, mu, sd = fit["clf"], fit["mu"], fit["sd"]

    # existing Epicurus full_model (DEV: Gartner-CV, optimistic) — Gartner comparator only
    from event_b.decision_benchmark import run_gartner_benchmark
    _, gb = run_gartner_benchmark()
    g = fit["gartner"].merge(
        gb[["candidate_id", "full_model"]].rename(columns={"full_model": "cmp_full_model_dev"}),
        on="candidate_id", how="left")

    gartner_eval = _evaluate(g, clf, mu, sd, extra_arms=["cmp_mhcflurry", "cmp_pres_ens", "cmp_full_model_dev"])
    improve_eval = _evaluate(fit["improve"], clf, mu, sd, extra_arms=["cmp_pres_ens"])

    def _headline(ev, name):
        p = ev["residual_vs_comparators_top20_paired_bootstrap"]["genuine_prime"]
        return (f"{name}: epicurus_residual vs genuine PRIME hits@20 Δ={p['residual_minus_arm_hits@20']:+.3f} "
                f"(verdict {p['verdict']})")

    return {
        "status": "executed",
        "frozen_model": fit["meta"],
        "design": "ONE residual trained on multimer, applied frozen to two untouched external cohorts.",
        "step2_transfer_test_gartner": gartner_eval,
        "step4_external_validation_improve": improve_eval,
        "headline": {
            "gartner": _headline(gartner_eval, "Gartner (transfer test)"),
            "improve": _headline(improve_eval, "IMPROVE (external validation)"),
            "north_star_status": "Epicurus does NOT yet beat genuine PRIME on untouched external patients "
                                 "unless BOTH cohorts show ACCEPT vs genuine_prime.",
        },
        "caveats": [
            "Gartner 46 pos/26 pts; IMPROVE 467 pos/70 pts. IMPROVE is far better powered.",
            "PRIME is genuine GfellerLab PRIME %rank (lower=better) computed by the installed tool on ALL "
            "three cohorts — no reliance on any precomputed/oriented PRIME column. A MixMHCpred/EL baseline "
            "is never labeled PRIME.",
            "cmp_full_model_dev was trained on Gartner (CV) -> optimistic DEVELOPMENT reference only.",
            "VAF not in the frozen formula (multimer training cohort lacks it).",
            "IMPROVE is a screened Event-A candidate set (~250 tested/patient), not the raw somatic universe; "
            "top-20 is over that tested universe (genuine tested negatives), a legitimate within-patient test.",
        ],
    }
