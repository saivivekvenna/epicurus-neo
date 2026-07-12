"""Rich-feature dynamic gate (Milestone 7) — IMPROVE 88-column table.

The v1/v2 gates used only the 7-column pool export (prime/el/expr) and hit a wall. The full IMPROVE table
(data/raw/improve/data.zip member .../03_3_final_peptide_features_Partition_for_CV.txt) carries the
orthogonal WES/RNA evidence that was missing: DNA VAF (VarAlFreq), mutant-allele RNA (rna_af/var/total,
ValMutRNACoef, rna_bin), stability, DAI / WT-rank agretopicity (RankEL_wt), foreignness/self-similarity,
physchem, mutation class. Patient-relative univariate AUC within the top-40 threat zone reaches ~0.64
(RankEL_wt) — real orthogonal signal exactly where reselection operates.

This module builds a NEGATIVE/POSITIVE-risk learner on the DECISION BOUNDARY (frozen-Epicurus top-20 +
ranks 21-60 + all positives; equal weight per patient) and a constrained COUNTERFACTUAL SWAP policy that
feeds SURVIVORS to the UNCHANGED frozen Epicurus ranker (it never replaces the ranker). Objective = NET
recognized hits@20. Leakage discipline: label/identity/TME/pipeline-score columns are excluded; missing
values get explicit indicators (never silently imputed as signal); features are within-patient percentiles
computed per row (fold-safe).
"""

from __future__ import annotations

import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from event_b.dynamic_gate import within_patient_percentile
from event_b.pool_size_sensitivity import score_arms

IMPROVE_ZIP = Path("data/raw/improve/data.zip")
IMPROVE_MEMBER = "data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt"
PRIME_CACHE = Path("data/raw/gartner_nci/_cache_improve_prime.tsv")
K = 20

# ------------------------------------------------------------------ feature families (safe-initial) ---
# Candidate-VARYING, deployable, pre-outcome. Grouped for ablation.
FAMILIES: dict[str, list[str]] = {
    "wt_rank": ["RankEL_wt", "RankEL_wt_4.1"],                 # WT presentation -> mut-vs-WT differential
    "dai": ["DAI", "DAI_4.1"],                                 # differential agretopicity index
    "rna": ["rna_af", "rna_var", "rna_total", "ValMutRNACoef", "rna_bin"],  # mutant-allele RNA evidence
    "dna_vaf": ["VarAlFreq"],                                  # tumor DNA VAF
    "stability": ["Stability"],                                # pMHC stability (half-life)
    "foreign": ["Foreigness", "SelfSim"],                      # foreignness / self-similarity
    "physchem": ["PropHydroAro", "HydroAll", "HydroCore", "PropSmall", "PropAro",
                 "PropBasic", "PropAcidic", "mw", "Aro", "Inst", "pI"],
    "mutclass": ["Cancer_Driver_Gene", "Misense_mutation", "Frameshift_mutation",
                 "Inframe_deletion_mutation", "Inframe_insertion", "PeptLen"],
}
SAFE_FAMILIES = list(FAMILIES.keys())

# EXCLUDED — label/identity/pipeline/TME. Documented so the audit is explicit.
EXCLUDED_LEAKAGE = {
    "PrioScore": "pipeline prioritization meta-score (selection/label-adjacent) — audit required",
    "CelPrev": "cellular prevalence/clonality — plausibly deployable but audit required before use",
    "IB_CB": "unknown semantics (IB/CB) — audit required", "IB_CB_cat": "categorical of IB_CB",
    "NetMHCExp": "NetMHC x expression composite — audit required",
    "validation": "assay validation status — LABEL-ADJACENT", "response": "OUTCOME LABEL",
    "pMHC": "identity string", "norm_pMHC": "identity string",
}
EXCLUDED_TME = ["Tcells", "TcellsCD8", "CytoxLympho", "Blinage", "NKcells", "Monocytes", "MyeloidDC",
                "Neutrophils", "Endothelial", "Fibroblasts", "MCPmean", "CYT"]  # patient-constant => useless within-patient
EXCLUDED_PRESENTATION = ["RankEL", "RankBA", "RankEL_4.1", "RankBA_4.1", "Prime", "Expression"]  # the ranker's own signals


def family_columns(families: list[str]) -> list[str]:
    return [c for fam in families for c in FAMILIES[fam]]


# ------------------------------------------------------------------------------------- data loading ---
def load_improve_rich() -> pd.DataFrame:
    with zipfile.ZipFile(IMPROVE_ZIP) as z, z.open(IMPROVE_MEMBER) as fh:
        df = pd.read_csv(fh, sep="\t")
    df = df[df["response"].isin([0, 1])].copy().reset_index(drop=True)
    prime = pd.read_csv(PRIME_CACHE, sep="\t")
    df = df.merge(prime.rename(columns={"mutant_peptide": "Mut_peptide", "hla_allele": "HLA_allele"}),
                  on=["Mut_peptide", "HLA_allele"], how="left")
    out = pd.DataFrame({
        "patient_id": "improve:" + df["Patient"].astype(str),
        "mutant_peptide": df["Mut_peptide"].astype(str),
        "hla_allele": df["HLA_allele"].astype(str),
        "label": np.where(df["response"] == 1, "POSITIVE", "TESTED_NEGATIVE"),
        "prime": pd.to_numeric(df["prime_rank"], errors="coerce"),
        "el": pd.to_numeric(df["RankEL"], errors="coerce"),
        "expr": pd.to_numeric(df["Expression"], errors="coerce"),
        "partition": df["Partition"].astype(int),
        "source": df["cohort"].astype(str),  # bladder / melanoma / Basket -> leave-source-out transport
    })
    for c in family_columns(SAFE_FAMILIES):
        out[c] = pd.to_numeric(df[c], errors="coerce") if c in df else np.nan
    return out


# ------------------------------------------------------------- feature matrix (within-patient pct) ---
def feature_matrix(frame: pd.DataFrame, cols: list[str]) -> tuple[np.ndarray, list[str]]:
    """Within-patient percentile of each feature (patient-relative; fold-safe) + explicit missing
    indicator per feature. Missing percentile imputed to 0.5 but the indicator preserves the fact."""
    blocks, names = [], []
    for c in cols:
        p = within_patient_percentile(frame, c, higher_better=True)  # sign learned by the model
        miss = np.isnan(p)
        blocks.append(np.where(miss, 0.5, p))
        blocks.append(miss.astype(float))
        names += [f"{c}__pct", f"{c}__missing"]
    return np.column_stack(blocks), names


# ------------------------------------------------------------------------- decision-boundary subset ---
def decision_boundary_mask(frame: pd.DataFrame, top: int = K, extra: int = 40) -> np.ndarray:
    """Per patient: frozen-Epicurus top-`top` + next `extra` ranks + ALL positives. Focuses the learner on
    hard boundary negatives instead of 17k easy ones."""
    keep = np.zeros(len(frame), bool)
    scored = score_arms(frame)
    for _, gp in scored.groupby("patient_id"):
        local = frame.index.get_indexer(gp.index)
        order = np.argsort(-gp["frozen_epicurus"].to_numpy(), kind="mergesort")
        keep[local[order[:top + extra]]] = True
        keep[local[(gp["label"] == "POSITIVE").to_numpy()]] = True
    return keep


def patient_equal_weights(frame: pd.DataFrame) -> np.ndarray:
    """Weight rows so each patient contributes equal TOTAL weight (no large-pool domination)."""
    w = np.ones(len(frame))
    for _, idx in frame.groupby("patient_id").groups.items():
        rows = frame.index.get_indexer(idx)
        w[rows] = 1.0 / len(rows)
    return w


def balanced_weights(frame: pd.DataFrame, balance_class: bool = True) -> np.ndarray:
    """Equal TOTAL weight per patient, and (if balance_class) equal TOTAL weight to positives vs
    tested-negatives — otherwise a ~3% positive rate lets HistGBT collapse toward all-negative."""
    w = patient_equal_weights(frame)
    if balance_class:
        y = (frame["label"].to_numpy() == "POSITIVE")
        wp, wn = w[y].sum(), w[~y].sum()
        if wp > 0 and wn > 0:
            w = np.where(y, w * (0.5 / wp), w * (0.5 / wn))
    return w


# ------------------------------------------------------------------------------------ risk learners ---
@dataclass
class RiskLearner:
    kind: str
    model: object
    cols: list[str]
    names: list[str]

    def q(self, frame: pd.DataFrame) -> np.ndarray:
        """q = P(POSITIVE) per candidate (higher = more likely recognized)."""
        X, _ = feature_matrix(frame, self.cols)
        if self.kind == "histgbt":
            return self.model.predict_proba(X)[:, 1]
        # pairwise-residual: decision_function on standardized features -> logistic
        mu, sd = self.model["mu"], self.model["sd"]
        return _sigmoid((X - mu) / sd @ self.model["w"])


def _sigmoid(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -30, 30)))


def fit_histgbt(train: pd.DataFrame, families: list[str], seed: int = 0,
                balance_class: bool = True) -> RiskLearner:
    from sklearn.ensemble import HistGradientBoostingClassifier

    cols = family_columns(families)
    X, names = feature_matrix(train, cols)
    y = (train["label"].to_numpy() == "POSITIVE").astype(int)
    w = balanced_weights(train, balance_class=balance_class)
    clf = HistGradientBoostingClassifier(max_depth=3, max_iter=200, learning_rate=0.05,
                                         l2_regularization=1.0, random_state=seed, early_stopping=False)
    clf.fit(X, y, sample_weight=w)
    return RiskLearner("histgbt", clf, cols, names)


def fit_pairwise(train: pd.DataFrame, families: list[str]) -> RiskLearner:
    """Within-patient positive-vs-blocking-negative feature-difference logistic (RankSVM-style residual).
    A 'blocking negative' is any tested-negative ranked at/above a positive by frozen Epicurus."""
    from sklearn.linear_model import LogisticRegression

    cols = family_columns(families)
    X, names = feature_matrix(train, cols)
    mu, sd = X.mean(0), X.std(0)
    sd = np.where(sd == 0, 1.0, sd)
    Z = (X - mu) / sd
    scored = score_arms(train)
    diffs, ys = [], []
    for _, gp in scored.groupby("patient_id"):
        frame_rows = train.index.get_indexer(gp.index)
        fe = gp["frozen_epicurus"].to_numpy()
        lab = gp["label"].to_numpy()
        pos = np.where(lab == "POSITIVE")[0]
        neg = np.where(lab == "TESTED_NEGATIVE")[0]
        for pi in pos:
            blockers = neg[fe[neg] >= fe[pi]]  # negatives ranked at/above this positive
            if len(blockers) > 40:
                blockers = blockers[np.argsort(-fe[blockers])[:40]]
            for bi in blockers:
                d = Z[frame_rows[pi]] - Z[frame_rows[bi]]
                diffs.append(d)
                ys.append(1)
                diffs.append(-d)  # symmetric
                ys.append(0)
    if not diffs:
        w = np.zeros(Z.shape[1])
    else:
        clf = LogisticRegression(max_iter=2000, C=1.0, fit_intercept=False).fit(np.array(diffs), np.array(ys))
        w = clf.coef_[0]
    return RiskLearner("pairwise", {"w": w, "mu": mu, "sd": sd}, cols, names)


# ------------------------------------------------------------ base-anchored residual utility (v2) ---
def frozen_base_score(frame: pd.DataFrame) -> np.ndarray:
    """Frozen-Epicurus score per candidate (the UNCHANGED ranker's own utility). Higher = better."""
    return score_arms(frame)["frozen_epicurus"].to_numpy()


def within_patient_z(frame: pd.DataFrame, values: np.ndarray) -> np.ndarray:
    """Standardize `values` within each patient (mean 0 / sd 1). Patient-relative, fold-safe."""
    z = np.zeros(len(values), float)
    for _, idx in frame.groupby("patient_id").groups.items():
        rows = frame.index.get_indexer(idx)
        v = values[rows]
        sd = v.std()
        z[rows] = (v - v.mean()) / (sd if sd > 0 else 1.0)
    return z


_HYDRO_ARO = set("AVILMFWYC")  # hydrophobic + aromatic residues (proxy for IMPROVE PropHydroAro)


def peptide_hydro_aro_fraction(peptides: pd.Series) -> np.ndarray:
    """Fraction of hydrophobic/aromatic residues per peptide — a peptide-DERIVED proxy for PropHydroAro,
    computable for ANY cohort (used to transfer the IMPROVE-frozen policy to external pools)."""
    out = []
    for p in peptides.astype(str):
        aa = [c for c in p.upper() if c.isalpha()]
        out.append(sum(c in _HYDRO_ARO for c in aa) / len(aa) if aa else np.nan)
    return np.array(out, float)


def feature_percentile(frame: pd.DataFrame, col: str, higher_better: bool = True) -> np.ndarray:
    """Within-patient percentile of one raw feature; missing -> 0.5 (neutral -> no base push -> abstains)."""
    p = within_patient_percentile(frame, col, higher_better=higher_better)
    return np.where(np.isnan(p), 0.5, p)


def base_anchored_hits(frame: pd.DataFrame, feat_pct: np.ndarray, alpha: float, threat: int = 60,
                       return_swaps: bool = False):
    """BASE-ANCHORED selection. Per patient: rank the top-`threat` candidates by frozen-Epicurus base
    percentile, re-order that set by U = base_pct + alpha*(feat_pct-0.5), take top-20. Epicurus is the
    UNCHANGED base anchor; the feature only NUDGES. alpha=0 => U==base => Epicurus top-20 (no-op).
    Returns {patient: hits@20}; with return_swaps also per-patient (pos/neg swapped in and out)."""
    base = frozen_base_score(frame)
    lab = frame["label"].to_numpy()
    hits, swaps = {}, {}
    for pid, gidx in frame.groupby("patient_id").groups.items():
        rows = frame.index.get_indexer(gidx)
        bp = pd.Series(base[rows]).rank(pct=True).to_numpy()
        fp = feat_pct[rows]
        thr = np.argsort(-bp, kind="mergesort")[:threat]                 # top-`threat` by base
        base_top = set(thr[np.argsort(-bp[thr], kind="mergesort")[:K]])  # Epicurus top-20
        U = bp + alpha * (fp - 0.5)
        new_top = set(thr[np.argsort(-U[thr], kind="mergesort")[:K]])    # base-anchored top-20
        labr = lab[rows]
        hits[str(pid)] = float(sum(labr[i] == "POSITIVE" for i in new_top))
        if return_swaps:
            swapped_in = new_top - base_top
            swapped_out = base_top - new_top
            swaps[str(pid)] = {
                "n_swaps": len(swapped_in),
                "pos_in": int(sum(labr[i] == "POSITIVE" for i in swapped_in)),
                "neg_in": int(sum(labr[i] == "TESTED_NEGATIVE" for i in swapped_in)),
                "pos_out": int(sum(labr[i] == "POSITIVE" for i in swapped_out)),
                "neg_out": int(sum(labr[i] == "TESTED_NEGATIVE" for i in swapped_out))}
    return (hits, swaps) if return_swaps else hits


def additive_utility(frame: pd.DataFrame, resid_q: np.ndarray, alpha: float) -> np.ndarray:
    """BASE-ANCHORED utility = patient-relative frozen-Epicurus base + alpha * patient-relative residual.
    alpha=0 => pure base => the swap policy makes ZERO swaps (removing the lowest-base top-20 and
    backfilling a lower-base rank-21 can never raise base utility). So any swap is driven by the residual,
    and only when it overcomes the base rank gap."""
    bz = within_patient_z(frame, frozen_base_score(frame))
    rz = within_patient_z(frame, np.asarray(resid_q, float))
    return bz + alpha * rz


# ---------------------------------------------------------------------- counterfactual swap policy ---
def counterfactual_swaps(frame: pd.DataFrame, q: np.ndarray, *, max_budget: int = 8, margin: float = 0.0,
                         fixed_m: int | None = None) -> np.ndarray:
    """Constrained top-20 replacement: per patient, remove the m lowest-q of the frozen top-20 while the
    replacements (ranks 21..20+m) have higher summed q than the removed (net expected-hit gain > margin).
    `fixed_m` forces a fixed budget (control). Returns a keep mask; survivors go to UNCHANGED Epicurus."""
    keep = np.ones(len(frame), bool)
    scored = score_arms(frame)
    for _, gp in scored.groupby("patient_id"):
        local = frame.index.get_indexer(gp.index)
        order = np.argsort(-gp["frozen_epicurus"].to_numpy(), kind="mergesort")
        ranked = local[order]
        top = ranked[:K]
        backfill = ranked[K:K + max_budget]
        if len(top) == 0 or len(backfill) == 0:
            continue
        rem_order = top[np.argsort(q[top], kind="mergesort")]  # lowest-q first
        if fixed_m is not None:
            best_m = min(fixed_m, len(backfill), len(top))
        else:
            best_m = 0
            for m in range(1, min(max_budget, len(backfill), len(top)) + 1):
                gain = float(q[backfill[:m]].sum() - q[rem_order[:m]].sum())
                if gain > margin:
                    best_m = m
                else:
                    break
        if best_m > 0:
            keep[rem_order[:best_m]] = False
    return keep


def patient_hits20(frame: pd.DataFrame, keep: np.ndarray) -> dict[str, float]:
    out = {}
    sub = frame[keep]
    for pid, gp in sub.groupby("patient_id"):
        r = score_arms(gp).sort_values("frozen_epicurus", ascending=False, kind="mergesort")
        out[str(pid)] = float((r["label"].to_numpy()[:K] == "POSITIVE").sum())
    for pid in frame["patient_id"].astype(str).unique():
        out.setdefault(pid, 0.0)
    return out
