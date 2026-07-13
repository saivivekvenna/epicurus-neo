"""Epicurus v0.3 DEVELOPMENT experiment — MIL-aware ranker on the frozen mil_dev_split_v1.

DEVELOPMENT ONLY. Follows `artifacts/milestone_7_decision/epicurus_v03/PREREGISTERED_PROTOCOL.md` exactly. The
Gartner TEST holdout is never loaded/scored here; the frozen split (patient->fold + recurrent-peptide
quarantine) is used verbatim. No external-superiority claim is produced.

Contents: feature assembly (§7), the model ladder (§6: presentation baseline / additive logistic / MIL
log-sum-exp ranker), and the nested out-of-fold evaluation + paired patient-level bootstrap (§4-5).
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize

from event_b.leakage_registry import canonical_peptide
from event_b.nci_crosswalk import build_train_crosswalk
from event_b.prime_training import prime_leakage_mask

FROZEN_SPLIT = Path("configs/frozen/mil_dev_split_v1.json")
G_PRIME = Path("data/raw/gartner_nci/_cache_gartner_muller_prime.tsv")
IMP_PRIME = Path("data/raw/gartner_nci/_cache_improve_primemix.tsv")
MM_PRIME = Path("data/raw/gartner_nci/_cache_multimer_primemix.tsv")
IMPROVE_ZIP = Path("data/raw/improve/data.zip")
IMPROVE_MEMBER = "data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt"

CORE = ["f_prime_pct", "f_mix_pct", "f_el_pct", "f_pres_abs"]     # present in every source (fail-closed)
ORTHO = ["f_expr", "f_agreto", "f_foreign", "f_bindstab", "f_proc"]  # source-specific, masked-neutral
FEATURES = CORE + ORTHO
_POS_BAG_INSTANCE = {"POSITIVE_EXACT", "AMBIGUOUS_POSITIVE_BAG"}
K_TOP = 20


# --------------------------------------------------------------------------------------------------
# frozen split
# --------------------------------------------------------------------------------------------------
def load_frozen() -> tuple[dict, set, int]:
    d = json.loads(FROZEN_SPLIT.read_text())
    return {p: int(f) for p, f in d["patient_fold"].items()}, set(d["quarantined_recurrent_peptides"]), int(d["k"])


# --------------------------------------------------------------------------------------------------
# per-source raw feature loaders -> a common partial frame
# --------------------------------------------------------------------------------------------------
def _gartner_rows() -> pd.DataFrame:
    cw = build_train_crosswalk().instances
    pr = pd.read_csv(G_PRIME, sep="\t")
    cw = cw.merge(pr, on=["peptide", "hla_allele"], how="left")
    bag_unit = cw["resolved_parent_key"].fillna(cw["candidate_parent_keys"]).astype(str)
    bag_label = np.where(cw["instance_label"].isin(_POS_BAG_INSTANCE), "POSITIVE",
                np.where(cw["instance_label"] == "NEGATIVE_BAG_CHILD", "NEGATIVE", "EXCLUDE"))
    return pd.DataFrame({
        "source": "gartner", "patient_id": "gartner:" + cw["patient_id"].astype(str),
        "bag_id": "gartner:" + cw["patient_id"].astype(str) + "#" + bag_unit,
        "peptide": cw["peptide"].astype(str), "hla_allele": cw["hla_allele"].astype(str),
        "bag_label": bag_label, "eval_positive": (cw["instance_label"] == "POSITIVE_EXACT").astype(int),
        "prime_rank": pd.to_numeric(cw["prime_rank"], errors="coerce"),
        "mix_rank": pd.to_numeric(cw["mix_rank"], errors="coerce"),
        "el_score": pd.to_numeric(cw["score_el"], errors="coerce"),   # 0-1 likelihood, HIGHER better
        "el_percent_rank": np.nan,
        "el_oriented": pd.to_numeric(cw["score_el"], errors="coerce"),  # higher=better, full resolution
        "el_strength": _strength_from_score(cw["score_el"]),
        "f_expr": np.nan, "f_agreto": pd.to_numeric(cw["agretopicity"], errors="coerce"),
        "f_foreign": np.nan, "f_bindstab": pd.to_numeric(cw["bind_stab"], errors="coerce"), "f_proc": np.nan,
    })


def _multimer_rows() -> pd.DataFrame:
    from event_b.cd8_multimer_corpus import load_cd8_multimer
    m = load_cd8_multimer().frame.copy()
    m = m[m["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    pr = pd.read_csv(MM_PRIME, sep="\t").rename(columns={"peptide": "mutant_peptide"})
    m = m.merge(pr, on=["mutant_peptide", "hla_allele"], how="left")
    pid = "multimer:" + m["patient_id"].astype(str)
    return pd.DataFrame({
        "source": "multimer", "patient_id": pid,
        "bag_id": pid + "#" + m["mutant_peptide"].astype(str) + "#" + m["hla_allele"].astype(str),
        "peptide": m["mutant_peptide"].astype(str), "hla_allele": m["hla_allele"].astype(str),
        "bag_label": np.where(m["label"] == "POSITIVE", "POSITIVE", "NEGATIVE"),
        "eval_positive": (m["label"] == "POSITIVE").astype(int),
        "prime_rank": pd.to_numeric(m["prime_rank"], errors="coerce"),
        "mix_rank": pd.to_numeric(m["mix_rank"], errors="coerce"),
        "el_score": np.nan,
        "el_percent_rank": pd.to_numeric(m["EL (%Rank score)"], errors="coerce"),  # %rank, LOWER better
        "el_oriented": -pd.to_numeric(m["EL (%Rank score)"], errors="coerce"),      # higher=better, full res
        "el_strength": _strength_from_rank(m["EL (%Rank score)"]),
        "f_expr": pd.to_numeric(m["RNA expression (TPM)"], errors="coerce"),
        "f_agreto": pd.to_numeric(m["Agretopicity"], errors="coerce"),
        "f_foreign": pd.to_numeric(m["Foreignness score"], errors="coerce"),
        "f_bindstab": np.nan, "f_proc": pd.to_numeric(m["Proteasomal processing score"], errors="coerce"),
    })


def _improve_rows() -> pd.DataFrame:
    from zipfile import ZipFile
    with ZipFile(IMPROVE_ZIP) as z, z.open(IMPROVE_MEMBER) as fh:
        imp = pd.read_csv(fh, sep="\t")
    imp = imp[imp["response"].isin([0, 1])].copy()
    pr = pd.read_csv(IMP_PRIME, sep="\t")
    imp = imp.merge(pr, left_on=["Mut_peptide", "HLA_allele"], right_on=["peptide", "hla_allele"], how="left")
    pid = "improve:" + imp["Patient"].astype(str)
    return pd.DataFrame({
        "source": "improve", "patient_id": pid,
        "bag_id": pid + "#" + imp["Mut_peptide"].astype(str) + "#" + imp["HLA_allele"].astype(str),
        "peptide": imp["Mut_peptide"].astype(str), "hla_allele": imp["HLA_allele"].astype(str),
        "bag_label": np.where(imp["response"] == 1, "POSITIVE", "NEGATIVE"),
        "eval_positive": (imp["response"] == 1).astype(int),
        "prime_rank": pd.to_numeric(imp["prime_rank"], errors="coerce"),
        "mix_rank": pd.to_numeric(imp["mix_rank"], errors="coerce"),
        "el_score": np.nan,
        "el_percent_rank": pd.to_numeric(imp["RankEL"], errors="coerce"),   # %rank, LOWER better
        "el_oriented": -pd.to_numeric(imp["RankEL"], errors="coerce"),       # higher=better, full res
        "el_strength": _strength_from_rank(imp["RankEL"]),
        "f_expr": pd.to_numeric(imp["Expression"], errors="coerce"),
        "f_agreto": np.nan, "f_foreign": np.nan, "f_bindstab": np.nan, "f_proc": np.nan,
    })


# --------------------------------------------------------------------------------------------------
# EL semantics differ by source and MUST NOT be conflated (audited):
#   * Gartner/Müller `Score_EL`  is a 0-1 EL LIKELIHOOD  -> HIGHER is better.
#   * IMPROVE `RankEL` / multimer `EL (%Rank score)` are %RANKS -> LOWER is better.
# We keep the raw source-specific fields explicit (`el_score` vs `el_percent_rank`) and derive ONE
# consistently-oriented, comparable presentation strength (`el_strength`, higher = better, ~[0,4]).
# --------------------------------------------------------------------------------------------------
def _strength_from_rank(rank: pd.Series) -> np.ndarray:
    """Lower-is-better %rank -> higher-is-better strength: -log10(rank/100), clipped to [0,4]."""
    r = pd.to_numeric(rank, errors="coerce") / 100.0
    return np.clip(-np.log10(np.clip(r, 1e-4, 1.0)), 0, 4)


def _strength_from_score(score: pd.Series) -> np.ndarray:
    """Higher-is-better 0-1 EL likelihood -> comparable strength: -log10(1 - score), clipped to [0,4]."""
    s = pd.to_numeric(score, errors="coerce")
    return np.clip(-np.log10(np.clip(1.0 - s, 1e-4, 1.0)), 0, 4)


# --------------------------------------------------------------------------------------------------
# feature engineering (§7): within-patient percentiles + absolute presentation + masked orthogonals
# --------------------------------------------------------------------------------------------------
def _wp_pct(value: pd.Series, patient: pd.Series, higher_better: bool) -> np.ndarray:
    """Within-patient percentile in [0,1], higher = better. Missing -> NaN (handled per feature)."""
    v = value if higher_better else -value
    return v.groupby(patient).rank(pct=True).to_numpy()


def _centered_ortho(value: pd.Series, patient: pd.Series, higher_better: bool) -> np.ndarray:
    """Orthogonal feature as a within-patient percentile centred to [-0.5, 0.5]; missing -> 0 (neutral abstain,
    inert in a within-patient ranking when the whole pool shares the absence). No present-indicator is exposed
    as a model feature (that would leak)."""
    pct = _wp_pct(value, patient, higher_better)
    out = np.where(np.isnan(pct), 0.0, pct - 0.5)
    return out


def assemble_frame() -> pd.DataFrame:
    patient_fold, quarantine, _ = load_frozen()
    df = pd.concat([_gartner_rows(), _multimer_rows(), _improve_rows()], ignore_index=True)
    df["canonical_peptide"] = [canonical_peptide(p) for p in df["peptide"].astype(str)]

    # fail-closed: a unit missing a core presentation score is dropped, never imputed
    core_ok = df[["prime_rank", "mix_rank", "el_strength"]].notna().all(axis=1)
    df = df[core_ok].copy()

    # attach frozen fold + quarantine (pure functions of patient / canonical peptide)
    df["fold"] = df["patient_id"].map(patient_fold)
    df = df[df["fold"].notna()].copy()
    df["fold"] = df["fold"].astype(int)
    df["quarantined"] = df["canonical_peptide"].isin(quarantine)

    # provenance-leakage guard: mask PRIME feature to neutral on near-PRIME-training peptides (fail-closed)
    pt_leak = np.asarray(prime_leakage_mask(df["peptide"].tolist(), near=True), dtype=bool)
    df["prime_masked"] = pt_leak

    p = df["patient_id"]
    # core within-patient percentiles (lower %rank = better -> higher_better=False on the rank)
    prime_for_model = df["prime_rank"].where(~df["prime_masked"], np.nan)
    df["f_prime_pct"] = _wp_pct(prime_for_model, p, higher_better=False)
    df["f_prime_pct"] = np.where(np.isnan(df["f_prime_pct"]), 0.5, df["f_prime_pct"])  # masked -> neutral
    df["f_mix_pct"] = _wp_pct(df["mix_rank"], p, higher_better=False)
    # el_oriented is source-correct higher=better at FULL resolution (unclipped) -> percentile ascends on it
    df["f_el_pct"] = _wp_pct(df["el_oriented"], p, higher_better=True)
    # absolute presentation strength (cross-patient calibration): mean of source-correct EL strength and the
    # uniform MixMHCpred %rank strength (both higher=better, ~[0,4]). Never applies a %rank transform to Score_EL.
    df["f_pres_abs"] = np.nanmean(np.vstack([df["el_strength"].to_numpy(),
                                             _strength_from_rank(df["mix_rank"])]), axis=0)

    # orthogonal features (masked-neutral, within-patient). orientation fixed a-priori:
    df["f_expr"] = _centered_ortho(df["f_expr"], p, higher_better=True)      # more expression -> better
    df["f_agreto"] = _centered_ortho(df["f_agreto"], p, higher_better=False)  # lower agretopicity -> better
    df["f_foreign"] = _centered_ortho(df["f_foreign"], p, higher_better=True)  # more foreign -> better
    df["f_bindstab"] = _centered_ortho(df["f_bindstab"], p, higher_better=True)  # more stable -> better
    df["f_proc"] = _centered_ortho(df["f_proc"], p, higher_better=True)       # better processing -> better

    return df.reset_index(drop=True)


# --------------------------------------------------------------------------------------------------
# fit weighting (§4): source-balanced x patient-balanced x bag-balanced
# --------------------------------------------------------------------------------------------------
def _fit_weights(train: pd.DataFrame) -> np.ndarray:
    """Positional weight array aligned to `train`'s row order: source-balanced x patient-balanced x
    bag-balanced, with a bag's instances sharing the bag weight (anti-inflation: never per-instance for
    Gartner). Reset internally so positional indexing is valid regardless of the caller's index."""
    t = train.reset_index(drop=True)
    w = np.ones(len(t), dtype=float)
    n_src = t["source"].nunique()
    for src, g in t.groupby("source"):
        n_pat = g["patient_id"].nunique()
        for pid, gp in g.groupby("patient_id"):
            nbag = gp["bag_id"].nunique()
            for bid, gb in gp.groupby("bag_id"):
                w[gb.index.to_numpy()] = 1.0 / (n_src * n_pat * max(nbag, 1) * len(gb))
    return w


# --------------------------------------------------------------------------------------------------
# model ladder (§6)
# --------------------------------------------------------------------------------------------------
@dataclass
class PresentationBaseline:
    """Rung 1: no fitting. Score = within-patient presentation strength (best available EL)."""
    def fit(self, train): return self
    def raw_score(self, df): return df["f_el_pct"].to_numpy()
    def to_dict(self): return {"rung": 1, "kind": "presentation_baseline", "feature": "f_el_pct"}


@dataclass
class AdditiveRanker:
    """Rung 2: L2 logistic on confirmed instances (Gartner POSITIVE_EXACT vs NEGATIVE_BAG_CHILD, bag-weighted;
    IMPROVE/multimer instance labels). Ignores the MIL bag structure (no aggregation)."""
    C: float = 0.3
    mean_: np.ndarray = None
    std_: np.ndarray = None
    coef_: np.ndarray = None
    intercept_: float = 0.0

    def _std(self, X, fit):
        if fit:
            self.mean_ = X.mean(0); self.std_ = X.std(0) + 1e-9
        return (X - self.mean_) / self.std_

    def fit(self, train):
        from sklearn.linear_model import LogisticRegression
        lab = train.copy()
        # confirmed-instance supervision: positives = eval_positive; negatives = Gartner NEG_BAG_CHILD + inst neg
        pos = lab["eval_positive"] == 1
        neg = ((lab["source"] == "gartner") & (lab["bag_label"] == "NEGATIVE")) | \
              ((lab["source"] != "gartner") & (lab["bag_label"] == "NEGATIVE"))
        use = lab[pos | neg].copy()
        y = (use["eval_positive"] == 1).astype(int).to_numpy()
        w = _fit_weights(use)
        X = self._std(use[FEATURES].to_numpy(float), fit=True)
        clf = LogisticRegression(C=self.C, max_iter=2000, class_weight=None)
        clf.fit(X, y, sample_weight=w)
        self.coef_ = clf.coef_.ravel(); self.intercept_ = float(clf.intercept_[0])
        return self

    def raw_score(self, df):
        X = self._std(df[FEATURES].to_numpy(float), fit=False)
        return X @ self.coef_ + self.intercept_

    def to_dict(self):
        return {"rung": 2, "kind": "additive_logistic", "C": self.C,
                "coefficients": dict(zip(FEATURES, [round(float(c), 4) for c in self.coef_]))}


@dataclass
class MILRanker:
    """Rung 3: linear instance scorer f(x)=w·x; Gartner bags aggregate instances by log-sum-exp (temperature
    tau) and are supervised by the BAG label; instance sources are singleton bags (reduces to logistic). L2,
    solved by L-BFGS. Anti-inflation is intrinsic: the loss is summed over BAGS, not instances."""
    C: float = 0.3
    tau: float = 1.0
    mean_: np.ndarray = None
    std_: np.ndarray = None
    coef_: np.ndarray = None
    intercept_: float = 0.0

    def _std(self, X, fit):
        if fit:
            self.mean_ = X.mean(0); self.std_ = X.std(0) + 1e-9
        return (X - self.mean_) / self.std_

    def fit(self, train):
        lab = train[train["bag_label"].isin(["POSITIVE", "NEGATIVE"])].copy().reset_index(drop=True)
        # sort instances so each bag occupies a contiguous segment (enables vectorized segmented log-sum-exp)
        codes = pd.factorize(lab["bag_id"])[0]
        order = np.argsort(codes, kind="stable")
        lab = lab.iloc[order].reset_index(drop=True)
        Xs = self._std(lab[FEATURES].to_numpy(float), fit=True)
        codes_s = codes[order]
        boundaries = np.concatenate([[0], np.flatnonzero(np.diff(codes_s)) + 1])
        counts = np.diff(np.append(boundaries, len(codes_s)))
        yb, wb = self._bag_targets(lab, boundaries)
        tau, n, nbag = self.tau, Xs.shape[1], len(boundaries)
        reg = 1.0 / (self.C * nbag)
        rep = lambda a: np.repeat(a, counts)  # per-bag -> per-instance

        def negll(theta):
            w, b = theta[:n], theta[n]
            s = Xs @ w + b
            seg_max = np.maximum.reduceat(s, boundaries)
            e = np.exp((s - rep(seg_max)) / tau)
            seg_sum = np.add.reduceat(e, boundaries)
            lse = seg_max + tau * np.log(seg_sum / counts)              # soft-max bag score
            p = 1.0 / (1.0 + np.exp(-lse))
            loss = float(np.sum(wb * -(yb * np.log(p + 1e-12) + (1 - yb) * np.log(1 - p + 1e-12))))
            dl_bag = wb * (p - yb)                                      # per bag
            sw = e / rep(seg_sum)                                       # softmax weight per instance
            ds = rep(dl_bag) * sw                                       # per instance
            grad = np.empty(n + 1)
            grad[:n] = Xs.T @ ds + reg * w
            grad[n] = ds.sum()
            return loss + 0.5 * reg * (w @ w), grad

        res = minimize(negll, np.zeros(n + 1), jac=True, method="L-BFGS-B",
                       options={"maxiter": 500, "ftol": 1e-9})
        self.coef_ = res.x[:n]; self.intercept_ = float(res.x[n])
        return self

    @staticmethod
    def _bag_targets(lab_sorted, boundaries):
        """Per-bag label and weight (source-balanced x patient-balanced x bag-balanced). lab_sorted is
        bag-contiguous; boundaries[j] indexes the first instance of bag j."""
        first = lab_sorted.iloc[boundaries]
        yb = (first["bag_label"].to_numpy() == "POSITIVE").astype(float)
        n_src = lab_sorted["source"].nunique()
        npat = lab_sorted.groupby("source")["patient_id"].nunique().to_dict()
        nbag = lab_sorted.groupby("patient_id")["bag_id"].nunique().to_dict()
        wb = np.array([1.0 / (n_src * npat[s] * max(nbag[p], 1))
                       for s, p in zip(first["source"], first["patient_id"])])
        return yb, wb

    def raw_score(self, df):
        X = self._std(df[FEATURES].to_numpy(float), fit=False)
        return X @ self.coef_ + self.intercept_

    def to_dict(self):
        return {"rung": 3, "kind": "MIL_logsumexp", "C": self.C, "tau": self.tau,
                "coefficients": dict(zip(FEATURES, [round(float(c), 4) for c in self.coef_]))}


# baselines as pure score functions (orientation fixed & SOURCE-CORRECT: higher score = better)
def baseline_score(df: pd.DataFrame, which: str) -> np.ndarray:
    if which == "prime":
        return -df["prime_rank"].to_numpy()      # %rank lower better -> negate
    if which == "mix":
        return -df["mix_rank"].to_numpy()        # %rank lower better -> negate
    if which == "presentation":
        # el_oriented is each source's native EL oriented higher=better at full resolution (Gartner Score_EL
        # higher-better; IMPROVE/multimer %rank lower-better) -> correct within-patient ranking, no saturation
        return df["el_oriented"].to_numpy()
    raise ValueError(which)


# --------------------------------------------------------------------------------------------------
# evaluation (§4-5)
# --------------------------------------------------------------------------------------------------
def hits_at_k(scores: np.ndarray, eval_pos: np.ndarray, k: int = K_TOP) -> int:
    order = np.argsort(-scores, kind="stable")
    return int(eval_pos[order[:k]].sum())


def per_patient_metrics(df: pd.DataFrame, score: np.ndarray, k: int = K_TOP) -> pd.DataFrame:
    df = df.assign(_s=score)
    rows = []
    for (src, pid), g in df.groupby(["source", "patient_id"]):
        npos = int(g["eval_positive"].sum())
        if npos == 0:
            continue
        s = g["_s"].to_numpy(); ep = g["eval_positive"].to_numpy()
        h = hits_at_k(s, ep, k)
        rows.append({"source": src, "patient_id": pid, "n_pos": npos, "pool": len(g),
                     "hits": h, "recall": h / npos})
    return pd.DataFrame(rows)


def _source_balanced_weight(mt: pd.DataFrame) -> np.ndarray:
    n_src = mt["source"].nunique()
    npat = mt.groupby("source")["patient_id"].nunique().to_dict()
    return (1.0 / (n_src * mt["source"].map(npat))).to_numpy()


def paired_bootstrap(mt: pd.DataFrame, a: str, b: str, col: str = "hits", n: int = 2000,
                     seed: int = 0) -> dict:
    """Source-balanced paired patient bootstrap of weighted-mean [a - b]."""
    d = (mt[f"{col}_{a}"] - mt[f"{col}_{b}"]).to_numpy()
    w = _source_balanced_weight(mt)
    idx = np.arange(len(mt))
    rng = np.random.default_rng(seed)
    point = float(np.average(d, weights=w))
    boots = np.empty(n)
    # resample patients within source strata (keeps source balance)
    strata = {s: idx[mt["source"].to_numpy() == s] for s in mt["source"].unique()}
    for i in range(n):
        pick = np.concatenate([rng.choice(ix, size=len(ix), replace=True) for ix in strata.values()])
        boots[i] = np.average(d[pick], weights=w[pick])
    lo, hi = np.percentile(boots, [2.5, 97.5])
    return {"delta": round(point, 4), "ci_lo": round(float(lo), 4), "ci_hi": round(float(hi), 4),
            "n_patients": int(len(mt))}


@dataclass
class OOFResult:
    metrics: pd.DataFrame
    spec: dict = field(default_factory=dict)
    models: list = field(default_factory=list)   # [(fold, fitted_model)]


def run_oof(frame: pd.DataFrame, make_model, *, nested_grid=None, k: int = K_TOP,
            restrict_train_source: str | None = None) -> OOFResult:
    """Outer OOF on the frozen folds. Trains on non-quarantined outer-train; scores non-quarantined held-out.
    Optional nested selection of a hyperparameter grid inside outer-train."""
    dev = frame[~frame["quarantined"]].copy()
    folds = sorted(dev["fold"].unique())
    parts, specs, models = [], [], []
    for f in folds:
        tr = dev[dev["fold"] != f]
        te = dev[dev["fold"] == f]
        if restrict_train_source is not None:
            tr = tr[tr["source"] == restrict_train_source]
        best = _nested_select(tr, make_model, nested_grid, k) if nested_grid else make_model()
        model = best.fit(tr)
        mt = per_patient_metrics(te, model.raw_score(te), k)
        parts.append(mt)
        specs.append({"fold": int(f), **(model.to_dict() if hasattr(model, "to_dict") else {})})
        models.append((int(f), model))
    metrics = pd.concat(parts, ignore_index=True)
    return OOFResult(metrics, {"folds": specs}, models)


def quarantine_stratum(frame: pd.DataFrame, oof: OOFResult, k: int = K_TOP) -> dict:
    """Robustness stratum (§8.8): score each held-out fold's FULL pool (kept + quarantined recurrent antigens)
    with that fold's OOF model, counting ONLY quarantined positives as retrieval targets. No selection here."""
    by_fold = dict(oof.models)
    rows_model, rows_prime = [], []
    for f, te_full in frame.groupby("fold"):
        model = by_fold.get(int(f))
        if model is None:
            continue
        for (src, pid), g in te_full.groupby(["source", "patient_id"]):
            qpos = g["quarantined"].to_numpy() & (g["eval_positive"].to_numpy() == 1)
            if qpos.sum() == 0:
                continue
            hm = hits_at_k(model.raw_score(g), qpos.astype(int), k)
            hp = hits_at_k(baseline_score(g, "prime"), qpos.astype(int), k)
            rows_model.append(hm); rows_prime.append(hp)
    if not rows_model:
        return {"n_patients_with_quarantined_pos": 0}
    return {"n_patients_with_quarantined_pos": len(rows_model),
            "mean_hits_model": round(float(np.mean(rows_model)), 3),
            "mean_hits_prime": round(float(np.mean(rows_prime)), 3),
            "note": "recurrent-antigen robustness stratum; reported only, never used for selection."}


def _nested_select(tr: pd.DataFrame, make_model, grid: list[dict], k: int):
    inner_folds = sorted(tr["fold"].unique())
    best_score, best_cfg = -1e9, grid[0]
    for cfg in grid:
        vals = []
        for vf in inner_folds:
            itr, ite = tr[tr["fold"] != vf], tr[tr["fold"] == vf]
            if ite.empty or itr.empty:
                continue
            m = make_model(**cfg).fit(itr)
            mt = per_patient_metrics(ite, m.raw_score(ite), k)
            if len(mt):
                vals.append(np.average(mt["hits"], weights=_source_balanced_weight(mt)))
        s = float(np.mean(vals)) if vals else -1e9
        if s > best_score:
            best_score, best_cfg = s, cfg
    return make_model(**best_cfg)


def evaluate_model(frame: pd.DataFrame, oof: OOFResult, baselines=("prime", "mix", "presentation")) -> dict:
    """Attach baseline hits per patient (on the same held-out pools) and run the paired bootstrap vs each."""
    dev = frame[~frame["quarantined"]].copy()
    base_mt = {}
    for b in baselines:
        base_mt[b] = per_patient_metrics(dev, baseline_score(dev, b)).set_index(["source", "patient_id"])
    m = oof.metrics.set_index(["source", "patient_id"])
    joined = m.rename(columns={"hits": "hits_model", "recall": "recall_model"})
    for b in baselines:
        joined[f"hits_{b}"] = base_mt[b]["hits"]
        joined[f"recall_{b}"] = base_mt[b]["recall"]
    joined = joined.reset_index()
    out = {"n_scored_patients": int(len(joined)),
           "overall_hits_model": round(float(np.average(joined["hits_model"],
                                                        weights=_source_balanced_weight(joined))), 4)}
    for b in baselines:
        out[f"vs_{b}"] = paired_bootstrap(joined, "model", b)
        out[f"overall_hits_{b}"] = round(float(np.average(joined[f"hits_{b}"],
                                                          weights=_source_balanced_weight(joined))), 4)
    out["per_source"] = _per_source(joined, baselines)
    return out, joined


def _per_source(joined: pd.DataFrame, baselines) -> dict:
    res = {}
    for src, g in joined.groupby("source"):
        d = {"patients": int(len(g)), "mean_hits_model": round(float(g["hits_model"].mean()), 3),
             "mean_recall_model": round(float(g["recall_model"].mean()), 3)}
        for b in baselines:
            d[f"mean_hits_{b}"] = round(float(g[f"hits_{b}"].mean()), 3)
            d[f"delta_vs_{b}"] = round(float((g["hits_model"] - g[f"hits_{b}"]).mean()), 3)
        res[src] = d
    return res
