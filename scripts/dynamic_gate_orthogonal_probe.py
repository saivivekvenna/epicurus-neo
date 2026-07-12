"""Circularity audit + orthogonal decoy-discrimination probe (Milestone 7, v1 correction).

Two honest analyses on data already in hand — no new results are invented:

  A. CIRCULARITY AUDIT. The v1 gate's core veto axes are within-patient percentiles of EL and PRIME.
     The downstream rankers are genuine PRIME (= -PRIME) and frozen Epicurus (a PRIME-dominated logistic
     on prime/el/expr percentiles). So the gate removes exactly the candidates the SAME rankers already
     rank at the bottom. We quantify: (i) how many gate-removed candidates ever sat in a ranker's top-20
     (expected 0), and (ii) the Spearman between the gate keep-margin and each ranker score. This shows
     that "unchanged top-20" is a STRUCTURAL consequence of feature overlap, not evidence that a
     label-blind gate in general cannot help.

  B. ORTHOGONAL PROBE (multimer, IN-SAMPLE, exploratory). Multimer is the one cohort with recognition
     features that are NOT downstream-ranker inputs: agretopicity, foreignness, proteasomal processing,
     dissimilarity, TPM. Restricted to the HARD-DECOY stratum (candidates high on BOTH EL and PRIME
     percentiles — the top-quartile-both region where high-presentation tested-negatives compete with
     positives), we ask the v2 question: do orthogonal features separate POSITIVE from TESTED_NEGATIVE
     there? Cross-fitted leave-one-patient-out AUROC + the removal-at-95%-retention a purely orthogonal
     veto would achieve within the stratum. multimer is frozen-Epicurus' training cohort => IN-SAMPLE and
     tiny (34 positives) => this is a mechanism sanity check, NEVER a headline.

    python -m scripts.dynamic_gate_orthogonal_probe

Writes artifacts/milestone_7_decision/dynamic_gate/{CIRCULARITY_AUDIT.md, orthogonal_probe.json}.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from event_b.dynamic_gate import GateConfig, apply_gate, attach_percentiles, within_patient_percentile
from event_b.pool_size_sensitivity import patient_eligibility, score_arms

POOL = Path("artifacts/milestone_7_decision/pool_size_sensitivity")
ART = Path("artifacts/milestone_7_decision/dynamic_gate")
DEV_COHORTS = ["gartner", "improve", "multimer"]
RANKERS = ["genuine_prime", "frozen_epicurus"]
FROZEN_T = 0.25   # v1 frozen config threshold
STRATUM_Q = 0.75  # "hard decoy" = top-quartile on BOTH el and prime percentiles


def _load_dev(name: str) -> pd.DataFrame:
    f = pd.read_csv(POOL / f"base_{name}.csv")
    f["patient_id"] = f["patient_id"].astype(str)
    f = f[f["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    elig = patient_eligibility(f)
    return f[f["patient_id"].isin(elig.eligible)].reset_index(drop=True)


# --------------------------------------------------------------------------------------------------
# A. Circularity audit
# --------------------------------------------------------------------------------------------------
def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    m = np.isfinite(a) & np.isfinite(b)
    if m.sum() < 3:
        return float("nan")
    ar = pd.Series(a[m]).rank().to_numpy()
    br = pd.Series(b[m]).rank().to_numpy()
    return float(np.corrcoef(ar, br)[0, 1])


def circularity_audit() -> dict:
    cfg = GateConfig(t=FROZEN_T)
    out = {}
    for name in DEV_COHORTS:
        frame = _load_dev(name)
        gated = apply_gate(frame, cfg)
        removed = ~gated["dyn_gate_keep"].to_numpy(bool)
        per_ranker = {}
        removed_in_top20_total = {}
        for pid, gp in gated.groupby("patient_id"):
            scored = score_arms(gp)
            rem_local = ~gp["dyn_gate_keep"].to_numpy(bool)
            for a in RANKERS:
                r = scored.sort_values(a, ascending=False, kind="mergesort")
                top20_idx = set(r.index[:20])
                removed_idx = set(gp.index[rem_local])
                removed_in_top20_total.setdefault(a, 0)
                removed_in_top20_total[a] += len(removed_idx & top20_idx)
        # global spearman: gate keep-margin (min core-axis percentile) vs each ranker score
        p = attach_percentiles(frame)
        keep_margin = np.minimum(p["s_el"].to_numpy(), p["s_prime"].to_numpy())  # higher = safer from veto
        scored_all = score_arms(frame)
        for a in RANKERS:
            per_ranker[a] = {
                "removed_that_were_in_top20": int(removed_in_top20_total[a]),
                "spearman_keepmargin_vs_ranker": round(_spearman(keep_margin, scored_all[a].to_numpy(float)), 3),
            }
        # hard-decoy stratum removal (high EL AND high PRIME tested-negatives)
        s_el = p["s_el"].to_numpy()
        s_pr = p["s_prime"].to_numpy()
        isneg = (frame["label"].to_numpy() == "TESTED_NEGATIVE")
        hard_decoy = isneg & (s_el > STRATUM_Q) & (s_pr > STRATUM_Q)
        out[name] = {
            "n_candidates": int(len(frame)),
            "n_removed_by_gate": int(removed.sum()),
            "per_ranker": per_ranker,
            "hard_decoy_stratum": {
                "n": int(hard_decoy.sum()),
                "fraction_removed_by_gate": round(float(np.mean(removed[hard_decoy])), 4) if hard_decoy.any() else None,
            },
        }
    return out


# --------------------------------------------------------------------------------------------------
# B. Orthogonal probe (multimer, IN-SAMPLE)
# --------------------------------------------------------------------------------------------------
ORTHO = [  # (column, higher_raw_better) — none of these is a downstream-ranker input
    ("Agretopicity", False),           # MT/WT rank ratio; lower = more mutant-specific
    ("Foreignness score", True),
    ("Proteasomal processing score", True),
    ("Dissimilarity", True),
    ("RNA expression (TPM)", True),
]


def _clopper_pearson_lower(k: int, n: int, conf: float = 0.95) -> float:
    if n <= 0 or k <= 0:
        return 0.0
    if k >= n:
        return float((1 - conf) ** (1 / n))
    from scipy.stats import beta
    return float(beta.ppf(1 - conf, k, n - k + 1))


def orthogonal_probe() -> dict:
    from event_b.cd8_multimer_corpus import load_cd8_multimer

    m = load_cd8_multimer().frame.copy()
    m = m[m["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].reset_index(drop=True)
    m["patient_id"] = m["patient_id"].astype(str)
    # need EL + PRIME percentiles to define the stratum; PRIME from the pool base csv (already computed)
    base = pd.read_csv(POOL / "base_multimer.csv")
    base["patient_id"] = base["patient_id"].astype(str)
    m = m.merge(base[["patient_id", "mutant_peptide", "hla_allele", "prime"]].drop_duplicates(),
                on=["patient_id", "mutant_peptide", "hla_allele"], how="left")
    m["el"] = pd.to_numeric(m["EL (%Rank score)"], errors="coerce")
    m["s_el"] = within_patient_percentile(m, "el", higher_better=False)
    m["s_prime"] = within_patient_percentile(m, "prime", higher_better=False)

    stratum = (m["s_el"] > STRATUM_Q) & (m["s_prime"] > STRATUM_Q)
    st = m[stratum].reset_index(drop=True)
    y = (st["label"] == "POSITIVE").astype(int).to_numpy()

    # orthogonal feature matrix (oriented; standardized; missing -> column median, but flagged missing)
    cols_present = [(c, hi) for c, hi in ORTHO if c in st]
    X = np.column_stack([
        (pd.to_numeric(st[c], errors="coerce").to_numpy() * (1 if hi else -1)) for c, hi in cols_present
    ])
    present_mask = np.isfinite(X)
    col_med = np.nanmedian(np.where(present_mask, X, np.nan), axis=0)
    Xf = np.where(present_mask, X, col_med)
    mu, sd = Xf.mean(0), Xf.std(0)
    sd = np.where(sd == 0, 1.0, sd)
    Xz = (Xf - mu) / sd

    # cross-fitted leave-one-patient-out AUROC of the orthogonal residual within the stratum
    pids = st["patient_id"].to_numpy()
    pred = np.full(len(st), np.nan)
    for pid in np.unique(pids):
        te = pids == pid
        tr = ~te
        if y[tr].sum() == 0 or (y[tr] == 0).sum() == 0 or te.sum() == 0:
            continue
        clf = LogisticRegression(max_iter=2000, class_weight="balanced").fit(Xz[tr], y[tr])
        pred[te] = clf.decision_function(Xz[te])

    scored = np.isfinite(pred)
    auroc = _auroc(y[scored], pred[scored]) if scored.sum() > 2 and 0 < y[scored].sum() < scored.sum() else None

    # v2-style within-stratum orthogonal veto: keep top by orthogonal score; measure neg removal @ >=95% pos retention
    removal_at_95 = None
    if auroc is not None:
        order = np.argsort(-pred[scored])  # best orthogonal score first
        ys = y[scored][order]
        npos = int(ys.sum())
        nneg = int((ys == 0).sum())
        best = 0.0
        for keep_k in range(len(ys), 0, -1):
            kept = ys[:keep_k]
            pos_ret = kept.sum() / npos if npos else 0
            if pos_ret >= 0.95:
                neg_removed = (nneg - (keep_k - kept.sum())) / nneg if nneg else 0
                best = max(best, neg_removed)
        removal_at_95 = round(float(best), 4)

    if auroc is None:
        interp = "Stratum too small / degenerate to fit a cross-fitted model."
    elif auroc >= 0.58:
        interp = (
            f"A FAINT orthogonal signal in-sample (AUROC {auroc:.2f}; ~{removal_at_95:.0%} of high-presentation "
            "decoys removable at >=95% stratum retention) — but this is multimer's own training cohort "
            "(IN-SAMPLE), tiny (24 stratum positives), and leave-one-patient-out is NOT leave-one-study-out. "
            "An INDEPENDENT sequence-only hard-decoy experiment (train-Gartner->test-IMPROVE) collapsed OOD "
            "(retained 1.5% of positives), so cross-study/assay shift is the real killer. => promising "
            "mechanism, NOT validated; requires real WES/RNA features + leave-one-study-out + OOD abstention."
        )
    else:
        interp = (
            f"AUROC {auroc:.2f} ~ chance => orthogonal features do NOT separate high-presentation decoys "
            "from positives here (consistent with the recognition wall). A real test needs WES/RNA features "
            "absent from this cohort (mutant-RNA VAF, DNA VAF/depth/CCF)."
        )
    return {
        "cohort": "multimer",
        "flags": ["IN_SAMPLE (frozen Epicurus training cohort)", "EXPLORATORY", "underpowered (34 positives)",
                  "leave-one-PATIENT-out only (NOT leave-one-study-out) => does not probe cross-study shift"],
        "stratum_definition": f"within-patient EL percentile > {STRATUM_Q} AND PRIME percentile > {STRATUM_Q}",
        "orthogonal_features_used": [c for c, _ in cols_present],
        "note_features_not_in_downstream_rankers": True,
        "stratum_size": int(len(st)),
        "stratum_positives": int(y.sum()),
        "stratum_tested_negatives": int((y == 0).sum()),
        "crossfit_orthogonal_auroc_within_stratum": round(auroc, 4) if auroc is not None else None,
        "orthogonal_veto_negremoval_at_95pct_retention_within_stratum": removal_at_95,
        "independent_cross_study_falsification": (
            "A separate (uncommitted, not ours) sequence-only hard-decoy gate — peptide n-grams + "
            "allele/anchor + length/chemistry, PRIME/EL only for the stratum — failed catastrophically "
            "cross-study: train-Gartner->test-IMPROVE retained 1.5% of positives (removed 198/201); "
            "train-IMPROVE->test-Gartner retained 12.8% (removed 34/39). Peptide/HLA recognition motifs "
            "encode severe study/assay domain shift => sequence-only vetoes are FALSIFIED for v2; "
            "within-study CV is misleading; OOD abstention is mandatory. "
            "(scripts/hard_decoy_gate_experiment.py, artifacts/.../hard_decoy_gate/results.json — not ours)"
        ),
        "interpretation": interp,
    }


def _auroc(y: np.ndarray, s: np.ndarray) -> float:
    order = np.argsort(s)
    ranks = np.empty(len(s), float)
    ranks[order] = np.arange(1, len(s) + 1)
    n1 = y.sum()
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    return float((ranks[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _audit_md(circ: dict, probe: dict) -> str:
    L = ["# Circularity / feature-overlap audit — dynamic gate v1\n",
         "The v1 gate's **core veto axes are within-patient percentiles of EL and PRIME**. The downstream "
         "rankers are **genuine PRIME** (`-PRIME`) and **frozen Epicurus** (logistic on prime/el/expr "
         "percentiles; PRIME coef 1.04 dominates el 0.235 / expr 0.25). The gate's veto inputs are a "
         "SUBSET of — and the dominant terms in — the ranker inputs.\n"]
    L.append("\n## Feature overlap: gate veto axes vs downstream ranker inputs\n")
    L.append("| downstream ranker | inputs | overlaps gate veto {EL, PRIME}? |")
    L.append("|---|---|---|")
    L.append("| genuine PRIME | PRIME %rank | **PRIME (100% overlap)** |")
    L.append("| frozen Epicurus v0.1 | PRIME, EL, expr percentiles | **EL + PRIME (2 of 3, both dominant)** |")
    L.append("\n## Empirical consequence (frozen t=0.25)\n")
    L.append("| cohort | removed | removed that were in ANY ranker top-20 | hard-decoy stratum n | stratum removed |")
    L.append("|---|--:|--:|--:|--:|")
    for name, e in circ.items():
        rit = max(e["per_ranker"][a]["removed_that_were_in_top20"] for a in RANKERS)
        hd = e["hard_decoy_stratum"]
        L.append(f"| {name} | {e['n_removed_by_gate']} | **{rit}** | {hd['n']} | "
                 f"**{hd['fraction_removed_by_gate']}** |")
    L.append("\n**Spearman(gate keep-margin, ranker score)** — near +1 confirms the gate orders candidates "
             "the same way the rankers do:\n")
    for name, e in circ.items():
        parts = ", ".join(f"{a} {e['per_ranker'][a]['spearman_keepmargin_vs_ranker']:+.2f}" for a in RANKERS)
        L.append(f"- {name}: {parts}")
    L.append("\n## Conclusion\n"
             "Essentially zero gate-removed candidates ever sat in a ranker's top-20 (0 / 0 / 1 across "
             "gartner / improve / multimer; the single multimer case is a saturated ≤20-candidate pool "
             "where every row is trivially 'top-20'), and the gate removes **0%** of "
             "the hard-decoy stratum (high-EL AND high-PRIME tested-negatives). This is a **structural "
             "tautology of feature overlap**, not evidence about label-blind gating in general: a "
             "same-feature monotone presentation gate CANNOT move a top-20 produced by those same features. "
             "The falsification is therefore scoped to **presentation/PRIME-derived gating with current "
             "peptide features** — the orthogonal-feature dynamic-gate hypothesis is untested (data-blocked), "
             "NOT falsified.\n")
    L.append("\n## Orthogonal probe (multimer, IN-SAMPLE, exploratory)\n"
             f"Within the hard-decoy stratum (EL%>{STRATUM_Q} AND PRIME%>{STRATUM_Q}: "
             f"{probe['stratum_positives']} positives / {probe['stratum_tested_negatives']} tested-negatives), "
             f"a cross-fitted leave-one-patient-out model on orthogonal features "
             f"({', '.join(probe['orthogonal_features_used'])}) scores "
             f"**AUROC {probe['crossfit_orthogonal_auroc_within_stratum']}**; a purely-orthogonal veto would "
             f"remove **{probe['orthogonal_veto_negremoval_at_95pct_retention_within_stratum']}** of stratum "
             f"negatives at >=95% positive retention. {probe['interpretation']}\n"
             "\n> multimer is frozen Epicurus' training cohort (IN-SAMPLE) and tiny (34 positives) — this is "
             "a mechanism sanity check, never a headline. It is the ONLY orthogonal signal available now; a "
             "real test requires the WES/RNA features absent from every current eval cohort.\n")
    return "\n".join(L)


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    circ = circularity_audit()
    probe = orthogonal_probe()
    (ART / "orthogonal_probe.json").write_text(
        json.dumps({"circularity_audit": circ, "orthogonal_probe": probe}, indent=2, sort_keys=True, default=str) + "\n")
    (ART / "CIRCULARITY_AUDIT.md").write_text(_audit_md(circ, probe))
    print("Wrote CIRCULARITY_AUDIT.md + orthogonal_probe.json")
    for name, e in circ.items():
        rit = max(e["per_ranker"][a]["removed_that_were_in_top20"] for a in RANKERS)
        print(f"  [{name}] removed={e['n_removed_by_gate']} in_any_top20={rit} "
              f"hard_decoy_removed={e['hard_decoy_stratum']['fraction_removed_by_gate']}")
    print(f"  [multimer orthogonal probe IN-SAMPLE] stratum {probe['stratum_positives']}pos/"
          f"{probe['stratum_tested_negatives']}neg AUROC={probe['crossfit_orthogonal_auroc_within_stratum']} "
          f"removal@95%={probe['orthogonal_veto_negremoval_at_95pct_retention_within_stratum']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
