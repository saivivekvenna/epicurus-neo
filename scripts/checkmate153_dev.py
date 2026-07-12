"""CheckMate 153 external-validation lane — genuine PRIME vs frozen Epicurus v0.1 on an UNTOUCHED cohort.

CheckMate 153 (Alban et al., Nat Med 2024, s41591-024-03240-y) is a combinatorial-tetramer neoantigen
screen in NSCLC (nivolumab). It is:
  * class-I, HLA-resolved (each candidate 9-mer carries its single restricting allele);
  * a genuine within-patient decision problem (14 patients, ~1,197 predicted-binder candidates, 162
    tetramer-POSITIVE vs ~1,035 tetramer-TESTED_NEGATIVE);
  * PRIME-untouched — the tetramer labels were published Oct-2024, long after PRIME 2.0's 2023 training
    set (Lausanne/NeoDisc), and the frozen Epicurus v0.1 residual was trained ONLY on the CD8 multimer
    cohort. Leakage vs PRIME's training set and vs the multimer training peptides is checked + reported.

We DISCARD the paper's own `scores` column (their trained model = circular) and compute every feature
ourselves: genuine GfellerLab PRIME 2.1 %rank (the incumbent), MHCflurry presentation percentile (the
`el` presentation feature — NetMHCpan-EL is not installed; MHCflurry is the repo's EL-type predictor and
is within-patient percentile-ranked, so the substitution is scale-free), and gene-level RNA expression
from the study's own RNA-seq count table (MOESM3). VAF is carried for a richer future model.

    python -m scripts.checkmate153_dev

The ONLY gate is the pre-registered one: source-balanced within-patient mean hits@20, frozen Epicurus
vs genuine PRIME, patient-paired bootstrap — ACCEPT only if the Δ CI lower bound > 0.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from epicurus_neo.mhcflurry_features import add_mhcflurry_predictions
from event_b.prime_adapter import score_prime
from event_b.prime_transfer import external_validate
from event_b.prime_training import prime_leakage_mask

RAW = Path("data/raw/checkmate153")
MOESM4 = RAW / "CM153_NatMed_MOESM4.xlsx"
MOESM3 = RAW / "CM153_NatMed_MOESM3.xlsx"
PRIME_CACHE = RAW / "_cache_prime.tsv"
NORM_OUT = Path("data/processed/checkmate153.normalized.csv")
ARTIFACTS = Path("artifacts/milestone_7_decision/checkmate153")

ALLELE_MAP = {
    "A0101": "A*01:01", "A0201": "A*02:01", "A0301": "A*03:01", "A1101": "A*11:01",
    "A2402": "A*24:02", "B0702": "B*07:02", "B0801": "B*08:01", "B2705": "B*27:05",
}


def build_cohort() -> pd.DataFrame:
    """Merge model_training + model_testing into one within-patient candidate universe with metadata."""
    xl = pd.ExcelFile(MOESM4)

    def sheet(name: str) -> pd.DataFrame:
        d = xl.parse(name)
        onehot = [c for c in d.columns if c.startswith("V1_") and c.replace("V1_", "") in ALLELE_MAP]
        allele = d[onehot].fillna(0).astype(float).idxmax(axis=1).str.replace("V1_", "")
        return pd.DataFrame({
            "split": name.replace("model_", ""),
            "patient_id": d["pt_pep"].astype(str).str.split("_", n=1).str[0],
            "mutant_peptide": d["pt_pep"].astype(str).str.split("_", n=1).str[1].str.upper(),
            "hla_allele": "HLA-" + allele.map(ALLELE_MAP),
            "y": d["tet_pos_or_neg"].astype(int),
            "screen_class": d["class"].astype(str) if "class" in d.columns else "NA",
        })

    cohort = pd.concat([sheet("model_training"), sheet("model_testing")], ignore_index=True)

    # maf join (Patient, MT.Peptide.x) -> gene, wildtype peptide, DNA VAF
    maf = xl.parse("maf")
    m = (maf[["Patient", "MT.Peptide.x", "Hugo_Symbol", "WT.Peptide", "vaf_maf.pre"]]
         .dropna(subset=["MT.Peptide.x"]).copy())
    m["Patient"] = m["Patient"].astype(str)
    m["MT.Peptide.x"] = m["MT.Peptide.x"].astype(str).str.upper()
    m = m.drop_duplicates(["Patient", "MT.Peptide.x"])
    cohort = cohort.merge(
        m, left_on=["patient_id", "mutant_peptide"], right_on=["Patient", "MT.Peptide.x"], how="left")
    cohort["wildtype_peptide"] = cohort["WT.Peptide"].astype(str).str.upper().where(cohort["WT.Peptide"].notna(), "")
    cohort["gene"] = cohort["Hugo_Symbol"]
    cohort["vaf"] = pd.to_numeric(cohort["vaf_maf.pre"], errors="coerce")

    # gene-level RNA expression from the study's own count table (MOESM3), baseline (_pre) sample.
    ct = pd.ExcelFile(MOESM3).parse("CountTableBMS153")
    ct = ct.rename(columns={ct.columns[0]: "gene"}).set_index("gene")
    def expr_for(row) -> float:
        pre = f"{row['patient_id']}_pre"
        on = f"{row['patient_id']}_on"
        g = row["gene"]
        if not isinstance(g, str) or g not in ct.index:
            return np.nan
        for col in (pre, on):
            if col in ct.columns:
                v = pd.to_numeric(ct.loc[g, col], errors="coerce")
                if isinstance(v, pd.Series):
                    v = v.mean()
                if pd.notna(v):
                    return float(v)
        return np.nan
    cohort["expr"] = cohort.apply(expr_for, axis=1)

    cohort["label"] = np.where(cohort["y"] == 1, "POSITIVE", "TESTED_NEGATIVE")
    cohort["candidate_id"] = ["cm153:" + str(i) for i in range(len(cohort))]
    cohort["source_dataset"] = "checkmate153"
    cohort["study_id"] = "checkmate153_alban_2024"
    return cohort


def featurize(cohort: pd.DataFrame) -> pd.DataFrame:
    """Attach genuine PRIME %rank (cached) + MHCflurry presentation percentile."""
    pairs = cohort[["mutant_peptide", "hla_allele"]].drop_duplicates()
    if PRIME_CACHE.exists():
        pr = pd.read_csv(PRIME_CACHE, sep="\t")
    else:
        res = score_prime(pairs, peptide_col="mutant_peptide", hla_col="hla_allele")
        pr = (res.scored.drop_duplicates(["mutant_peptide", "hla_allele"])
              [["mutant_peptide", "hla_allele", "prime_rank", "prime_score", "mixmhcpred_rank", "prime_status"]])
        pr.to_csv(PRIME_CACHE, sep="\t", index=False)
        print("PRIME status:", res.provenance["status_counts"])
    cohort = cohort.merge(pr, on=["mutant_peptide", "hla_allele"], how="left")

    mhc = add_mhcflurry_predictions(cohort, peptide_col="mutant_peptide", allele_col="hla_allele")
    cohort["el"] = pd.to_numeric(mhc["mhcflurry_presentation_percentile"], errors="coerce")
    cohort["prime"] = pd.to_numeric(cohort["prime_rank"], errors="coerce")
    return cohort


def main() -> int:
    ARTIFACTS.mkdir(parents=True, exist_ok=True)
    cohort = build_cohort()
    cohort = featurize(cohort)

    # persist normalized cohort (repo schema-ish)
    NORM_OUT.parent.mkdir(parents=True, exist_ok=True)
    keep = ["candidate_id", "source_dataset", "study_id", "patient_id", "split", "hla_allele",
            "mutant_peptide", "wildtype_peptide", "gene", "label", "y", "prime", "el", "expr", "vaf",
            "prime_score", "mixmhcpred_rank", "prime_status"]
    cohort[[c for c in keep if c in cohort.columns]].to_csv(NORM_OUT, index=False)

    # leakage audits (report-only; nothing is dropped from an external TEST cohort)
    peps = cohort["mutant_peptide"].astype(str).tolist()
    prime_leak = int(np.asarray(prime_leakage_mask(peps, near=True), dtype=bool).sum())

    # frame for the frozen external-validation entry point
    frame = pd.DataFrame({
        "patient_id": cohort["patient_id"].astype(str),
        "hla_allele": cohort["hla_allele"].astype(str),
        "mutant_peptide": cohort["mutant_peptide"].astype(str),
        "label": cohort["label"].to_numpy(),
        "prime": cohort["prime"].to_numpy(float),
        "el": cohort["el"].to_numpy(float),
        "expr": cohort["expr"].to_numpy(float),
    })
    # extra comparator arm: MixMHCpred %rank (companion presentation predictor; higher=better -> negate)
    frame["cmp_mixmhcpred"] = -pd.to_numeric(cohort["mixmhcpred_rank"], errors="coerce")

    frame["screen_class"] = cohort["screen_class"].to_numpy()

    n_scored = int(cohort["prime"].notna().sum())
    result = external_validate(frame.drop(columns=["screen_class"]), extra_arms=["cmp_mixmhcpred"])

    # Sensitivity: recognition-only problem among PREDICTED BINDERS (drop NonBinder easy negatives),
    # i.e. TETRAMER+ vs TETRAMER-. Keep only patients that still have >=1 positive.
    binders = frame[frame["screen_class"].isin(["TETRAMER+", "TETRAMER-"])].copy()
    keep_pts = binders.groupby("patient_id")["label"].apply(lambda s: (s == "POSITIVE").any())
    binders = binders[binders["patient_id"].isin(keep_pts[keep_pts].index)]
    binders_result = external_validate(binders.drop(columns=["screen_class"]), extra_arms=["cmp_mixmhcpred"])

    report = {
        "cohort": "CheckMate 153 (Alban et al., Nat Med 2024) — combinatorial-tetramer NSCLC neoantigen screen",
        "status": "executed",
        "design": "Frozen Epicurus v0.1 (multimer-trained residual) + genuine PRIME 2.1 applied ONCE to an "
                  "untouched external cohort. el=MHCflurry presentation percentile (NetMHCpan-EL not installed). "
                  "expr=study RNA-seq gene counts (MOESM3). The paper's own model score is discarded.",
        "n_candidates": int(len(cohort)),
        "n_patients": int(cohort["patient_id"].nunique()),
        "label_counts": cohort["label"].value_counts().to_dict(),
        "prime_scored": f"{n_scored}/{len(cohort)}",
        "expr_available": int(cohort["expr"].notna().sum()),
        "leakage_audit": {
            "prime_training_near_or_exact": prime_leak,
            "note": "CheckMate 153 tetramer labels post-date PRIME 2.0 training (2023). Frozen v0.1 was "
                    "trained on the CD8 multimer cohort only; CheckMate was never used to fit it.",
        },
        "evaluation_all_negatives": result,
        "evaluation_binders_only": {
            "note": "Recognition-only sensitivity: negatives restricted to tetramer-tested predicted BINDERS "
                    "(TETRAMER-), dropping predicted-non-binder easy negatives. The hard recognition test.",
            "screen_class_counts": cohort["screen_class"].value_counts().to_dict(),
            **binders_result,
        },
    }
    (ARTIFACTS / "CM153_EXTERNAL.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")

    def _print_block(title, ev):
        print(f"\n---- {title}: n={ev['n_candidates']} candidates, {ev['n_patients']} patients, "
              f"labels={ev['label_counts']} ----")
        for arm, m in ev["per_arm_metrics"].items():
            print(f"  {arm:18s} hits@20={m['mean_hits@20']:.3f}  recall@20={m['recall_top20']}  "
                  f"AUROC={m['tested_auroc']:.3f}" if m['tested_auroc'] is not None else
                  f"  {arm:18s} hits@20={m['mean_hits@20']:.3f}")
        gp = ev["residual_vs_comparators_top20_paired_bootstrap"]["genuine_prime"]
        print(f"  GATE frozen-Epicurus-v0.1 vs genuine PRIME: Δhits@20={gp['residual_minus_arm_hits@20']:+.4f}  "
              f"CI={[round(x,3) for x in gp['delta_ci']]}  verdict={gp['verdict']}")

    print("\n================ CheckMate 153 external validation (UNTOUCHED, PRIME-untouched) ================")
    print(f"PRIME scored: {report['prime_scored']}  expr available: {report['expr_available']}  "
          f"PRIME-training leakage flag: {report['leakage_audit']['prime_training_near_or_exact']}/{len(cohort)}")
    _print_block("ALL negatives (predicted-non-binders + tetramer-neg)", result)
    _print_block("BINDERS-ONLY (TETRAMER+ vs TETRAMER-; hard recognition test)", binders_result)
    print(f"\nartifact: {ARTIFACTS/'CM153_EXTERNAL.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
