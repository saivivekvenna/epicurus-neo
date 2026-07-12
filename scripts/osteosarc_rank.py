"""Recognition DIAGNOSTIC on the osteosarc.com patient (Sid Sijbrandij) — the first deployment-grade
patient for whom we hold BOTH the full North-Star input AND a *measured* T-cell recognition label.

Distinct from RTTP SR24-58221 (a different patient; DUA-gated Personalis report, NO labels). This patient
is PUBLIC (openly published on osteosarc.com; Backblaze b2://osteosarc-data, no DUA), so results are
committable.

INPUTS (all public, staged in data/raw/osteosarc/):
  * pvactools_all_epitopes.tsv     — pVACtools candidate universe for the 2025.01 (T2) tumor WGS:
                                      21 curated somatic mutations x ~700 peptide x HLA each = 14,780 rows,
                                      carrying a full ensemble already scored on these exact peptides
                                      (NetMHCpan, NetMHCpanEL, MHCflurry EL, BigMHC_EL/IM, DeepImmuno, ...).
  * pvactools_curated_aggregated.tsv — per-mutation summary (best peptide, tier, DNA VAF).
  * rsem.2025.01.genes.results     — matched RSEM gene TPM (SAME tumor/timepoint) -> the expr feature.
  * {May,Aug}_all_expanders.tsv    — Hudson-Lab IFNg peptide-expansion assay: PBMCs stimulated with tumor
                                      mutation-derived peptides, sorted IFNg+/-, TCR-seq'd (MiXCR); the
                                      `mutation`/`is_mutation_specific` columns are the MEASURED label.

MEASURED recognized mutations (mutation-specific IFNg+ clonal expansion):
    ASPM p.G2179R (May+Aug), DYNC1H1 p.V314I (May+Aug), MAP2 p.*868fs (May).
Of these, ONLY DYNC1H1 is in the pVACtools candidate universe. ASPM & MAP2 never entered the shortlist.

WHAT THIS IS: a method-shortcoming DIAGNOSTIC, NOT a tuned model and NOT an ACCEPT-gate win. n=3 measured
positives (1 in-universe) is far too small to fit anything. Frozen Epicurus v0.1 is applied OUT-OF-SAMPLE,
unchanged (configs/frozen/epicurus_v0_1.json). We report, side by side, where every method places the one
recognized in-universe mutation among the 21 curated candidates, plus the upstream recall gap. Any
improvement idea must be validated on the independent cohorts (multimer/Gartner/IMPROVE/CheckMate), never
on these 3 points.

    .venv/bin/python -m scripts.osteosarc_rank
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from event_b.prime_adapter import score_prime            # noqa: E402
from event_b.prime_transfer import score_with_frozen     # noqa: E402

RAW = ROOT / "data/raw/osteosarc"
ART = ROOT / "artifacts/milestone_7_decision/osteosarc_sid"   # COMMITTABLE — public data (no DUA)
ALLEP = RAW / "pvactools_all_epitopes.tsv"
AGG = RAW / "pvactools_curated_aggregated.tsv"
RSEM = RAW / "rsem.2025.01.genes.results"
EXPANDERS = [RAW / "May_all_expanders.tsv", RAW / "Aug_all_expanders.tsv"]
PRIME_CACHE = RAW / "_cache_prime.tsv"
PATIENT = "osteosarc_sid"

# Stable Ensembl gene IDs for the measured-recognized genes, so we can report their expression from the
# matched RSEM quant even when they are ABSENT from the pVACtools candidate table (that absence is the point).
RECOGNIZED_ENSG = {"ASPM": "ENSG00000066279", "DYNC1H1": "ENSG00000197102", "MAP2": "ENSG00000078018"}

# Comparators already present in the pVACtools table. (column, higher_is_better)
# Presentation / binding models and the two immunogenicity models that directly compete with Epicurus.
COMPARATORS = {
    "pres_bestmt_pctile": ("Best MT Percentile", False),          # pVACtools headline presentation %ile
    "netmhcpan_el_pctile": ("NetMHCpanEL MT Percentile", False),  # NetMHCpan-EL presentation
    "mhcflurry_el_pctile": ("MHCflurryEL Presentation MT Percentile", False),
    "bigmhc_im": ("BigMHC_IM MT Score", True),                    # immunogenicity model
    "deepimmuno": ("DeepImmuno MT Score", True),                  # immunogenicity model
}


def _num(s: pd.Series) -> pd.Series:
    return pd.to_numeric(s, errors="coerce")


def load_labels() -> dict:
    """Measured recognized mutations -> {gene: {'protein': .., 'timepoints': [..]}} from the IFNg assay."""
    recog: dict[str, dict] = {}
    for fp in EXPANDERS:
        tp = fp.name.split("_")[0]
        df = pd.read_csv(fp, sep="\t")
        ms = df[df["is_mutation_specific"].astype(str).str.upper() == "TRUE"]
        for mut in sorted(set(ms["mutation"]) - {"NA", ""}):
            gene = str(mut).split(".")[0]
            rec = recog.setdefault(gene, {"mutation_label": mut, "timepoints": []})
            if tp not in rec["timepoints"]:
                rec["timepoints"].append(tp)
    return recog


def load_expression() -> dict:
    """Ensembl gene id (version-stripped) -> TPM from the matched T2 RSEM quant."""
    r = pd.read_csv(RSEM, sep="\t", usecols=["gene_id", "TPM"])
    r["ensg"] = r["gene_id"].astype(str).str.split(".").str[0]
    return dict(zip(r["ensg"], _num(r["TPM"])))


def load_candidates(expr_map: dict) -> pd.DataFrame:
    df = pd.read_csv(ALLEP, sep="\t", dtype=str)
    df = df.rename(columns={
        "Gene Name": "gene", "Ensembl Gene ID": "ensg", "HGVSp": "hgvsp",
        "HLA Allele": "hla_allele", "MT Epitope Seq": "mutant_peptide", "Peptide Length": "pep_len",
        "Best MT IC50 Score": "best_ic50",
    })
    df["ensg_base"] = df["ensg"].astype(str).str.split(".").str[0]
    df["expr"] = df["ensg_base"].map(expr_map)           # gene-level TPM (matched timepoint)
    for _, (col, _hb) in COMPARATORS.items():
        df[col] = _num(df[col])
    df["best_ic50"] = _num(df["best_ic50"])
    df["patient_id"] = PATIENT
    return df


def attach_prime(df: pd.DataFrame) -> pd.DataFrame:
    pairs = df[["mutant_peptide", "hla_allele"]].drop_duplicates()
    if PRIME_CACHE.exists():
        pr = pd.read_csv(PRIME_CACHE, sep="\t")
    else:
        print(f"[prime] scoring {len(pairs)} unique peptide x HLA with genuine PRIME 2.1 ...", flush=True)
        pr = score_prime(pairs, peptide_col="mutant_peptide", hla_col="hla_allele").scored
        pr.to_csv(PRIME_CACHE, sep="\t", index=False)
        print(f"[prime] cached -> {PRIME_CACHE}", flush=True)
    pr = pr[["mutant_peptide", "hla_allele", "prime_rank"]].drop_duplicates(["mutant_peptide", "hla_allele"])
    df = df.merge(pr, on=["mutant_peptide", "hla_allele"], how="left")
    df["prime"] = _num(df["prime_rank"])
    # Epicurus `el` = MHCflurry presentation %ile (same model family as v0.1 training)
    df["el"] = df[COMPARATORS["mhcflurry_el_pctile"][0]]
    return df


def _oriented(v: pd.Series, higher_better: bool) -> pd.Series:
    """Orient so larger = better candidate (NaN preserved)."""
    return v if higher_better else -v


def per_mutation(df: pd.DataFrame, methods: dict) -> pd.DataFrame:
    """One row per mutation; each method's score = its BEST peptide's oriented score (larger = better)."""
    rows = []
    for gene, g in df.groupby("gene"):
        row = {"gene": gene, "mutation": g["hgvsp"].iloc[0], "n_candidates": len(g),
               "expr_tpm": float(g["expr"].iloc[0]) if pd.notna(g["expr"].iloc[0]) else None}
        for name, (col, hb) in methods.items():
            s = _oriented(g[col], hb)
            row[name] = float(s.max()) if s.notna().any() else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def method_scorecard(pm: pd.DataFrame, methods: dict, positive_genes: set) -> dict:
    """For each method: rank of the recognized in-universe mutation among all curated mutations,
    1-positive AUROC, and hits@{1,3,5}. Larger method score = better (already oriented)."""
    n = len(pm)
    pos = pm[pm["gene"].isin(positive_genes)]["gene"].tolist()
    out = {}
    for name in methods:
        order = pm.sort_values(name, ascending=False, kind="mergesort").reset_index(drop=True)
        ranks = {gene: int(order.index[order["gene"] == gene][0]) + 1 for gene in pos}
        best_rank = min(ranks.values()) if ranks else None
        # single-positive AUROC (1 pos vs n-1 neg): (n - rank)/(n - 1)
        auroc = round((n - best_rank) / (n - 1), 3) if best_rank else None
        out[name] = {
            "rank_of_recognized": ranks,            # gene -> 1-based rank (1 = top pick)
            "best_recognized_rank": best_rank,
            "n_candidates": n,
            "single_pos_auroc": auroc,
            "hits@1": int(best_rank <= 1) if best_rank else 0,
            "hits@3": int(best_rank <= 3) if best_rank else 0,
            "hits@5": int(best_rank <= 5) if best_rank else 0,
        }
    return out


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    recog = load_labels()
    expr_map = load_expression()
    df = load_candidates(expr_map)
    df = attach_prime(df)

    # Epicurus v0.1 residual (frozen, out-of-sample). Within-patient percentiles of prime/el/expr.
    df["epicurus"] = score_with_frozen(df[["patient_id", "prime", "el", "expr"]])

    methods = {
        "epicurus_v0_1": ("epicurus", True),
        "genuine_prime": ("prime", False),
        **COMPARATORS,
    }
    pm = per_mutation(df, methods)

    universe_genes = set(pm["gene"])
    recog_genes = set(recog)
    in_universe = sorted(recog_genes & universe_genes)
    missing = sorted(recog_genes - universe_genes)

    scorecard = method_scorecard(pm, methods, set(in_universe))

    # ---- coverage funnel: measured recognized -> candidate universe (recall), with expression ----
    coverage = {"measured_recognized": {}, "recall_into_candidate_universe": f"{len(in_universe)}/{len(recog_genes)}"}
    for gene, rec in sorted(recog.items()):
        ensg = RECOGNIZED_ENSG.get(gene)
        tpm = expr_map.get(ensg) if ensg else None
        tpm = float(tpm) if tpm is not None and pd.notna(tpm) else None
        coverage["measured_recognized"][gene] = {
            "mutation": rec["mutation_label"], "timepoints": rec["timepoints"],
            "in_pvactools_candidate_universe": gene in universe_genes,
            "gene_tpm_matched_timepoint": tpm,
        }

    # ---- peptide-level: where does the recognized peptide land among all 14,780 candidates? ----
    peptide_level = {}
    for gene in in_universe:
        # where the recognized gene's best peptide lands among ALL peptide x HLA candidates
        for name, (col, hb) in methods.items():
            s = _oriented(df[col], hb)
            order = df.assign(_s=s).sort_values("_s", ascending=False, kind="mergesort").reset_index(drop=True)
            hit_rows = order.index[order["gene"] == gene].tolist()
            peptide_level.setdefault(gene, {})[name] = {
                "best_candidate_rank_of_gene": int(hit_rows[0]) + 1 if hit_rows else None,
                "total_candidates": int(len(order)),
            }

    report = {
        "patient": PATIENT,
        "source": "osteosarc.com / Research to the People (Sid Sijbrandij) — PUBLIC, no DUA",
        "role": "RECOGNITION DIAGNOSTIC — measured IFNg label present; n=3 positives (1 in-universe) => "
                "descriptive method comparison, NOT a tuned model and NOT an ACCEPT-gate result",
        "hla": sorted(df["hla_allele"].unique().tolist()),
        "candidate_universe": {
            "curated_mutations": int(pm.shape[0]),
            "peptide_x_hla_candidates": int(len(df)),
            "prime_scored": int(df["prime"].notna().sum()),
            "expr_mapped": int(df["expr"].notna().sum()),
        },
        "measured_label": {
            "assay": "Hudson-Lab IFNg peptide-expansion (PBMC stim -> IFNg+/- sort -> MiXCR TCR-seq)",
            "recognized_genes": sorted(recog_genes),
            "in_universe": in_universe, "missing_from_candidate_universe": missing,
        },
        "coverage_funnel": coverage,
        "per_method_scorecard": scorecard,
        "peptide_level_rank_of_recognized_gene": peptide_level,
        "caveats": [
            "n=3 measured positives (only DYNC1H1 is in the candidate universe). Single-positive AUROC and "
            "hits@k are DESCRIPTIVE for one patient; not powered, not an ACCEPT/REJECT verdict.",
            "Negatives assumed = the other curated mutations (tested pool ~= curated set). The exact "
            "stimulation-pool composition (true denominator) should be requested from Hudson Lab / RTTP.",
            "ASPM/MAP2 absent from the 2025.01 pVACtools output; whether upstream expression/binding "
            "filtering or a callset/timepoint mismatch dropped them is flagged for follow-up.",
            "Frozen Epicurus v0.1 is applied unchanged (configs/frozen/epicurus_v0_1.json). Nothing is fit "
            "to this patient. Any improvement must be validated on the independent cohorts.",
        ],
    }
    (ART / "report.json").write_text(json.dumps(report, indent=2, default=str) + "\n")

    # ---- ordered per-mutation table under Epicurus (with the label + every method's score) ----
    pm_out = pm.copy()
    pm_out["recognized"] = pm_out["gene"].isin(in_universe)
    pm_out = pm_out.sort_values("epicurus_v0_1", ascending=False, kind="mergesort")
    pm_out.to_csv(ART / "per_mutation_scores.csv", index=False)

    # ---- console summary ----
    print("\n================ osteosarc.com (Sid) — recognition diagnostic ================")
    print(f"HLA: {report['hla']}")
    print(f"candidate universe: {pm.shape[0]} curated mutations, {len(df)} peptide x HLA; "
          f"PRIME scored {report['candidate_universe']['prime_scored']}/{len(df)}")
    print(f"\nMEASURED recognized: {sorted(recog_genes)}")
    print(f"  in candidate universe: {in_universe}   MISSING (recall gap): {missing}")
    print("\ncoverage funnel (recognized -> in shortlist? + expression):")
    for gene, c in coverage["measured_recognized"].items():
        print(f"  {gene:9s} {c['mutation']:16s} tps={c['timepoints']} in_universe={c['in_pvactools_candidate_universe']} "
              f"TPM={c['gene_tpm_matched_timepoint']}")
    print(f"\nRECALL into candidate shortlist: {coverage['recall_into_candidate_universe']}")
    print("\n--- per-method placement of the recognized in-universe mutation (DYNC1H1) among 21 curated ---")
    print(f"{'method':22s} {'rank':>5s} {'1pos_AUROC':>11s} {'hit@1':>6s} {'hit@3':>6s} {'hit@5':>6s}")
    for name, sc in scorecard.items():
        print(f"{name:22s} {str(sc['best_recognized_rank']):>5s} {str(sc['single_pos_auroc']):>11s} "
              f"{sc['hits@1']:>6d} {sc['hits@3']:>6d} {sc['hits@5']:>6d}")
    print("\n--- 21 curated mutations ranked by frozen Epicurus v0.1 (recognized flagged) ---")
    show = pm_out[["gene", "mutation", "epicurus_v0_1", "genuine_prime", "bigmhc_im", "deepimmuno",
                   "pres_bestmt_pctile", "expr_tpm", "recognized"]]
    with pd.option_context("display.width", 200, "display.max_columns", 20,
                           "display.float_format", lambda x: f"{x:.3f}"):
        print(show.to_string(index=False))
    print(f"\nartifacts -> {ART}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
