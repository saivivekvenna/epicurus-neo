"""First end-to-end Epicurus deployment on a REAL clinical patient (RTTP / Research to the People).

Patient SR24-58221 (Personalis ImmunoID NeXT-style report, DUA-gated Google Drive). This is the literal
North Star input: WES somatic mutations + RNA expression (TPM) + full class-I HLA typing (+ HLA LOH) ->
a RANKED neoantigen list out. We run the genuine tool stack:

    prime = genuine GfellerLab PRIME 2.1 %rank   (event_b.prime_adapter.score_prime)
    el    = MHCflurry 2.2.1 presentation percentile
    expr  = Gene Level Expression TPM (from the RNA pipeline, embedded in the neoantigen report)
    epicurus = frozen Epicurus v0.1 within-patient-percentile logistic residual (score_with_frozen)

IMPORTANT — this is a DEPLOYMENT DEMO, not a benchmark. The RTTP report carries NO experimental T-cell
recognition label (the "Immunogenicity Score" is a *predicted* value). So we emit a ranked list and report
method CONCORDANCE (Epicurus vs genuine PRIME vs the report's SHERPA presentation rank vs the report's
predicted immunogenicity). We do NOT compute hits@20 / AUROC / any ACCEPT-gate: there is no answer key.

HLA LOH: this tumor has lost one class-I allele (HIGH confidence). Candidates presentable ONLY by the lost
allele cannot be displayed by the tumor, so we also produce a deployable LOH-aware ranking (drop LOH-lost
presentations, re-percentile) and quantify how much the shortlist changes.

    .venv/bin/python -m scripts.rttp_rank
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from event_b.prime_adapter import score_prime                     # noqa: E402
from epicurus_neo.mhcflurry_features import add_mhcflurry_predictions  # noqa: E402
from event_b.prime_transfer import score_with_frozen              # noqa: E402

RAW = ROOT / "data/raw/rttp"
SNV = RAW / "Neoantigen/Neoantigen/tsv/DNA_SR24-58221_C1_neoantigen_class_I_report_SNV_Indel.tsv"
ART = ROOT / "artifacts/milestone_7_decision/rttp_sr24"
PRIME_CACHE = RAW / "_cache_prime.tsv"
MHC_CACHE = RAW / "_cache_mhcflurry.tsv"
PATIENT = "RTTP_SR24-58221"


def to_mhcflurry(a: str) -> str:
    a = str(a)
    return f"HLA-{a[0]}*{a[1:3]}:{a[3:5]}"


def _spearman(a: pd.Series, b: pd.Series) -> float:
    m = a.notna() & b.notna()
    if m.sum() < 3:
        return float("nan")
    return float(a[m].rank().corr(b[m].rank()))


def load_candidates() -> pd.DataFrame:
    df = pd.read_csv(SNV, sep="\t")
    df = df.rename(columns={
        "Peptide": "mutant_peptide", "HLA": "hla_prime", "Gene Level Expression TPM": "expr",
        "Variant": "variant", "Gene Symbol": "gene", "Protein Variant": "protein_variant",
        "Source Variant Type": "src_type", "HLA LOH": "hla_loh",
        "SHERPA Presentation Rank": "sherpa_rank", "NetMHCpan Binding Rank": "netmhcpan_rank",
        "Immunogenicity Score": "rttp_immuno", "Expressed": "expressed",
    })
    df["expr"] = pd.to_numeric(df["expr"], errors="coerce")
    df["sherpa_rank"] = pd.to_numeric(df["sherpa_rank"], errors="coerce")
    df["rttp_immuno"] = pd.to_numeric(df["rttp_immuno"], errors="coerce")
    df["patient_id"] = PATIENT
    return df


def attach_prime(df: pd.DataFrame) -> pd.DataFrame:
    pairs = (df[["mutant_peptide", "hla_prime"]].drop_duplicates()
             .rename(columns={"hla_prime": "hla_allele"}))
    if PRIME_CACHE.exists():
        pr = pd.read_csv(PRIME_CACHE, sep="\t")
    else:
        print(f"[prime] scoring {len(pairs)} unique peptide x HLA pairs with genuine PRIME 2.1 ...", flush=True)
        res = score_prime(pairs, peptide_col="mutant_peptide", hla_col="hla_allele")
        pr = res.scored
        pr.to_csv(PRIME_CACHE, sep="\t", index=False)
        print(f"[prime] cached -> {PRIME_CACHE}", flush=True)
    key_p = "mutant_peptide" if "mutant_peptide" in pr.columns else ("peptide" if "peptide" in pr.columns else None)
    key_h = "hla_allele" if "hla_allele" in pr.columns else ("hla" if "hla" in pr.columns else None)
    pr = pr.rename(columns={key_p: "mutant_peptide", key_h: "hla_prime"})
    df = df.merge(pr[["mutant_peptide", "hla_prime", "prime_rank"]], on=["mutant_peptide", "hla_prime"], how="left")
    df["prime"] = pd.to_numeric(df["prime_rank"], errors="coerce")
    return df


def attach_mhcflurry(df: pd.DataFrame) -> pd.DataFrame:
    if MHC_CACHE.exists():
        mh = pd.read_csv(MHC_CACHE, sep="\t")
    else:
        u = df[["mutant_peptide", "hla_prime"]].drop_duplicates().copy()
        u["hla_allele"] = u["hla_prime"].map(to_mhcflurry)
        u["wildtype_peptide"] = ""   # required by schema.add_normalized_columns; unused for MHCflurry EL
        print(f"[mhcflurry] scoring {len(u)} pairs ...", flush=True)
        scored = add_mhcflurry_predictions(u, peptide_col="mutant_peptide", allele_col="hla_allele")
        mh = scored[["mutant_peptide", "hla_prime", "mhcflurry_presentation_percentile"]]
        mh.to_csv(MHC_CACHE, sep="\t", index=False)
        print(f"[mhcflurry] cached -> {MHC_CACHE}", flush=True)
    df = df.merge(mh, on=["mutant_peptide", "hla_prime"], how="left")
    df["el"] = pd.to_numeric(df["mhcflurry_presentation_percentile"], errors="coerce")
    return df


def rank_arms(df: pd.DataFrame) -> pd.DataFrame:
    """Score the candidate table with the frozen formula (within-patient percentiles recomputed on THIS set)."""
    df = df.copy()
    df["epicurus"] = score_with_frozen(df[["patient_id", "prime", "el", "expr"]])
    df["genuine_prime"] = -df["prime"]          # higher = better
    df["sherpa"] = -df["sherpa_rank"]           # lower rank = better -> negate
    return df


def mutation_shortlist(df: pd.DataFrame, score_col: str) -> pd.DataFrame:
    """Best candidate per mutation by `score_col`; ranked descending."""
    idx = df.groupby("variant")[score_col].idxmax()
    best = df.loc[idx].copy()
    return best.sort_values(score_col, ascending=False, kind="mergesort").reset_index(drop=True)


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    df = load_candidates()
    df = attach_prime(df)
    df = attach_mhcflurry(df)

    cov = {
        "candidates_peptide_x_hla": int(len(df)),
        "unique_peptides": int(df["mutant_peptide"].nunique()),
        "unique_mutations": int(df["variant"].nunique()),
        "hla_alleles": sorted(df["hla_prime"].unique().tolist()),
        "prime_scored": int(df["prime"].notna().sum()),
        "mhcflurry_scored": int(df["el"].notna().sum()),
        "expr_present": int(df["expr"].notna().sum()),
        "loh_lost_rows_on_lost_allele": int((df["hla_loh"].astype(str).str.upper() == "Y").sum()),
    }

    # ---- naive ranking (all candidates) ----
    scored = rank_arms(df)
    epi = mutation_shortlist(scored, "epicurus")
    prime_sl = mutation_shortlist(scored, "genuine_prime")
    sherpa_sl = mutation_shortlist(scored, "sherpa")
    immuno_sl = mutation_shortlist(scored, "rttp_immuno")

    # ---- deployable LOH-aware ranking: drop presentations on the lost allele, re-percentile ----
    dep = df[df["hla_loh"].astype(str).str.upper() != "Y"].copy()
    dep = rank_arms(dep)
    epi_dep = mutation_shortlist(dep, "epicurus")

    # ---- method concordance over per-mutation best Epicurus candidates ----
    m = mutation_shortlist(scored, "epicurus")[["variant", "epicurus", "genuine_prime", "sherpa", "rttp_immuno", "expr"]]
    concord = {
        "per_mutation_best": {
            "epicurus_vs_genuine_prime": round(_spearman(m["epicurus"], m["genuine_prime"]), 3),
            "epicurus_vs_sherpa_presentation": round(_spearman(m["epicurus"], m["sherpa"]), 3),
            "epicurus_vs_rttp_immunogenicity": round(_spearman(m["epicurus"], m["rttp_immuno"]), 3),
            "genuine_prime_vs_sherpa": round(_spearman(m["genuine_prime"], m["sherpa"]), 3),
        },
        "per_candidate_all_rows": {
            "epicurus_vs_genuine_prime": round(_spearman(scored["epicurus"], scored["genuine_prime"]), 3),
            "epicurus_vs_sherpa_presentation": round(_spearman(scored["epicurus"], scored["sherpa"]), 3),
            "epicurus_vs_rttp_immunogenicity": round(_spearman(scored["epicurus"], scored["rttp_immuno"]), 3),
            "genuine_prime_vs_sherpa": round(_spearman(scored["genuine_prime"], scored["sherpa"]), 3),
        },
    }

    def topset(sl, n=20):
        return set(sl["variant"].head(n))
    overlap = {
        "epicurus_vs_prime_top20": len(topset(epi) & topset(prime_sl)),
        "epicurus_vs_sherpa_top20": len(topset(epi) & topset(sherpa_sl)),
        "epicurus_naive_vs_LOHaware_top20": len(topset(epi) & topset(epi_dep)),
        "top20_epicurus_on_lost_allele": int(epi.head(20)["hla_loh"].astype(str).str.upper().eq("Y").sum()),
    }

    # ---- write full ranked candidate table + shortlists ----
    outcols = ["variant", "gene", "protein_variant", "src_type", "mutant_peptide", "hla_prime", "hla_loh",
               "prime", "el", "expr", "expressed", "sherpa_rank", "rttp_immuno", "epicurus"]
    scored.sort_values("epicurus", ascending=False, kind="mergesort")[outcols].to_csv(
        ART / "RTTP_SR24_ranked_candidates.csv", index=False)
    sl_cols = ["variant", "gene", "protein_variant", "src_type", "mutant_peptide", "hla_prime", "hla_loh",
               "prime", "el", "expr", "sherpa_rank", "rttp_immuno", "epicurus"]
    epi.head(50)[sl_cols].to_csv(ART / "RTTP_SR24_shortlist_epicurus.csv", index=False)
    epi_dep.head(50)[sl_cols].to_csv(ART / "RTTP_SR24_shortlist_epicurus_LOHaware.csv", index=False)

    summary = {
        "patient": PATIENT, "source": "Research to the People (DUA-gated); Personalis ImmunoID NeXT-style report",
        "role": "DEPLOYMENT DEMO — ranked neoantigen list; NO experimental recognition label -> NOT a benchmark",
        "coverage": cov, "method_concordance_spearman": concord, "top20_overlap": overlap,
        "loh": {"lost_allele": "one class-I allele (HIGH confidence)", "note": "LOH-aware ranking drops lost-allele-only presentations"},
        "artifacts": ["RTTP_SR24_ranked_candidates.csv", "RTTP_SR24_shortlist_epicurus.csv",
                      "RTTP_SR24_shortlist_epicurus_LOHaware.csv"],
    }
    (ART / "RTTP_SR24.json").write_text(json.dumps(summary, indent=2, default=str) + "\n")

    # ---- console ----
    pd.set_option("display.width", 200); pd.set_option("display.max_columns", 30)
    print("\n================ RTTP SR24-58221 — Epicurus end-to-end deployment ================")
    print(json.dumps(cov, indent=2))
    print("\nmethod concordance (Spearman, per-mutation best):", json.dumps(concord))
    print("top-20 overlap:", json.dumps(overlap))
    print("\n--- TOP 15 neoantigens by frozen Epicurus v0.1 (naive) ---")
    show = ["gene", "protein_variant", "mutant_peptide", "hla_prime", "hla_loh", "prime", "el", "expr", "rttp_immuno", "epicurus"]
    with pd.option_context("display.float_format", lambda x: f"{x:.3f}"):
        print(epi.head(15)[show].to_string(index=False))
    print("\n--- TOP 15 by deployable LOH-aware Epicurus (lost-allele presentations dropped) ---")
    with pd.option_context("display.float_format", lambda x: f"{x:.3f}"):
        print(epi_dep.head(15)[show].to_string(index=False))
    print(f"\nartifacts -> {ART}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
