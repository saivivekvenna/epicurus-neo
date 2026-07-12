"""Sid (osteosarc) identical-input end-to-end benchmark — universe + leakage guard + metrics.

Falsification-first. The prior "lossless 1/3 -> 3/3" claim is INVALID as an end-to-end benchmark: the
generator (scripts/osteosarc_peptide_recovery.py) hard-codes TARGETS = {ASPM, MAP2, DYNC1H1} — the exact
recognized positives — and generates candidates only for those. Target selection leaks the answer even
though Hudson labels are joined later. This module enforces the correction:

  * the candidate-generation input is the COMPLETE, LABEL-BLIND variant universe (all 200 public variants
    from data/raw/osteosarc/site_cache/variant_vafs_long.tsv), with consequence eligibility declared BEFORE
    any label join — never a TARGETS list;
  * a hard leakage guard (assert_generation_label_blind) fails if the generated variant set is a subset
    selected from the positives, or does not cover the declared eligible universe;
  * evaluation is mutation-level recognized hits@k on stable mutation identity (dedup peptide/HLA rows so
    they cannot inflate mutation hits). Hudson labels (ASPM/DYNC1H1/MAP2) are EVALUATION-ONLY, joined only
    after every pipeline output is frozen. Post-hoc n=1 patient / 3 positives — always labelled as such.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

VAF_TABLE = Path("data/raw/osteosarc/site_cache/variant_vafs_long.tsv")

# Class-I generation eligibility — declared BEFORE label join (a variant is eligible for class-I peptide
# generation iff it has a protein-coding consequence AND somatic tumor read support).
CLASS_I_ELIGIBLE_CONSEQUENCES = frozenset({
    "missense_variant", "frameshift_variant", "stop_gained",
    "inframe_deletion", "inframe_insertion", "inframe_indel", "indel",
})

# EVALUATION-ONLY Hudson-recognized mutations (IFNy/TCR). Exact mutation identity matters: the public
# universe contains additional, unrecognized MAP2 and DYNC1H1 variants. These labels must never gate
# which variants are generated.
HUDSON_RECOGNIZED_VARIANT_IDS = frozenset({
    "ASPM-chr1-197102716",
    "DYNC1H1-chr14-101980529",
    "MAP2-chr2-209694772",
})


class GenerationLeakageError(AssertionError):
    """Raised when candidate generation was conditioned on the positive labels."""


def load_variant_universe() -> pd.DataFrame:
    """Collapse the longitudinal multi-pipeline VAF table to one row per variant_id (LABEL-BLIND).
    Eligibility uses only consequence + somatic tumor read support — no recognition label."""
    d = pd.read_csv(VAF_TABLE, sep="\t")
    tumor = d["tissue"].astype(str).str.lower().eq("tumor")
    g = d.assign(
        _tumor_alt=d["alt_reads"].where(tumor, 0),
        _tumor_vaf=d["vaf"].where(tumor),
    ).groupby("variant_id").agg(
        gene=("gene", "first"), chrom=("chrom", "first"), pos=("pos", "first"),
        ref=("ref", "first"), alt=("alt", "first"),
        consequence=("consequence", lambda s: s.dropna().iloc[0] if s.notna().any() else ""),
        max_tumor_alt_reads=("_tumor_alt", "max"), max_tumor_vaf=("_tumor_vaf", "max"),
        on_variants_page=("on_variants_page", "any"),
    ).reset_index()
    g["class_i_eligible"] = (g["consequence"].isin(CLASS_I_ELIGIBLE_CONSEQUENCES)
                             & (g["max_tumor_alt_reads"] > 0))
    return g


def eligible_universe_ids(universe: pd.DataFrame | None = None) -> set[str]:
    u = universe if universe is not None else load_variant_universe()
    return set(u.loc[u["class_i_eligible"], "variant_id"])


def hudson_positive_variant_ids(universe: pd.DataFrame | None = None) -> set[str]:
    """EVAL-ONLY: exact recognized variant IDs present in the universe. Not a generation input."""
    u = universe if universe is not None else load_variant_universe()
    return set(HUDSON_RECOGNIZED_VARIANT_IDS & set(u["variant_id"]))


def assert_generation_label_blind(generated_variant_ids: set[str], *, universe: pd.DataFrame | None = None,
                                  min_coverage: float = 0.95) -> dict:
    """HARD GUARD. Fails if the generated variant set looks label-conditioned:
      (1) it is a subset selected from the positives (covers positives but misses most of the universe), or
      (2) it does not cover the declared eligible universe (>= min_coverage of eligible variants).
    Returns coverage stats on success."""
    u = universe if universe is not None else load_variant_universe()
    elig = eligible_universe_ids(u)
    pos = hudson_positive_variant_ids(u)
    gen = set(generated_variant_ids)
    covered = gen & elig
    coverage = len(covered) / len(elig) if elig else 0.0
    pos_only = pos & gen
    # (1) positive-selected subset: generated hits the positives but covers almost none of the universe
    if pos_only and coverage < 0.5:
        raise GenerationLeakageError(
            f"generation is label-conditioned: covers {len(pos_only)} positive variant(s) but only "
            f"{coverage:.1%} of the {len(elig)} eligible universe variants (target-conditioned leakage).")
    # (2) incomplete universe coverage
    if coverage < min_coverage:
        raise GenerationLeakageError(
            f"generation covers only {coverage:.1%} of the eligible universe (< {min_coverage:.0%}); "
            f"end-to-end generation must be over the COMPLETE label-blind universe, not a subset.")
    return {"eligible_universe": len(elig), "generated_eligible_covered": len(covered),
            "coverage": round(coverage, 4), "positives_in_universe": len(pos),
            "positives_covered": len(pos & gen)}


# ------------------------------------------------------------------- mutation-level scoring ---
def mutation_hits_at_k(ranked: pd.DataFrame, positive_variant_ids: set[str], k: int = 20, *,
                       variant_col: str = "variant_id", score_col: str = "score",
                       ascending: bool = False) -> dict:
    """Mutation-level recognized hits@k. Collapse peptide/HLA rows to their best-scoring row per
    variant_id BEFORE ranking, so duplicate peptides cannot inflate mutation hits."""
    df = ranked.dropna(subset=[score_col]).sort_values(score_col, ascending=ascending, kind="mergesort")
    per_variant = df.drop_duplicates(variant_col, keep="first")
    topk = per_variant.head(k)
    hit_ids = set(topk[variant_col]) & positive_variant_ids
    return {"k": k, "n_variants_ranked": int(per_variant[variant_col].nunique()),
            "hits_at_k": len(hit_ids), "hit_variant_ids": sorted(hit_ids),
            "recall_at_k": round(len(hit_ids) / len(positive_variant_ids), 4) if positive_variant_ids else None,
            "positive_variant_ids_present": sorted(set(per_variant[variant_col]) & positive_variant_ids)}
