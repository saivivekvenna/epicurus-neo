"""Minimal-epitope expansion for scoring genuine PRIME on mutation-level (25mer) candidates.

Gartner NCI Nmers are 25mers with the SNV at the centre (position 13, 1-based; verified: 1996/2000
single-substitution 25mers differ at 0-based index 12). PRIME scores 8-14mers, so a mutation-level
candidate is expanded into all sub-peptides that CONTAIN the mutated residue, each is scored by
genuine PRIME for the patient's restricting allele(s), and the results are aggregated back to the
mutation level (best %rank) — mirroring how Gartner's "top ranked minimal" columns were built.

Mutation-level labels are preserved: one label per 25mer; the expansion only affects the score, not
the label. Multiple expansion rules (window length range, aggregation) are supported so the
head-to-head can report sensitivity rather than depend on a single arbitrary rule.

BLOCKER for a FAIR Gartner head-to-head: Gartner Nmers omits the per-patient restricting HLA
allele (its predictor columns are "best minimal across the patient's UNLISTED alleles"), and the
Gartner Testing patients are disjoint from the Müller min file that does carry HLA. So a fair PRIME
score needs an external NCI patient-HLA table. This module runs the moment that mapping is supplied;
using a broad allele panel instead is ILLUSTRATIVE ONLY (it gives PRIME an unfair best-of-many-
alleles advantage) and must never be reported as a fair PRIME-vs-us benchmark.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from event_b.prime_adapter import score_prime

DEFAULT_MUT_POS0 = 12  # 0-based centre of a 25mer (Gartner)
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class ExpansionRule:
    name: str
    min_len: int
    max_len: int
    require_mutation: bool = True   # only keep windows that include the mutated residue
    aggregate: str = "best"          # 'best' = min PRIME %rank across windows x alleles


def mutation_windows(peptide: str, mut_pos0: int, min_len: int, max_len: int,
                     require_mutation: bool = True) -> list[str]:
    """All min_len..max_len sub-peptides of `peptide`; if require_mutation, only those covering
    the mutated residue at 0-based `mut_pos0`."""
    peptide = str(peptide).strip().upper()
    n = len(peptide)
    windows = []
    for length in range(min_len, max_len + 1):
        for start in range(0, n - length + 1):
            end = start + length  # exclusive
            if require_mutation and not (start <= mut_pos0 < end):
                continue
            sub = peptide[start:end]
            if set(sub).issubset(STD_AA):
                windows.append(sub)
    return sorted(set(windows))


def score_mutations_with_prime(
    mutations: pd.DataFrame,
    *,
    peptide_col: str = "mutant_peptide",
    allele_col: str = "hla_allele",
    mut_pos_col: str | None = None,
    default_mut_pos0: int = DEFAULT_MUT_POS0,
    rule: ExpansionRule,
) -> pd.DataFrame:
    """Score genuine PRIME at the mutation level under one expansion rule.

    `mutations` needs a mutant peptide, a restricting allele (or comma-list), and optionally a
    per-row 0-based mutation position (else `default_mut_pos0`). Returns one row per input mutation
    with the aggregated PRIME %rank/score and the number of windows scored.
    """
    # Build the (window peptide, allele) pairs, remembering which mutation row each came from.
    pair_rows = []
    for idx, row in mutations.reset_index(drop=True).iterrows():
        pep = str(row[peptide_col]).strip().upper()
        mut_pos0 = int(row[mut_pos_col]) if mut_pos_col and pd.notna(row.get(mut_pos_col)) else default_mut_pos0
        windows = mutation_windows(pep, mut_pos0, rule.min_len, rule.max_len, rule.require_mutation)
        alleles = [a.strip() for a in str(row[allele_col]).replace(";", ",").split(",") if a.strip()]
        for w in windows:
            for a in alleles:
                pair_rows.append({"_mut_idx": idx, "peptide": w, "hla_allele": a})
    pairs = pd.DataFrame(pair_rows)
    if pairs.empty:
        out = mutations.reset_index(drop=True).copy()
        out["prime_rank"] = np.nan
        out["prime_windows_scored"] = 0
        out["prime_status"] = "NO_WINDOWS"
        return out

    scored = score_prime(pairs[["peptide", "hla_allele"]].drop_duplicates(),
                         peptide_col="peptide", hla_col="hla_allele").scored
    pairs = pairs.merge(scored[["peptide", "hla_allele", "prime_rank", "prime_score"]],
                        on=["peptide", "hla_allele"], how="left")

    # Aggregate to mutation level: best = min PRIME %rank across all windows x alleles.
    agg = pairs.groupby("_mut_idx").agg(
        prime_rank=("prime_rank", "min"),
        prime_score=("prime_score", "max"),
        prime_windows_scored=("prime_rank", lambda s: int(s.notna().sum())),
    )
    out = mutations.reset_index(drop=True).copy()
    out = out.join(agg, how="left")
    out["prime_windows_scored"] = out["prime_windows_scored"].fillna(0).astype(int)
    out["prime_status"] = np.where(out["prime_rank"].notna(), "SCORED", "PRIME_UNSCORABLE")
    out["expansion_rule"] = rule.name
    return out


DEFAULT_RULES = (
    ExpansionRule("classI_8_11_mut", 8, 11, True, "best"),
    ExpansionRule("classI_8_14_mut", 8, 14, True, "best"),
    ExpansionRule("classI_9_10_mut", 9, 10, True, "best"),
)
