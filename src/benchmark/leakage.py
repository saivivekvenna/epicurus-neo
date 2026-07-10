"""Deterministic guards for extraction and temporal leakage."""

from __future__ import annotations

import pandas as pd


def candidate_keys(
    df: pd.DataFrame,
    *,
    peptide_col: str = "mutant_peptide",
    hla_col: str = "hla_allele",
) -> set[tuple[str, str]]:
    return set(zip(df[peptide_col].astype(str), df[hla_col].astype(str), strict=True))


def assert_no_candidate_overlap(
    curated: pd.DataFrame,
    benchmark: pd.DataFrame,
    *,
    peptide_col: str = "mutant_peptide",
    hla_col: str = "hla_allele",
) -> None:
    overlap = candidate_keys(curated, peptide_col=peptide_col, hla_col=hla_col) & candidate_keys(
        benchmark, peptide_col=peptide_col, hla_col=hla_col
    )
    if overlap:
        examples = sorted(overlap)[:5]
        raise ValueError(
            f"Extraction leakage: {len(overlap)} candidate overlaps; examples={examples}"
        )


def assert_temporal_cutoff(
    df: pd.DataFrame,
    *,
    decision_col: str = "decision_date",
    max_sample_col: str = "max_sample_date",
) -> None:
    missing = {decision_col, max_sample_col}.difference(df.columns)
    if missing:
        raise ValueError(f"Temporal replay requires date fields: {sorted(missing)}")
    decision = pd.to_datetime(df[decision_col], errors="raise", utc=True)
    sample = pd.to_datetime(df[max_sample_col], errors="raise", utc=True)
    leaked = sample > decision
    if leaked.any():
        raise ValueError(
            f"Temporal leakage in {int(leaked.sum())} rows: sample data postdates decision"
        )
