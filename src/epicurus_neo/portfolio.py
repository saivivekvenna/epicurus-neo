from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class PortfolioConstraints:
    k: int = 20
    max_per_hla: int | None = None
    max_per_gene: int | None = None
    min_score: float | None = None


def _under_limit(selected: list[pd.Series], candidate: pd.Series, column: str, limit: int | None) -> bool:
    if limit is None or column not in candidate.index:
        return True
    value = candidate[column]
    if pd.isna(value):
        return True
    count = sum(1 for row in selected if row.get(column) == value)
    return count < limit


def select_portfolio(
    frame: pd.DataFrame,
    *,
    score_col: str = "epicurus_neo_score",
    constraints: PortfolioConstraints = PortfolioConstraints(),
    hla_col: str = "hla_allele",
    gene_col: str = "gene_symbol",
) -> pd.DataFrame:
    """Greedily select a top-k submitted set with optional diversity constraints."""
    ranked = frame.sort_values(score_col, ascending=False, kind="mergesort")
    selected: list[pd.Series] = []

    for _, candidate in ranked.iterrows():
        if constraints.min_score is not None and candidate[score_col] < constraints.min_score:
            continue
        if not _under_limit(selected, candidate, hla_col, constraints.max_per_hla):
            continue
        if not _under_limit(selected, candidate, gene_col, constraints.max_per_gene):
            continue
        selected.append(candidate)
        if len(selected) >= constraints.k:
            break

    if not selected:
        return ranked.iloc[0:0].copy()
    return pd.DataFrame(selected).reset_index(drop=True)

