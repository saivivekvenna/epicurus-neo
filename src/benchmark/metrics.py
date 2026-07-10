"""Per-group ranking metrics with a mandatory identity-hash tie-break."""

from __future__ import annotations

from collections.abc import Callable
from hashlib import md5

import numpy as np
import pandas as pd


_PEPTIDE_COLUMNS = ("mutant_peptide", "Mut_peptide", "mutant_peptide_norm")
_HLA_COLUMNS = ("hla_allele", "HLA_allele", "hla_allele_norm")


def _resolve_column(df: pd.DataFrame, requested: str, aliases: tuple[str, ...]) -> str:
    if requested in df.columns:
        return requested
    match = next((column for column in aliases if column in df.columns), None)
    if match is None:
        raise ValueError(
            f"Cannot construct the mandatory tie-break: missing {requested!r}. "
            f"Accepted columns: {aliases}"
        )
    return match


def identity_tiebreak(
    df: pd.DataFrame,
    *,
    peptide_col: str = "mutant_peptide",
    hla_col: str = "hla_allele",
) -> pd.Series:
    """Return the exact plan §9.1 tie-break for every candidate."""
    peptide_col = _resolve_column(df, peptide_col, _PEPTIDE_COLUMNS)
    hla_col = _resolve_column(df, hla_col, _HLA_COLUMNS)
    peptides = df[peptide_col].fillna("").astype(str)
    alleles = df[hla_col].fillna("").astype(str)
    values = (
        int(md5(f"{peptide}|{hla}".encode()).hexdigest()[:8], 16)
        for peptide, hla in zip(peptides, alleles, strict=True)
    )
    return pd.Series(values, index=df.index, dtype="uint64")


def _positive_mask(labels: pd.Series) -> np.ndarray:
    """Map supported label spellings to positive/non-positive without binarizing storage."""
    raw = labels.map(lambda value: getattr(value, "name", value))
    text = raw.astype(str).str.strip().str.upper()
    numeric = pd.to_numeric(raw, errors="coerce")
    known = text.isin(
        {
            "POSITIVE",
            "TESTED_NEGATIVE",
            "NEGATIVE",
            "UNTESTED",
            "UNKNOWN",
            "1",
            "0",
            "-1",
        }
    ) | numeric.isin([1, 0, -1])
    if not bool(known.all()):
        bad = sorted(text.loc[~known].unique().tolist())
        raise ValueError(f"Unsupported labels: {bad}")
    return ((text == "POSITIVE") | (numeric == 1)).to_numpy(dtype=bool)


def _ranked_groups(
    df: pd.DataFrame,
    *,
    group_col: str,
    score_col: str,
    label_col: str,
    ascending: bool,
) -> list[tuple[object, np.ndarray]]:
    required = {group_col, score_col, label_col}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Missing metric columns: {sorted(missing)}")
    if df[group_col].isna().any():
        raise ValueError(f"{group_col} must be non-null")
    scores = pd.to_numeric(df[score_col], errors="coerce")
    invalid = df[score_col].notna() & scores.isna()
    if invalid.any() or np.isinf(scores.dropna()).any():
        raise ValueError(f"{score_col} must contain numeric scores or missing values")

    work = df.copy()
    work[score_col] = scores
    work["_tiebreak"] = identity_tiebreak(work)
    groups: list[tuple[object, np.ndarray]] = []
    for group_id, group in work.groupby(group_col, sort=True, dropna=False):
        ranked = group.sort_values(
            [score_col, "_tiebreak"],
            ascending=[ascending, True],
            kind="mergesort",
            na_position="last",
        )
        groups.append((group_id, _positive_mask(ranked[label_col])))
    return groups


def _metric(
    df: pd.DataFrame,
    function: Callable[[np.ndarray], float],
    *,
    group_col: str,
    score_col: str,
    label_col: str,
    ascending: bool,
) -> np.ndarray:
    groups = _ranked_groups(
        df,
        group_col=group_col,
        score_col=score_col,
        label_col=label_col,
        ascending=ascending,
    )
    return np.asarray([function(labels) for _, labels in groups], dtype=float)


def hits_at_k(
    df: pd.DataFrame,
    group_col: str = "patient_id",
    score_col: str = "score",
    k: int = 20,
    ascending: bool = False,
    label_col: str = "label",
) -> np.ndarray:
    """Return true-positive counts in the top k, one value per sorted group."""
    if k <= 0:
        raise ValueError("k must be positive")
    return _metric(
        df,
        lambda labels: float(labels[:k].sum()),
        group_col=group_col,
        score_col=score_col,
        label_col=label_col,
        ascending=ascending,
    )


def capture_fraction(
    df: pd.DataFrame,
    group_col: str = "patient_id",
    score_col: str = "score",
    k: int = 20,
    ascending: bool = False,
    label_col: str = "label",
) -> np.ndarray:
    """Return hits / min(number of positives, k); NaN for zero-positive groups."""
    if k <= 0:
        raise ValueError("k must be positive")

    def calculate(labels: np.ndarray) -> float:
        positives = int(labels.sum())
        return float(labels[:k].sum() / min(positives, k)) if positives else float("nan")

    return _metric(
        df,
        calculate,
        group_col=group_col,
        score_col=score_col,
        label_col=label_col,
        ascending=ascending,
    )


def p_at_least_one(
    df: pd.DataFrame,
    group_col: str = "patient_id",
    score_col: str = "score",
    k: int = 20,
    ascending: bool = False,
    label_col: str = "label",
) -> np.ndarray:
    """Return whether top k contains a hit; NaN for zero-positive groups."""
    if k <= 0:
        raise ValueError("k must be positive")

    def calculate(labels: np.ndarray) -> float:
        if not labels.any():
            return float("nan")
        return float(labels[:k].any())

    return _metric(
        df,
        calculate,
        group_col=group_col,
        score_col=score_col,
        label_col=label_col,
        ascending=ascending,
    )


def precision_at_k(
    df: pd.DataFrame,
    group_col: str = "patient_id",
    score_col: str = "score",
    k: int = 20,
    ascending: bool = False,
    label_col: str = "label",
) -> np.ndarray:
    """Return hits / k, as pre-registered (including lists shorter than k)."""
    return (
        hits_at_k(
            df,
            group_col=group_col,
            score_col=score_col,
            k=k,
            ascending=ascending,
            label_col=label_col,
        )
        / k
    )


def mrr(
    df: pd.DataFrame,
    group_col: str = "patient_id",
    score_col: str = "score",
    k: int = 20,
    ascending: bool = False,
    label_col: str = "label",
) -> np.ndarray:
    """Return reciprocal rank of the first positive (not truncated at k)."""
    del k  # Kept in the common metric signature; MRR is not top-k truncated.

    def calculate(labels: np.ndarray) -> float:
        positive_positions = np.flatnonzero(labels)
        return 1.0 / float(positive_positions[0] + 1) if len(positive_positions) else 0.0

    return _metric(
        df,
        calculate,
        group_col=group_col,
        score_col=score_col,
        label_col=label_col,
        ascending=ascending,
    )


def group_positive_counts(
    df: pd.DataFrame,
    *,
    group_col: str = "patient_id",
    label_col: str = "label",
) -> np.ndarray:
    """Return total positives per group in the same sorted order as every metric."""
    if group_col not in df or label_col not in df:
        raise ValueError(f"Missing {group_col!r} or {label_col!r}")
    return np.asarray(
        [
            float(_positive_mask(group[label_col]).sum())
            for _, group in df.groupby(group_col, sort=True, dropna=False)
        ]
    )
