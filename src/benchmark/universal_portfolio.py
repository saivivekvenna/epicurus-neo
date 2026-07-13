"""Patient-agnostic, label-free mutation portfolio candidates for Miller calibration."""

from __future__ import annotations

import hashlib

import pandas as pd


def _tie(row: pd.Series) -> str:
    raw = "|".join(str(row.get(c, "")) for c in ("mutation_id", "mutant_peptide", "hla_allele"))
    return hashlib.sha256(raw.encode()).hexdigest()


def _rank_utility(values: pd.Series) -> pd.Series:
    """Map higher-is-better values to [0,1]; missing evidence receives zero utility."""
    numeric = pd.to_numeric(values, errors="coerce")
    present = numeric.notna()
    out = pd.Series(0.0, index=values.index)
    if present.any():
        ranks = numeric[present].rank(method="average", ascending=False)
        n = int(present.sum())
        out.loc[present] = 1.0 if n == 1 else 1.0 - (ranks - 1.0) / (n - 1.0)
    return out


def mutation_representatives(frame: pd.DataFrame, score_col: str) -> pd.DataFrame:
    """Best scored peptide/HLA route per mutation, deterministically."""
    required = {"mutation_id", "mutant_peptide", "hla_allele", score_col}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"portfolio input missing columns: {sorted(missing)}")
    work = frame.copy()
    work[score_col] = pd.to_numeric(work[score_col], errors="coerce")
    work = work[work[score_col].notna()].copy()
    work["_tie"] = work.apply(_tie, axis=1)
    return (
        work.sort_values([score_col, "_tie"], ascending=[False, True], kind="mergesort")
        .drop_duplicates("mutation_id", keep="first")
        .drop(columns="_tie")
        .reset_index(drop=True)
    )


def select_rank_fusion_cap1(frame: pd.DataFrame, score_cols: tuple[str, ...], *, k: int = 20) -> pd.DataFrame:
    """Fuse within-patient mutation ranks, then select one best fused route per mutation."""
    if not score_cols:
        raise ValueError("rank fusion requires at least one score")
    work = frame.copy()
    utilities = []
    for col in score_cols:
        if col not in work:
            raise ValueError(f"portfolio input missing score column: {col}")
        utility = _rank_utility(work[col])
        work[f"_utility_{col}"] = utility
        utilities.append(utility)
    work["portfolio_score"] = pd.concat(utilities, axis=1).mean(axis=1)
    reps = mutation_representatives(work, "portfolio_score")
    return reps.head(k).reset_index(drop=True)


def select_evidence_lane_portfolio(frame: pd.DataFrame, score_cols: tuple[str, ...], *, k: int = 20) -> pd.DataFrame:
    """Round-robin independent evidence lanes, skipping unavailable lanes and duplicate mutations."""
    lanes: list[tuple[str, list[pd.Series]]] = []
    for col in score_cols:
        if col not in frame:
            continue
        reps = mutation_representatives(frame, col)
        if len(reps):
            lanes.append((col, [row for _, row in reps.iterrows()]))
    if not lanes:
        return frame.iloc[0:0].copy()
    cursors = [0] * len(lanes)
    selected: list[pd.Series] = []
    seen: set[str] = set()
    while len(selected) < k:
        progress = False
        for lane_i, (lane_name, lane) in enumerate(lanes):
            while cursors[lane_i] < len(lane):
                row = lane[cursors[lane_i]]
                cursors[lane_i] += 1
                mutation = str(row["mutation_id"])
                if mutation in seen:
                    continue
                picked = row.copy()
                picked["portfolio_lane"] = lane_name
                selected.append(picked)
                seen.add(mutation)
                progress = True
                break
            if len(selected) >= k:
                break
        if not progress:
            break
    return pd.DataFrame(selected).reset_index(drop=True) if selected else frame.iloc[0:0].copy()
