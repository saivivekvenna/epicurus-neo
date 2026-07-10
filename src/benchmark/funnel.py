"""Candidate-reachability accounting from validated target to selected portfolio."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from statistics import NormalDist
from typing import Any

import numpy as np
import pandas as pd


class ReachabilityStatus(str, Enum):
    REACHED = "reached"
    LOST = "lost"
    NOT_ASSESSED = "not_assessed"


class FunnelStage(str, Enum):
    MUTATION_CALLED = "mutation_called"
    TRANSCRIPT_REPRESENTED = "transcript_represented"
    PEPTIDE_GENERATED = "peptide_generated"
    SURVIVES_GATING = "survives_gating"
    HLA_INCLUDED = "hla_included"
    PRESENTATION_CANDIDATE = "presentation_candidate"
    RANKING_STAGE = "ranking_stage"
    TOP_K = "top_k"


STAGES = tuple(stage.value for stage in FunnelStage)


@dataclass(frozen=True)
class FunnelStageSummary:
    stage: str
    reached: int
    lost_here: int
    cumulative_lost: int
    not_assessed: int
    total_validated_positives: int
    recall: float
    ci_lo: float
    ci_hi: float
    reachability_lower_bound: float
    reachability_upper_bound: float


def _has_status(series: pd.Series, status: ReachabilityStatus) -> pd.Series:
    """Compare str-backed enums without pandas coercing the scalar to text."""
    return series.map(lambda value: value is status)


def _coerce_status(value: Any) -> ReachabilityStatus:
    if isinstance(value, ReachabilityStatus):
        return value
    if value is None or pd.isna(value):
        return ReachabilityStatus.NOT_ASSESSED
    if isinstance(value, (bool, np.bool_)):
        return ReachabilityStatus.REACHED if value else ReachabilityStatus.LOST
    if isinstance(value, (int, np.integer)) and value in {0, 1}:
        return ReachabilityStatus.REACHED if value else ReachabilityStatus.LOST
    text = str(value).strip().lower()
    aliases = {
        "reached": ReachabilityStatus.REACHED,
        "present": ReachabilityStatus.REACHED,
        "true": ReachabilityStatus.REACHED,
        "1": ReachabilityStatus.REACHED,
        "lost": ReachabilityStatus.LOST,
        "absent": ReachabilityStatus.LOST,
        "false": ReachabilityStatus.LOST,
        "0": ReachabilityStatus.LOST,
        "not_assessed": ReachabilityStatus.NOT_ASSESSED,
        "unknown": ReachabilityStatus.NOT_ASSESSED,
        "na": ReachabilityStatus.NOT_ASSESSED,
    }
    try:
        return aliases[text]
    except KeyError as error:
        raise ValueError(f"Unknown reachability status: {value!r}") from error


def validate_reachability_ledger(
    frame: pd.DataFrame,
    *,
    positive_id_col: str = "positive_id",
    patient_col: str = "patient_id",
) -> pd.DataFrame:
    """Validate stage order and return a canonical status ledger.

    A downstream stage may never be reached when an upstream stage is lost or
    unassessed. Once a positive is confirmed lost, downstream unassessed values
    are canonically marked lost because the target is no longer reachable.
    """
    required = {positive_id_col, patient_col, *STAGES}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Reachability ledger missing columns: {sorted(missing)}")
    if frame[positive_id_col].isna().any() or frame[patient_col].isna().any():
        raise ValueError("positive and patient identifiers must be non-null")
    if frame[positive_id_col].duplicated().any():
        duplicates = frame.loc[frame[positive_id_col].duplicated(), positive_id_col].tolist()
        raise ValueError(f"Duplicate validated positive identifiers: {duplicates[:5]}")

    out = frame.copy()
    for stage in STAGES:
        out[stage] = out[stage].map(_coerce_status)

    for index, row in out.iterrows():
        upstream_lost = False
        upstream_unknown = False
        for stage in STAGES:
            status = row[stage]
            if status is ReachabilityStatus.REACHED and (upstream_lost or upstream_unknown):
                reason = "lost" if upstream_lost else "not assessed"
                raise ValueError(
                    f"{row[positive_id_col]!r} reaches {stage} after an upstream stage was {reason}"
                )
            if upstream_lost and status is ReachabilityStatus.NOT_ASSESSED:
                out.at[index, stage] = ReachabilityStatus.LOST
                status = ReachabilityStatus.LOST
            upstream_lost = upstream_lost or status is ReachabilityStatus.LOST
            upstream_unknown = upstream_unknown or status is ReachabilityStatus.NOT_ASSESSED
    return out


def annotate_reachability(
    positives: pd.DataFrame,
    stage_tables: dict[str, pd.DataFrame | None],
    stage_keys: dict[str, tuple[str, ...]],
    stage_complete: dict[str, bool],
    *,
    positive_id_col: str = "positive_id",
    patient_col: str = "patient_id",
) -> pd.DataFrame:
    """Build a ledger by testing validated-positive identities against each stage table.

    A missing stage table means ``not_assessed``. Presence in an incomplete table
    proves reachability, but absence remains ``not_assessed``. Only absence from
    a stage explicitly declared complete means ``lost``. Keys are explicit per
    stage so mutation, transcript, peptide, HLA, and selected-candidate universes
    are never conflated.
    """
    required = {positive_id_col, patient_col}
    missing = required.difference(positives.columns)
    if missing:
        raise ValueError(f"Validated-positive table missing columns: {sorted(missing)}")
    ledger = positives.copy()
    for stage in STAGES:
        table = stage_tables.get(stage)
        keys = stage_keys.get(stage)
        if table is None:
            ledger[stage] = ReachabilityStatus.NOT_ASSESSED
            continue
        if not keys:
            raise ValueError(f"Explicit identity keys are required for supplied stage {stage}")
        missing_positive = set(keys).difference(positives.columns)
        missing_stage = set(keys).difference(table.columns)
        if missing_positive or missing_stage:
            raise ValueError(
                f"Stage {stage} key mismatch: positives missing {sorted(missing_positive)}, "
                f"stage table missing {sorted(missing_stage)}"
            )
        stage_identities = set(
            table.loc[:, list(keys)].astype(str).itertuples(index=False, name=None)
        )
        positive_identities = (
            positives.loc[:, list(keys)].astype(str).itertuples(index=False, name=None)
        )
        absent_status = (
            ReachabilityStatus.LOST
            if stage_complete.get(stage, False)
            else ReachabilityStatus.NOT_ASSESSED
        )
        ledger[stage] = [
            ReachabilityStatus.REACHED if identity in stage_identities else absent_status
            for identity in positive_identities
        ]
    return validate_reachability_ledger(
        ledger, positive_id_col=positive_id_col, patient_col=patient_col
    )


def wilson_ci(successes: int, total: int, *, alpha: float = 0.05) -> tuple[float, float]:
    """Two-sided Wilson interval for peptide-level candidate recall."""
    if total <= 0 or not 0 <= successes <= total:
        raise ValueError("Wilson interval requires 0 <= successes <= total and total > 0")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    z = NormalDist().inv_cdf(1 - alpha / 2)
    proportion = successes / total
    denominator = 1 + z**2 / total
    center = (proportion + z**2 / (2 * total)) / denominator
    half_width = (
        z * np.sqrt(proportion * (1 - proportion) / total + z**2 / (4 * total**2)) / denominator
    )
    return float(center - half_width), float(center + half_width)


def candidate_reachability_funnel(
    frame: pd.DataFrame,
    *,
    positive_id_col: str = "positive_id",
    patient_col: str = "patient_id",
) -> dict[str, Any]:
    """Summarize cumulative reachability with CIs and missing-evidence bounds."""
    ledger = validate_reachability_ledger(
        frame, positive_id_col=positive_id_col, patient_col=patient_col
    )
    total = len(ledger)
    if total == 0:
        raise ValueError("candidate reachability requires at least one validated positive")

    summaries: list[FunnelStageSummary] = []
    previous_reached = total
    for stage in STAGES:
        statuses = ledger[stage]
        reached_mask = _has_status(statuses, ReachabilityStatus.REACHED)
        lost_mask = _has_status(statuses, ReachabilityStatus.LOST)
        not_assessed_mask = _has_status(statuses, ReachabilityStatus.NOT_ASSESSED)
        reached = int(reached_mask.sum())
        cumulative_lost = int(lost_mask.sum())
        not_assessed = int(not_assessed_mask.sum())
        lost_here = (
            previous_reached - reached
            if not_assessed == 0
            else int(
                (
                    lost_mask
                    & (
                        True
                        if stage == STAGES[0]
                        else _has_status(
                            ledger[STAGES[STAGES.index(stage) - 1]],
                            ReachabilityStatus.REACHED,
                        )
                    )
                ).sum()
            )
        )
        lo, hi = wilson_ci(reached, total)
        summaries.append(
            FunnelStageSummary(
                stage=stage,
                reached=reached,
                lost_here=lost_here,
                cumulative_lost=cumulative_lost,
                not_assessed=not_assessed,
                total_validated_positives=total,
                recall=reached / total,
                ci_lo=lo,
                ci_hi=hi,
                reachability_lower_bound=reached / total,
                reachability_upper_bound=(reached + not_assessed) / total,
            )
        )
        previous_reached = reached
    return {
        "unit": "validated_positive",
        "total_validated_positives": total,
        "patients": int(ledger[patient_col].nunique()),
        "stages": [asdict(summary) for summary in summaries],
    }
