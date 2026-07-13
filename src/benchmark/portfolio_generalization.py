"""Scorer × portfolio crossed benchmark for mutation-resolved candidate universes.

Selection is label-free. Recognition labels are accepted only by ``evaluate_frozen``
after selections have been materialized by the caller.
"""

from __future__ import annotations

from dataclasses import replace
import hashlib
from typing import Iterable

import pandas as pd

from epicurus_neo.evidence_router import (
    DEFAULT_ROUTER_POLICY,
    RouterPolicy,
    route_candidates,
    select_route_aware_topk,
)


def _text(value: object) -> str:
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


def candidate_identity(row: pd.Series) -> str:
    """Stable route identity used in the pre-label freeze payload."""
    fields = (
        _text(row.get("patient_id")),
        _text(row.get("mutation_id")),
        _text(row.get("mutant_peptide")),
        _text(row.get("hla_allele")),
    )
    return "|".join(fields)


def _ordered_rankable(frame: pd.DataFrame, score_col: str) -> pd.DataFrame:
    routed = route_candidates(frame)
    pool = routed[
        routed["router_eligible"].astype(bool) & routed["rankable"].astype(bool)
    ].copy()
    pool[score_col] = pd.to_numeric(pool[score_col], errors="coerce")
    pool = pool[pool[score_col].notna()].copy()
    pool["_tie_key"] = [
        hashlib.md5(
            f"{_text(peptide)}|{_text(hla)}".encode()
        ).hexdigest()
        for peptide, hla in zip(pool["mutant_peptide"], pool["hla_allele"])
    ]
    return pool.sort_values(
        [score_col, "_tie_key"], ascending=[False, True], kind="mergesort"
    )


def select_plain(frame: pd.DataFrame, score_col: str, *, k: int) -> pd.DataFrame:
    """Ordinary top-k over the same router-valid/rankable rows as the selector."""
    return _ordered_rankable(frame, score_col).head(k).copy()


def select_route_aware(
    frame: pd.DataFrame,
    score_col: str,
    *,
    k: int,
    max_per_mutation: int | None = 2,
    include_route_reserves: bool = True,
    base_policy: RouterPolicy = DEFAULT_ROUTER_POLICY,
) -> pd.DataFrame:
    """Apply the frozen selector, optionally disabling reserves for attribution."""
    policy = replace(
        base_policy,
        k=k,
        max_per_mutation=max_per_mutation,
        reserve_routes=(base_policy.reserve_routes if include_route_reserves else ()),
        reserve_per_route=(base_policy.reserve_per_route if include_route_reserves else 0),
        max_reserve=(base_policy.max_reserve if include_route_reserves else 0),
    )
    routed = route_candidates(frame, policy)
    selected = select_route_aware_topk(
        routed, score_column=score_col, policy=policy
    )
    return selected[selected["route_selected"].astype(bool)].sort_values(
        "route_rank", kind="mergesort"
    ).copy()


def freeze_selection(selection: pd.DataFrame, *, k: int) -> dict:
    """Serialize a selection without consulting recognition labels."""
    mutation_ids = selection["mutation_id"].astype(str).tolist()
    identities = [candidate_identity(row) for _, row in selection.iterrows()]
    return {
        "requested_k": int(k),
        "n_selected": int(len(selection)),
        "saturated": bool(len(selection) == k),
        "selected_candidate_ids": identities,
        "selected_mutation_ids": mutation_ids,
        "n_unique_selected_mutations": int(len(set(mutation_ids))),
        "duplicate_slot_burden": int(len(mutation_ids) - len(set(mutation_ids))),
    }


def crossed_selections(
    frame: pd.DataFrame,
    *,
    prime_col: str,
    epicurus_col: str,
    k: int = 20,
    max_per_mutation: int | None = 2,
) -> dict[str, dict]:
    """Freeze the primary 2×2 scorer/selector comparison plus cap-only controls."""
    required = {
        "patient_id",
        "mutation_id",
        "mutant_peptide",
        "hla_allele",
        prime_col,
        epicurus_col,
    }
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(f"mutation-resolved portfolio benchmark missing columns: {missing}")

    selections = {
        "prime_plain": select_plain(frame, prime_col, k=k),
        "prime_route_aware": select_route_aware(
            frame, prime_col, k=k, max_per_mutation=max_per_mutation
        ),
        "epicurus_plain": select_plain(frame, epicurus_col, k=k),
        "epicurus_route_aware": select_route_aware(
            frame, epicurus_col, k=k, max_per_mutation=max_per_mutation
        ),
        "prime_cap_only": select_route_aware(
            frame,
            prime_col,
            k=k,
            max_per_mutation=max_per_mutation,
            include_route_reserves=False,
        ),
        "epicurus_cap_only": select_route_aware(
            frame,
            epicurus_col,
            k=k,
            max_per_mutation=max_per_mutation,
            include_route_reserves=False,
        ),
    }
    return {name: freeze_selection(value, k=k) for name, value in selections.items()}


def evaluate_frozen(selection: dict, positives: Iterable[str]) -> dict:
    """Join recognition outcomes only after a label-free selection has been frozen."""
    positive_set = {str(value) for value in positives}
    selected = {str(value) for value in selection["selected_mutation_ids"]}
    hits = sorted(selected & positive_set)
    out = dict(selection)
    out.update(
        {
            "n_recognized_mutations": len(positive_set),
            "hits_at_k_unique_mutations": len(hits),
            "recall_at_k": round(len(hits) / len(positive_set), 4)
            if positive_set
            else None,
            "hit_mutation_ids": hits,
        }
    )
    return out


def paired_deltas(evaluated: dict[str, dict]) -> dict[str, int]:
    h = {k: v["hits_at_k_unique_mutations"] for k, v in evaluated.items()}
    return {
        "selector_delta_on_prime": h["prime_route_aware"] - h["prime_plain"],
        "selector_delta_on_epicurus": h["epicurus_route_aware"] - h["epicurus_plain"],
        "score_delta_under_plain": h["epicurus_plain"] - h["prime_plain"],
        "score_delta_under_route_aware": h["epicurus_route_aware"]
        - h["prime_route_aware"],
        "reserve_delta_on_prime": h["prime_route_aware"] - h["prime_cap_only"],
        "reserve_delta_on_epicurus": h["epicurus_route_aware"]
        - h["epicurus_cap_only"],
    }
