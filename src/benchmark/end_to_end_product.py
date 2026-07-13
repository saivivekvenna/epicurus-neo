"""Label-isolated evaluation helpers for the canonical Epicurus product path."""

from __future__ import annotations

import hashlib
from typing import Iterable

import pandas as pd

from epicurus_neo.product import InferenceConfig, score_product_candidates


FORBIDDEN_LABEL_COLUMNS = frozenset(
    {"label", "y", "is_recognized", "recognized_label", "hudson_label"}
)


def assert_label_free(frame: pd.DataFrame) -> None:
    present = FORBIDDEN_LABEL_COLUMNS & set(frame.columns)
    if present:
        raise ValueError(f"evaluation labels reached product inference: {sorted(present)}")


def _tie_key(row: pd.Series) -> str:
    return hashlib.md5(
        f"{row.get('mutant_peptide', '')}|{row.get('hla_allele', '')}".encode()
    ).hexdigest()


def _prime_selection(
    scored: pd.DataFrame, *, k: int, max_per_mutation: int | None
) -> list[str]:
    pool = scored[scored["deterministic_gate_pass"].astype(bool)].copy()
    pool["recognition_score"] = pd.to_numeric(pool["recognition_score"], errors="coerce")
    pool = pool[pool["recognition_score"].notna()].copy()
    pool["_tie"] = pool.apply(_tie_key, axis=1)
    ordered = pool.sort_values(
        ["recognition_score", "_tie"], ascending=[False, True], kind="mergesort"
    )
    selected: list[str] = []
    counts: dict[str, int] = {}
    for row in ordered.itertuples():
        mutation = str(row.mutation_id)
        if max_per_mutation is not None and counts.get(mutation, 0) >= max_per_mutation:
            continue
        selected.append(mutation)
        counts[mutation] = counts.get(mutation, 0) + 1
        if len(selected) == k:
            break
    return selected


def freeze_product_pipeline(
    frame: pd.DataFrame, config: InferenceConfig = InferenceConfig()
) -> dict:
    """Run the shipped product path and serialize stage membership without labels."""
    assert_label_free(frame)
    scored = score_product_candidates(frame, config)
    assert_label_free(scored)
    selected = scored[scored["selected"].astype(bool)].sort_values("rank")

    def mutations(mask: pd.Series | None = None) -> list[str]:
        rows = scored if mask is None else scored[mask]
        return sorted(set(rows["mutation_id"].astype(str)))

    product_selected = selected["mutation_id"].astype(str).tolist()
    return {
        "config": {
            "k": config.k,
            "max_per_mutation": config.max_per_mutation,
            "max_per_gene": config.max_per_gene,
            "max_per_hla": config.max_per_hla,
            "core_threshold": config.core_threshold,
            "supporting_threshold": config.supporting_threshold,
            "apply_validity_gate": config.apply_validity_gate,
        },
        "stages": {
            "generated": mutations(),
            "deterministic_valid": mutations(scored["deterministic_gate_pass"].astype(bool)),
            "product_eligible": mutations(scored["eligible"].astype(bool)),
            "selected": sorted(set(product_selected)),
        },
        "counts": {
            "candidate_rows": int(len(scored)),
            "generated_mutations": int(scored["mutation_id"].nunique()),
            "deterministic_valid_rows": int(scored["deterministic_gate_pass"].sum()),
            "product_eligible_rows": int(scored["eligible"].sum()),
            "selected_routes": int(len(selected)),
            "selected_unique_mutations": int(selected["mutation_id"].nunique()),
            "duplicate_slot_burden": int(len(selected) - selected["mutation_id"].nunique()),
        },
        "removal_reasons": {
            str(reason): int(count)
            for reason, count in scored.loc[
                ~scored["eligible"].astype(bool), "exclusion_reason"
            ].value_counts().items()
        },
        "product_selected_mutation_ids": product_selected,
        "product_selected_candidate_ids": selected["candidate_id"].astype(str).tolist(),
        "prime_plain_selected_mutation_ids": _prime_selection(
            scored, k=config.k, max_per_mutation=None
        ),
        "prime_cap2_selected_mutation_ids": _prime_selection(
            scored, k=config.k, max_per_mutation=config.max_per_mutation
        ),
    }


def evaluate_frozen_pipeline(frozen: dict, positives: Iterable[str]) -> dict:
    """Join recognized mutation IDs after the full stage/selection freeze."""
    pos = {str(value) for value in positives}
    stage_hits = {}
    for stage, mutation_ids in frozen["stages"].items():
        hit = sorted(pos & set(mutation_ids))
        stage_hits[stage] = {"n": len(hit), "of": len(pos), "ids": hit}

    losses = {}
    stage_order = ["generated", "deterministic_valid", "product_eligible", "selected"]
    for mutation in sorted(pos):
        last = "not_generated"
        for stage in stage_order:
            if mutation in set(frozen["stages"][stage]):
                last = stage
            else:
                break
        losses[mutation] = last

    def hits(key: str) -> dict:
        ids = sorted(pos & set(frozen[key]))
        return {"n": len(ids), "of": len(pos), "ids": ids}

    return {
        "stage_reachability": stage_hits,
        "last_reached_stage_by_positive": losses,
        "product_hits_at_20": hits("product_selected_mutation_ids"),
        "prime_plain_hits_at_20": hits("prime_plain_selected_mutation_ids"),
        "prime_cap2_hits_at_20": hits("prime_cap2_selected_mutation_ids"),
    }
