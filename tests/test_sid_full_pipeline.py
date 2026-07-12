from __future__ import annotations

import pandas as pd
import pytest

from event_b.sid_full_pipeline import assert_label_blind, evaluate_frozen, freeze_mutation_topk, freeze_portfolio


def _frame() -> pd.DataFrame:
    return pd.DataFrame({
        "candidate_id": ["a1", "a2", "b1", "c1"],
        "mutation_id": ["A", "A", "B", "C"],
        "score": [0.9, 0.8, 0.7, 0.6],
        "chosen": [True, False, True, False],
        "route_rank": [1, 4, 2, 3],
    })


def test_label_column_fails_closed():
    with pytest.raises(ValueError, match="evaluation label"):
        assert_label_blind(_frame().assign(label="POSITIVE"))


def test_mutation_topk_deduplicates_routes_before_selection():
    frozen = freeze_mutation_topk(_frame(), "score", k=2)
    assert frozen["selected_mutation_ids"] == ["A", "B"]
    assert frozen["n_ranked_mutations"] == 3
    assert evaluate_frozen(frozen, {"A", "C"})["hits_at_20"] == 1


def test_portfolio_preserves_policy_selected_routes():
    frozen = freeze_portfolio(_frame(), "chosen", "route_rank", k=20)
    assert frozen["selected_mutation_ids"] == ["A", "B"]
    assert frozen["n_unique_selected_mutations"] == 2
