"""Tests for the Sid identical-input benchmark leakage guard + universe (src/event_b/sid_benchmark.py).

The critical test: candidate generation restricted to the positive TARGETS must FAIL the guard.
"""

from __future__ import annotations

import pandas as pd
import pytest

from event_b.sid_benchmark import (
    GenerationLeakageError,
    assert_generation_label_blind,
    eligible_universe_ids,
    hudson_positive_variant_ids,
    load_variant_universe,
    mutation_hits_at_k,
)


def test_universe_is_label_blind_and_contains_targets():
    u = load_variant_universe()
    assert len(u) == 200
    elig = eligible_universe_ids(u)
    assert 100 < len(elig) < 200  # ~147 eligible
    pos = hudson_positive_variant_ids(u)
    assert pos == {
        "ASPM-chr1-197102716",
        "DYNC1H1-chr14-101980529",
        "MAP2-chr2-209694772",
    }
    # eligibility never references the recognized-gene list
    assert pos & elig  # positives happen to be eligible, but not BECAUSE they are positive


def test_target_conditioned_generation_fails_guard():
    """Generating only the 3 positive TARGETS (the old leakage) must raise."""
    pos = hudson_positive_variant_ids()
    with pytest.raises(GenerationLeakageError):
        assert_generation_label_blind(pos)  # positives only -> label-conditioned


def test_incomplete_universe_fails_guard():
    elig = eligible_universe_ids()
    half = set(list(elig)[: len(elig) // 2])
    with pytest.raises(GenerationLeakageError):
        assert_generation_label_blind(half)  # < 95% coverage


def test_complete_universe_passes_guard():
    elig = eligible_universe_ids()
    pos = hudson_positive_variant_ids()
    stats = assert_generation_label_blind(elig)
    assert stats["coverage"] == 1.0
    assert stats["positives_in_universe"] == 3
    assert stats["positives_covered"] == 3
    assert pos <= elig


def test_mutation_hits_dedup_peptides():
    # two peptides for the same positive variant must count as ONE mutation hit
    ranked = pd.DataFrame({
        "variant_id": ["POS-1", "POS-1", "NEG-1", "NEG-2"],
        "score": [0.9, 0.8, 0.5, 0.4],
    })
    out = mutation_hits_at_k(ranked, {"POS-1"}, k=20, ascending=False)
    assert out["hits_at_k"] == 1
    assert out["n_variants_ranked"] == 3  # deduped to 3 variants, not 4 rows
