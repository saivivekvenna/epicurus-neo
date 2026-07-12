"""Tests for the biology-first WES/RNA gate investigation core.

The gate is a pure SELECTION change on top of the frozen Epicurus base order: a
label-blind demotion predicate removes top-k candidates that fail a biological
prerequisite, and freed slots are backfilled by base order (the reranker is
never touched). These tests pin the swap accounting, within-patient binning, the
matched-random control, and the paired bootstrap — the places bugs hide.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from benchmark.improve_wes_rna_gate import (
    apply_demotion_gate,
    matched_random_removal,
    paired_bootstrap,
    partial_effect,
    within_patient_bin,
)


def _toy():
    # one patient, ranks by score desc: neg,neg,pos,neg,pos
    return pd.DataFrame(
        {
            "patient_id": ["A"] * 5,
            "score": [10.0, 9.0, 8.0, 7.0, 6.0],
            "pos": [0, 0, 1, 0, 1],
        }
    )


def test_demotion_backfills_positive_from_challengers():
    f = _toy()
    f["demote"] = [True, False, False, False, False]  # demote rank-1 negative
    out = apply_demotion_gate(f, k=2)
    # base top2 = [neg,neg] = 0 hits; after removing rank1, top2 = [9(neg),8(pos)] = 1
    assert out["deltas"]["A"] == 1
    assert out["mean_delta"] == 1.0


def test_demotion_of_a_positive_is_counted_as_harm():
    f = _toy()
    f["demote"] = [False, False, True, False, False]  # demote the rank-3 positive
    out = apply_demotion_gate(f, k=3)
    # base top3 = [neg,neg,pos] = 1 hit; removing the pos -> top3 = [neg,neg,neg] = 0
    assert out["deltas"]["A"] == -1


def test_demoting_a_challenger_outside_topk_is_a_noop():
    f = _toy()
    f["demote"] = [False, False, False, False, True]  # demote rank-5 (already excluded)
    out = apply_demotion_gate(f, k=2)
    assert out["deltas"]["A"] == 0


def test_demoted_challenger_is_never_backfilled():
    # If a demoted candidate sits at rank>k it must not be pulled in to refill.
    f = _toy()
    # demote rank1 (neg, top-k) AND rank3 (pos, challenger). The freed slot must skip
    # the demoted positive and take rank2.. -> top2 = [9(neg),7(neg)] = 0 hits.
    f["demote"] = [True, False, True, False, False]
    out = apply_demotion_gate(f, k=2)
    assert out["deltas"]["A"] == 0  # not +1, because the demoted pos cannot backfill


def test_within_patient_bin_is_relative_and_missing_is_own_bin():
    f = pd.DataFrame(
        {
            "patient_id": ["A", "A", "A", "A", "B", "B"],
            "x": [1.0, 2.0, 3.0, np.nan, 100.0, 200.0],
        }
    )
    b = within_patient_bin(f, "x", n_bins=2)
    # within A: 1,2 -> low bin 0 ; 3 -> high bin 1 ; nan -> -1 ; B: 100->0, 200->1
    assert b.tolist() == [0, 0, 1, -1, 0, 1]


def test_partial_effect_reports_per_bin_rate_and_support():
    f = pd.DataFrame(
        {
            "patient_id": ["A"] * 4 + ["B"] * 4,
            "x": [1, 2, 3, 4, 1, 2, 3, 4],
            "pos": [0, 0, 1, 1, 0, 0, 1, 1],
        }
    )
    pe = partial_effect(f, "x", n_bins=2)
    lo = [r for r in pe if r["bin"] == 0][0]
    hi = [r for r in pe if r["bin"] == 1][0]
    assert lo["pos_rate"] == 0.0 and lo["n"] == 4
    assert hi["pos_rate"] == 1.0 and hi["n"] == 4
    assert hi["n_patients"] == 2


def test_matched_random_removal_matches_demote_count_and_is_seed_deterministic():
    f = _toy()
    f["demote"] = [True, False, False, False, False]  # 1 removed from top-k
    a = matched_random_removal(f, k=2, seed=0, reps=50)
    b = matched_random_removal(f, k=2, seed=0, reps=50)
    assert a["mean_delta"] == b["mean_delta"]  # deterministic under fixed seed
    assert a["n_removed_per_patient"]["A"] == 1


def test_paired_bootstrap_ci_brackets_mean():
    deltas = {f"p{i}": v for i, v in enumerate([1, 1, 1, 0, 0, -1, 1, 0, 1, 1])}
    out = paired_bootstrap(deltas, reps=500, seed=0)
    assert out["lo"] <= out["mean"] <= out["hi"]
    assert abs(out["mean"] - np.mean(list(deltas.values()))) < 1e-9
