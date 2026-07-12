"""Tests for the RNA-expression ranking-policy forms and the portfolio reserve selector.

Three candidate roles for RNA expression, evaluated label-blind on development cohorts:
  (a) rank penalty        — expression co-drives the score (demotes low-expression candidates)
  (b) confidence-only     — expression annotates only; ranking stays = genuine PRIME (protected incumbent)
  (c) soft-saturating     — expression demotes ONLY candidates that are BOTH low-presentation AND
                            low-expression (strong presenters are protected)

Plus a Pareto/portfolio selector that reserves top-presentation candidates across expression strata and
predictor-disagreement to preserve reachability, without tuning any constant to a specific patient.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from benchmark.expression_policy import (
    EXPRESSION_STRATA,
    expression_confidence_annotation,
    expr_penalty_score,
    no_regression_verdict,
    prime_only_score,
    select_portfolio_reserved,
    soft_saturating_score,
    within_patient_percentile,
)


def _frame(prime, expr, patient="p"):
    return pd.DataFrame({
        "patient_id": [patient] * len(prime),
        "mutant_peptide": [f"PEP{i:04d}" for i in range(len(prime))],
        "hla_allele": ["HLA-A*02:01"] * len(prime),
        "prime": prime,
        "expr": expr,
    })


# ---------------------------------------------------------------------------
# within-patient percentile
# ---------------------------------------------------------------------------
def test_within_patient_percentile_orientation():
    f = _frame(prime=[0.1, 0.5, 0.9], expr=[1, 10, 100])
    pp = within_patient_percentile(f, "prime", higher_better=False)  # lower prime = better
    assert pp.iloc[0] > pp.iloc[2]  # 0.1 is the best presenter


# ---------------------------------------------------------------------------
# (b) confidence-only == protected PRIME incumbent
# ---------------------------------------------------------------------------
def test_prime_only_score_ignores_expression():
    a = prime_only_score(_frame(prime=[0.1, 0.5, 0.9], expr=[1, 2, 3]))
    b = prime_only_score(_frame(prime=[0.1, 0.5, 0.9], expr=[300, 200, 100]))
    assert list(a) == list(b)  # expression permutation does not change the score


def test_confidence_annotation_labels_strata_without_touching_rank():
    f = _frame(prime=[0.1, 0.2, 0.3, 0.4], expr=[0.1, 5, 50, 5000])
    ann = expression_confidence_annotation(f)
    assert set(ann.unique()) <= set(EXPRESSION_STRATA)
    # highest-expression candidate is in the top stratum
    assert ann.iloc[3] == EXPRESSION_STRATA[-1]


# ---------------------------------------------------------------------------
# (a) rank penalty demotes low-expression candidates
# ---------------------------------------------------------------------------
def test_expr_penalty_can_demote_a_low_expression_strong_presenter():
    # candidate 0 is the best presenter but lowest expression; candidate 2 is a weak presenter, high expr
    f = _frame(prime=[0.01, 0.5, 0.9], expr=[1, 50, 5000])
    s = expr_penalty_score(f)
    # under the penalty the high-expression weak presenter is pulled up relative to prime-only
    prime = prime_only_score(f)
    assert (s.iloc[2] - s.iloc[0]) > (prime.iloc[2] - prime.iloc[0])


# ---------------------------------------------------------------------------
# (c) soft-saturating PROTECTS strong presenters
# ---------------------------------------------------------------------------
def test_soft_saturating_does_not_demote_a_strong_presenter():
    # candidate 0: strong presenter (top), bottom expression -> must NOT be demoted
    f = _frame(prime=[0.01, 0.4, 0.5, 0.6, 0.9], expr=[1, 50, 60, 70, 5000])
    s = soft_saturating_score(f)
    pp = prime_only_score(f)
    assert s.iloc[0] == pp.iloc[0]  # strong presenter protected


def test_soft_saturating_demotes_a_weak_presenter_that_is_also_low_expression():
    # candidate 0 is in the bottom presentation AND bottom expression stratum -> demoted
    f = _frame(prime=[0.9, 0.4, 0.3, 0.2, 0.1], expr=[1, 50, 60, 70, 5000])
    s = soft_saturating_score(f)
    pp = prime_only_score(f)
    worst = pp.idxmin()
    assert s.loc[worst] < pp.loc[worst]


# ---------------------------------------------------------------------------
# portfolio reserve selector
# ---------------------------------------------------------------------------
def test_portfolio_reserves_a_low_expression_high_presenter_slot():
    # 22 candidates; the low-expression high-presenter sits just outside PRIME top-20 by construction.
    prime = list(np.linspace(0.02, 0.9, 22))  # ascending prime %rank = descending presentation
    expr = [100.0] * 22
    expr[20] = 0.5  # a bottom-expression candidate at PRIME rank 21
    f = _frame(prime=prime, expr=expr)
    # make candidate 20 a strong presenter so it is a reachability target
    f.loc[20, "prime"] = 0.03
    reserved = select_portfolio_reserved(f, k=20, reserve=2)
    assert len(reserved) == 20
    assert 20 in reserved.index  # the reserved low-expression high-presenter is retained


def test_portfolio_is_deterministic_under_row_permutation():
    prime = list(np.linspace(0.02, 0.9, 25))
    f = _frame(prime=prime, expr=[10.0] * 25)
    a = set(select_portfolio_reserved(f, k=20).index)
    b = set(select_portfolio_reserved(f.iloc[::-1].reset_index(drop=True), k=20).index)
    # same peptides selected regardless of input order
    fa = set(f.loc[sorted(a), "mutant_peptide"])
    fr = f.iloc[::-1].reset_index(drop=True)
    fb = set(fr.loc[sorted(b), "mutant_peptide"])
    assert fa == fb


# ---------------------------------------------------------------------------
# no-regression verdict
# ---------------------------------------------------------------------------
def test_no_regression_verdict_flags_a_drop():
    # per-patient hits: policy strictly below incumbent
    delta = np.array([-1, -2, 0, -1], dtype=float)
    v = no_regression_verdict(delta)
    assert v["regresses"] is True
    assert v["mean_delta"] < 0


def test_no_regression_verdict_passes_when_equal():
    v = no_regression_verdict(np.zeros(5))
    assert v["regresses"] is False
    assert v["mean_delta"] == 0.0


# ---------------------------------------------------------------------------
# Frozen-decision guard (locks the label-blind policy decision against silent drift)
# ---------------------------------------------------------------------------
def test_frozen_expression_policy_config_locks_confidence_only():
    import json
    from pathlib import Path

    cfg_path = Path(__file__).resolve().parents[1] / "configs" / "frozen" / "expression_policy_v1.json"
    if not cfg_path.exists():
        import pytest
        pytest.skip("frozen config not generated yet (run scripts.expression_policy_analysis)")
    cfg = json.loads(cfg_path.read_text())
    assert cfg["decision"]["chosen"] == "confidence_only"
    assert "PRIME" in cfg["protected_incumbent"]
    assert cfg["constants"]["tuned_to_sid"] is False
    assert cfg["constants"]["tuned_to_any_eval_cohort"] is False
    # the rank-penalty form MUST be recorded as regressing at least one development cohort
    regresses = {c: ev["policy_regresses"]["expr_rank_penalty"]
                 for c, ev in cfg["development_evidence"].items()}
    assert any(regresses.values()), "expr_rank_penalty should regress >=1 dev cohort in the frozen record"
