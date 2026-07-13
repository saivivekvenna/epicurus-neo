from __future__ import annotations

import pandas as pd

from benchmark.portfolio_generalization import (
    crossed_selections,
    evaluate_frozen,
    paired_deltas,
)


def _frame() -> pd.DataFrame:
    rows = []
    for mutation, prime, epicurus in (
        ("m1", 100.0, 10.0),
        ("m2", 50.0, 30.0),
        ("m3", 40.0, 20.0),
    ):
        for i in range(5):
            rows.append(
                {
                    "patient_id": "p1",
                    "mutation_id": mutation,
                    "gene_symbol": mutation,
                    "mutant_peptide": f"ACDEFGH{chr(73 + i)}",
                    "hla_allele": f"HLA-A*0{i + 1}:01",
                    "source_variant_type": "SNV",
                    "mhc_class": "I",
                    "prime": prime - i,
                    "epicurus": epicurus - i,
                    "expression_tpm": 1.0,
                    "binding_percentile_rank": 1.0,
                }
            )
    return pd.DataFrame(rows)


def test_crossed_selector_is_label_free_and_limits_duplicate_slots():
    frozen = crossed_selections(
        _frame(), prime_col="prime", epicurus_col="epicurus", k=6, max_per_mutation=2
    )
    assert frozen["prime_plain"]["n_unique_selected_mutations"] == 2
    assert frozen["prime_route_aware"]["n_unique_selected_mutations"] == 3
    assert frozen["prime_route_aware"]["duplicate_slot_burden"] == 3
    assert all("hit" not in key for arm in frozen.values() for key in arm)


def test_labels_join_only_in_evaluation_and_deltas_attribute_selector():
    frozen = crossed_selections(
        _frame(), prime_col="prime", epicurus_col="epicurus", k=6, max_per_mutation=2
    )
    evaluated = {arm: evaluate_frozen(value, {"m3"}) for arm, value in frozen.items()}
    deltas = paired_deltas(evaluated)
    assert evaluated["prime_plain"]["hits_at_k_unique_mutations"] == 0
    assert evaluated["prime_route_aware"]["hits_at_k_unique_mutations"] == 1
    assert deltas["selector_delta_on_prime"] == 1


def test_cap_only_control_excludes_route_reserve_effects():
    frozen = crossed_selections(
        _frame(), prime_col="prime", epicurus_col="epicurus", k=6, max_per_mutation=2
    )
    assert frozen["prime_cap_only"]["selected_mutation_ids"] == frozen["prime_route_aware"]["selected_mutation_ids"]
