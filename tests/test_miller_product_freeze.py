"""Universal product portfolio freeze: label isolation, PRIME direction, and Hu parity."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from benchmark.miller_product_freeze import ARM_IDS, build_selections, freeze_product_selections


ROOT = Path(__file__).resolve().parents[1]


def _tiny() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "patient_id": ["P", "P"],
            "mutation_id": ["1:1:A:T", "1:2:C:G"],
            "mutant_peptide": ["AAAAAAAAA", "CCCCCCCCC"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:01"],
            "source_variant_type": ["MISSENSE", "MISSENSE"],
            "expression_tpm": [10.0, 10.0],
            "tumor_vaf": [0.3, 0.3],
            "prime_rank": [1.0, 90.0],
            "epicurus": [0.2, 0.9],
            "mixmhcpred_rank": [2.0, 2.0],
        }
    )


def test_prime_percentile_direction_is_lower_rank_first():
    _, selections = build_selections(_tiny(), k=2)
    assert selections["prime_plain"]["mutation_id"].tolist() == ["1:1:A:T", "1:2:C:G"]
    assert selections["prime_mutation_cap1"]["mutation_id"].tolist() == ["1:1:A:T", "1:2:C:G"]


def test_label_column_is_rejected_before_scoring():
    raw = _tiny().assign(label=[1, 0])
    with pytest.raises(ValueError, match="labels reached product inference"):
        build_selections(raw)


def test_freeze_writes_all_ordered_arms_and_hashes(tmp_path):
    meta = freeze_product_selections(_tiny(), tmp_path, k=2)
    assert meta["labels_opened"] is False
    assert set(meta["arms"]) == set(ARM_IDS)
    assert set(meta["sha256"]) == {f"product_select_{arm}.csv" for arm in ARM_IDS}
    for arm in ARM_IDS:
        path = tmp_path / meta["arms"][arm]["selection_file"]
        assert path.exists()
        assert pd.read_csv(path)["selection_rank"].tolist() == [1, 2]


def test_hu287_shipped_product_exact_ordered_parity():
    raw = pd.read_csv(ROOT / "data/raw/miller_ipv/hu_287/freeze/universe.csv", low_memory=False)
    expected = json.loads(
        (ROOT / "artifacts/milestone_7_decision/end_to_end_product/FROZEN_PIPELINE.json").read_text()
    )["patients"]["Hu_287"]
    _, selections = build_selections(raw)
    shipped = selections["shipped_epicurus_product"]
    assert shipped["candidate_id"].astype(str).tolist() == expected["product_selected_candidate_ids"]
    assert shipped["mutation_id"].astype(str).tolist() == expected["product_selected_mutation_ids"]
