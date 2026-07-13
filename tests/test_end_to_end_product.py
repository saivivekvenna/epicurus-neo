from __future__ import annotations

import pandas as pd
import pytest

from benchmark.end_to_end_product import (
    evaluate_frozen_pipeline,
    freeze_product_pipeline,
)
from epicurus_neo.product import normalize_product_candidates


def _frame() -> pd.DataFrame:
    rows = []
    for mutation, rank, reads in (("m1", 1.0, 5), ("m2", 2.0, 5), ("m3", 3.0, 0)):
        for i in range(3):
            rows.append(
                {
                    "patient_id": "p1",
                    "mutation_id": mutation,
                    "gene_symbol": mutation,
                    "mutant_peptide": ("ACDEFGHIK", "CDEFGHIKL", "DEFGHIKLM")[i],
                    "hla_allele": f"HLA-A*0{i + 1}:01",
                    "source_variant_type": "SNV",
                    "mhc_class": "I",
                    "expression_tpm": 10.0,
                    "rna_depth": 20,
                    "rna_mutant_reads": reads,
                    "tumor_vaf": 0.2,
                    "mixmhcpred_rank": rank,
                    "prime_rank": rank,
                }
            )
    return normalize_product_candidates(pd.DataFrame(rows))


def test_full_product_freeze_has_no_label_fields_and_respects_caps():
    frozen = freeze_product_pipeline(_frame())
    assert frozen["counts"]["selected_routes"] == 4
    assert frozen["counts"]["selected_unique_mutations"] == 2
    assert frozen["removal_reasons"] == {"NO_MUTANT_RNA_SUPPORT": 3}
    assert "hits" not in str(frozen)


def test_evaluation_identifies_product_eligibility_loss():
    frozen = freeze_product_pipeline(_frame())
    result = evaluate_frozen_pipeline(frozen, {"m2", "m3"})
    assert result["stage_reachability"]["generated"]["n"] == 2
    assert result["stage_reachability"]["product_eligible"]["n"] == 1
    assert result["last_reached_stage_by_positive"]["m3"] == "deterministic_valid"


def test_pipeline_rejects_research_labels_at_inference_boundary():
    frame = _frame()
    frame["label"] = "POSITIVE"
    with pytest.raises(ValueError, match="labels reached"):
        freeze_product_pipeline(frame)
