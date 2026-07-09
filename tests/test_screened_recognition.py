import numpy as np
import pandas as pd
import pytest

from epicurus_neo.screened_recognition import (
    add_screened_recognition_features,
    aggregate_screened_reference,
)


def _unit(*values: float) -> np.ndarray:
    vector = np.asarray(values, dtype=np.float32)
    return vector / np.linalg.norm(vector)


def test_aggregate_screened_reference_preserves_repeated_evidence():
    reference = pd.DataFrame(
        {
            "mutant_peptide": ["AAAAAAAAA", "AAAAAAAAA", "AAAAAAAAA", "CCCCCCCCC"],
            "hla_allele": ["A0201", "A0201", "A0201", "HLA-A*03:01"],
            "label": ["positive", "negative", "negative", "positive"],
        }
    )

    aggregated = aggregate_screened_reference(reference)

    first = aggregated.loc[aggregated["mutant_peptide"] == "AAAAAAAAA"].iloc[0]
    assert first["hla_allele"] == "HLA-A*02:01"
    assert first["positive_count"] == 1
    assert first["negative_count"] == 2
    assert first["response_rate"] == pytest.approx(1 / 3)
    assert first["hla_family"] == "HLA-A*02"


def test_screened_recognition_features_distinguish_hla_family_and_global_neighbors():
    reference = pd.DataFrame(
        {
            "mutant_peptide": ["AAAAAAAAA", "CCCCCCCCC", "DDDDDDDDD"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:02", "HLA-B*07:02"],
            "label": ["positive", "negative", "negative"],
        }
    )
    query = pd.DataFrame(
        {
            "candidate_id": ["q"],
            "mutant_peptide": ["QQQQQQQQQ"],
            "hla_allele": ["A0201"],
            "label": ["unknown"],
        }
    )
    embeddings = {
        "AAAAAAAAA": _unit(1.0, 0.0, 0.0),
        "CCCCCCCCC": _unit(0.8, 0.6, 0.0),
        "DDDDDDDDD": _unit(0.0, 1.0, 0.0),
        "QQQQQQQQQ": _unit(1.0, 0.1, 0.0),
    }

    scored = add_screened_recognition_features(
        query,
        reference,
        embeddings,
        top_ks=(1, 2),
    )

    assert scored.loc[0, "screened_recognition_hla_reference_count"] == 1
    assert scored.loc[0, "screened_recognition_family_reference_count"] == 2
    assert scored.loc[0, "screened_recognition_global_reference_count"] == 3
    assert scored.loc[0, "screened_recognition_hla_top1_response_rate"] == 1.0
    assert scored.loc[0, "screened_recognition_family_top2_response_rate"] == 0.5
    assert (
        scored.loc[0, "screened_recognition_global_max_positive_similarity"]
        > scored.loc[0, "screened_recognition_global_max_negative_similarity"]
    )


def test_screened_recognition_features_reject_invalid_top_k():
    with pytest.raises(ValueError, match="positive integers"):
        add_screened_recognition_features(
            pd.DataFrame(),
            pd.DataFrame(),
            {},
            top_ks=(0,),
        )
