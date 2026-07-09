import numpy as np
import pandas as pd
import pytest

from epicurus_neo.transfer_ranker import (
    TransferRankerConfig,
    refit_and_score_transfer_ranker,
    select_transfer_ranker,
)


def _frame(source: str, groups: list[str], peptides: list[str], labels: list[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": [f"{source}:{idx}" for idx in range(len(peptides))],
            "source_dataset": source,
            "study_id": source,
            "patient_id": groups,
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:01", "HLA-B*07:02", "HLA-B*07:02"],
            "mutant_peptide": peptides,
            "wildtype_peptide": [""] * len(peptides),
            "label": labels,
            "label_weight": [1.0] * len(peptides),
            "assay_type": ["tcell"] * len(peptides),
            "mhcflurry_presentation_score": [0.9, 0.1, 0.8, 0.2],
        }
    )


def test_transfer_ranker_selects_and_refits_without_external_overlap():
    pytest.importorskip("xgboost")
    external = _frame(
        "external",
        ["p1", "p1", "p2", "p2"],
        ["AAAAAAAAA", "CCCCCCCCC", "DDDDDDDDD", "EEEEEEEEE"],
        ["positive", "negative", "positive", "negative"],
    )
    train = _frame(
        "train",
        ["a", "a", "b", "b"],
        ["FFFFFFFFF", "GGGGGGGGG", "HHHHHHHHH", "IIIIIIIII"],
        ["positive", "negative", "positive", "negative"],
    )
    validation = _frame(
        "validation",
        ["v1", "v1", "v2", "v2"],
        ["KKKKKKKKK", "LLLLLLLLL", "MMMMMMMMM", "NNNNNNNNN"],
        ["positive", "negative", "positive", "negative"],
    )
    target = _frame(
        "target",
        ["t1", "t1", "t2", "t2"],
        ["PPPPPPPPP", "QQQQQQQQQ", "RRRRRRRRR", "SSSSSSSSS"],
        ["unknown"] * 4,
    )
    peptides = (
        external["mutant_peptide"].tolist()
        + train["mutant_peptide"].tolist()
        + validation["mutant_peptide"].tolist()
        + target["mutant_peptide"].tolist()
    )
    embeddings = {
        peptide: np.eye(len(peptides), dtype=np.float32)[idx]
        for idx, peptide in enumerate(peptides)
    }
    configs = (
        TransferRankerConfig(
            strategy="external_only",
            pretrain_estimators=2,
            target_estimators=2,
        ),
        TransferRankerConfig(
            strategy="teacher",
            pretrain_estimators=2,
            target_estimators=2,
        ),
    )

    selection, scored = select_transfer_ranker(
        external,
        train,
        validation,
        embeddings,
        model_name="toy",
        configs=configs,
    )
    target_scored = refit_and_score_transfer_ranker(
        external,
        pd.concat([train, validation], ignore_index=True),
        target,
        embeddings,
        selection,
    )

    assert selection.config in configs
    assert "epicurus_transfer_ranker_score" in scored
    assert target_scored["epicurus_transfer_ranker_score"].notna().all()


def test_transfer_ranker_rejects_external_validation_overlap():
    external = _frame(
        "external",
        ["p1", "p1", "p2", "p2"],
        ["AAAAAAAAA", "CCCCCCCCC", "DDDDDDDDD", "EEEEEEEEE"],
        ["positive", "negative", "positive", "negative"],
    )
    validation = external.copy()

    with pytest.raises(ValueError, match="leakage"):
        select_transfer_ranker(
            external,
            external,
            validation,
            {},
            model_name="toy",
        )
