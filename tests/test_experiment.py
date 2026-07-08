import pandas as pd

from epicurus_neo.experiment import grouped_cross_validate, summarize_cross_validation


def _rows(study: str, patient: str, offset: int) -> list[dict]:
    rows = []
    for idx in range(10):
        positive = idx in {1, 6}
        rows.append(
            {
                "candidate_id": f"{study}-{patient}-{idx}",
                "source_dataset": "toy",
                "study_id": study,
                "patient_id": patient,
                "hla_allele": "HLA-A*02:01",
                "mutant_peptide": f"MUT{offset}{idx}AAAA",
                "wildtype_peptide": f"WTX{offset}{idx}AAAA",
                "label": "positive" if positive else "negative",
                "label_weight": 1.0,
                "assay_type": "synthetic",
                "binding_affinity_nm": 20.0 if positive else 500.0,
                "presentation_score": 0.9 if positive else 0.1,
                "expression_tpm": 40.0 if positive else 2.0,
                "foreignness_score": 0.8 if positive else 0.2,
            }
        )
    return rows


def test_grouped_cross_validate_runs_leave_study_out():
    frame = pd.DataFrame(
        _rows("s1", "p1", 1)
        + _rows("s2", "p2", 2)
        + _rows("s3", "p3", 3)
    )
    folds = grouped_cross_validate(frame, group_col="study_id", k=5)
    assert len(folds) == 3
    assert all(fold.status == "ok" for fold in folds)
    summary = summarize_cross_validation(folds)
    assert summary["ok_folds"] == 3
    assert "epicurus_score" in summary["aggregate"]


def test_grouped_cross_validate_blocks_leaky_patient_fold():
    frame = pd.DataFrame(_rows("s1", "p1", 1) + _rows("s1", "p2", 2))
    folds = grouped_cross_validate(frame, group_col="patient_id", k=5)
    assert any(fold.status == "leakage_blocked" for fold in folds)

