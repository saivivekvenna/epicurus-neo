import pandas as pd

from epicurus_neo.data_manifest import downloadable_sources, load_dataset_manifest
from epicurus_neo.splits import assign_holdout_split, leave_group_out_splits, split_from_column


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": ["a", "b", "c", "d"],
            "source_dataset": ["toy"] * 4,
            "study_id": ["s1", "s1", "s2", "s2"],
            "patient_id": ["p1", "p1", "p2", "p3"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:01", "HLA-B*07:02", "HLA-C*07:01"],
            "mutant_peptide": ["AAA", "BBB", "CCC", "DDD"],
            "wildtype_peptide": ["AAT", "BBT", "CCT", "DDT"],
            "label": ["positive", "negative", "positive", "negative"],
            "label_weight": [1.0] * 4,
            "assay_type": ["synthetic"] * 4,
        }
    )


def test_dataset_manifest_loads_sources():
    sources = load_dataset_manifest("configs/datasets.yml")
    keys = {source.key for source in sources}
    assert "tesla" in keys
    assert any(source.locked_test for source in sources)
    assert any(source.key == "cedar" for source in downloadable_sources(sources))


def test_assign_holdout_split():
    assigned = assign_holdout_split(_frame(), group_col="patient_id", holdout_values=["p2"])
    assert assigned.loc[assigned["patient_id"] == "p2", "split"].unique().tolist() == ["test"]
    assert set(assigned["split"]) == {"train", "test"}


def test_split_from_column_reports_no_exact_leakage_for_disjoint_patient_and_study():
    assigned = assign_holdout_split(_frame(), group_col="patient_id", holdout_values=["p1"])
    split = split_from_column(assigned)
    assert split.name == "split"
    assert not split.leakage.has_leakage


def test_patient_holdout_can_still_leak_study():
    assigned = assign_holdout_split(_frame(), group_col="patient_id", holdout_values=["p3"])
    split = split_from_column(assigned)
    assert split.leakage.has_leakage
    assert split.leakage.shared_studies == ("s2",)


def test_leave_group_out_splits():
    splits = leave_group_out_splits(_frame(), group_col="patient_id", max_splits=2)
    assert len(splits) == 2
    assert all(split.test["patient_id"].nunique() == 1 for split in splits)
