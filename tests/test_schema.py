import pandas as pd

from epicurus_neo.schema import add_normalized_columns, normalize_hla, supervised_rows, validate_schema


def _frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": ["a", "b", "c"],
            "source_dataset": ["toy", "toy", "toy"],
            "study_id": ["s1", "s1", "s2"],
            "patient_id": ["p1", "p1", "p2"],
            "hla_allele": ["A*02:01", "HLA-A*02:01", "B*07:02"],
            "mutant_peptide": [" siinfekla ", "AAAA", "BBBB"],
            "wildtype_peptide": ["SIINFEKLT", "AAAT", "BBBT"],
            "label": ["positive", "negative", "unknown"],
            "label_weight": [1.0, 1.0, 0.0],
            "assay_type": ["elispot", "elispot", "not_tested"],
        }
    )


def test_schema_validation_and_supervised_filter():
    frame = _frame()
    report = validate_schema(frame)
    assert report.ok
    assert len(supervised_rows(frame)) == 2


def test_normalized_columns():
    norm = add_normalized_columns(_frame())
    assert norm.loc[0, "mutant_peptide_norm"] == "SIINFEKLA"
    assert norm.loc[0, "hla_allele_norm"] == "HLA-A*02:01"
    assert norm.loc[0, "mutant_hla_key"] == "SIINFEKLA|HLA-A*02:01"


def test_normalize_hla_canonicalizes_compact_four_digit_alleles():
    assert normalize_hla("A0201") == "HLA-A*02:01"
    assert normalize_hla("HLA-B0702") == "HLA-B*07:02"
