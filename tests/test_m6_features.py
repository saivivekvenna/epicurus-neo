from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.features import (
    assert_no_banned_features,
    build_feature_matrix,
    feature_columns,
)


def test_core_features_are_universal_and_leakage_free():
    frame = load_label_frame()
    matrix = build_feature_matrix(frame, "core")
    columns = feature_columns(matrix)
    assert "seq_len" in columns
    assert "peptide_length" in columns
    assert "is_class_i" in columns
    # Banned identifiers/labels never survive as features.
    for banned in ("label", "study_id", "patient_id", "mhc_class", "hla_allele", "mutant_peptide"):
        assert banned not in columns
    assert_no_banned_features(columns)  # raises if violated
    # Core features are fully populated for all rows (no study-correlated missingness).
    assert matrix[columns].notna().all().all()


def test_contrastive_anchor_features_are_class_gated():
    frame = load_label_frame()
    matrix = build_feature_matrix(frame, "contrastive")
    class_i = matrix.mhc_class == "CLASS_I"
    # Anchor/TCR-face counts are defined only where the class-I register applies...
    assert matrix.loc[~class_i, "mutation_anchor_count"].isna().all()
    # ...and are present for at least some class-I rows that have a paired wildtype.
    assert matrix.loc[class_i, "mutation_anchor_count"].notna().any()
