from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.features import (
    assert_no_banned_features,
    build_feature_matrix,
    feature_columns,
)
from epicurus_neo.m6.loso import loso_folds


def test_no_banned_feature_reaches_any_tier():
    frame = load_label_frame()
    for tier in ("core", "contrastive", "presentation"):
        assert_no_banned_features(feature_columns(build_feature_matrix(frame, tier)))


def test_loso_training_never_sees_the_held_out_study():
    frame = load_label_frame()
    for fold in loso_folds(frame):
        assert fold.held_out_study not in set(fold.train.study_id)
        # The core feature matrix carries no study_id/patient_id/label as a usable feature.
        cols = feature_columns(build_feature_matrix(fold.train, "core"))
        assert not ({"study_id", "patient_id", "label", "cancer_type"} & set(cols))
