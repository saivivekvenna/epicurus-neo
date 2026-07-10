from epicurus_neo.m6.confounds import prevalence_by_study, study_only_classifier
from epicurus_neo.m6.dataset import load_label_frame


def test_prevalence_varies_sharply_by_study():
    table = prevalence_by_study(load_label_frame()).set_index("study_id")
    assert int(table.n.sum()) == 965
    # Study prevalence is heterogeneous — the core confound M6A must expose.
    assert table.positive_rate.max() - table.positive_rate.min() > 0.1


def test_study_only_classifier_reports_confound_strength():
    result = study_only_classifier(load_label_frame())
    assert 0.0 <= result["accuracy"] <= 1.0
    assert 0.0 <= result["majority_rate"] <= 1.0
    assert "per_study" in result
    # Features carry substantial study identity (the confound the design warns about).
    assert result["accuracy"] > result["majority_rate"]
