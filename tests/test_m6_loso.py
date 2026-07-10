from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.loso import loso_folds


def test_loso_holds_out_each_study_exactly_once():
    folds = loso_folds(load_label_frame())
    assert [f.held_out_study for f in folds] == [
        "braun_rcc_2025",
        "hu_neovax_2021",
        "mkras_vax_2026",
        "pdac_neovax_2023",
    ]
    for fold in folds:
        # Evaluation is exactly the held-out study; training excludes it.
        assert set(fold.evaluation.study_id.unique()) == {fold.held_out_study}
        assert fold.held_out_study not in set(fold.train.study_id.unique())
        # No patient or candidate leaks across the split.
        assert not set(fold.train.patient_id) & set(fold.evaluation.patient_id)
        assert not set(fold.train.candidate_id) & set(fold.evaluation.candidate_id)
        assert len(fold.train) + len(fold.evaluation) == 965


def test_loso_on_two_study_subset_yields_two_folds():
    frame = load_label_frame()
    subset = frame[frame.study_id.isin(["hu_neovax_2021", "pdac_neovax_2023"])]
    assert len(loso_folds(subset)) == 2
