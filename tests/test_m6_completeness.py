from epicurus_neo.m6.dataset import completeness_report, load_label_frame


def test_completeness_gate_grounded_counts():
    report = completeness_report(load_label_frame())
    assert len(report) == 45
    # Only hu patients clear n_eligible > 20 (selection is non-degenerate).
    assert int(report.ranking_informative.sum()) == 8
    assert report.loc[report.ranking_informative, "study_id"].unique().tolist() == [
        "hu_neovax_2021"
    ]
    # mKRAS is a fixed 6-peptide shared panel: k_patient == 6, never informative.
    mkras = report[report.study_id == "mkras_vax_2026"]
    assert (mkras.k_patient == 6).all()
    assert (~mkras.ranking_informative).all()
    # 38 patients carry a tested negative (HAS_TESTED_NEGATIVE); 7 mKRAS 6/6-responders
    # carry none (NO_TESTED_NEGATIVE). This is a rankability flag, not a denominator-bias
    # claim -- the mKRAS six-peptide panel is a complete shared denominator; those 7 simply
    # responded to all of it, so they are barred from primary top-k (no negative to rank).
    counts = report.denominator_type.value_counts()
    assert counts["HAS_TESTED_NEGATIVE"] == 38
    assert counts["NO_TESTED_NEGATIVE"] == 7
    no_negative = report[report.denominator_type == "NO_TESTED_NEGATIVE"]
    assert (no_negative.study_id == "mkras_vax_2026").all()
    assert (no_negative.n_positive == no_negative.n_candidates).all()
    # k_patient never exceeds 20 and never exceeds the patient's candidate count.
    assert (report.k_patient <= 20).all()
    assert (report.k_patient <= report.n_candidates).all()
