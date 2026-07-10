from epicurus_neo.m6.event_a import load_event_a_frame


def test_event_a_frame_grounded_counts():
    frame = load_event_a_frame()
    # IMPROVE Event-A: 17,082 pre-existing-reactivity observations, all class-I short epitopes.
    assert len(frame) == 17082
    assert (frame.study_id == "improve").all()
    assert (frame.mhc_class == "CLASS_I").all()
    # Binary recognition label: 458 POSITIVE vs 16,624 TESTED_NEGATIVE.
    assert int((frame.label == 1).sum()) == 458
    assert int((frame.label == 0).sum()) == 16624
    # Short minimal epitopes (8-11mer) — this is the teacher's training domain.
    lengths = frame.mutant_peptide.map(len)
    assert lengths.min() == 8
    assert lengths.max() == 11
    # Feature-bearing columns present; never-merge invariant means only Event-A rows are here.
    for col in ("mutant_peptide", "wildtype_peptide", "mhc_class", "hla_allele", "label"):
        assert col in frame.columns
