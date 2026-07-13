from scripts.miller_generalization_split import build_split


def test_split_is_complete_disjoint_and_label_independent():
    split = build_split()
    calibration = {row["patient_id"] for row in split["calibration"]}
    final = {row["patient_id"] for row in split["final_held_out"]}
    assert len(calibration) == len(final) == 6
    assert not calibration & final
    assert "Hu_287" not in calibration | final
    assert split["label_columns_read"] == []
    assert all(row["complete_raw_input_crosswalk"] for row in split["calibration"] + split["final_held_out"])


def test_split_membership_is_locked():
    split = build_split()
    assert [row["patient_id"] for row in split["calibration"]] == [
        "Hu_182", "Hu_315", "Hu_277", "Hu_268", "Hu_254", "Hu_343"
    ]
    assert [row["patient_id"] for row in split["final_held_out"]] == [
        "Hu_333", "Hu_159", "Hu_344", "Hu_048", "Hu_293", "Hu_250"
    ]
