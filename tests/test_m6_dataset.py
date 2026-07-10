from epicurus_neo.m6.dataset import load_label_frame, parse_alleles


def test_label_frame_matches_registered_population():
    frame = load_label_frame()
    assert len(frame) == 965
    assert int(frame.label.sum()) == 272
    assert int((frame.label == 0).sum()) == 693
    assert frame.patient_id.nunique() == 45
    assert sorted(frame.study_id.unique()) == [
        "braun_rcc_2025",
        "hu_neovax_2021",
        "mkras_vax_2026",
        "pdac_neovax_2023",
    ]
    # UNTESTED dropped; labels are strictly binary.
    assert set(frame.label.unique()) == {0, 1}
    # Every row carries a scalar tie-break allele (possibly empty), never a list.
    assert frame.hla_allele.map(lambda v: isinstance(v, str)).all()
    # gene/protein_change ride along for the pdac presentation join.
    assert {"gene", "protein_change"}.issubset(frame.columns)


def test_parse_alleles_handles_json_list_and_scalar():
    assert parse_alleles('["HLA-A*02:01", "HLA-B*07:02"]') == ["HLA-A*02:01", "HLA-B*07:02"]
    assert parse_alleles("HLA-A*02:01") == ["HLA-A*02:01"]
    assert parse_alleles(None) == []
    assert parse_alleles("[]") == []
