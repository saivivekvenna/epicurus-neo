from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.presentation import (
    presentation_availability,
    resolve_class_i_alleles,
)


def test_hla_resolution_covers_hu_and_pdac_only():
    frame = resolve_class_i_alleles(load_label_frame())
    availability = presentation_availability(frame).set_index("study_id")
    # hu carries candidate-level HLA; pdac gets predicted best-binder alleles by join.
    assert availability.loc["hu_neovax_2021", "resolved"] > 0
    assert availability.loc["pdac_neovax_2023", "resolved"] > 0
    # braun HLA is not public; mKRAS long peptides are NOT_ASSESSED -> zero resolved.
    assert availability.loc["braun_rcc_2025", "resolved"] == 0
    assert availability.loc["mkras_vax_2026", "resolved"] == 0
