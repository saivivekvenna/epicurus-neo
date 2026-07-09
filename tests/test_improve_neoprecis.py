from epicurus_neo.neoprecis_adapter import (
    normalize_neoprecis_allele,
    wildtype_pseudo_core,
)


def test_normalize_allele_converts_improve_notation():
    assert normalize_neoprecis_allele("HLA-A02:01") == "A*02:01"
    assert normalize_neoprecis_allele("B*07:02") == "B*07:02"


def test_wildtype_pseudo_core_preserves_corresponding_residues():
    assert wildtype_pseudo_core("CIDFQPEIY", "CIDFQPDIY", "CIDFQPDIY") == "CIDFQPEIY"
    assert wildtype_pseudo_core("SCIDFQPEIY", "SCIDFQPDIY", "SIDFQPDIY") == "SIDFQPEIY"
