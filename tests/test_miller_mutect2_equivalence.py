from pathlib import Path

import pytest

from scripts import miller_mutect2_equivalence as equivalence


LEFT = "1\t100\t.\tA\tT\t10\tPASS\tDP=20\tGT:AD\t0/0:20,0\t0/1:10,10"
RIGHT = "1\t200\t.\tG\tC\t20\tPASS\tDP=30\tGT:AD\t0/0:30,0\t0/1:15,15"


def test_exact_records_are_deployable_basis():
    result = equivalence.compare_records([LEFT, RIGHT], [LEFT, RIGHT])
    assert result["records_exact"] is True
    assert result["site_keys_equal_in_order"] is True
    assert result["missing_from_scattered"] == [] and result["extra_in_scattered"] == []
    assert result["record_stream_sha256"]["serial"] == result["record_stream_sha256"]["scattered"]


def test_annotation_or_genotype_difference_fails_even_at_same_site():
    changed = LEFT.replace("DP=20", "DP=19")
    result = equivalence.compare_records([LEFT], [changed])
    assert result["site_keys_equal_in_order"] is True
    assert result["records_exact"] is False
    assert result["first_record_mismatches"][0]["serial"] == LEFT


def test_missing_and_extra_sites_are_reported():
    result = equivalence.compare_records([LEFT], [RIGHT])
    assert result["missing_from_scattered"] == ["1:100:A:T"]
    assert result["extra_in_scattered"] == ["1:200:G:C"]


def test_malformed_record_fails_closed():
    with pytest.raises(ValueError, match="malformed"):
        equivalence.variant_key("1\t2\tA")
