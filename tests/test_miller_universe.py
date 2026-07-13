"""Pure-helper tests for the Miller Hu_287 universe/scoring script (no VCF/network/labels)."""

from __future__ import annotations

import importlib

u = importlib.import_module("scripts.miller_hu287_universe")


def test_norm_chrom_and_variant_key_are_join_safe():
    # label uses 'chrX'; Ensembl VCF uses 'X' -> the position key must match either way
    assert u.norm_chrom("chrX") == "X" and u.norm_chrom("6") == "6" and u.norm_chrom("chr17") == "17"
    assert u.variant_key("chrX", 77618864, "a", "g") == u.variant_key("X", 77618864, "A", "G")
    assert u.variant_key("6", "31000000", "C", "T") == "6:31000000:C:T"


def test_base_filters_frozen_thresholds():
    # tumor VAF>=0.05, normal VAF<=0.05, depth>=10 both, T/N ratio>=1
    assert u.passes_base_filters(0.30, 0.00, 40, 38) is True
    assert u.passes_base_filters(0.04, 0.00, 40, 38) is False        # tumor VAF too low
    assert u.passes_base_filters(0.30, 0.10, 40, 38) is False        # germline leakage (normal VAF high)
    assert u.passes_base_filters(0.30, 0.00, 8, 38) is False         # tumor depth too low
    assert u.passes_base_filters(0.30, 0.00, 40, 6) is False         # normal depth too low
    assert u.passes_base_filters(0.05, 0.05, 10, 10) is True         # exact boundary (ratio=1, depth=10)


def test_vaf_depth_prefers_AD_then_falls_back():
    assert u._vaf_depth({"AD": (30, 10)}) == (0.25, 40)
    v, d = u._vaf_depth({"DP": 50, "AF": (0.2,)})
    assert d == 50 and abs(v - 0.2) < 1e-9
    assert u._vaf_depth({"AD": (0, 0), "DP": 0}) == (0.0, 0)          # degenerate -> zero (filtered out)
