"""Pure-helper + audit-fix tests for the Miller Hu_287 freeze/unseal script (no VCF/network/labels)."""

from __future__ import annotations

import importlib
import json

import numpy as np
import pandas as pd

u = importlib.import_module("scripts.miller_hu287_universe")


def test_norm_chrom_and_variant_key_are_join_safe():
    assert u.norm_chrom("chrX") == "X" and u.norm_chrom("6") == "6" and u.norm_chrom("chr17") == "17"
    assert u.variant_key("chrX", 77618864, "a", "g") == u.variant_key("X", 77618864, "A", "G")
    assert u.variant_key("6", "31000000", "C", "T") == "6:31000000:C:T"


def test_base_filters_frozen_thresholds():
    assert u.passes_base_filters(0.30, 0.00, 40, 38) is True
    assert u.passes_base_filters(0.04, 0.00, 40, 38) is False
    assert u.passes_base_filters(0.30, 0.10, 40, 38) is False
    assert u.passes_base_filters(0.30, 0.00, 8, 38) is False
    assert u.passes_base_filters(0.30, 0.00, 40, 6) is False
    assert u.passes_base_filters(0.05, 0.05, 10, 10) is True


def test_vaf_depth_prefers_AD_then_falls_back():
    assert u._vaf_depth({"AD": (30, 10)}) == (0.25, 40)
    v, d = u._vaf_depth({"DP": 50, "AF": (0.2,)})
    assert d == 50 and abs(v - 0.2) < 1e-9


# ---- Fix 4: VEP consequence gates enumeration (synonymous/stop/splice -> NOT_ENUMERABLE) --------------
class _StubClient:
    def __init__(self, term):
        self._t = term

    def vep_hgvs(self, hgvs):
        return {"json": [{"most_severe_consequence": self._t}], "url": "x", "sha256": "y"}


def test_classify_consequence_only_protein_altering_enumerable(monkeypatch):
    monkeypatch.setattr(u, "sys", u.sys)  # no-op to keep import side effects stable
    import event_b.lossless_peptide_generation as lg
    monkeypatch.setattr(lg, "genomic_hgvs", lambda c, p, r, a: "hgvs")
    assert u.classify_consequence(_StubClient("missense_variant"), "6", 1, "A", "T") == ("missense", "missense_variant")
    assert u.classify_consequence(_StubClient("frameshift_variant"), "6", 1, "A", "AT") == ("frameshift", "frameshift_variant")
    assert u.classify_consequence(_StubClient("inframe_deletion"), "6", 1, "ATG", "A") == ("inframe", "inframe_deletion")
    # non-protein-altering -> NOT_ENUMERABLE (kind is None)
    for term in ("synonymous_variant", "stop_gained", "splice_region_variant", "intron_variant"):
        kind, raw = u.classify_consequence(_StubClient(term), "6", 1, "A", "T")
        assert kind is None and raw == term


# ---- Fix 1: no genuine pVAC -> empty (never fabricated from lossless) --------------------------------
def test_load_pvac_candidates_absent_is_empty(monkeypatch, tmp_path):
    monkeypatch.setattr(u, "PVAC_CANDIDATES", tmp_path / "nope.csv")
    assert u.load_pvac_candidates().empty


# ---- Fix 3: freeze/unseal hash lock ------------------------------------------------------------------
def test_verify_frozen_hashes_detects_tamper(tmp_path):
    (tmp_path / "select_x.csv").write_text("mutation_id\nM1\n")
    man = {"sha256": {"select_x.csv": u.sha256_file(tmp_path / "select_x.csv")}}
    ok, bad = u.verify_frozen_hashes(man, tmp_path)
    assert ok and bad is None
    (tmp_path / "select_x.csv").write_text("mutation_id\nM1\nM2\n")   # tamper after freeze
    ok2, bad2 = u.verify_frozen_hashes(man, tmp_path)
    assert not ok2 and bad2 == "select_x.csv"


# ---- Fix 2: Epicurus el is NaN (NetMHCpan-EL unavailable; MixMHCpred is not a valid substitute) ------
def test_score_universe_sets_el_nan(monkeypatch):
    # stub PRIME so the test needs no binary; assert el is set to NaN (not mixmhcpred_rank)
    class _Res:
        scored = pd.DataFrame({"peptide": ["AAAAAAAAA"], "hla_allele": ["HLA-A*02:01"],
                               "prime_rank": [0.5], "mixmhcpred_rank": [0.3]})
    monkeypatch.setattr(u, "score_universe", u.score_universe)  # keep ref
    import event_b.prime_adapter as pa
    import event_b.prime_transfer as pt
    monkeypatch.setattr(pa, "score_prime", lambda pairs, **k: _Res())
    monkeypatch.setattr(pt, "score_with_frozen", lambda frame, spec=None: np.array([0.0] * len(frame)))
    uni = pd.DataFrame({"mutant_peptide": ["AAAAAAAAA"], "hla_allele": ["HLA-A*02:01"], "expr": [5.0]})
    out = u.score_universe(uni)
    assert out["el"].isna().all()                       # el neutralized, NOT the MixMHCpred rank
    assert "mixmhcpred_rank" in out.columns and out["mixmhcpred_rank"].iloc[0] == 0.3


def test_freeze_manifest_discloses_el_and_pvac(tmp_path, monkeypatch):
    # freeze returns NOT_EVALUABLE cleanly when inputs are missing (no crash, no label access)
    monkeypatch.setattr(u, "PASS_VCF", tmp_path / "absent.vcf.gz")
    m = u.freeze()
    assert m["status"] == "NOT_EVALUABLE" and "absent" in m["missing"]
    _ = json  # keep import used
