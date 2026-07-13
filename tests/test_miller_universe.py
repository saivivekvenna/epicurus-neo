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


# ---- Fix 1: genuine pVAC only — a CSV alone must NOT establish a pVAC lane ---------------------------
def test_load_pvac_candidates_requires_genuine_provenance(monkeypatch, tmp_path):
    csv, prov = tmp_path / "pvac.csv", tmp_path / "prov.json"
    monkeypatch.setattr(u, "PVAC_CANDIDATES", csv)
    monkeypatch.setattr(u, "PVAC_PROVENANCE", prov)
    assert u.load_pvac_candidates().empty                         # both absent
    # a valid-looking CSV with NO provenance is rejected (no file alone establishes genuine pVAC)
    csv.write_text("mutation_id,mutant_peptide,hla_allele\n6:31000000:C:T,AAAAAAAAA,A*02:01\n")
    assert u.load_pvac_candidates().empty
    # provenance with wrong tool -> rejected
    prov.write_text(json.dumps({"tool": "homemade", "version": "1"}))
    assert u.load_pvac_candidates().empty
    # genuine provenance + schema + position-based id -> accepted, HLA normalized, source=pvac
    prov.write_text(json.dumps({"tool": "pvacseq", "version": "4.0.1"}))
    got = u.load_pvac_candidates()
    assert len(got) == 1 and got["candidate_source"].iloc[0] == "pvac"
    assert got["hla_allele"].iloc[0] == "HLA-A*02:01"
    # a non-position-based mutation_id is rejected even with provenance
    csv.write_text("mutation_id,mutant_peptide,hla_allele\nATRX-chrX-1,AAAAAAAAA,A*02:01\n")
    assert u.load_pvac_candidates().empty


def test_normalize_hla_form():
    assert u.normalize_hla("A*02:01") == "HLA-A*02:01" and u.normalize_hla("HLA-B*07:02") == "HLA-B*07:02"


# ---- Fix 3 (strengthened): freeze NEVER opens LABELS on a fully-mocked SUCCESS path ------------------
def test_freeze_never_reads_labels_on_success(monkeypatch, tmp_path):
    import builtins
    for name in ("PASS_VCF", "QUANT", "HLA_JSON", "REF"):
        p = tmp_path / name
        p.write_text("x")
        monkeypatch.setattr(u, name, p)
    (tmp_path / "HLA_JSON").write_text(json.dumps({"class_i_alleles": ["HLA-A*02:01"]}))
    monkeypatch.setattr(u, "FREEZE_DIR", tmp_path / "freeze")
    monkeypatch.setattr(u, "ART", tmp_path / "art")
    monkeypatch.setattr(u, "normalize_pass_vcf", lambda pv, ref, out: pv)
    monkeypatch.setattr(u, "load_filtered_variants",
                        lambda p: pd.DataFrame({"key": ["6:1:C:T"], "pass_filters": [True]}))
    uni = pd.DataFrame({"patient_id": ["Hu_287"], "mutation_id": ["6:1:C:T"], "candidate_source": ["lossless_recovery"],
                        "mutant_peptide": ["AAAAAAAAA"], "hla_allele": ["HLA-A*02:01"],
                        "genuine_prime": [0.9], "epicurus": [0.5]})
    monkeypatch.setattr(u, "build_universe", lambda v, h, t: (uni, [], False))
    monkeypatch.setattr(u, "score_universe", lambda x: x)
    monkeypatch.setattr(u, "gene_tpm_by_ensg", lambda q: {})
    monkeypatch.setattr(u, "arm_selection", lambda uni, arm: pd.DataFrame(
        {"mutation_id": ["6:1:C:T"], "mutant_peptide": ["AAAAAAAAA"], "hla_allele": ["HLA-A*02:01"]}))
    # tripwire: opening the sealed label path during freeze is a hard failure
    real_open = builtins.open
    monkeypatch.setattr(builtins, "open",
                        lambda f, *a, **k: (_ for _ in ()).throw(AssertionError("freeze opened LABELS"))
                        if str(f) == str(u.LABELS) else real_open(f, *a, **k))
    m = u.freeze()
    assert m["LOCK"] == "FROZEN_NO_LABELS" and m["labels_opened"] is False


# ---- Fix A: unseal reports null (not 0) for NOT_EVALUABLE arms ---------------------------------------
def test_unseal_reports_null_for_non_evaluable_arm(monkeypatch, tmp_path):
    # craft a frozen dir with one evaluable + one non-evaluable arm + a tiny labels file
    fd = tmp_path / "freeze"
    fd.mkdir()
    monkeypatch.setattr(u, "FREEZE_DIR", fd)
    monkeypatch.setattr(u, "ART", tmp_path / "art")
    (tmp_path / "art").mkdir()
    (fd / "variants.csv").write_text("key,pass_filters\n6:1:C:T,True\n")
    (fd / "select_lossless_prime.csv").write_text("mutation_id\n6:1:C:T\n")
    (fd / "select_pvac_prime.csv").write_text("mutation_id\n")
    data_files = ("variants.csv", "select_lossless_prime.csv", "select_pvac_prime.csv")
    man = {"sha256": {f: u.sha256_file(fd / f) for f in data_files}, "arms": {
        "lossless_prime": {"evaluable": True, "n_selected": 1, "saturated": False, "selection_file": "select_lossless_prime.csv"},
        "pvac_prime": {"evaluable": False, "missing": ["pvac_candidates"], "n_selected": 0, "saturated": False, "selection_file": "select_pvac_prime.csv"}}}
    (fd / "FREEZE_MANIFEST.json").write_text(json.dumps(man))
    labels = tmp_path / "labels.csv"
    labels.write_text("patient_id,gene_symbol,chrom,pos,ref,alt,label\nHu_287,G,6,1,C,T,POSITIVE\n")
    monkeypatch.setattr(u, "LABELS", labels)
    out = u.unseal()
    assert out["endpoint_b_class_i_four_arm"]["pvac_prime"]["hits_at_20_unique_mutations"] is None
    assert out["endpoint_b_class_i_four_arm"]["pvac_prime"]["n_selected"] is None
    assert out["endpoint_b_class_i_four_arm"]["lossless_prime"]["hits_at_20_unique_mutations"] == 1


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
