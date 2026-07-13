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
    monkeypatch.setattr(u, "rna_alt_evidence", lambda v, rna_bam=None: ({}, "COMPUTED"))
    monkeypatch.setattr(u, "build_universe", lambda v, h, t, r: (uni, [], False))
    monkeypatch.setattr(u, "score_universe", lambda x: x)
    monkeypatch.setattr(u, "gene_tpm_by_ensg", lambda q: {})
    # provenance inputs mocked to a single existing file so freeze can proceed on the success path
    sfile = tmp_path / "srcfile"
    sfile.write_text("code")
    monkeypatch.setattr(u, "_source_inputs", lambda: (sfile,))
    monkeypatch.setattr(u, "_all_inputs", lambda has_pvac: (sfile,))
    monkeypatch.setattr(u, "verify_tool_commits", lambda: (True, {"PRIME": {"match": True}}))
    monkeypatch.setattr(u, "verify_git_tracked_clean", lambda: (True, {"x": "CLEAN"}))
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
    monkeypatch.setattr(u, "verify_input_hashes", lambda man: (True, None, None))   # not under test here
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
def test_score_universe_el_nan_and_presentation_from_mixmhcpred(monkeypatch):
    # stub PRIME so the test needs no binary; el must be NaN, presentation must come from MixMHCpred rank
    class _Res:
        scored = pd.DataFrame({"peptide": ["AAAAAAAAA"], "hla_allele": ["HLA-A*02:01"],
                               "prime_rank": [0.5], "mixmhcpred_rank": [0.3]})
    import event_b.prime_adapter as pa
    import event_b.prime_transfer as pt
    monkeypatch.setattr(pa, "score_prime", lambda pairs, **k: _Res())
    monkeypatch.setattr(pt, "score_with_frozen", lambda frame, spec=None: np.array([0.0] * len(frame)))
    uni = pd.DataFrame({"mutant_peptide": ["AAAAAAAAA"], "hla_allele": ["HLA-A*02:01"], "expr": [5.0]})
    out = u.score_universe(uni)
    assert out["el"].isna().all()                       # el neutralized, NOT the MixMHCpred rank
    # Fix 5: router presentation evidence wired from the real MixMHCpred %rank (not NetMHCpan-EL, not absent)
    assert out["binding_percentile_rank"].iloc[0] == 0.3
    assert "MixMHCpred" in out["binding_rank_provenance"].iloc[0]


def test_class_i_length_filter_excludes_12_to_14mers():
    uni = pd.DataFrame({"mutant_peptide": ["A" * n for n in (8, 9, 11, 12, 13, 14)],
                        "mutation_id": ["m"] * 6, "candidate_source": ["lossless_recovery"] * 6})
    filt, counts = u.filter_class_i_lengths(uni)
    lens = sorted(filt["mutant_peptide"].str.len().unique())
    assert lens == [8, 9, 11] and counts["post"] == 3 and counts["dropped_len_12_14"] == 3
    assert (filt["mutant_peptide"].str.len() <= 11).all()   # 12-14mers never enter any arm


def test_build_universe_true_gene_type_and_evidence(monkeypatch):
    # stub VEP consequence + generator so no network is needed; assert schema fixes 6/7
    monkeypatch.setattr(u, "classify_consequence", lambda c, ch, p, r, a: ("inframe", "inframe_deletion"))

    class _Client:
        def __init__(self, *a, **k):
            pass

    def _gen(variant, client, panel, **k):
        cand = pd.DataFrame({"mutant_peptide": ["AAAAAAAAA"], "hla_allele": panel,
                             "source_variant_type": ["SNV"]})      # generator mislabels non-frameshift SNV
        return {"candidates": cand, "provenance": {"gene_symbol": "BRAF", "gene_id": "ENSG0.1"}, "windows": []}
    import event_b.lossless_peptide_generation as lg
    monkeypatch.setattr(lg, "EnsemblClient", _Client)
    monkeypatch.setattr(lg, "generate_variant_candidates", _gen)
    monkeypatch.setattr(u, "load_pvac_candidates", lambda: pd.DataFrame())
    variants = pd.DataFrame({"key": ["7:1:A:ATG"], "chrom": ["7"], "pos": [1], "ref": ["A"], "alt": ["ATG"],
                             "tumor_vaf": [0.4], "normal_vaf": [0.0], "tumor_dp": [50], "normal_dp": [40],
                             "pass_filters": [True]})
    rna = {"7:1:A:ATG": {"rna_alt_obs": 7, "rna_depth": 20, "rna_vaf": 0.35}}
    uni, notes, has_pvac = u.build_universe(variants, ["HLA-A*02:01"], {"ENSG0": 12.3}, rna)
    row = uni.iloc[0]
    assert row["gene_symbol"] == "BRAF"                   # Fix 6: real gene, not blank
    assert row["source_variant_type"] == "INFRAME"        # Fix 7: not mislabeled SNV
    assert row["tumor_vaf"] == 0.4 and row["normal_dp"] == 40   # Fix 6: WES evidence propagated
    assert row["rna_vaf"] == 0.35 and row["rna_alt_obs"] == 7   # RNA evidence attached
    assert row["expr"] == 12.3 and not has_pvac


def test_freeze_is_immutable_and_fails_closed(monkeypatch, tmp_path):
    fd = tmp_path / "freeze"
    fd.mkdir()
    monkeypatch.setattr(u, "FREEZE_DIR", fd)
    monkeypatch.setattr(u, "verify_input_hashes", lambda man: (True, None, None))   # input gate not under test here
    # valid LOCK with matching derived hash -> ALREADY_FROZEN (never overwrite)
    (fd / "variants.csv").write_text("k\n1\n")
    man = {"LOCK": "FROZEN_NO_LABELS", "sha256": {"variants.csv": u.sha256_file(fd / "variants.csv")}, "arms": {}}
    (fd / "FREEZE_MANIFEST.json").write_text(json.dumps(man))
    assert u.freeze()["status"] == "ALREADY_FROZEN"
    # corrupt lock -> fail closed, never overwrite
    (fd / "FREEZE_MANIFEST.json").write_text("{not json")
    assert u.freeze()["status"] == "FROZEN_CORRUPT"
    # tampered derived file under a valid lock -> FROZEN_HASH_MISMATCH (not ALREADY_FROZEN, not overwrite)
    (fd / "FREEZE_MANIFEST.json").write_text(json.dumps(man))
    (fd / "variants.csv").write_text("k\n1\n2\n")
    assert u.freeze()["status"] == "FROZEN_HASH_MISMATCH"


def test_freeze_not_evaluable_when_rna_missing_no_manifest(monkeypatch, tmp_path):
    fd = tmp_path / "freeze"
    monkeypatch.setattr(u, "FREEZE_DIR", fd)
    present = tmp_path / "present"
    present.write_text("x")
    monkeypatch.setattr(u, "_source_inputs", lambda: (present, tmp_path / "rna.bam"))   # RNA BAM absent
    out = u.freeze()
    assert out["status"] == "NOT_EVALUABLE" and any("rna.bam" in m for m in out["missing_inputs"])
    assert not (fd / "FREEZE_MANIFEST.json").exists()      # NEVER writes a lock when an input is missing


def test_freeze_not_evaluable_when_pvac_inputs_missing_no_manifest(monkeypatch, tmp_path):
    fd = tmp_path / "freeze"
    monkeypatch.setattr(u, "FREEZE_DIR", fd)
    present = tmp_path / "present"
    present.write_text("x")
    monkeypatch.setattr(u, "_source_inputs", lambda: (present,))
    monkeypatch.setattr(u, "verify_tool_commits", lambda: (True, {}))
    monkeypatch.setattr(u, "verify_git_tracked_clean", lambda: (True, {}))
    # genuine_pvac_lane True but the pVAC CSV/provenance are absent -> _all_inputs includes them -> refuse
    monkeypatch.setattr(u, "_all_inputs", lambda has_pvac: (present, tmp_path / "pvac.csv") if has_pvac else (present,))
    monkeypatch.setattr(u, "normalize_pass_vcf", lambda pv, ref, out: present)
    monkeypatch.setattr(u, "load_filtered_variants", lambda p: pd.DataFrame({"key": ["6:1:C:T"], "pass_filters": [True]}))
    monkeypatch.setattr(u, "rna_alt_evidence", lambda v, rna_bam=None: ({}, "COMPUTED"))
    uni = pd.DataFrame({"patient_id": ["Hu_287"], "mutation_id": ["6:1:C:T"], "candidate_source": ["pvac"],
                        "mutant_peptide": ["AAAAAAAAA"], "hla_allele": ["HLA-A*02:01"], "genuine_prime": [0.9], "epicurus": [0.5]})
    monkeypatch.setattr(u, "build_universe", lambda v, h, t, r: (uni, [], True))    # has_pvac=True
    monkeypatch.setattr(u, "score_universe", lambda x: x)
    monkeypatch.setattr(u, "gene_tpm_by_ensg", lambda q: {})
    monkeypatch.setattr(u, "HLA_JSON", tmp_path / "hla.json")
    (tmp_path / "hla.json").write_text(json.dumps({"class_i_alleles": ["HLA-A*02:01"]}))
    out = u.freeze()
    assert out["status"] == "NOT_EVALUABLE" and "pvac.csv" in out.get("missing_input", "")
    assert not (fd / "FREEZE_MANIFEST.json").exists()


def test_code_files_cover_all_six_semantics_files():
    names = {p.name for p in u.CODE_FILES}
    assert names == {"miller_hu287_universe.py", "four_arm.py", "lossless_peptide_generation.py",
                     "prime_adapter.py", "prime_transfer.py", "evidence_router.py"}
    assert len(u.CODE_FILES) == 6


def test_verify_tool_commits_matches_adapter_constants():
    # the on-disk PRIME/MixMHCpred repos must match the adapter's pinned commits AND be tracked-clean
    ok, info = u.verify_tool_commits()
    assert ok is True
    assert info["PRIME"]["match"] and info["MixMHCpred"]["match"]
    assert info["PRIME"]["tracked_clean"] and info["MixMHCpred"]["tracked_clean"]
    assert info["PRIME"]["dir_head"] == info["PRIME"]["adapter_constant"]


def test_verify_tool_commits_fails_closed_on_dirty_tracked(monkeypatch):
    from event_b.prime_adapter import MIX_COMMIT, PRIME_COMMIT
    # HEAD matches, but tracked files are dirty (git diff --quiet exits non-zero) -> fail closed
    monkeypatch.setattr(u.subprocess, "check_output",
                        lambda cmd, **k: (PRIME_COMMIT if "MixMHCpred" not in cmd[2] else MIX_COMMIT) + "\n")

    def _dirty(cmd, **k):
        raise u.subprocess.CalledProcessError(1, cmd)
    monkeypatch.setattr(u.subprocess, "check_call", _dirty)
    ok, info = u.verify_tool_commits()
    assert ok is False
    assert info["PRIME"]["match"] and info["PRIME"]["tracked_clean"] is False


def test_freeze_not_evaluable_on_tool_commit_mismatch_no_manifest(monkeypatch, tmp_path):
    fd = tmp_path / "freeze"
    monkeypatch.setattr(u, "FREEZE_DIR", fd)
    present = tmp_path / "present"
    present.write_text("x")
    monkeypatch.setattr(u, "_source_inputs", lambda: (present,))
    monkeypatch.setattr(u, "verify_tool_commits",
                        lambda: (False, {"PRIME": {"dir_head": "abc", "adapter_constant": "def", "match": False}}))
    out = u.freeze()
    assert out["status"] == "NOT_EVALUABLE" and "tool_commit_issue" in out
    assert not (fd / "FREEZE_MANIFEST.json").exists()


def test_verify_git_tracked_clean_states(tmp_path, monkeypatch):
    import subprocess as sp
    repo = tmp_path / "repo"
    repo.mkdir()
    env = {"GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@t", "GIT_COMMITTER_NAME": "t",
           "GIT_COMMITTER_EMAIL": "t@t", "PATH": __import__("os").environ["PATH"]}
    sp.check_call(["git", "init", "-q", str(repo)])
    clean, staged, unstaged = repo / "clean.py", repo / "staged.py", repo / "unstaged.py"
    for f, t in ((clean, "a=1\n"), (staged, "b=1\n"), (unstaged, "c=1\n")):
        f.write_text(t)
        sp.check_call(["git", "-C", str(repo), "add", f.name])
    sp.check_call(["git", "-C", str(repo), "commit", "-qm", "init"], env=env)
    untracked = repo / "untracked.py"
    untracked.write_text("d=1\n")
    staged.write_text("b=2\n")
    sp.check_call(["git", "-C", str(repo), "add", "staged.py"])       # staged modification
    unstaged.write_text("c=2\n")                                       # unstaged modification
    monkeypatch.setattr(u, "ROOT", repo)
    monkeypatch.setattr(u, "_repo_semantic_files", lambda: (clean, staged, unstaged, untracked))
    ok, info = u.verify_git_tracked_clean()
    assert not ok
    assert info["clean.py"] == "CLEAN"
    assert info["staged.py"] == "STAGED_MODIFIED"
    assert info["unstaged.py"] == "UNSTAGED_MODIFIED"
    assert info["untracked.py"] == "UNTRACKED"


def test_freeze_not_evaluable_on_untracked_semantic_file_no_manifest(monkeypatch, tmp_path):
    fd = tmp_path / "freeze"
    monkeypatch.setattr(u, "FREEZE_DIR", fd)
    present = tmp_path / "present"
    present.write_text("x")
    monkeypatch.setattr(u, "_source_inputs", lambda: (present,))
    monkeypatch.setattr(u, "verify_tool_commits", lambda: (True, {}))
    monkeypatch.setattr(u, "verify_git_tracked_clean",
                        lambda: (False, {"src/event_b/prime_adapter.py": "UNTRACKED"}))
    out = u.freeze()
    assert out["status"] == "NOT_EVALUABLE" and "git_tracking_issue" in out
    assert not (fd / "FREEZE_MANIFEST.json").exists()


def test_unseal_is_once_only(monkeypatch, tmp_path):
    fd = tmp_path / "freeze"
    fd.mkdir()
    monkeypatch.setattr(u, "FREEZE_DIR", fd)
    (fd / "UNSEALED.json").write_text(json.dumps({"unsealed": True, "manifest_file_sha256": "abc"}))
    # tripwire: if unseal reads LABELS after already-unsealed, fail
    monkeypatch.setattr(u.pd, "read_csv",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError("re-unseal read LABELS")))
    out = u.unseal()
    assert out["status"] == "ALREADY_UNSEALED"


def test_verify_frozen_hashes_fails_closed_on_missing_file(tmp_path):
    man = {"sha256": {"gone.csv": "deadbeef"}}                    # file does not exist
    ok, bad = u.verify_frozen_hashes(man, tmp_path)
    assert not ok and bad == "gone.csv"
    assert u.verify_frozen_hashes({"sha256": "notadict"}, tmp_path)[0] is False   # malformed manifest


def test_unseal_missing_frozen_file_is_hash_mismatch_no_label_read(monkeypatch, tmp_path):
    fd = tmp_path / "freeze"
    fd.mkdir()
    monkeypatch.setattr(u, "FREEZE_DIR", fd)
    # manifest references a derived file that is absent -> HASH_MISMATCH, and LABELS is never read
    man = {"sha256": {"variants.csv": "deadbeef"}, "arms": {}, "input_sha256": {}}
    (fd / "FREEZE_MANIFEST.json").write_text(json.dumps(man))
    monkeypatch.setattr(u.pd, "read_csv",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError("unseal read LABELS on bad hash")))
    out = u.unseal()
    assert out["status"] == "HASH_MISMATCH" and out["file"] == "variants.csv"


def test_sha256_file_streams_match(tmp_path):
    import hashlib as _h
    blob = b"hello-world" * 300000                        # ~3.3MB, forces multiple chunks
    p = tmp_path / "big.bin"
    p.write_bytes(blob)
    assert u.sha256_file(p, chunk=4096) == _h.sha256(blob).hexdigest()


def test_verify_input_hashes_requires_complete_valid_set(monkeypatch, tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("A")
    b.write_text("B")
    monkeypatch.setattr(u, "_all_inputs", lambda has_pvac: (a, b))
    complete = {u._rel(a): u.sha256_file(a), u._rel(b): u.sha256_file(b)}
    assert u.verify_input_hashes({"input_sha256": complete})[0] is True
    assert u.verify_input_hashes({"input_sha256": {u._rel(a): u.sha256_file(a)}})[0] is False   # key missing
    assert u.verify_input_hashes({"input_sha256": {u._rel(a): "MISSING", u._rel(b): u.sha256_file(b)}})[0] is False
    assert u.verify_input_hashes({"input_sha256": {}})[0] is False                               # empty fails closed
    assert u.verify_input_hashes({})[0] is False                                                 # absent fails closed


def test_unseal_fails_closed_on_incomplete_input_hashes_no_label_read(monkeypatch, tmp_path):
    fd = tmp_path / "freeze"
    fd.mkdir()
    monkeypatch.setattr(u, "FREEZE_DIR", fd)
    (fd / "variants.csv").write_text("k\n1\n")
    man = {"sha256": {"variants.csv": u.sha256_file(fd / "variants.csv")}, "arms": {}, "input_sha256": {}}
    (fd / "FREEZE_MANIFEST.json").write_text(json.dumps(man))
    monkeypatch.setattr(u.pd, "read_csv",
                        lambda *a, **k: (_ for _ in ()).throw(AssertionError("unseal read LABELS on incomplete inputs")))
    out = u.unseal()
    assert out["status"] == "INPUT_HASH_INCOMPLETE_OR_MISMATCH"    # derived hashes OK, input gate fails closed


def test_rna_alt_snv_counts_indel_not_assessed(monkeypatch, tmp_path):
    # no BAM -> NOT_ASSESSED status (never fabricated)
    monkeypatch.setattr(u, "RNA_BAM", tmp_path / "absent.bam")
    got, status = u.rna_alt_evidence(pd.DataFrame({"key": [], "chrom": [], "pos": [], "ref": [], "alt": [],
                                                   "pass_filters": []}))
    assert got == {} and "NOT_ASSESSED" in status


def test_freeze_not_evaluable_on_missing_inputs(tmp_path, monkeypatch):
    # freeze returns NOT_EVALUABLE cleanly when a source input is missing (no crash, no label access)
    monkeypatch.setattr(u, "FREEZE_DIR", tmp_path / "freeze")
    monkeypatch.setattr(u, "PASS_VCF", tmp_path / "absent.vcf.gz")
    m = u.freeze()
    assert m["status"] == "NOT_EVALUABLE" and any("absent.vcf.gz" in x for x in m["missing_inputs"])
