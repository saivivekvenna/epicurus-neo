"""Generic label-blind FREEZE core (Miller arbitrary patient). Freeze-ONLY: no unseal, no label path.

These tests never touch reconstructed data, the network, PRIME binaries, or the recognition-label table.
"""

from __future__ import annotations

import dataclasses
import json
import sys

import pandas as pd

import scripts.miller_hu287_universe as u
from benchmark import miller_universe_core as core
from benchmark.miller_patient import load_patient


def _tmp_config(tmp_path, patient_id="Hu_315", **overrides):
    """A real resolved config with all write targets redirected under tmp_path (no repo side effects)."""
    cfg = core.UniverseConfig.for_patient(load_patient(patient_id))
    base = {"freeze_dir": tmp_path / "freeze", "artifact_dir": tmp_path / "art"}
    base.update(overrides)
    return dataclasses.replace(cfg, **base)


# ---------------------------------------------------------------------------
# Config resolution: Hu_287 equivalence + Hu_315 dynamic
# ---------------------------------------------------------------------------
def test_config_for_hu287_matches_frozen_script_constants():
    """The generic resolver reproduces the frozen Hu_287 script's constants EXACTLY (faithful param, not a
    re-derivation), so the two lanes would consume byte-identical inputs for Hu_287."""
    c = core.UniverseConfig.for_patient(load_patient("Hu_287"))
    assert c.pass_vcf == u.PASS_VCF and c.norm_vcf == u.NORM_VCF and c.quant == u.QUANT
    assert c.hla_json == u.HLA_JSON and c.rna_bam == u.RNA_BAM and c.ens_cache == u.ENS_CACHE
    assert c.pvac_candidates == u.PVAC_CANDIDATES and c.pvac_provenance == u.PVAC_PROVENANCE
    assert c.artifact_dir == u.ART and c.freeze_dir == u.FREEZE_DIR and c.ref == u.REF
    assert (c.sample_normal, c.sample_tumor) == ("Hu_287_N", "Hu_287_T")


def test_config_for_hu315_is_dynamic_and_milestone8():
    c = core.UniverseConfig.for_patient(load_patient("Hu_315"))
    assert c.raw_dir.name == "hu_315"
    assert c.pass_vcf.name == "Hu_315.somatic.pass.vcf.gz"
    assert c.rna_bam.name == "Hu_315_tumor_rna.sorted.bam"
    assert (c.sample_normal, c.sample_tumor) == ("Hu_315_N", "Hu_315_T")
    assert c.artifact_dir.as_posix().endswith("artifacts/milestone_8_generalization/patients/Hu_315")
    assert c.freeze_dir.as_posix().endswith("data/raw/miller_ipv/hu_315/freeze")


def test_code_files_pin_the_frozen_script_and_patient_resolver():
    names = {p.name for p in core.CODE_FILES}
    # the generic core executes miller_hu287_universe helpers -> its bytes shape output -> must be pinned
    assert "miller_hu287_universe.py" in names
    # miller_patient resolves every path/sample id -> must be pinned
    assert "miller_patient.py" in names
    assert {"miller_universe_core.py", "miller_patient_universe.py", "four_arm.py",
            "lossless_peptide_generation.py", "prime_adapter.py", "prime_transfer.py",
            "evidence_router.py"} <= names
    # semantic_files (via for_patient) = code files + the three frozen configs
    c = core.UniverseConfig.for_patient(load_patient("Hu_315"))
    sem = {p.name for p in c.semantic_files}
    assert names <= sem and "epicurus_v0_1.json" in sem


# ---------------------------------------------------------------------------
# Freeze-only surface: NO unseal, NO label path anywhere in this module
# ---------------------------------------------------------------------------
def test_module_is_freeze_only_and_has_no_label_surface():
    assert not hasattr(core, "unseal")
    assert not hasattr(core, "LABELS")
    src = (u.ROOT / "src/benchmark/miller_universe_core.py").read_text()
    assert "recognition_labels" not in src and "miller_recognition_labels" not in src


# ---------------------------------------------------------------------------
# Preservation of the Hu_287 frozen provenance
# ---------------------------------------------------------------------------
def test_freeze_refuses_hu287_and_writes_nothing(tmp_path):
    c = _tmp_config(tmp_path, patient_id="Hu_287")
    out = core.freeze(c)
    assert out["status"] == "REFUSED_HU287"
    assert not (c.freeze_dir / "FREEZE_MANIFEST.json").exists()
    assert not c.freeze_dir.exists()


# ---------------------------------------------------------------------------
# Fail-closed: missing reconstructed inputs -> NOT_EVALUABLE, no manifest
# ---------------------------------------------------------------------------
def test_freeze_not_evaluable_when_inputs_missing_no_manifest(tmp_path):
    # real Hu_315 config: none of its reconstructed inputs exist on disk yet
    c = _tmp_config(tmp_path)
    out = core.freeze(c)
    assert out["status"] == "NOT_EVALUABLE"
    assert any("hu_315" in m for m in out["missing_inputs"])
    assert not (c.freeze_dir / "FREEZE_MANIFEST.json").exists()


# ---------------------------------------------------------------------------
# The ONLY generalization of the variant loader: patient-specific sample columns
# ---------------------------------------------------------------------------
class _Rec:
    def __init__(self, chrom, pos, ref, alt, t_smp, n_ad, sample_t, sample_n):
        self.chrom, self.pos, self.ref, self.alts = chrom, pos, ref, (alt,)
        self.samples = {sample_t: t_smp, sample_n: {"AD": n_ad}}


def test_load_filtered_variants_reads_patient_sample_columns(monkeypatch, tmp_path):
    c = _tmp_config(tmp_path)  # Hu_315 -> expects Hu_315_T / Hu_315_N
    good = _Rec("6", 100, "C", "T", {"AD": (200, 100)}, (300, 0), "Hu_315_T", "Hu_315_N")
    wrong = _Rec("6", 200, "A", "G", {"AD": (200, 100)}, (300, 0), "Hu_287_T", "Hu_287_N")  # wrong sample -> skip

    class _VF:
        def __init__(self, path):
            pass

        def __iter__(self):
            return iter([good, wrong])

    monkeypatch.setitem(sys.modules, "pysam", type("pysam", (), {"VariantFile": _VF}))
    df = core.load_filtered_variants(c, tmp_path / "unused.vcf")
    # only the row with THIS patient's sample columns survives; the Hu_287-named record is silently skipped
    assert list(df["pos"]) == [100]
    assert df["tumor_alt_reads"].iloc[0] == 100 and bool(df["pass_filters"].iloc[0]) is True
    assert "strict5_pass" in df.columns


def test_build_universe_stamps_real_patient_id_and_uses_patient_cache(monkeypatch, tmp_path):
    c = _tmp_config(tmp_path)  # Hu_315
    monkeypatch.setattr(core.u, "classify_consequence", lambda cl, ch, p, r, a: ("missense", "missense_variant"))

    seen = {}

    class _Client:
        def __init__(self, cache_dir, *a, **k):
            seen["cache"] = cache_dir
            self.accessed = {}
            self.cache_dir = cache_dir

    def _gen(variant, client, panel, **k):
        cand = pd.DataFrame({"mutant_peptide": ["AAAAAAAAA"], "hla_allele": panel,
                             "source_variant_type": ["SNV"]})
        return {"candidates": cand, "provenance": {"gene_symbol": "KRAS", "gene_id": "ENSG42.3"}, "windows": []}

    import event_b.lossless_peptide_generation as lg
    monkeypatch.setattr(lg, "EnsemblClient", _Client)
    monkeypatch.setattr(lg, "generate_variant_candidates", _gen)
    variants = pd.DataFrame({"key": ["12:1:C:T"], "chrom": ["12"], "pos": [1], "ref": ["C"], "alt": ["T"],
                             "tumor_vaf": [0.4], "normal_vaf": [0.0], "tumor_dp": [50], "normal_dp": [40],
                             "tumor_alt_reads": [20], "pass_filters": [True]})
    uni, notes, has_pvac, used = core.build_universe(c, variants, ["HLA-A*02:01"], {"ENSG42": 9.9},
                                                     {"12:1:C:T": {"rna_alt_obs": 5, "rna_depth": 15, "rna_vaf": 0.33}})
    assert uni["patient_id"].iloc[0] == "Hu_315"        # real patient id, NOT hard-coded Hu_287
    assert uni["gene_symbol"].iloc[0] == "KRAS"
    assert uni["source_variant_type"].iloc[0] == "MISSENSE"
    assert uni["expr"].iloc[0] == 9.9 and uni["rna_vaf"].iloc[0] == 0.33
    assert seen["cache"] == c.ens_cache and not has_pvac  # patient-specific Ensembl cache


# ---------------------------------------------------------------------------
# Parameterized provenance gates
# ---------------------------------------------------------------------------
def test_verify_input_hashes_requires_complete_set(monkeypatch, tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    a.write_text("A")
    b.write_text("B")
    c = _tmp_config(tmp_path)
    monkeypatch.setattr(core, "_all_inputs", lambda config, has_pvac: (a, b))
    complete = {u._rel(a): u.sha256_file(a), u._rel(b): u.sha256_file(b)}
    assert core.verify_input_hashes(c, {"input_sha256": complete})[0] is True
    assert core.verify_input_hashes(c, {"input_sha256": {u._rel(a): u.sha256_file(a)}})[0] is False
    assert core.verify_input_hashes(c, {"input_sha256": {}})[0] is False


def test_verify_ensembl_used_resolves_against_patient_cache(tmp_path):
    c = _tmp_config(tmp_path, ens_cache=tmp_path / "ens")
    c.ens_cache.mkdir()
    f = c.ens_cache / "resp.json"
    f.write_text('{"x":1}')
    rec = {"url": "https://rest.ensembl.org/vep/x", "cache_path": u._rel(f), "sha256": u.sha256_file(f)}
    assert core.verify_ensembl_used(c, [rec], require_nonempty=True)[0] is True
    assert core.verify_ensembl_used(c, None, require_nonempty=False)[0] is False
    # a cache_path OUTSIDE this patient's cache is rejected (traversal guard)
    outside = tmp_path / "secret.json"
    outside.write_text("s")
    bad = {"url": "u2", "cache_path": u._rel(outside), "sha256": u.sha256_file(outside)}
    ok, reason = core.verify_ensembl_used(c, [bad], require_nonempty=False)
    assert ok is False and "path_escapes_cache" in reason


# ---------------------------------------------------------------------------
# One-shot immutability + fully-mocked success path (never opens a label)
# ---------------------------------------------------------------------------
def test_freeze_is_immutable_already_frozen(monkeypatch, tmp_path):
    c = _tmp_config(tmp_path)
    fd = c.freeze_dir
    fd.mkdir(parents=True)
    (fd / "variants.csv").write_text("k\n1\n")
    man = {"LOCK": "FROZEN_NO_LABELS", "sha256": {"variants.csv": u.sha256_file(fd / "variants.csv")},
           "arms": {}, "n_variants_pass": 0, "ensembl_used_responses": []}
    (fd / "FREEZE_MANIFEST.json").write_text(json.dumps(man))
    monkeypatch.setattr(core, "verify_input_hashes", lambda config, m: (True, None, None))
    assert core.freeze(c)["status"] == "ALREADY_FROZEN"
    # tampered derived file under a valid lock -> hash mismatch, never overwrites
    (fd / "variants.csv").write_text("k\n1\n2\n")
    assert core.freeze(c)["status"] == "FROZEN_HASH_MISMATCH"


def test_freeze_success_writes_lock_with_real_patient_id_no_label_read(monkeypatch, tmp_path):
    c = _tmp_config(tmp_path, hla_json=tmp_path / "hla.json")
    (tmp_path / "hla.json").write_text(json.dumps({"class_i_alleles": ["HLA-A*02:01"]}))
    sfile = tmp_path / "src"
    sfile.write_text("code")
    monkeypatch.setattr(core, "_source_inputs", lambda config: (sfile,))
    monkeypatch.setattr(core, "_all_inputs", lambda config, has_pvac: (sfile,))
    monkeypatch.setattr(core, "verify_git_tracked_clean", lambda config: (True, {"x": "CLEAN"}))
    monkeypatch.setattr(core, "verify_ensembl_used", lambda config, records, require_nonempty: (True, None))
    monkeypatch.setattr(core.u, "verify_tool_commits", lambda: (True, {"PRIME": {"match": True}}))
    monkeypatch.setattr(core.u, "verify_frozen_module_integrity", lambda: (True, {"module": "m"}))
    monkeypatch.setattr(core.u, "normalize_pass_vcf", lambda pv, ref, out: pv)
    monkeypatch.setattr(core.u, "rna_alt_evidence", lambda v, rna_bam=None: ({}, "COMPUTED"))
    monkeypatch.setattr(core.u, "gene_tpm_by_ensg", lambda q: {})
    monkeypatch.setattr(core.u, "score_universe", lambda x: x)
    monkeypatch.setattr(core, "load_filtered_variants",
                        lambda config, p: pd.DataFrame({"key": ["6:1:C:T"], "pass_filters": [True],
                                                        "strict5_pass": [True]}))
    uni = pd.DataFrame({"patient_id": ["Hu_315"], "mutation_id": ["6:1:C:T"],
                        "candidate_source": ["lossless_recovery"], "mutant_peptide": ["AAAAAAAAA"],
                        "hla_allele": ["HLA-A*02:01"], "genuine_prime": [0.9], "epicurus": [0.5]})
    monkeypatch.setattr(core, "build_universe", lambda config, v, h, t, r: (uni, [], False, [{"used": 1}]))
    monkeypatch.setattr(core.u, "arm_selection", lambda uni, arm: pd.DataFrame(
        {"mutation_id": ["6:1:C:T"], "mutant_peptide": ["AAAAAAAAA"], "hla_allele": ["HLA-A*02:01"]}))
    # tripwire: freeze must open NO recognition-label file
    import builtins
    real_open = builtins.open
    labels_path = str(u.LABELS)
    monkeypatch.setattr(builtins, "open", lambda f, *a, **k: (_ for _ in ()).throw(
        AssertionError("freeze opened LABELS")) if str(f) == labels_path else real_open(f, *a, **k))

    m = core.freeze(c)
    assert m["LOCK"] == "FROZEN_NO_LABELS" and m["labels_opened"] is False
    assert m["patient_id"] == "Hu_315"                          # real id in the manifest
    assert m["sample_names"] == {"normal": "Hu_315_N", "tumor": "Hu_315_T"}
    assert m["provenance_lane"] == "generic-miller-universe-core"
    assert (c.freeze_dir / "FREEZE_MANIFEST.json").exists()
    assert (c.artifact_dir / "FREEZE_MANIFEST.json").exists()   # provenance copied to milestone-8 art dir
    # the pinned code_files list carries the frozen script + the patient resolver
    assert "scripts/miller_hu287_universe.py" in m["code_files"]
    assert "src/benchmark/miller_patient.py" in m["code_files"]
