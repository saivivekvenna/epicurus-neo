"""Fail-closed post-freeze storage lifecycle: hardened verification (path/symlink/traversal safety,
code-file pinning, frozen-module integrity, commit existence), dry-run safety, containment, and guarded
execution against a hermetic fixture. No real cohort files are asserted."""

from __future__ import annotations

import json
import os
from pathlib import Path

from benchmark import miller_storage_lifecycle as life
from benchmark.miller_storage_lifecycle import PatientPaths, classify_entries, execute_cleanup, plan_cleanup


def CLEAN(rel, root):                         # inject: pinned provenance is git-clean in the fixture
    return "CLEAN"


def COMMIT_OK(sha, root):                     # inject: the recorded commit resolves to an object
    return True


def BLOB(sha, rel, root):                     # inject: fixture's tracked file is its historical blob
    path = root / rel
    return life._sha256(path) if path.is_file() else None


def _snapshot(root: Path) -> set[str]:
    return {os.path.relpath(os.path.join(d, f), root)
            for d, _, fs in os.walk(root) for f in fs}


def _fixture(tmp_path: Path, *, valid_manifest: bool = True, mutate=None) -> PatientPaths:
    raw = tmp_path / "data/raw/miller_ipv/hu_test"
    freeze = raw / "freeze"
    for rel, content in {
        "SRR1.sra": b"\x00sra-archive",
        "DOWNLOAD_MANIFEST.json": b"{}",
        "CONVERT_PROVENANCE.json": b"{}",
        "fastq/SRR1_1.fastq": b"@r1\nACGT\n+\nFFFF\n",
        "fastq/SRR1_2.fastq": b"@r2\nTGCA\n+\nFFFF\n",
        "hla/hla_1.fq": b"@h1\nAC\n+\nFF\n",
        "hla/hla_2.fq": b"@h2\nGT\n+\nFF\n",
        "hla/optitype/hu_test_result.tsv": b"A1\tA2\n",
        "somatic/hu_test.somatic.pass.vcf.gz": b"\x1f\x8bvcf",
        "somatic/hu_test_N.md.bam": b"BAMbam",
        "somatic/hu_test_N.md.bam.bai": b"bai",
        "salmon_quant/quant.sf": b"Name\tTPM\n",
        "ensembl_cache/abc.json": b'{"cached":true}',
        "freeze/universe.csv": b"mutation_id,epicurus\n1:1:A:T,0.2\n",
        "freeze/variants.csv": b"key\n1:1:A:T\n",
    }.items():
        p = raw / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_bytes(content)
    # repo-root code/config inputs the freeze pins (relative to the fixture root=tmp_path)
    (tmp_path / "code").mkdir(exist_ok=True)
    (tmp_path / "code/mod.py").write_bytes(b"# pinned code file\n")
    (tmp_path / "code/fmod.py").write_bytes(b"# frozen scoring module\n")

    vcf_rel = "data/raw/miller_ipv/hu_test/somatic/hu_test.somatic.pass.vcf.gz"
    man = {
        "patient_id": "Hu_test",
        "LOCK": "FROZEN_NO_LABELS" if valid_manifest else "OPEN",
        "labels_opened": False,
        "git_commit": "a" * 40,
        "sha256": {"universe.csv": life._sha256(freeze / "universe.csv"),
                   "variants.csv": life._sha256(freeze / "variants.csv")},
        "input_sha256": {
            vcf_rel: life._sha256(raw / "somatic/hu_test.somatic.pass.vcf.gz"),
            "code/mod.py": life._sha256(tmp_path / "code/mod.py"),
        },
        "code_files": ["code/mod.py"],
        "git_tracked_clean": {"code/mod.py": "CLEAN"},
        "frozen_module_integrity": {"module": "code/fmod.py",
                                    "module_sha256": life._sha256(tmp_path / "code/fmod.py")},
    }
    if mutate is not None:
        mutate(man, raw, freeze, tmp_path)
    (freeze / "FREEZE_MANIFEST.json").write_text(json.dumps(man))
    return PatientPaths(patient_id="Hu_test", raw_dir=raw, freeze_dir=freeze, root=tmp_path)


def _verify(paths):
    return life.verify_frozen_no_labels(
        paths.freeze_dir, paths.root, git_status=CLEAN, commit_exists=COMMIT_OK,
        git_blob_sha256=BLOB,
    )


# ---- classification ---------------------------------------------------------------------------------
def test_classification_reclaims_only_documented_regenerable(tmp_path):
    paths = _fixture(tmp_path)
    by_rel = {e["rel"]: e for e in classify_entries(paths)}
    for rel in ("fastq/SRR1_1.fastq", "fastq/SRR1_2.fastq", "hla/hla_1.fq", "hla/hla_2.fq"):
        assert by_rel[rel]["category"] == "REMOVE"
        assert by_rel[rel]["regeneration"]
    for rel in ("SRR1.sra", "DOWNLOAD_MANIFEST.json", "somatic/hu_test_N.md.bam",
                "somatic/hu_test.somatic.pass.vcf.gz", "salmon_quant/quant.sf",
                "hla/optitype/hu_test_result.tsv", "ensembl_cache/abc.json",
                "freeze/universe.csv", "freeze/FREEZE_MANIFEST.json"):
        assert by_rel[rel]["category"] == "PRESERVE", rel


def test_unclassified_file_is_preserved_fail_closed(tmp_path):
    paths = _fixture(tmp_path)
    (paths.raw_dir / "mystery.dat").write_bytes(b"???")
    e = {x["rel"]: x for x in classify_entries(paths)}["mystery.dat"]
    assert e["category"] == "PRESERVE" and "fail-closed" in e["reason"]


# ---- dry-run is the default and never mutates -------------------------------------------------------
def test_plan_is_dry_run_and_never_mutates(tmp_path):
    paths = _fixture(tmp_path)
    before = _snapshot(tmp_path)
    report = plan_cleanup(paths, git_status=CLEAN, commit_exists=COMMIT_OK, git_blob_sha256=BLOB)
    assert report["mode"] == "dry_run"
    assert _snapshot(tmp_path) == before


def test_reclaim_withheld_until_freeze_verifies(tmp_path):
    paths = _fixture(tmp_path)
    (paths.freeze_dir / "FREEZE_MANIFEST.json").unlink()
    report = plan_cleanup(paths, git_status=CLEAN, commit_exists=COMMIT_OK, git_blob_sha256=BLOB)
    assert report["frozen_verified"] is False
    assert report["reclaimable"] == [] and report["reclaimable_bytes"] == 0
    assert {e["rel"] for e in report["withheld_pending_freeze"]} >= {"fastq/SRR1_1.fastq", "hla/hla_1.fq"}


# ---- successful planning against the fixture --------------------------------------------------------
def test_successful_planning_lists_reclaimable_with_bytes_and_reasons(tmp_path):
    paths = _fixture(tmp_path)
    report = plan_cleanup(paths, git_status=CLEAN, commit_exists=COMMIT_OK, git_blob_sha256=BLOB)
    assert report["frozen_verified"] is True
    v = report["verification"]["checks"]
    assert v["outputs_verified"] == 2 and v["inputs_verified"] == 2 and v["frozen_module_intact"] is True
    assert report["manifest_sha256"] and len(report["manifest_sha256"]) == 64
    rels = {e["rel"] for e in report["reclaimable"]}
    assert rels == {"fastq/SRR1_1.fastq", "fastq/SRR1_2.fastq", "hla/hla_1.fq", "hla/hla_2.fq"}
    assert report["reclaimable_bytes"] == sum(e["size_bytes"] for e in report["reclaimable"]) > 0
    assert report["summary"]["REMOVE"]["n"] == 4 and report["summary"]["PRESERVE"]["n"] >= 9


# ---- fail-closed verification gate: basic -----------------------------------------------------------
def test_verify_fails_on_missing_manifest(tmp_path):
    paths = _fixture(tmp_path)
    (paths.freeze_dir / "FREEZE_MANIFEST.json").unlink()
    ok, detail = _verify(paths)
    assert ok is False and "no FROZEN_NO_LABELS manifest" in detail["reason"]


def test_verify_fails_on_wrong_lock_or_opened_labels(tmp_path):
    ok, _ = _verify(_fixture(tmp_path, valid_manifest=False))
    assert ok is False
    paths2 = _fixture(tmp_path / "b", mutate=lambda m, *a: m.__setitem__("labels_opened", True))
    ok2, d2 = _verify(paths2)
    assert ok2 is False and "labels_opened" in d2["reason"]


def test_verify_fails_on_tampered_output_or_input(tmp_path):
    paths = _fixture(tmp_path)
    (paths.freeze_dir / "universe.csv").write_bytes(b"TAMPERED")
    ok, d = _verify(paths)
    assert ok is False and "output hash mismatch" in d["reason"]
    paths2 = _fixture(tmp_path / "b")
    (paths2.raw_dir / "somatic/hu_test.somatic.pass.vcf.gz").write_bytes(b"TAMPERED")
    ok2, d2 = _verify(paths2)
    assert ok2 is False and "input hash mismatch" in d2["reason"]


def test_verify_fails_on_dirty_provenance_recorded_but_uses_historical_blob_not_live_tree(tmp_path):
    # recorded non-CLEAN in the manifest
    paths = _fixture(tmp_path, mutate=lambda m, *a: m["git_tracked_clean"].__setitem__("code/mod.py", "UNSTAGED_MODIFIED"))
    ok, d = _verify(paths)
    assert ok is False and "recorded non-CLEAN" in d["reason"]
    # Recorded CLEAN and historical blob intact: later worktree evolution must not invalidate the freeze.
    paths2 = _fixture(tmp_path / "b")
    historical = {
        rel: life._sha256(paths2.root / rel) for rel in ("code/mod.py", "code/fmod.py")
    }
    (paths2.root / "code/mod.py").write_bytes(b"# evolved after freeze\n")
    ok2, d2 = life.verify_frozen_no_labels(paths2.freeze_dir, paths2.root,
                                           git_status=lambda rel, root: "STAGED_MODIFIED",
                                           commit_exists=COMMIT_OK,
                                           git_blob_sha256=lambda sha, rel, root: historical.get(rel))
    assert ok2 is True, d2


# ---- adversarial path safety (traversal / symlink / absolute) --------------------------------------
def test_verify_rejects_output_key_traversal_absolute_and_symlink(tmp_path):
    hexsha = "b" * 64
    p1 = _fixture(tmp_path / "a", mutate=lambda m, *a: m["sha256"].__setitem__("../escape.txt", hexsha))
    assert "unsafe output key" in _verify(p1)[1]["reason"]
    p2 = _fixture(tmp_path / "b", mutate=lambda m, *a: m["sha256"].__setitem__("/abs.csv", hexsha))
    assert "unsafe output key" in _verify(p2)[1]["reason"]

    def add_symlink_output(m, raw, freeze, root):
        (freeze / "evil.csv").symlink_to(freeze / "universe.csv")
        m["sha256"]["evil.csv"] = life._sha256(freeze / "universe.csv")
    p3 = _fixture(tmp_path / "c", mutate=add_symlink_output)
    assert "unsafe output key" in _verify(p3)[1]["reason"]


def test_verify_rejects_input_key_traversal_and_symlink(tmp_path):
    hexsha = "c" * 64
    p1 = _fixture(tmp_path / "a", mutate=lambda m, *a: m["input_sha256"].__setitem__("../../etc/passwd", hexsha))
    assert "unsafe input key" in _verify(p1)[1]["reason"]

    def add_symlink_input(m, raw, freeze, root):
        (root / "linky").symlink_to(raw / "somatic/hu_test.somatic.pass.vcf.gz")
        m["input_sha256"]["linky"] = life._sha256(raw / "somatic/hu_test.somatic.pass.vcf.gz")
    p2 = _fixture(tmp_path / "b", mutate=add_symlink_input)
    assert "unsafe input key" in _verify(p2)[1]["reason"]


# ---- adversarial: symlinked PARENT component whose target stays inside base ------------------------
def test_verify_rejects_symlinked_parent_component_inside_base(tmp_path):
    def out_parent(m, raw, freeze, root):
        (freeze / "realsub").mkdir()
        (freeze / "realsub" / "u2.csv").write_bytes(b"x")
        (freeze / "linksub").symlink_to(freeze / "realsub", target_is_directory=True)
        m["sha256"]["linksub/u2.csv"] = life._sha256(freeze / "realsub" / "u2.csv")
    assert "unsafe output key" in _verify(_fixture(tmp_path / "a", mutate=out_parent))[1]["reason"]

    def in_parent(m, raw, freeze, root):
        (root / "realdir").mkdir()
        (root / "realdir" / "x").write_bytes(b"y")
        (root / "linkdir").symlink_to(root / "realdir", target_is_directory=True)
        m["input_sha256"]["linkdir/x"] = life._sha256(root / "realdir" / "x")
    assert "unsafe input key" in _verify(_fixture(tmp_path / "b", mutate=in_parent))[1]["reason"]

    def mod_parent(m, raw, freeze, root):
        (root / "realmod").mkdir()
        (root / "realmod" / "fmod.py").write_bytes(b"# f\n")
        (root / "linkmod").symlink_to(root / "realmod", target_is_directory=True)
        m["frozen_module_integrity"] = {"module": "linkmod/fmod.py",
                                        "module_sha256": life._sha256(root / "realmod" / "fmod.py")}
    assert "unsafe frozen module path" in _verify(_fixture(tmp_path / "c", mutate=mod_parent))[1]["reason"]


def test_verify_rejects_symlinked_manifest(tmp_path):
    paths = _fixture(tmp_path)
    mp = paths.freeze_dir / "FREEZE_MANIFEST.json"
    real = paths.freeze_dir / "real_manifest.json"
    mp.rename(real)
    mp.symlink_to(real)
    ok, d = _verify(paths)
    assert ok is False and "is a symlink" in d["reason"]


# ---- classification proves the regeneration SOURCE exists before REMOVE ----------------------------
def test_missing_sra_withholds_fastq_fail_closed(tmp_path):
    paths = _fixture(tmp_path)
    (paths.raw_dir / "SRR1.sra").unlink()                     # regeneration source gone
    by = {e["rel"]: e for e in classify_entries(paths)}
    for rel in ("fastq/SRR1_1.fastq", "fastq/SRR1_2.fastq"):
        assert by[rel]["category"] == "PRESERVE" and "regeneration source" in by[rel]["reason"]
    rep = plan_cleanup(paths, git_status=CLEAN, commit_exists=COMMIT_OK, git_blob_sha256=BLOB)
    assert not any(e["rel"].startswith("fastq/") for e in rep["reclaimable"])


def test_missing_normal_bam_withholds_hla_reads_fail_closed(tmp_path):
    paths = _fixture(tmp_path)
    (paths.raw_dir / "somatic/hu_test_N.md.bam").unlink()     # no normal alignment to re-extract from
    by = {e["rel"]: e for e in classify_entries(paths)}
    for rel in ("hla/hla_1.fq", "hla/hla_2.fq"):
        assert by[rel]["category"] == "PRESERVE" and "no preserved normal alignment" in by[rel]["reason"]


def test_hla_recipe_is_generic_not_hu287_specific(tmp_path):
    paths = _fixture(tmp_path)
    by = {e["rel"]: e for e in classify_entries(paths)}
    recipe = by["hla/hla_1.fq"]["regeneration"]
    assert "miller_hu287_hla.sh" not in recipe
    assert "miller_patient_reconstruct.py" in recipe and "hla" in recipe


# ---- code_files pinning + frozen module + commit ---------------------------------------------------
def test_verify_requires_nonempty_code_files_pinned_in_both_maps(tmp_path):
    p_empty = _fixture(tmp_path / "a", mutate=lambda m, *a: m.__setitem__("code_files", []))
    assert "code_files missing/empty" in _verify(p_empty)[1]["reason"]
    p_untracked = _fixture(tmp_path / "b", mutate=lambda m, *a: m["code_files"].append("code/extra.py"))
    assert "not fully pinned" in _verify(p_untracked)[1]["reason"]

    def add_to_tracked_only(m, raw, freeze, root):
        m["code_files"].append("code/extra.py")
        m["git_tracked_clean"]["code/extra.py"] = "CLEAN"        # present in tracked but NOT in input_sha256
    p_half = _fixture(tmp_path / "c", mutate=add_to_tracked_only)
    d = _verify(p_half)[1]
    assert "not fully pinned" in d["reason"] and "input_sha256=['code/extra.py']" in d["reason"]


def test_verify_requires_complete_matching_frozen_module(tmp_path):
    p_missing = _fixture(tmp_path / "a", mutate=lambda m, *a: m.pop("frozen_module_integrity"))
    assert "frozen_module_integrity missing/incomplete" in _verify(p_missing)[1]["reason"]
    p_bad = _fixture(tmp_path / "b",
                     mutate=lambda m, *a: m["frozen_module_integrity"].__setitem__("module_sha256", "d" * 64))
    assert "frozen module hash mismatch" in _verify(p_bad)[1]["reason"]
    p_trav = _fixture(tmp_path / "c",
                      mutate=lambda m, *a: m["frozen_module_integrity"].__setitem__("module", "../fmod.py"))
    assert "unsafe frozen module path" in _verify(p_trav)[1]["reason"]


def test_verify_requires_valid_existing_commit(tmp_path):
    p_short = _fixture(tmp_path / "a", mutate=lambda m, *a: m.__setitem__("git_commit", "abc123"))
    assert "not a full 40-hex sha" in _verify(p_short)[1]["reason"]
    paths = _fixture(tmp_path / "b")
    ok, d = life.verify_frozen_no_labels(paths.freeze_dir, paths.root, git_status=CLEAN,
                                         commit_exists=lambda sha, root: False,
                                         git_blob_sha256=BLOB)
    assert ok is False and "does not resolve to an existing commit" in d["reason"]


# ---- execute refuses unless verified AND confirmed --------------------------------------------------
def test_execute_refuses_without_confirm_even_when_verified(tmp_path):
    paths = _fixture(tmp_path)
    before = _snapshot(tmp_path)
    report = execute_cleanup(paths, confirm=False, git_status=CLEAN, commit_exists=COMMIT_OK,
                             git_blob_sha256=BLOB)
    assert report["status"] == "REFUSED_NO_CONFIRM" and report["deleted"] == []
    assert _snapshot(tmp_path) == before


def test_execute_refuses_when_freeze_unverified(tmp_path):
    paths = _fixture(tmp_path)
    (paths.freeze_dir / "universe.csv").write_bytes(b"TAMPERED")
    before = _snapshot(tmp_path)
    report = execute_cleanup(paths, confirm=True, git_status=CLEAN, commit_exists=COMMIT_OK,
                             git_blob_sha256=BLOB)
    assert report["status"] == "REFUSED_UNVERIFIED_FREEZE" and report["deleted"] == []
    assert _snapshot(tmp_path) == before


def test_execute_deletes_only_regenerable_when_verified_and_confirmed(tmp_path):
    paths = _fixture(tmp_path)
    report = execute_cleanup(paths, confirm=True, git_status=CLEAN, commit_exists=COMMIT_OK,
                             git_blob_sha256=BLOB)
    assert report["status"] == "EXECUTED"
    survivors = _snapshot(tmp_path)
    raw = "data/raw/miller_ipv/hu_test"
    for rel in ("fastq/SRR1_1.fastq", "hla/hla_1.fq", "hla/hla_2.fq"):
        assert f"{raw}/{rel}" not in survivors
    for rel in ("SRR1.sra", "somatic/hu_test_N.md.bam", "somatic/hu_test.somatic.pass.vcf.gz",
                "salmon_quant/quant.sf", "freeze/universe.csv", "freeze/FREEZE_MANIFEST.json"):
        assert f"{raw}/{rel}" in survivors
    assert report["deleted_bytes"] > 0
    assert {d["rel"] for d in report["deleted"]} == {
        "fastq/SRR1_1.fastq", "fastq/SRR1_2.fastq", "hla/hla_1.fq", "hla/hla_2.fq"}


# ---- symlink + label-table safety (hermetic: label path monkeypatched to a tmp file) ---------------
def test_symlinks_are_skipped_and_labels_are_never_touched(tmp_path, monkeypatch):
    fake_labels = tmp_path / "fake_cohort_labels.csv"
    fake_labels.write_bytes(b"peptide,recognized\nX,1\n")
    monkeypatch.setattr(life, "LABELS_PATH", fake_labels)
    paths = _fixture(tmp_path)
    link = paths.raw_dir / "fastq" / "SRR1_3.fastq"          # looks reclaimable, points at the label table
    link.symlink_to(fake_labels)
    entries = {e["rel"]: e for e in classify_entries(paths)}
    assert entries["fastq/SRR1_3.fastq"]["category"] == "SKIP_SYMLINK"
    report = execute_cleanup(paths, confirm=True, git_status=CLEAN, commit_exists=COMMIT_OK,
                             git_blob_sha256=BLOB)
    assert report["status"] == "EXECUTED"
    assert link.is_symlink()                                  # symlink untouched
    assert fake_labels.exists() and fake_labels.read_bytes() == b"peptide,recognized\nX,1\n"  # target untouched
    assert "fastq/SRR1_3.fastq" not in {d["rel"] for d in report["deleted"]}


def test_missing_raw_dir_yields_empty_classification(tmp_path):
    paths = PatientPaths(patient_id="Hu_absent", raw_dir=tmp_path / "nope",
                         freeze_dir=tmp_path / "nope/freeze", root=tmp_path)
    assert classify_entries(paths) == []
    report = plan_cleanup(paths, git_status=CLEAN, commit_exists=COMMIT_OK, git_blob_sha256=BLOB)
    assert report["frozen_verified"] is False and report["reclaimable"] == []
