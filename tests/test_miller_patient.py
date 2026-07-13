"""Label-blind patient resolution + reconstruction-driver env contract (no labels, no network, no tools)."""

from __future__ import annotations

import scripts.miller_patient_reconstruct as mpr
from benchmark.miller_patient import load_patient
from scripts.miller_patient_reconstruct import patient_manifest, script_env, script_for


# ---------------------------------------------------------------------------
# Exact Hu_287 backward-compatible resolution (the frozen north-star test case)
# ---------------------------------------------------------------------------
def test_resolves_hu287_runs_without_labels():
    p = load_patient("Hu_287")
    assert (p.normal_exome_run, p.tumor_exome_run, p.tumor_rna_run) == (
        "SRR24836184", "SRR24836169", "SRR24836183"
    )
    assert (p.normal_sample, p.tumor_sample) == ("Hu_287_N", "Hu_287_T")
    assert p.raw_dir.name == "hu_287"
    # Hu_287 keeps its milestone-7 reconstruction provenance dir (backward compatible)
    assert p.artifact_dir.as_posix().endswith(
        "artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction")


def test_hu287_manifest_declares_no_label_columns():
    manifest = patient_manifest("Hu_287")
    assert manifest["patient_id"] == "Hu_287"
    assert manifest["label_columns_read"] == []
    assert manifest["sample_names"] == {"normal": "Hu_287_N", "tumor": "Hu_287_T"}
    assert manifest["runs"] == {"normal_exome": "SRR24836184", "tumor_exome": "SRR24836169",
                                "tumor_rna": "SRR24836183"}


# ---------------------------------------------------------------------------
# Dynamic new-patient resolution (Hu_315 = a calibration-split patient)
# ---------------------------------------------------------------------------
def test_resolves_hu315_dynamically_without_labels():
    p = load_patient("Hu_315")
    # a real trio, distinct from Hu_287, resolved purely from SRA metadata
    assert len({p.normal_exome_run, p.tumor_exome_run, p.tumor_rna_run}) == 3
    assert {p.normal_exome_run, p.tumor_exome_run, p.tumor_rna_run}.isdisjoint(
        {"SRR24836184", "SRR24836169", "SRR24836183"})
    assert (p.normal_sample, p.tumor_sample) == ("Hu_315_N", "Hu_315_T")
    assert p.raw_dir.name == "hu_315"
    # non-Hu_287 patients route to the milestone-8 generalization tree
    assert p.artifact_dir.as_posix().endswith("artifacts/milestone_8_generalization/patients/Hu_315")


def test_hu315_manifest_declares_no_label_columns():
    p = load_patient("Hu_315")
    manifest = patient_manifest("Hu_315")
    assert manifest["patient_id"] == "Hu_315"
    assert manifest["label_columns_read"] == []
    assert manifest["runs"]["tumor_rna"] == p.tumor_rna_run
    assert "milestone_8_generalization" in manifest["artifact_dir"]


def test_unknown_patient_raises():
    import pytest
    with pytest.raises(ValueError, match="unknown Miller patient"):
        load_patient("Hu_000000")


# ---------------------------------------------------------------------------
# Reconstruction shell-script env contract (the four label-blind parameters)
# ---------------------------------------------------------------------------
def test_script_for_maps_stages_to_shell_scripts():
    assert script_for("wes").name == "miller_hu287_somatic.sh"
    assert script_for("hla").name == "miller_hu287_hla.sh"
    assert script_for("rna").name == "miller_hu287_rna.sh"


def test_script_env_is_label_blind_and_exact_for_hu287():
    env = script_env(load_patient("Hu_287"))
    # default-compatible: exactly the accessions the frozen scripts fall back to
    assert env == {"PATIENT_ID": "Hu_287", "NORMAL_EXOME_RUN": "SRR24836184",
                   "TUMOR_EXOME_RUN": "SRR24836169", "TUMOR_RNA_RUN": "SRR24836183"}
    # no recognition-label column is ever surfaced into the shell env
    assert "label" not in {k.lower() for k in env}


def test_script_env_is_dynamic_for_new_patient():
    p = load_patient("Hu_315")
    env = script_env(p)
    assert env["PATIENT_ID"] == "Hu_315"
    assert env["NORMAL_EXOME_RUN"] == p.normal_exome_run
    assert env["TUMOR_EXOME_RUN"] == p.tumor_exome_run
    assert env["TUMOR_RNA_RUN"] == p.tumor_rna_run
    # distinct from the Hu_287 defaults so the parameterized scripts do not silently reuse them
    assert env["NORMAL_EXOME_RUN"] != "SRR24836184"


def test_driver_dispatches_stage_with_parameterized_env(monkeypatch, capsys):
    """main('Hu_287', 'wes') runs the somatic script with the four env vars set (subprocess mocked)."""
    calls = {}

    def _fake_run(cmd, cwd=None, env=None, check=None):
        calls["cmd"] = cmd
        calls["cwd"] = cwd
        calls["env"] = env
        return None

    monkeypatch.setattr(mpr.subprocess, "run", _fake_run)
    rc = mpr.main(["Hu_287", "wes"])
    assert rc == 0
    assert calls["cmd"] == [str(script_for("wes"))]
    # the four label-blind parameters are threaded into the child environment
    for k, v in script_env(load_patient("Hu_287")).items():
        assert calls["env"][k] == v
    assert calls["env"]["THREADS"] == "4"
    assert str(calls["cwd"]) == str(mpr.ROOT)


def test_driver_threads_are_forwarded_to_reconstruction_script(monkeypatch):
    calls = {}

    def _fake_run(cmd, cwd=None, env=None, check=None):
        calls["env"] = env

    monkeypatch.setattr(mpr.subprocess, "run", _fake_run)
    assert mpr.main(["Hu_315", "rna", "--threads", "3"]) == 0
    assert calls["env"]["THREADS"] == "3"
