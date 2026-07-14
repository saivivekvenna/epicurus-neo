from __future__ import annotations

import json
from argparse import Namespace
from dataclasses import dataclass
from pathlib import Path
from types import SimpleNamespace

import pytest

from scripts import miller_patient_pipeline as pipeline


def test_stage_order_is_complete_end_to_end_and_sliceable():
    assert pipeline.STAGES == (
        "metadata", "download", "convert", "quant", "wes", "hla", "rna", "freeze", "verify"
    )
    assert pipeline.selected_stages("wes", "rna") == ("wes", "hla", "rna")
    with pytest.raises(ValueError, match="precedes"):
        pipeline.selected_stages("freeze", "convert")


def test_commands_use_modules_without_shell_and_keep_patient_identity():
    command = pipeline.stage_command("Hu_277", "wes", 12)
    assert command[1:5] == ["-m", "scripts.miller_patient_reconstruct", "Hu_277", "wes"]
    assert command[-2:] == ["--threads", "12"]
    freeze = pipeline.stage_command("Hu_277", "freeze", 12)
    assert freeze[1:] == ["-m", "scripts.miller_patient_universe", "Hu_277", "freeze"]


def test_atomic_state_write(tmp_path):
    target = tmp_path / "nested" / "state.json"
    pipeline.atomic_json(target, {"status": "COMPLETE"})
    assert json.loads(target.read_text()) == {"status": "COMPLETE"}
    assert not list(target.parent.glob(f".{target.name}.*"))


def test_sleep_guard_is_recursive_safe(monkeypatch):
    monkeypatch.setattr(pipeline.sys, "platform", "darwin")
    monkeypatch.setattr(pipeline.shutil, "which", lambda name: "/usr/bin/caffeinate")
    monkeypatch.delenv("EPICURUS_CAFFEINATED", raising=False)
    command = pipeline.sleep_guard_command(["Hu_277", "--threads", "12"])
    assert command[:2] == ["/usr/bin/caffeinate", "-dimsu"]
    monkeypatch.setenv("EPICURUS_CAFFEINATED", "1")
    assert pipeline.sleep_guard_command(["Hu_277"]) is None


def test_pipeline_mandates_dry_run_and_executes_cleanup_only_when_explicit(tmp_path, monkeypatch):
    @dataclass
    class Patient:
        patient_id: str = "Hu_test"
        artifact_dir: Path = tmp_path / "artifacts"

    monkeypatch.setattr(pipeline, "load_patient", lambda patient_id: Patient())
    commands = []

    def fake_run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(pipeline.subprocess, "run", fake_run)
    args = Namespace(
        patient_id="Hu_test",
        threads=8,
        from_stage="freeze",
        through_stage="verify",
        execute_cleanup=False,
    )
    result = pipeline.run_pipeline(args)
    assert result["labels_opened"] is False and result["status"] == "COMPLETE"
    cleanup = [c for c in commands if "scripts.miller_storage_cleanup" in c]
    assert len(cleanup) == 1 and "--execute" not in cleanup[0]

    commands.clear()
    args.execute_cleanup = True
    pipeline.run_pipeline(args)
    cleanup = [c for c in commands if "scripts.miller_storage_cleanup" in c]
    assert len(cleanup) == 2
    assert "--execute" not in cleanup[0]
    assert cleanup[1][-4:] == ["--confirm-patient", "Hu_test", "--report", str(tmp_path / "artifacts" / "STORAGE_CLEANUP_EXECUTED.json")]
