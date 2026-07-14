"""Resumable, label-blind end-to-end Miller patient reconstruction driver.

The driver composes the already-audited stage CLIs; it does not implement a
second biological path.  Every stage is idempotent/fail-closed, state is written
outside the raw patient directory, and destructive cleanup remains impossible
unless ``--execute-cleanup`` is explicit and the independent freeze verifier
passes.  No recognition or outcome source is imported or opened here.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from benchmark.miller_patient import load_patient


ROOT = Path(__file__).resolve().parents[1]
STAGES = ("metadata", "download", "convert", "quant", "wes", "hla", "rna", "freeze", "verify")


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def stage_command(patient_id: str, stage: str, threads: int) -> list[str]:
    python = sys.executable
    if stage in {"metadata", "download", "convert", "attest-fastq", "quant", "wes", "hla", "rna"}:
        return [
            python,
            "-m",
            "scripts.miller_patient_reconstruct",
            patient_id,
            stage,
            "--threads",
            str(threads),
        ]
    if stage == "freeze":
        return [python, "-m", "scripts.miller_patient_universe", patient_id, "freeze"]
    if stage == "verify":
        return [python, "-m", "scripts.miller_patient_pipeline", patient_id, "--verify-only"]
    raise ValueError(f"unknown stage: {stage}")


def selected_stages(from_stage: str, through_stage: str) -> tuple[str, ...]:
    start = STAGES.index(from_stage)
    end = STAGES.index(through_stage)
    if end < start:
        raise ValueError("--through-stage precedes --from-stage")
    return STAGES[start : end + 1]


def sleep_guard_command(argv: list[str]) -> list[str] | None:
    """Return the macOS guard command, or None when unavailable/already guarded."""
    if os.environ.get("EPICURUS_CAFFEINATED") == "1" or sys.platform != "darwin":
        return None
    caffeinate = shutil.which("caffeinate")
    if caffeinate is None:
        return None
    return [caffeinate, "-dimsu", sys.executable, "-m", "scripts.miller_patient_pipeline", *argv]


def verify_only(patient_id: str) -> dict:
    # Local import keeps the normal orchestration surface small and makes the
    # independent verifier, rather than pipeline state, authoritative.
    from benchmark.miller_storage_lifecycle import PatientPaths, verify_frozen_no_labels

    paths = PatientPaths.from_patient_id(patient_id)
    ok, detail = verify_frozen_no_labels(paths.freeze_dir, paths.root)
    return {"patient_id": patient_id, "frozen_verified": ok, "verification": detail}


def run_pipeline(args: argparse.Namespace) -> dict:
    patient = load_patient(args.patient_id)
    state_path = patient.artifact_dir / "PIPELINE_STATE.json"
    state = {
        "patient_id": patient.patient_id,
        "isolation": "LOCKED_TEST: no recognition/outcome data read",
        "labels_opened": False,
        "started_at": utc_now(),
        "threads": args.threads,
        "requested_stages": list(selected_stages(args.from_stage, args.through_stage)),
        "stages": {},
        "status": "RUNNING",
    }
    atomic_json(state_path, state)

    for stage in state["requested_stages"]:
        command = stage_command(patient.patient_id, stage, args.threads)
        record = {"status": "RUNNING", "started_at": utc_now(), "command": command}
        state["stages"][stage] = record
        atomic_json(state_path, state)
        completed = subprocess.run(command, cwd=ROOT)
        record.update(
            status="COMPLETE" if completed.returncode == 0 else "FAILED",
            completed_at=utc_now(),
            returncode=completed.returncode,
        )
        if completed.returncode != 0:
            state["status"] = "FAILED"
            state["failed_stage"] = stage
            atomic_json(state_path, state)
            raise subprocess.CalledProcessError(completed.returncode, command)
        atomic_json(state_path, state)

    # A dry run is mandatory after a verified freeze.  Execution is optional and
    # remains guarded again inside the lifecycle CLI by typed patient confirmation.
    if args.through_stage == "verify":
        # Before reclaiming FASTQs, pin the exact raw bytes Epicurus consumed.
        # The nextNEOpi Track-A bundle later rejects regenerated FASTQs unless
        # they match this attestation, preserving identical-input semantics.
        attest = stage_command(patient.patient_id, "attest-fastq", args.threads)
        subprocess.run(attest, cwd=ROOT, check=True)
        state["raw_input_attestation"] = str(patient.artifact_dir / "CONVERT_PROVENANCE.json")
        atomic_json(state_path, state)
        dry_report = patient.artifact_dir / "STORAGE_CLEANUP_DRY_RUN.json"
        dry = [
            sys.executable,
            "-m",
            "scripts.miller_storage_cleanup",
            patient.patient_id,
            "--report",
            str(dry_report),
        ]
        subprocess.run(dry, cwd=ROOT, check=True)
        state["cleanup_dry_run"] = str(dry_report)
        if args.execute_cleanup:
            executed_report = patient.artifact_dir / "STORAGE_CLEANUP_EXECUTED.json"
            execute = [
                sys.executable,
                "-m",
                "scripts.miller_storage_cleanup",
                patient.patient_id,
                "--execute",
                "--confirm-patient",
                patient.patient_id,
                "--report",
                str(executed_report),
            ]
            subprocess.run(execute, cwd=ROOT, check=True)
            state["cleanup_executed"] = str(executed_report)

    state["status"] = "COMPLETE"
    state["completed_at"] = utc_now()
    atomic_json(state_path, state)
    return state


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("patient_id")
    p.add_argument("--threads", type=int, default=12)
    p.add_argument("--from-stage", choices=STAGES, default="metadata")
    p.add_argument("--through-stage", choices=STAGES, default="verify")
    p.add_argument("--execute-cleanup", action="store_true")
    p.add_argument("--verify-only", action="store_true", help=argparse.SUPPRESS)
    p.add_argument("--no-sleep-guard", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    raw_argv = list(sys.argv[1:] if argv is None else argv)
    args = parser().parse_args(raw_argv)
    if args.threads < 1:
        raise SystemExit("--threads must be positive")
    if args.execute_cleanup and args.through_stage != "verify":
        raise SystemExit("--execute-cleanup requires --through-stage verify")
    if args.verify_only:
        result = verify_only(args.patient_id)
        print(json.dumps(result, indent=2, sort_keys=True))
        return 0 if result["frozen_verified"] else 1
    if not args.no_sleep_guard:
        guard = sleep_guard_command(raw_argv)
        if guard is not None:
            env = os.environ.copy()
            env["EPICURUS_CAFFEINATED"] = "1"
            os.execvpe(guard[0], guard, env)
    result = run_pipeline(args)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
