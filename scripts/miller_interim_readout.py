"""Once-only immutable interim readout for one frozen Miller calibration patient.

This is intentionally not a policy-selection surface.  It reuses the hardened
multi-patient evaluator, snapshots every label-blind byte before outcome access,
creates a durable exclusive claim, reads the outcome table exactly once through
the sole audited reader, and persists aggregate per-arm metrics only.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from benchmark import miller_storage_lifecycle as life
from benchmark.miller_generalization_eval import (
    REGISTERED_ARMS,
    _atomic_write_json,
    _exclusive_claim,
    _read_labels_once,
    _snapshot_stage_inputs,
    _validate_phase_label_support,
    default_config,
    evaluate_stage,
    preflight_patient,
)


ROOT = Path(__file__).resolve().parents[1]
ADDENDUM = ROOT / "artifacts/milestone_8_generalization/HU_315_INTERIM_READOUT_ADDENDUM_2026-07-13.md"
OUTPUT_ROOT = ROOT / "artifacts/milestone_8_generalization/interim_readouts"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_commit() -> str | None:
    try:
        return subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
            capture_output=True,
            check=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def run_interim(patient_id: str) -> dict:
    output_dir = OUTPUT_ROOT / patient_id
    result_path = output_dir / "INTERIM_RESULT.json"
    claim_path = output_dir / "INTERIM_UNSEAL_STARTED.json"
    if result_path.is_file():
        return json.loads(result_path.read_text())
    if claim_path.exists():
        return {
            "status": "INTERIM_UNSEAL_INCOMPLETE",
            "patient_id": patient_id,
            "reason": "a prior attempt crossed the one-time outcome boundary; outcomes will not be reopened",
        }
    if not ADDENDUM.is_file():
        return {"status": "MISSING_PREREGISTERED_ADDENDUM", "patient_id": patient_id}

    config = default_config(output_dir=output_dir)
    by_id = {patient.patient_id: patient for patient in config.calibration_patients}
    if patient_id not in by_id:
        return {
            "status": "REFUSED_NOT_CALIBRATION",
            "patient_id": patient_id,
            "reason": "interim readout is forbidden for final-held-out or non-calibration patients",
        }
    patient = by_id[patient_id]
    ok, preflight = preflight_patient(config, patient.patient_id, patient.freeze_dir)
    if not ok:
        return {"status": "PREFLIGHT_FAILED", "patient_id": patient_id, "preflight": preflight}
    manifest_hash = preflight["manifest_sha256"]
    try:
        snapshots = _snapshot_stage_inputs(
            config,
            (patient,),
            REGISTERED_ARMS,
            {patient_id: manifest_hash},
        )
    except (OSError, KeyError, ValueError) as exc:
        return {"status": "SNAPSHOT_FAILED", "patient_id": patient_id, "reason": str(exc)}

    runtime_commit = git_commit()
    claim = {
        "CLAIM": "IMMUTABLE_CALIBRATION_INTERIM_UNSEAL_STARTED",
        "patient_id": patient_id,
        "labels_opened": False,
        "policy_selection_allowed": False,
        "future_status": "OBSERVED_DEVELOPMENT_ONLY",
        "freeze_manifest_sha256": manifest_hash,
        "registered_arms": list(REGISTERED_ARMS),
        "addendum": {"path": str(ADDENDUM.relative_to(ROOT)), "sha256": sha256_file(ADDENDUM)},
        "evaluator": {
            "path": str(Path(__file__).resolve().relative_to(ROOT)),
            "sha256": sha256_file(Path(__file__).resolve()),
            "git_commit": runtime_commit,
        },
    }
    if not _exclusive_claim(claim_path, claim):
        return {"status": "INTERIM_UNSEAL_INCOMPLETE", "patient_id": patient_id}

    try:
        labels = _read_labels_once(config.labels_path)  # sole audited outcome read
        _validate_phase_label_support(labels, (patient,))
    except (OSError, ValueError) as exc:
        return {"status": "LABELS_INVALID", "patient_id": patient_id, "reason": str(exc)}

    evaluation = evaluate_stage(
        config,
        (patient,),
        labels,
        snapshots=snapshots,
        arm_ids=REGISTERED_ARMS,
    )
    result = {
        "status": "INTERIM_EVALUATED",
        "patient_id": patient_id,
        "stage": "CALIBRATION_DEVELOPMENT_INTERIM",
        "labels_opened": True,
        "label_file_reads_this_run": 1,
        "immutable": True,
        "policy_selection_allowed": False,
        "may_inform_model_changes": False,
        "final_held_out_affected": False,
        "freeze_manifest_sha256": manifest_hash,
        "registered_arms": list(REGISTERED_ARMS),
        "evaluation": evaluation,
        "provenance": claim,
    }
    _atomic_write_json(result_path, result)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patient_id", choices=("Hu_315",))
    args = parser.parse_args(argv)
    result = run_interim(args.patient_id)
    print(json.dumps(result, indent=2, default=str))
    return 0 if result.get("status") == "INTERIM_EVALUATED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
