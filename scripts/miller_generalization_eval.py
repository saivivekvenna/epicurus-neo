"""Miller multi-patient generalization evaluator — CLI.

    PYTHONPATH=src .venv/bin/python -m scripts.miller_generalization_eval preflight-calibration
    PYTHONPATH=src .venv/bin/python -m scripts.miller_generalization_eval preflight-final
    PYTHONPATH=src .venv/bin/python -m scripts.miller_generalization_eval calibrate
    PYTHONPATH=src .venv/bin/python -m scripts.miller_generalization_eval finalize

See ``src/benchmark/miller_generalization_eval.py`` for the gated once-only phase contracts.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from benchmark.miller_generalization_eval import (
    default_config,
    preflight_stage,
    run_calibration,
    run_final,
)

_FAILURE_STATUSES = {
    "PREFLIGHT_FAILED", "NO_UNIVERSAL_POLICY_LOCK", "LOCK_INVALID", "LOCK_UNREADABLE",
    "CALIBRATION_UNSEAL_INCOMPLETE", "FINAL_UNSEAL_INCOMPLETE", "LABELS_INVALID",
    "CALIBRATION_NOT_EVALUABLE",
    "FINAL_RESULT_UNREADABLE",
    "SNAPSHOT_FAILED",
}


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    cmd = argv[0] if argv else "preflight-calibration"
    valid_commands = {
        "preflight-calibration", "preflight-final", "calibrate", "finalize",
    }
    if cmd not in valid_commands:
        print("usage: miller_generalization_eval.py [preflight-calibration|preflight-final|calibrate|finalize]")
        return 2
    config = default_config()
    if cmd == "preflight-calibration":
        out = preflight_stage(config, "calibration")
        ok = bool(out.get("ok"))
    elif cmd == "preflight-final":
        out = preflight_stage(config, "final")
        ok = bool(out.get("ok"))
    elif cmd == "calibrate":
        out = run_calibration(config)
        ok = out.get("status") not in _FAILURE_STATUSES
    elif cmd == "finalize":
        out = run_final(config)
        ok = out.get("status") not in _FAILURE_STATUSES
    print(json.dumps(out, indent=2, default=str))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
