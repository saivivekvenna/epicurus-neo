"""CLI for the fail-closed post-freeze Miller storage lifecycle.

Default is a DRY-RUN plan that prints a JSON reclamation report to stdout — it never mutates anything.
Destructive reclamation requires BOTH ``--execute`` and a typed ``--confirm-patient`` that exactly matches
the positional patient id, AND an independently re-verified FROZEN_NO_LABELS freeze. The recognition-label
table is never opened; candidate discovery/deletion is confined to the patient's raw dir.

Examples:
  # dry-run plan to stdout
  python scripts/miller_storage_cleanup.py Hu_315
  # dry-run plan, also written atomically to a report file OUTSIDE the patient dir
  python scripts/miller_storage_cleanup.py Hu_315 --report artifacts/.../Hu_315_cleanup.json
  # destructive (only if the freeze verifies): both flags required, must match the patient id
  python scripts/miller_storage_cleanup.py Hu_315 --execute --confirm-patient Hu_315
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from benchmark.miller_storage_lifecycle import PatientPaths, execute_cleanup, plan_cleanup


def _atomic_write_json(path: Path, obj: dict) -> None:
    """Write JSON atomically (temp file in the same dir + os.replace) so a partial file is never observed."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=str(path.parent), prefix=".cleanup_report_", suffix=".json")
    try:
        with os.fdopen(fd, "w") as fh:
            json.dump(obj, fh, indent=2, default=str)
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.remove(tmp)
        raise


def _reject_report_inside_raw(report_path: Path, raw_dir: Path) -> str | None:
    """A report written under the patient raw dir could itself be a cleanup candidate; require it outside."""
    try:
        report_path.resolve().relative_to(raw_dir.resolve())
    except ValueError:
        return None
    return (f"--report path {report_path} is inside the patient raw dir {raw_dir}; choose a path outside the "
            "cleanup candidate set")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Fail-closed post-freeze Miller storage lifecycle")
    ap.add_argument("patient_id")
    ap.add_argument("--execute", action="store_true",
                    help="perform destructive reclamation (requires --confirm-patient and a verified freeze)")
    ap.add_argument("--confirm-patient", default=None,
                    help="typed confirmation; must exactly equal patient_id to permit deletion")
    ap.add_argument("--report", default=None, help="also write the JSON report to this path (atomic; must be "
                    "outside the patient raw dir)")
    args = ap.parse_args(argv)

    paths = PatientPaths.from_patient_id(args.patient_id)

    report_path = Path(args.report) if args.report else None
    if report_path is not None:
        err = _reject_report_inside_raw(report_path, paths.raw_dir)
        if err:
            print(json.dumps({"status": "REFUSED_BAD_REPORT_PATH", "reason": err}, indent=2))
            return 2

    if args.execute:
        if args.confirm_patient != args.patient_id:
            report = {"status": "REFUSED_CONFIRMATION_MISMATCH", "patient_id": args.patient_id,
                      "reason": "--execute requires --confirm-patient to exactly match the patient id",
                      "confirm_patient": args.confirm_patient}
            print(json.dumps(report, indent=2, default=str))
            if report_path is not None:
                _atomic_write_json(report_path, report)
            return 2
        report = execute_cleanup(paths, confirm=True)
    else:
        report = plan_cleanup(paths)

    print(json.dumps(report, indent=2, default=str))
    if report_path is not None:                 # written only after any cleanup has completed, atomically
        _atomic_write_json(report_path, report)
    ok_status = report.get("status") in (None, "EXECUTED")
    return 0 if (report.get("frozen_verified") or not args.execute) and ok_status else 1


if __name__ == "__main__":
    raise SystemExit(main())
