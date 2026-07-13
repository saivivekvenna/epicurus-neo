"""Create the locked, label-independent Miller calibration/final split."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CROSSWALK = ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv/INPUT_CROSSWALK.csv"
OUT = ROOT / "artifacts/milestone_8_generalization/SPLIT.json"
SALT = "epicurus-generalization-v1|"


def build_split() -> dict:
    with CROSSWALK.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    candidates = [row for row in rows if row["patient_id"] != "Hu_287"]
    if len(candidates) != 12:
        raise ValueError(f"expected 12 non-Hu_287 Miller patients, found {len(candidates)}")
    for row in candidates:
        row["split_hash"] = hashlib.sha256(
            f"{SALT}{row['patient_id']}".encode()
        ).hexdigest()
    ordered = sorted(candidates, key=lambda row: row["split_hash"])

    def compact(row: dict) -> dict:
        return {
            "patient_id": row["patient_id"],
            "split_hash": row["split_hash"],
            "size_gb": float(row["size_gb"]),
            "complete_raw_input_crosswalk": row["complete"].lower() == "true",
        }

    return {
        "split_id": "miller-generalization-v1",
        "method": "sort sha256('epicurus-generalization-v1|' + patient_id); first 6 calibration, last 6 final",
        "label_columns_read": [],
        "known_development": ["Hu_287", "Sid"],
        "calibration": [compact(row) for row in ordered[:6]],
        "final_held_out": [compact(row) for row in ordered[6:]],
        "crosswalk_sha256": hashlib.sha256(CROSSWALK.read_bytes()).hexdigest(),
    }


def main() -> int:
    payload = build_split()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
