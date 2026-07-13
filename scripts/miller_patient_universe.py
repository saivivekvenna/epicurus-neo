"""Generic label-blind candidate-universe FREEZE driver for an arbitrary Miller IPV patient.

    PYTHONPATH=src python -m scripts.miller_patient_universe Hu_315 freeze

This is the going-forward twin of ``scripts.miller_hu287_universe`` for the calibration/final Miller patients
(the label-blind generalization split in artifacts/milestone_8_generalization/SPLIT.json). It resolves every
per-patient path from public SRA metadata (never a recognition label), then runs the FREEZE in
``benchmark.miller_universe_core``. Hu_287 is intentionally REFUSED — its dedicated frozen script owns its
provenance and must not be re-frozen through this lane. There is deliberately no ``unseal`` command here.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from benchmark.miller_patient import load_patient
from benchmark.miller_universe_core import UniverseConfig, freeze


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("patient_id")
    parser.add_argument("command", choices=("freeze",), default="freeze", nargs="?")
    args = parser.parse_args(argv)

    if args.patient_id == "Hu_287":
        print("REFUSED: Hu_287 is owned by scripts/miller_hu287_universe.py (frozen provenance); "
              "use that script, not the generic lane.")
        return 1

    config = UniverseConfig.for_patient(load_patient(args.patient_id))
    m = freeze(config)
    print("FREEZE:", m.get("LOCK", m.get("status")), "| patient", m.get("patient_id", args.patient_id),
          "| universe rows", m.get("n_universe_rows"), "| pvac", m.get("genuine_pvac_lane"), "| arms",
          {a: (v.get("n_selected"), v.get("evaluable")) for a, v in (m.get("arms") or {}).items()})
    if m.get("LOCK") != "FROZEN_NO_LABELS" and m.get("status") != "ALREADY_FROZEN":
        print(json.dumps({k: v for k, v in m.items() if k != "arms"}, indent=2, default=str))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
