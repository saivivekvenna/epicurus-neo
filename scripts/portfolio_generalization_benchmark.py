"""Run the locked Hu_287/Sid portfolio generalization stress test.

Selections are written to ``FROZEN_SELECTIONS.json`` before either recognition
label source is opened. Results are descriptive: both patients were previously
inspected, and Hu_287 is the discovery patient.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from benchmark.portfolio_generalization import (  # noqa: E402
    crossed_selections,
    evaluate_frozen,
    paired_deltas,
)


OUT = ROOT / "artifacts/milestone_7_decision/portfolio_generalization"
POLICY = ROOT / "configs/frozen/evidence_router_v1.json"
HU_UNIVERSE = ROOT / "data/raw/miller_ipv/hu_287/freeze/universe.csv"
HU_LABELS = ROOT / "data/raw/miller_ipv/miller_recognition_labels.csv"
SID_UNIVERSE = ROOT / "artifacts/milestone_7_decision/sid_benchmark/scored_candidates.csv.gz"
RUNNER_CODE = ROOT / "scripts/portfolio_generalization_benchmark.py"
SELECTION_CODE = ROOT / "src/benchmark/portfolio_generalization.py"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_universes() -> dict[str, dict]:
    hu = pd.read_csv(HU_UNIVERSE, low_memory=False)
    sid = pd.read_csv(SID_UNIVERSE, low_memory=False)
    return {
        "Hu_287": {
            "frame": hu,
            "prime_col": "genuine_prime",
            "epicurus_col": "epicurus",
            "role": "DISCOVERY_REPLAY_NOT_INDEPENDENT",
            "source": HU_UNIVERSE,
        },
        "Sid": {
            "frame": sid,
            "prime_col": "arm_genuine_prime",
            "epicurus_col": "arm_frozen_epicurus_v0_1",
            "role": "POST_HOC_STRESS_TEST_NOT_BLIND",
            "source": SID_UNIVERSE,
        },
    }


def _hu_positives() -> set[str]:
    labels = pd.read_csv(HU_LABELS)
    pos = labels[(labels["patient_id"] == "Hu_287") & (labels["label"] == "POSITIVE")]
    return {
        f"{str(row.chrom).removeprefix('chr')}:{int(row.pos)}:{row.ref}:{row.alt}"
        for row in pos.itertuples()
    }


def _sid_positives() -> set[str]:
    from event_b.sid_benchmark import hudson_positive_variant_ids

    return set(hudson_positive_variant_ids())


def _eligibility_audit() -> list[dict]:
    return [
        {"cohort": "Hu_287", "evaluable": True, "reason": "mutation-resolved lossless universe"},
        {"cohort": "Sid", "evaluable": True, "reason": "mutation-resolved lossless universe"},
        {"cohort": "Gartner", "evaluable": False, "reason": "scored table lacks underlying mutation identity"},
        {"cohort": "IMPROVE", "evaluable": False, "reason": "one tested peptide row is not a multi-route generated mutation universe"},
        {"cohort": "multimer", "evaluable": False, "reason": "scored table lacks underlying mutation identity"},
        {"cohort": "Zhao", "evaluable": False, "reason": "tested peptide table lacks a reconstructed multi-route candidate universe"},
        {"cohort": "CEDAR", "evaluable": False, "reason": "assay corpus lacks patient-level generated mutation universe"},
        {"cohort": "Event-B backbone", "evaluable": False, "reason": "label corpus lacks reconstructed denominators and multi-route scores"},
        {"cohort": "remaining Miller patients", "evaluable": False, "reason": "labels exist; raw WES/RNA/HLA reconstruction not yet run locally"},
    ]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    universes = _load_universes()

    # LABEL-FREE phase. Freeze every primary and sensitivity selection first.
    frozen: dict = {
        "status": "FROZEN_BEFORE_LABEL_JOIN",
        "protocol": "PROTOCOL.md",
        "policy_sha256": _sha256(POLICY),
        "code_sha256": {
            str(RUNNER_CODE.relative_to(ROOT)): _sha256(RUNNER_CODE),
            str(SELECTION_CODE.relative_to(ROOT)): _sha256(SELECTION_CODE),
        },
        "inputs": {},
        "primary_k20_cap2": {},
        "k_sensitivity_cap2": {},
        "cap_sensitivity_k20": {},
        "eligibility_audit": _eligibility_audit(),
    }
    for patient, spec in universes.items():
        frame = spec["frame"]
        frozen["inputs"][patient] = {
            "role": spec["role"],
            "path": str(spec["source"].relative_to(ROOT)),
            "sha256": _sha256(spec["source"]),
            "n_candidate_rows": int(len(frame)),
            "n_mutations": int(frame["mutation_id"].nunique()),
            "prime_col": spec["prime_col"],
            "epicurus_col": spec["epicurus_col"],
        }
        frozen["primary_k20_cap2"][patient] = crossed_selections(
            frame,
            prime_col=spec["prime_col"],
            epicurus_col=spec["epicurus_col"],
            k=20,
            max_per_mutation=2,
        )
        frozen["k_sensitivity_cap2"][patient] = {
            str(k): crossed_selections(
                frame,
                prime_col=spec["prime_col"],
                epicurus_col=spec["epicurus_col"],
                k=k,
                max_per_mutation=2,
            )
            for k in (10, 30)
        }
        frozen["cap_sensitivity_k20"][patient] = {
            ("none" if cap is None else str(cap)): crossed_selections(
                frame,
                prime_col=spec["prime_col"],
                epicurus_col=spec["epicurus_col"],
                k=20,
                max_per_mutation=cap,
            )
            for cap in (1, 2, 3, 5, None)
        }

    freeze_path = OUT / "FROZEN_SELECTIONS.json"
    freeze_path.write_text(json.dumps(frozen, indent=2, sort_keys=True) + "\n")

    # EVALUATION-ONLY phase. Label sources are opened only after the freeze exists.
    positives = {"Hu_287": _hu_positives(), "Sid": _sid_positives()}

    def evaluate_block(block: dict) -> dict:
        return {
            patient: {
                "arms": {
                    arm: evaluate_frozen(selection, positives[patient])
                    for arm, selection in arms.items()
                },
            }
            for patient, arms in block.items()
        }

    primary = evaluate_block(frozen["primary_k20_cap2"])
    for patient in primary:
        primary[patient]["paired_deltas"] = paired_deltas(primary[patient]["arms"])

    k_sensitivity = {
        patient: {
            k: {
                "arms": {
                    arm: evaluate_frozen(selection, positives[patient])
                    for arm, selection in arms.items()
                }
            }
            for k, arms in by_k.items()
        }
        for patient, by_k in frozen["k_sensitivity_cap2"].items()
    }
    cap_sensitivity = {
        patient: {
            cap: {
                "arms": {
                    arm: evaluate_frozen(selection, positives[patient])
                    for arm, selection in arms.items()
                }
            }
            for cap, arms in by_cap.items()
        }
        for patient, by_cap in frozen["cap_sensitivity_k20"].items()
    }

    aggregate = {}
    for arm in next(iter(primary.values()))["arms"]:
        aggregate[arm] = {
            "total_hits": sum(primary[p]["arms"][arm]["hits_at_k_unique_mutations"] for p in primary),
            "total_recognized": sum(len(positives[p]) for p in primary),
            "patients_helped_vs_plain_same_score": None,
        }
    aggregate["prime_route_aware"]["patients_helped_vs_plain_same_score"] = sum(
        primary[p]["paired_deltas"]["selector_delta_on_prime"] > 0 for p in primary
    )
    aggregate["epicurus_route_aware"]["patients_helped_vs_plain_same_score"] = sum(
        primary[p]["paired_deltas"]["selector_delta_on_epicurus"] > 0 for p in primary
    )

    result = {
        "status": "POST_HOC_TWO_PATIENT_STRESS_TEST_NOT_INDEPENDENT_VALIDATION",
        "freeze_sha256": _sha256(freeze_path),
        "labels_opened_after_freeze": True,
        "primary_k20_cap2": primary,
        "aggregate_descriptive_only": aggregate,
        "k_sensitivity_cap2": k_sensitivity,
        "cap_sensitivity_k20": cap_sensitivity,
        "eligibility_audit": frozen["eligibility_audit"],
    }
    result_path = OUT / "RESULT.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    (OUT / "REPORT.md").write_text(_markdown(result))
    print(json.dumps({p: primary[p]["paired_deltas"] for p in primary}, indent=2))
    return 0


def _markdown(result: dict) -> str:
    lines = [
        "# Frozen portfolio generalization stress test",
        "",
        "> Two previously inspected patients; descriptive and post-hoc. Hu_287 is the discovery replay, Sid is a stress test. This is not independent validation.",
        "",
        "## Primary result — k=20, frozen cap=2",
        "",
        "| Patient | PRIME plain | PRIME + selector | Epicurus plain | Epicurus + selector | Selector delta on PRIME | Selector delta on Epicurus |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for patient, row in result["primary_k20_cap2"].items():
        a, d = row["arms"], row["paired_deltas"]
        lines.append(
            f"| {patient} | {a['prime_plain']['hits_at_k_unique_mutations']}/3 | "
            f"{a['prime_route_aware']['hits_at_k_unique_mutations']}/3 | "
            f"{a['epicurus_plain']['hits_at_k_unique_mutations']}/3 | "
            f"{a['epicurus_route_aware']['hits_at_k_unique_mutations']}/3 | "
            f"{d['selector_delta_on_prime']:+d} | {d['selector_delta_on_epicurus']:+d} |"
        )
    lines += [
        "",
        "## Slot-use diagnostics",
        "",
        "| Patient | Arm | Slots | Unique mutations | Duplicate burden |",
        "|---|---|---:|---:|---:|",
    ]
    for patient, row in result["primary_k20_cap2"].items():
        for arm, value in row["arms"].items():
            if arm not in {"prime_plain", "prime_route_aware", "epicurus_plain", "epicurus_route_aware"}:
                continue
            lines.append(
                f"| {patient} | `{arm}` | {value['n_selected']} | "
                f"{value['n_unique_selected_mutations']} | {value['duplicate_slot_burden']} |"
            )
    lines += [
        "",
        "## Interpretation guardrail",
        "",
        "The crossed control is decisive for attribution: if PRIME gains similarly from the selector, the supported mechanism is scorer-agnostic portfolio diversification. Sid determines whether the Hu_287 mechanism survives a second mutation-resolved patient; a tie or regression must be reported as such.",
        "",
        "Sensitivity results and the full eligibility audit are in `RESULT.json`. Caps other than 2 and k values other than 20 are post-hoc diagnostics, not alternative headline endpoints.",
    ]
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
