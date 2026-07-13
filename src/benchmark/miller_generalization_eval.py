"""Label-isolated Miller multi-patient generalization evaluator (Milestone 8).

Two gated, once-only phases over the six-patient calibration/final-held-out split locked in
``artifacts/milestone_8_generalization/SPLIT.json`` (see ``PROTOCOL.md`` and
``PRECALIBRATION_ARMS.md``):

  calibrate — hard-fails BEFORE any label read unless every one of the six calibration patients has
              a valid, independently re-verified ``FROZEN_NO_LABELS`` freeze manifest. Opens the
              recognition-label CSV exactly once, computes every registered arm's per-patient metrics,
              selects the universal policy by the locked lexicographic objective (excluding the PRIME
              controls), and atomically writes a once-only ``FROZEN_NO_FINAL_LABELS`` policy lock.
  finalize  — requires the policy lock plus all six final-held-out patients' freeze manifests to
              independently re-verify, then opens the label CSV exactly once and evaluates ONLY the
              selected arm against ``prime_plain`` on the untouched final cohort, producing one of
              ``GENERALIZES`` / ``TIES_PRIME`` / ``DOES_NOT_GENERALIZE`` / ``NOT_EVALUABLE``.

Historical provenance: repo-tracked semantic inputs (anything recorded in a manifest's
``git_tracked_clean``/``code_files``) are verified against the git blob at the manifest's own
``git_commit`` — NOT the current working tree — so an older calibration freeze survives later,
unrelated code evolution. Non-repo data inputs (raw sequencing/derived patient files) are verified
against current disk, as they carry no git history of their own.

This module never reads the real recognition-label CSV outside of the two gated label-open call
sites (``_read_labels_once``), and never mutates a freeze; it only reads already-frozen artifacts.
"""

from __future__ import annotations

import hashlib
import importlib
import io
import json
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import pandas as pd

from benchmark import miller_storage_lifecycle as life
from benchmark.miller_product_freeze import ARM_IDS, POLICY_ID as PRODUCT_POLICY_ID

# Reused pure helpers/constants — see [[miller_universe_core.py]] for the same reuse pattern. Importing
# this frozen script does NOT mutate its bytes (its own provenance is untouched).
u = importlib.import_module("scripts.miller_hu287_universe")

ROOT = Path(__file__).resolve().parents[2]
SPLIT_PATH = ROOT / "artifacts/milestone_8_generalization/SPLIT.json"
SPLIT_SHA256 = "a3b344f79f92f36e97b4218673588bf7061fc4662d085618628133089a4a996f"
EXPECTED_CALIBRATION_IDS: tuple[str, ...] = (
    "Hu_182", "Hu_315", "Hu_277", "Hu_268", "Hu_254", "Hu_343",
)
EXPECTED_FINAL_IDS: tuple[str, ...] = (
    "Hu_333", "Hu_159", "Hu_344", "Hu_048", "Hu_293", "Hu_250",
)

EVAL_POLICY_ID = "miller-generalization-eval-v1"
REGISTERED_ARMS: tuple[str, ...] = ARM_IDS
EXCLUDED_CONTROL_ARMS: tuple[str, ...] = ("prime_plain", "prime_mutation_cap1")
# Fixed simplicity order (simplest first) — also the deterministic tie-break in universal policy selection.
SELECTABLE_ARMS: tuple[str, ...] = (
    "shipped_epicurus_product",
    "epicurus_plain",
    "epicurus_mutation_cap1",
    "rank_fusion_cap1",
    "evidence_lane_portfolio",
)
COMPARATOR_ARM = "prime_plain"
_FULL_SHA = re.compile(r"[0-9a-f]{40}")

CODE_FILES: tuple[Path, ...] = (ROOT / "src/benchmark/miller_generalization_eval.py",)

VERDICTS = ("GENERALIZES", "TIES_PRIME", "DOES_NOT_GENERALIZE", "NOT_EVALUABLE")

LOCK_FILENAME = "UNIVERSAL_POLICY_LOCK.json"
CALIBRATION_RESULT_FILENAME = "CALIBRATION_RESULT.json"
FINAL_RESULT_FILENAME = "FINAL_RESULT.json"
CALIBRATION_CLAIM_FILENAME = "CALIBRATION_UNSEAL_STARTED.json"
FINAL_CLAIM_FILENAME = "FINAL_UNSEAL_STARTED.json"

_EXPECTED_AVAILABILITY_KEYS = {
    "translated", "presented", "recognized", "coverage", "genuine_prime_available",
    "frozen_epicurus_available", "shipped_epicurus_score_available",
}
_ALLOWED_LABELS = {"POSITIVE", "TESTED_NEGATIVE"}
_EXPECTED_PRIME_PROVENANCE = (
    "genuine PRIME percentile rank from frozen universe prime_rank; converted once to higher-is-better "
    "recognition_score"
)
_EXPECTED_EPICURUS_ARM_SCORE = (
    "shipped epicurus_lower_evidence_score; legacy frozen_epicurus_score is diagnostic only"
)
_EXPECTED_TOOL_HEADS = {
    "PRIME": "7b18d4e11042141e7102f7c69be2b0e03d138dab",
    "MixMHCpred": "0a7f9b9e20d1cf02236f4a0a90d16735be879b38",
}


def _default_git_blob_sha256(commit: str, rel: str, root: Path) -> str | None:
    """SHA-256 of the blob recorded for ``rel`` at ``commit`` — NOT the working tree. Returns None if the
    commit lacks that path (renamed/deleted/never existed at that point)."""
    try:
        proc = subprocess.run(
            ["git", "-C", str(root), "cat-file", "-p", f"{commit}:{rel}"],
            check=True,
            capture_output=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return hashlib.sha256(proc.stdout).hexdigest()


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PatientRef:
    patient_id: str
    freeze_dir: Path


@dataclass(frozen=True)
class EvalConfig:
    root: Path
    output_dir: Path
    labels_path: Path
    calibration_patients: tuple[PatientRef, ...]
    final_patients: tuple[PatientRef, ...]
    expected_calibration_ids: tuple[str, ...] = EXPECTED_CALIBRATION_IDS
    expected_final_ids: tuple[str, ...] = EXPECTED_FINAL_IDS
    split_path: Path | None = SPLIT_PATH
    expected_split_sha256: str | None = SPLIT_SHA256
    commit_exists: Callable[[str, Path], bool] = life._commit_exists
    git_blob_sha256: Callable[[str, str, Path], str | None] = _default_git_blob_sha256


def _refs_from_split(split: dict, key: str) -> tuple[PatientRef, ...]:
    from benchmark.miller_patient import load_patient

    return tuple(
        PatientRef(
            patient_id=row["patient_id"],
            freeze_dir=load_patient(row["patient_id"]).raw_dir / "freeze",
        )
        for row in split[key]
    )


def default_config(output_dir: Path | None = None) -> EvalConfig:
    """Real-repo configuration. Reads only ``SPLIT.json`` (no labels)."""
    if life._sha256(SPLIT_PATH) != SPLIT_SHA256:
        raise RuntimeError("locked SPLIT.json hash mismatch")
    split = json.loads(SPLIT_PATH.read_text())
    return EvalConfig(
        root=ROOT,
        output_dir=output_dir or (ROOT / "artifacts/milestone_8_generalization/generalization_eval"),
        labels_path=u.LABELS,
        calibration_patients=_refs_from_split(split, "calibration"),
        final_patients=_refs_from_split(split, "final_held_out"),
    )


# ---------------------------------------------------------------------------
# Phase preflight (label-free)
# ---------------------------------------------------------------------------
def _load_manifest(freeze_dir: Path) -> tuple[dict | None, str | None]:
    manifest_path = freeze_dir / "FREEZE_MANIFEST.json"
    if manifest_path.is_symlink():
        return None, "FREEZE_MANIFEST.json is a symlink; the manifest must be a regular file in freeze_dir"
    if not manifest_path.is_file():
        return None, "no FROZEN_NO_LABELS manifest (freeze absent or incomplete)"
    try:
        return json.loads(manifest_path.read_text()), None
    except (OSError, ValueError) as exc:
        return None, f"manifest unreadable: {exc}"


def _verify_path_hash(
    root: Path,
    rel: str,
    want: str,
    *,
    historical: bool,
    commit: str,
    git_blob_sha256: Callable[[str, str, Path], str | None],
) -> str | None:
    """Return None if ``rel`` hashes to ``want``, else a failure reason. Historical inputs (repo-tracked
    semantic files) are verified against the git blob at ``commit``; everything else against current disk."""
    if not life._is_hex(want):
        return f"hash for {rel!r} is not a valid sha256"
    rel_path = Path(rel)
    if rel_path.is_absolute() or ".." in rel_path.parts:
        return f"unsafe input key (absolute/traversal): {rel!r}"
    if historical:
        got = git_blob_sha256(commit, rel, root)
        if got is None:
            return f"historical blob unreadable for {rel!r} at commit {commit}"
        if got != want:
            return f"historical blob hash mismatch: {rel!r}"
        return None
    target = life._safe_member(root, rel)
    if target is None:
        return f"unsafe input key (absolute/traversal/symlink/escape): {rel!r}"
    if not target.is_file() or life._sha256(target) != want:
        return f"input hash mismatch or file absent: {rel!r}"
    return None


def preflight_patient(config: EvalConfig, expected_patient_id: str, freeze_dir: Path) -> tuple[bool, dict]:
    """Independently re-verify one patient's freeze without opening any label file. Fail-CLOSED on the
    first problem; ``ok`` is True only if every declared output hash, every declared input hash (historical
    for repo-tracked semantic files, current-disk otherwise), the pinned code/frozen-module provenance, the
    git commit, and the product-portfolio arm/selection metadata all verify."""
    detail: dict = {"patient_id": expected_patient_id, "freeze_dir": life._rel(freeze_dir, config.root)}
    man, err = _load_manifest(freeze_dir)
    if man is None:
        detail["reason"] = err
        return False, detail

    if man.get("LOCK") != "FROZEN_NO_LABELS":
        detail["reason"] = f"LOCK marker is {man.get('LOCK')!r}, expected 'FROZEN_NO_LABELS'"
        return False, detail
    if man.get("labels_opened") is not False:
        detail["reason"] = f"labels_opened is {man.get('labels_opened')!r}, expected False"
        return False, detail
    if man.get("patient_id") != expected_patient_id:
        detail["reason"] = f"patient_id mismatch: manifest={man.get('patient_id')!r} expected={expected_patient_id!r}"
        return False, detail

    outputs = man.get("sha256")
    if not isinstance(outputs, dict) or not outputs:
        detail["reason"] = "sha256 (output hashes) missing/empty/malformed"
        return False, detail
    if not {"variants.csv", "universe.csv"}.issubset(outputs):
        detail["reason"] = "sha256 must pin variants.csv and universe.csv"
        return False, detail
    for name, want in outputs.items():
        target = life._safe_member(freeze_dir, name)
        if target is None:
            detail["reason"] = f"unsafe output key (absolute/traversal/symlink/escape): {name!r}"
            return False, detail
        if not life._is_hex(want) or not target.is_file() or life._sha256(target) != want:
            detail["reason"] = f"output hash mismatch or file absent: {name}"
            return False, detail

    code_files = man.get("code_files")
    tracked = man.get("git_tracked_clean")
    inputs = man.get("input_sha256")
    if not isinstance(code_files, list) or not code_files:
        detail["reason"] = "code_files missing/empty"
        return False, detail
    if not isinstance(tracked, dict) or not tracked:
        detail["reason"] = "git_tracked_clean provenance missing/empty/malformed"
        return False, detail
    if not isinstance(inputs, dict) or not inputs:
        detail["reason"] = "input_sha256 missing/empty/malformed"
        return False, detail
    missing_tracked = [c for c in code_files if c not in tracked]
    missing_input = [c for c in code_files if c not in inputs]
    if missing_tracked or missing_input:
        detail["reason"] = (
            f"code_files not fully pinned: absent from git_tracked_clean={missing_tracked}, "
            f"from input_sha256={missing_input}"
        )
        return False, detail
    recorded_bad = {rel: st for rel, st in tracked.items() if st != "CLEAN"}
    if recorded_bad:
        detail["reason"] = f"git_tracked_clean recorded non-CLEAN statuses: {recorded_bad}"
        return False, detail

    commit = man.get("git_commit")
    if not (isinstance(commit, str) and _FULL_SHA.fullmatch(commit)):
        detail["reason"] = f"git_commit is not a full 40-hex sha: {commit!r}"
        return False, detail
    if not config.commit_exists(commit, config.root):
        detail["reason"] = f"git_commit does not resolve to an existing commit object: {commit}"
        return False, detail

    for rel, want in inputs.items():
        reason = _verify_path_hash(
            config.root, rel, want, historical=rel in tracked, commit=commit,
            git_blob_sha256=config.git_blob_sha256,
        )
        if reason is not None:
            detail["reason"] = reason
            return False, detail

    integ = man.get("frozen_module_integrity")
    if not isinstance(integ, dict) or not integ.get("module") or not integ.get("module_sha256"):
        detail["reason"] = "frozen_module_integrity missing/incomplete (need module + module_sha256)"
        return False, detail
    mod_rel = integ["module"]
    reason = _verify_path_hash(
        config.root, mod_rel, integ["module_sha256"], historical=mod_rel in tracked, commit=commit,
        git_blob_sha256=config.git_blob_sha256,
    )
    if reason is not None:
        detail["reason"] = f"frozen_module_integrity: {reason}"
        return False, detail

    product = man.get("product_portfolios")
    if not isinstance(product, dict):
        detail["reason"] = "product_portfolios missing/malformed"
        return False, detail
    if product.get("policy_id") != PRODUCT_POLICY_ID:
        detail["reason"] = f"product_portfolios.policy_id mismatch: {product.get('policy_id')!r}"
        return False, detail
    if product.get("k") != 20:
        detail["reason"] = f"product_portfolios.k is {product.get('k')!r}, expected 20"
        return False, detail
    if product.get("labels_opened") is not False:
        detail["reason"] = "product_portfolios.labels_opened must be False"
        return False, detail
    if product.get("prime_provenance") != _EXPECTED_PRIME_PROVENANCE:
        detail["reason"] = "product_portfolios.prime_provenance mismatch"
        return False, detail
    if product.get("epicurus_arm_score") != _EXPECTED_EPICURUS_ARM_SCORE:
        detail["reason"] = "product_portfolios.epicurus_arm_score mismatch"
        return False, detail
    if product.get("shipped_product_arm") != "shipped_epicurus_product":
        detail["reason"] = "product_portfolios.shipped_product_arm mismatch"
        return False, detail
    expected_preregistered = [a for a in REGISTERED_ARMS if a != "shipped_epicurus_product"]
    if product.get("preregistered_arm_ids") != expected_preregistered:
        detail["reason"] = "product_portfolios.preregistered_arm_ids mismatch"
        return False, detail
    availability = product.get("feature_availability_rows")
    if not isinstance(availability, dict) or set(availability) != _EXPECTED_AVAILABILITY_KEYS:
        detail["reason"] = "product_portfolios.feature_availability_rows keys malformed"
        return False, detail
    if any(type(value) is not int or value < 0 for value in availability.values()):
        detail["reason"] = "product_portfolios.feature_availability_rows values must be non-negative integers"
        return False, detail
    arms = product.get("arms")
    if not isinstance(arms, dict):
        detail["reason"] = "product_portfolios.arms missing/malformed"
        return False, detail
    if set(arms) != set(REGISTERED_ARMS):
        detail["reason"] = (
            "product_portfolios.arms must exactly equal registered arms: "
            f"got={sorted(arms)} expected={sorted(REGISTERED_ARMS)}"
        )
        return False, detail
    product_hashes = product.get("sha256")
    if not isinstance(product_hashes, dict):
        detail["reason"] = "product_portfolios.sha256 missing/malformed"
        return False, detail
    selection_files: set[str] = set()
    for arm_id in REGISTERED_ARMS:
        meta = arms[arm_id]
        if not isinstance(meta, dict):
            detail["reason"] = f"arm metadata malformed: {arm_id}"
            return False, detail
        sel_file = meta.get("selection_file")
        if not sel_file or sel_file not in outputs:
            detail["reason"] = f"arm {arm_id} selection_file missing/unpinned"
            return False, detail
        selection_files.add(sel_file)
        if product_hashes.get(sel_file) != outputs.get(sel_file):
            detail["reason"] = f"arm {arm_id} product/top-level selection hash mismatch"
            return False, detail
        ordered_c = meta.get("ordered_candidate_ids")
        ordered_m = meta.get("ordered_mutation_ids")
        n_selected = meta.get("n_selected")
        if (
            not isinstance(ordered_c, list)
            or not isinstance(ordered_m, list)
            or len(ordered_c) != len(ordered_m)
            or type(n_selected) is not int
            or n_selected < 0
            or n_selected > product["k"]
            or len(ordered_c) != n_selected
        ):
            detail["reason"] = f"arm {arm_id} ordered id metadata malformed"
            return False, detail
        sel_path = life._safe_member(freeze_dir, sel_file)
        try:
            selected = pd.read_csv(sel_path)
        except (OSError, ValueError) as exc:
            detail["reason"] = f"arm {arm_id} selection CSV unreadable: {exc}"
            return False, detail
        required_selection = {"selection_rank", "candidate_id", "mutation_id"}
        if not required_selection.issubset(selected.columns):
            detail["reason"] = f"arm {arm_id} selection CSV missing required columns"
            return False, detail
        if selected[["selection_rank", "candidate_id", "mutation_id"]].isna().any().any():
            detail["reason"] = f"arm {arm_id} selection CSV contains null required values"
            return False, detail
        if len(selected) != n_selected:
            detail["reason"] = f"arm {arm_id} selection row count mismatch"
            return False, detail
        ranks = pd.to_numeric(selected["selection_rank"], errors="coerce").tolist()
        if ranks != list(range(1, n_selected + 1)):
            detail["reason"] = f"arm {arm_id} selection_rank is not contiguous 1..n"
            return False, detail
        actual_c = selected["candidate_id"].astype(str).tolist()
        actual_m = selected["mutation_id"].astype(str).tolist()
        if len(actual_c) != len(set(actual_c)):
            detail["reason"] = f"arm {arm_id} contains duplicate candidate IDs"
            return False, detail
        if actual_c != [str(v) for v in ordered_c] or actual_m != [str(v) for v in ordered_m]:
            detail["reason"] = f"arm {arm_id} ordered IDs disagree with selection CSV"
            return False, detail
        n_unique = meta.get("n_unique_mutations")
        if type(n_unique) is not int or n_unique != len(set(actual_m)):
            detail["reason"] = f"arm {arm_id} n_unique_mutations mismatch"
            return False, detail
        if type(meta.get("saturated")) is not bool or meta.get("saturated") != (
            n_selected >= product["k"]
        ):
            detail["reason"] = f"arm {arm_id} saturated flag mismatch"
            return False, detail
    if len(selection_files) != len(REGISTERED_ARMS):
        detail["reason"] = "each registered arm must have a distinct selection file"
        return False, detail
    if set(product_hashes) != selection_files:
        detail["reason"] = "product_portfolios.sha256 must exactly cover selection files"
        return False, detail

    try:
        universe = pd.read_csv(freeze_dir / "universe.csv")
    except pd.errors.EmptyDataError:
        universe = pd.DataFrame()
    except (OSError, ValueError) as exc:
        detail["reason"] = f"universe.csv unreadable: {exc}"
        return False, detail
    if len(universe):
        if not {"candidate_id", "mutation_id", "prime_rank"}.issubset(universe.columns):
            detail["reason"] = "universe.csv missing candidate_id/mutation_id/prime_rank"
            return False, detail
        if universe[["candidate_id", "mutation_id"]].isna().any().any():
            detail["reason"] = "universe.csv contains null candidate or mutation IDs"
            return False, detail
        universe_pairs = set(zip(universe["candidate_id"].astype(str), universe["mutation_id"].astype(str)))
        genuine_prime_pairs = set(zip(
            universe.loc[pd.to_numeric(universe["prime_rank"], errors="coerce").notna(), "candidate_id"].astype(str),
            universe.loc[pd.to_numeric(universe["prime_rank"], errors="coerce").notna(), "mutation_id"].astype(str),
        ))
    else:
        universe_pairs = set()
        genuine_prime_pairs = set()
    for arm_id, meta in arms.items():
        selected_pairs = set(zip(
            [str(v) for v in meta["ordered_candidate_ids"]],
            [str(v) for v in meta["ordered_mutation_ids"]],
        ))
        if not selected_pairs.issubset(universe_pairs):
            detail["reason"] = f"arm {arm_id} contains candidate/mutation pairs absent from universe.csv"
            return False, detail
        if arm_id in {"prime_plain", "prime_mutation_cap1"} and not selected_pairs.issubset(
            genuine_prime_pairs
        ):
            detail["reason"] = f"arm {arm_id} contains rows without genuine PRIME scores"
            return False, detail
    if "n_universe_rows" in man and man.get("n_universe_rows") != len(universe):
        detail["reason"] = "n_universe_rows disagrees with universe.csv"
        return False, detail
    if any(value > len(universe) for value in availability.values()):
        detail["reason"] = "feature availability count exceeds frozen universe rows"
        return False, detail
    if availability["genuine_prime_available"] != len(genuine_prime_pairs):
        detail["reason"] = "genuine_prime_available disagrees with non-null universe prime_rank rows"
        return False, detail

    tool_commits = man.get("tool_commits")
    if not isinstance(tool_commits, dict) or set(tool_commits) != set(_EXPECTED_TOOL_HEADS):
        detail["reason"] = "tool_commits missing/malformed"
        return False, detail
    for tool, expected_head in _EXPECTED_TOOL_HEADS.items():
        recorded = tool_commits.get(tool)
        if not isinstance(recorded, dict) or recorded != {
            "dir_head": expected_head,
            "adapter_constant": expected_head,
            "match": True,
            "tracked_clean": True,
        }:
            detail["reason"] = f"tool commit provenance mismatch for {tool}"
            return False, detail

    detail["manifest_sha256"] = life._sha256(freeze_dir / "FREEZE_MANIFEST.json")
    detail["git_commit"] = commit
    return True, detail


def preflight_stage(config: EvalConfig, stage: str) -> dict:
    """Preflight ALL SIX locked IDs for ``stage`` ('calibration' or 'final'). ``ok`` is True only if every
    one of the six independently re-verifies; a single failure fails the whole stage closed."""
    if stage not in {"calibration", "final"}:
        return {"ok": False, "stage": stage, "reason": "stage must be 'calibration' or 'final'"}
    patients = config.calibration_patients if stage == "calibration" else config.final_patients
    expected = config.expected_calibration_ids if stage == "calibration" else config.expected_final_ids
    actual = tuple(p.patient_id for p in patients)
    if actual != expected or len(set(actual)) != 6:
        return {
            "ok": False, "stage": stage,
            "reason": f"patient IDs/order differ from locked {stage} split: got={actual} expected={expected}",
        }
    if config.split_path is not None and config.expected_split_sha256 is not None:
        if not config.split_path.is_file() or life._sha256(config.split_path) != config.expected_split_sha256:
            return {"ok": False, "stage": stage, "reason": "locked split file hash mismatch or missing"}
    per_patient: dict[str, dict] = {}
    ok = True
    for p in patients:
        pok, detail = preflight_patient(config, p.patient_id, p.freeze_dir)
        per_patient[p.patient_id] = detail
        ok = ok and pok
    return {"ok": ok, "stage": stage, "patients": per_patient}


# ---------------------------------------------------------------------------
# Label-bound evaluation (only reachable after a stage preflight passes)
# ---------------------------------------------------------------------------
def _read_labels_once(labels_path: Path) -> pd.DataFrame:
    """The SOLE label-CSV read call site. Callers must invoke this at most once per gated phase."""
    df = pd.read_csv(labels_path)
    required = {"patient_id", "chrom", "pos", "ref", "alt", "label"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"labels file missing required columns: {sorted(missing)}")
    if df[list(required)].isna().any().any():
        raise ValueError("labels file contains nulls in required columns")
    labels = set(df["label"].astype(str))
    if not labels.issubset(_ALLOWED_LABELS):
        raise ValueError(f"labels file contains unsupported labels: {sorted(labels - _ALLOWED_LABELS)}")
    positions = pd.to_numeric(df["pos"], errors="coerce")
    if positions.isna().any() or (positions <= 0).any() or (positions % 1 != 0).any():
        raise ValueError("labels file contains invalid genomic positions")
    clean = df.copy()
    clean["patient_id"] = clean["patient_id"].astype(str)
    clean["pos"] = positions.astype(int)
    clean["label"] = clean["label"].astype(str)
    clean["_key"] = [
        u.variant_key(c, p, r, a)
        for c, p, r, a in zip(clean["chrom"], clean["pos"], clean["ref"], clean["alt"])
    ]
    conflicts = clean.groupby(["patient_id", "_key"])["label"].nunique()
    if (conflicts > 1).any():
        bad = [f"{pid}:{key}" for pid, key in conflicts[conflicts > 1].index[:5]]
        raise ValueError(f"labels file contains conflicting duplicate labels: {bad}")
    return clean.drop_duplicates(["patient_id", "_key", "label"]).drop(columns="_key").reset_index(drop=True)


def _validate_phase_label_support(labels: pd.DataFrame, patients: tuple[PatientRef, ...]) -> None:
    expected = {p.patient_id for p in patients}
    present = set(labels["patient_id"].astype(str))
    missing = sorted(expected - present)
    if missing:
        raise ValueError(f"labels file has no tested rows for phase patients: {missing}")


@dataclass(frozen=True)
class PatientLabels:
    recognized: frozenset
    tested: frozenset


@dataclass(frozen=True)
class FrozenPatientSnapshot:
    """Byte-verified, label-blind evaluation inputs retained in memory across label unseal."""

    manifest: dict
    variants: pd.DataFrame
    universe: pd.DataFrame
    selections: dict[str, pd.DataFrame]


def _read_verified_csv_snapshot(freeze_dir: Path, rel: str, expected_sha256: str) -> pd.DataFrame:
    target = life._safe_member(freeze_dir, rel)
    if target is None or not target.is_file():
        raise ValueError(f"unsafe or missing frozen CSV: {rel!r}")
    payload = target.read_bytes()
    if hashlib.sha256(payload).hexdigest() != expected_sha256:
        raise ValueError(f"frozen CSV changed after preflight: {rel!r}")
    try:
        return pd.read_csv(io.BytesIO(payload))
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def _snapshot_stage_inputs(
    config: EvalConfig,
    patients: tuple[PatientRef, ...],
    arm_ids: tuple[str, ...],
    expected_manifest_hashes: dict[str, str] | None = None,
) -> dict[str, FrozenPatientSnapshot]:
    """Snapshot all label-blind bytes after preflight and before label access, closing the TOCTOU window."""
    snapshots: dict[str, FrozenPatientSnapshot] = {}
    for patient in patients:
        manifest_path = patient.freeze_dir / "FREEZE_MANIFEST.json"
        if manifest_path.is_symlink() or not manifest_path.is_file():
            raise ValueError(f"missing/unsafe freeze manifest for {patient.patient_id}")
        manifest_bytes = manifest_path.read_bytes()
        manifest_hash = hashlib.sha256(manifest_bytes).hexdigest()
        if expected_manifest_hashes is not None and manifest_hash != expected_manifest_hashes.get(
            patient.patient_id
        ):
            raise ValueError(f"freeze manifest changed after preflight: {patient.patient_id}")
        try:
            manifest = json.loads(manifest_bytes)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"freeze manifest became unreadable: {patient.patient_id}") from exc
        if manifest.get("patient_id") != patient.patient_id:
            raise ValueError(f"freeze manifest patient changed after preflight: {patient.patient_id}")
        outputs = manifest["sha256"]
        product = manifest["product_portfolios"]
        selections = {
            arm_id: _read_verified_csv_snapshot(
                patient.freeze_dir,
                product["arms"][arm_id]["selection_file"],
                outputs[product["arms"][arm_id]["selection_file"]],
            )
            for arm_id in arm_ids
        }
        snapshots[patient.patient_id] = FrozenPatientSnapshot(
            manifest=manifest,
            variants=_read_verified_csv_snapshot(
                patient.freeze_dir, "variants.csv", outputs["variants.csv"]
            ),
            universe=_read_verified_csv_snapshot(
                patient.freeze_dir, "universe.csv", outputs["universe.csv"]
            ),
            selections=selections,
        )
    return snapshots


def _patient_labels(labels: pd.DataFrame, patient_id: str) -> PatientLabels:
    lab = labels[labels["patient_id"] == patient_id].copy()
    lab["key"] = [
        u.variant_key(c, p, r, a) for c, p, r, a in zip(lab["chrom"], lab["pos"], lab["ref"], lab["alt"])
    ]
    recognized = frozenset(lab.loc[lab["label"] == "POSITIVE", "key"])
    tested = frozenset(lab["key"])
    return PatientLabels(recognized=recognized, tested=tested)


def _reachability(
    snapshot: FrozenPatientSnapshot,
    manifest: dict,
    recognized: frozenset,
    arm_ids: tuple[str, ...],
) -> dict:
    variants = snapshot.variants
    called = set(variants["key"]) if "key" in variants else set()
    called_pass = set(variants.loc[variants["pass_filters"], "key"]) if "pass_filters" in variants else set()
    universe = snapshot.universe
    generated = set(universe["mutation_id"].astype(str)) if "mutation_id" in universe else set()
    selected: dict[str, int] = {}
    for arm_id in arm_ids:
        meta = manifest["product_portfolios"]["arms"][arm_id]
        chosen = set(meta["ordered_mutation_ids"])
        selected[arm_id] = len(recognized & chosen)
    return {
        "n_recognized_mutations": len(recognized),
        "reachability_called": len(recognized & called),
        "reachability_called_and_passed": len(recognized & called_pass),
        "reachability_generated": len(recognized & generated),
        "reachability_valid": len(recognized & generated),
        "reachability_eligible": len(recognized & generated),
        "reachability_selected_by_arm": selected,
        "stage_semantics": {
            "called": "mutation key persisted in variants.csv",
            "called_and_passed": "mutation key with pass_filters=True in variants.csv",
            "generated": "mutation represented by >=1 peptide/HLA row in frozen universe.csv",
            "valid": "same persisted lossless-universe boundary as generated; no separate full scored table was frozen",
            "eligible": "same persisted lossless-universe boundary as generated; product-gate intermediate was not frozen",
            "selected": "unique recognized mutation IDs in each frozen ordered top-k selection",
        },
    }


def _arm_hits(
    selected: pd.DataFrame,
    product: dict,
    arm_id: str,
    meta: dict,
    recognized: frozenset,
) -> dict:
    """A verified empty portfolio is a valid zero-hit pipeline output. Genuine PRIME is the exception:
    its arms are not evaluable when the frozen scorer reports zero rows with genuine PRIME evidence."""
    genuine_prime_rows = product["feature_availability_rows"]["genuine_prime_available"]
    if arm_id in {"prime_plain", "prime_mutation_cap1"} and genuine_prime_rows == 0:
        return {
            "evaluable": False, "missing": ["no_genuine_prime_rows"], "hits": None,
            "n_selected": meta["n_selected"], "n_unique_mutations": meta["n_unique_mutations"],
            "duplicate_burden": None,
        }
    mutation_ids = selected["mutation_id"].astype(str).tolist() if len(selected) else []
    unique = set(mutation_ids)
    return {
        "evaluable": True,
        "hits": len(unique & recognized),
        "n_selected": len(mutation_ids),
        "n_unique_mutations": len(unique),
        "duplicate_burden": len(mutation_ids) - len(unique),
    }


def _aggregate_arm(per_patient: dict[str, dict], recognized_by_patient: dict[str, frozenset]) -> dict:
    hits_list, p_at_least_one, recall_list, duplicate_list, missing_patients = [], [], [], [], []
    for pid, m in per_patient.items():
        if not m["evaluable"]:
            missing_patients.append(pid)
            continue
        hits_list.append(m["hits"])
        duplicate_list.append(m["duplicate_burden"])
        p_at_least_one.append(1 if m["hits"] >= 1 else 0)
        n_recognized = len(recognized_by_patient[pid])
        if n_recognized > 0:
            recall_list.append(m["hits"] / n_recognized)
    return {
        "macro_hits_at_20": (sum(hits_list) / len(hits_list)) if hits_list else None,
        "worst_patient_hits_at_20": min(hits_list) if hits_list else None,
        "p_at_least_one_hit": (sum(p_at_least_one) / len(p_at_least_one)) if p_at_least_one else None,
        "macro_recall_at_20": (sum(recall_list) / len(recall_list)) if recall_list else None,
        "mean_duplicate_slot_burden": (
            sum(duplicate_list) / len(duplicate_list) if duplicate_list else None
        ),
        "worst_patient_duplicate_slot_burden": max(duplicate_list) if duplicate_list else None,
        "n_evaluable_patients": len(hits_list),
        "all_patients_evaluable": len(hits_list) == len(per_patient),
        "n_recall_eligible_patients": len(recall_list),
        "missing_patients": missing_patients,
    }


def _helped_tied_harmed(arm_per_patient: dict[str, dict], prime_per_patient: dict[str, dict]) -> dict:
    helped = tied = harmed = 0
    excluded = []
    per_patient_delta: dict[str, int] = {}
    for pid, m in arm_per_patient.items():
        pm = prime_per_patient[pid]
        if not m["evaluable"] or not pm["evaluable"]:
            excluded.append(pid)
            continue
        per_patient_delta[pid] = m["hits"] - pm["hits"]
        if m["hits"] > pm["hits"]:
            helped += 1
        elif m["hits"] == pm["hits"]:
            tied += 1
        else:
            harmed += 1
    return {
        "helped": helped,
        "tied": tied,
        "harmed": harmed,
        "per_patient_delta_hits_at_20": per_patient_delta,
        "macro_mean_delta_hits_at_20": (
            sum(per_patient_delta.values()) / len(per_patient_delta) if per_patient_delta else None
        ),
        "excluded_missing": excluded,
    }


def evaluate_stage(
    config: EvalConfig,
    patients: tuple[PatientRef, ...],
    labels: pd.DataFrame,
    snapshots: dict[str, FrozenPatientSnapshot],
    arm_ids: tuple[str, ...] = REGISTERED_ARMS,
) -> dict:
    """Compute per-arm, per-patient hits@20 and macro diagnostics for ``patients`` against ``labels``.
    ``arm_ids`` must include ``COMPARATOR_ARM`` if helped/tied/harmed comparisons are needed."""
    manifests = {patient_id: snapshot.manifest for patient_id, snapshot in snapshots.items()}
    per_patient_labels = {p.patient_id: _patient_labels(labels, p.patient_id) for p in patients}
    reachability = {
        p.patient_id: _reachability(
            snapshots[p.patient_id],
            manifests[p.patient_id],
            per_patient_labels[p.patient_id].recognized,
            arm_ids,
        )
        for p in patients
    }
    recognized_by_patient = {pid: pl.recognized for pid, pl in per_patient_labels.items()}

    per_arm_per_patient: dict[str, dict] = {}
    for arm_id in arm_ids:
        per_patient_metrics = {}
        for p in patients:
            meta = manifests[p.patient_id]["product_portfolios"]["arms"][arm_id]
            per_patient_metrics[p.patient_id] = _arm_hits(
                snapshots[p.patient_id].selections[arm_id],
                manifests[p.patient_id]["product_portfolios"],
                arm_id,
                meta,
                per_patient_labels[p.patient_id].recognized,
            )
        per_arm_per_patient[arm_id] = per_patient_metrics

    arm_results: dict[str, dict] = {}
    prime_per_patient = per_arm_per_patient.get(COMPARATOR_ARM)
    for arm_id in arm_ids:
        result = _aggregate_arm(per_arm_per_patient[arm_id], recognized_by_patient)
        if prime_per_patient is not None:
            result["vs_prime_plain"] = _helped_tied_harmed(per_arm_per_patient[arm_id], prime_per_patient)
        arm_results[arm_id] = result

    return {"arms": arm_results, "per_arm_per_patient": per_arm_per_patient, "reachability": reachability}


# ---------------------------------------------------------------------------
# Universal policy selection (calibration only)
# ---------------------------------------------------------------------------
def select_universal_policy(arm_results: dict) -> dict:
    """Lexicographic: maximize worst-patient hits@20, then macro hits@20, then P(>=1 hit); minimize
    harmed-vs-``prime_plain``; break remaining ties by the fixed ``SELECTABLE_ARMS`` simplicity order.
    Never selects a PRIME control arm."""

    def sort_key(arm_id: str):
        r = arm_results[arm_id]
        return (
            -(r["worst_patient_hits_at_20"] if r["worst_patient_hits_at_20"] is not None else -1),
            -(r["macro_hits_at_20"] if r["macro_hits_at_20"] is not None else -1),
            -(r["p_at_least_one_hit"] if r["p_at_least_one_hit"] is not None else -1),
            r["vs_prime_plain"]["harmed"],
            SELECTABLE_ARMS.index(arm_id),
        )

    if not arm_results[COMPARATOR_ARM]["all_patients_evaluable"]:
        raise ValueError("genuine PRIME comparator is not evaluable on all six calibration patients")
    eligible = [arm for arm in SELECTABLE_ARMS if arm_results[arm]["all_patients_evaluable"]]
    if not eligible:
        raise ValueError("no selectable arm is evaluable on all six calibration patients")
    ranking = sorted(eligible, key=sort_key)
    return {"selected_arm": ranking[0], "ranking": ranking}


def _verdict(selected_metrics: dict, prime_metrics: dict, *, expected_patients: int = 6) -> str:
    if (
        selected_metrics["n_evaluable_patients"] != expected_patients
        or prime_metrics["n_evaluable_patients"] != expected_patients
        or not selected_metrics["all_patients_evaluable"]
        or not prime_metrics["all_patients_evaluable"]
    ):
        return "NOT_EVALUABLE"
    if selected_metrics["n_recall_eligible_patients"] == 0:
        return "NOT_EVALUABLE"
    sel, prime = selected_metrics["macro_hits_at_20"], prime_metrics["macro_hits_at_20"]
    if sel is None or prime is None:
        return "NOT_EVALUABLE"
    if sel > prime:
        return "GENERALIZES"
    if sel == prime:
        return "TIES_PRIME"
    return "DOES_NOT_GENERALIZE"


# ---------------------------------------------------------------------------
# Atomic writes
# ---------------------------------------------------------------------------
def _atomic_write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w") as handle:
            handle.write(json.dumps(payload, indent=2, default=str) + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        tmp.replace(path)
        dir_fd = os.open(path.parent, os.O_RDONLY)
        try:
            os.fsync(dir_fd)
        finally:
            os.close(dir_fd)
    finally:
        if tmp.exists():
            tmp.unlink()


def _exclusive_claim(path: Path, payload: dict) -> bool:
    """Create the durable pre-label-read claim. False means a prior attempt already crossed the gate."""
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o444)
    except FileExistsError:
        return False
    with os.fdopen(fd, "w") as handle:
        handle.write(json.dumps(payload, indent=2, default=str) + "\n")
        handle.flush()
        os.fsync(handle.fileno())
    dir_fd = os.open(path.parent, os.O_RDONLY)
    try:
        os.fsync(dir_fd)
    finally:
        os.close(dir_fd)
    return True


def _runtime_contract(config: EvalConfig) -> dict:
    module_path = Path(__file__).resolve()
    cli_path = ROOT / "scripts/miller_generalization_eval.py"
    split_hash = (
        life._sha256(config.split_path)
        if config.split_path is not None and config.split_path.is_file()
        else None
    )
    try:
        git_commit = subprocess.run(
            ["git", "-C", str(ROOT), "rev-parse", "HEAD"], check=True, capture_output=True, text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        git_commit = None
    return {
        "eval_policy_id": EVAL_POLICY_ID,
        "evaluator_sha256": life._sha256(module_path),
        "cli_sha256": life._sha256(cli_path) if cli_path.is_file() else None,
        "runtime_git_commit": git_commit,
        "split_sha256": split_hash,
        "expected_calibration_ids": list(config.expected_calibration_ids),
        "expected_final_ids": list(config.expected_final_ids),
        "registered_arms": list(REGISTERED_ARMS),
        "selectable_arms": list(SELECTABLE_ARMS),
        "excluded_control_arms": list(EXCLUDED_CONTROL_ARMS),
        "comparator_arm": COMPARATOR_ARM,
        "selection_objective": [
            "max_worst_patient_hits_at_20", "max_macro_hits_at_20", "max_p_at_least_one_hit",
            "min_patients_harmed_vs_prime_plain", "simplicity_order",
        ],
    }


def _lock_problem(config: EvalConfig, lock: dict, result_path: Path) -> str | None:
    if lock.get("LOCK") != "FROZEN_NO_FINAL_LABELS":
        return "missing FROZEN_NO_FINAL_LABELS marker"
    locked_contract = lock.get("contract")
    if not isinstance(locked_contract, dict):
        return "calibration lock has no runtime contract"
    locked_semantics = {k: v for k, v in locked_contract.items() if k != "runtime_git_commit"}
    current_contract = _runtime_contract(config)
    current_semantics = {k: v for k, v in current_contract.items() if k != "runtime_git_commit"}
    if locked_semantics != current_semantics:
        return "runtime evaluation contract differs from calibration lock"
    if not result_path.is_file() or life._sha256(result_path) != lock.get("calibration_result_sha256"):
        return "calibration_result hash mismatch or missing"
    try:
        result = json.loads(result_path.read_text())
    except (OSError, ValueError):
        return "calibration_result unreadable"
    if result.get("status") != "CALIBRATED" or result.get("policy_id") != EVAL_POLICY_ID:
        return "calibration_result identity/status mismatch"
    try:
        recomputed = select_universal_policy(result["evaluation"]["arms"])
    except (KeyError, TypeError, ValueError) as exc:
        return f"cannot recompute universal policy: {exc}"
    if lock.get("selected_arm") != recomputed["selected_arm"] or lock.get("ranking") != recomputed["ranking"]:
        return "selected arm/ranking does not match calibration result"
    if lock.get("selectable_arms") != list(SELECTABLE_ARMS):
        return "selectable arm registry mismatch"
    if lock.get("excluded_control_arms") != list(EXCLUDED_CONTROL_ARMS):
        return "excluded control registry mismatch"
    expected_metrics = {arm: result["evaluation"]["arms"][arm] for arm in SELECTABLE_ARMS}
    if lock.get("metrics") != expected_metrics:
        return "locked calibration metrics do not match calibration result"
    expected_patients = [p.patient_id for p in config.calibration_patients]
    if result.get("patients") != expected_patients or lock.get("calibration_patients") != expected_patients:
        return "calibration patient IDs mismatch"
    for p in config.calibration_patients:
        manifest_path = p.freeze_dir / "FREEZE_MANIFEST.json"
        if not manifest_path.is_file() or life._sha256(manifest_path) != lock.get(
            "calibration_manifest_sha256", {}
        ).get(p.patient_id):
            return f"calibration freeze manifest changed or missing: {p.patient_id}"
    return None


# ---------------------------------------------------------------------------
# Calibration (once-only)
# ---------------------------------------------------------------------------
def run_calibration(config: EvalConfig) -> dict:
    lock_path = config.output_dir / LOCK_FILENAME
    result_path = config.output_dir / CALIBRATION_RESULT_FILENAME
    claim_path = config.output_dir / CALIBRATION_CLAIM_FILENAME
    if lock_path.exists():                                 # once-only: never re-open labels
        try:
            return {"status": "ALREADY_CALIBRATED", "lock": json.loads(lock_path.read_text())}
        except (OSError, ValueError):
            return {"status": "LOCK_UNREADABLE"}
    if claim_path.exists():
        return {
            "status": "CALIBRATION_UNSEAL_INCOMPLETE",
            "reason": "a prior attempt crossed the label-unseal boundary; labels will not be reopened",
        }

    pre = preflight_stage(config, "calibration")
    if not pre["ok"]:
        return {"status": "PREFLIGHT_FAILED", "preflight": pre}

    manifest_hashes = {pid: detail["manifest_sha256"] for pid, detail in pre["patients"].items()}
    try:
        snapshots = _snapshot_stage_inputs(
            config, config.calibration_patients, REGISTERED_ARMS, manifest_hashes
        )
    except (OSError, KeyError, ValueError) as exc:
        return {"status": "SNAPSHOT_FAILED", "phase": "calibration", "reason": str(exc)}
    contract = _runtime_contract(config)
    claim = {
        "CLAIM": "CALIBRATION_UNSEAL_STARTED",
        "labels_opened": False,
        "patients": [p.patient_id for p in config.calibration_patients],
        "calibration_manifest_sha256": manifest_hashes,
        "contract": contract,
    }
    if not _exclusive_claim(claim_path, claim):
        return {"status": "CALIBRATION_UNSEAL_INCOMPLETE"}

    try:
        labels = _read_labels_once(config.labels_path)     # <-- the ONLY label read in this call
        _validate_phase_label_support(labels, config.calibration_patients)
    except (OSError, ValueError) as exc:
        return {"status": "LABELS_INVALID", "phase": "calibration", "reason": str(exc)}

    evaluation = evaluate_stage(
        config,
        config.calibration_patients,
        labels,
        arm_ids=REGISTERED_ARMS,
        snapshots=snapshots,
    )
    try:
        policy = select_universal_policy(evaluation["arms"])
    except ValueError as exc:
        calibration_result = {
            "status": "CALIBRATION_NOT_EVALUABLE", "policy_id": EVAL_POLICY_ID,
            "stage": "calibration", "patients": [p.patient_id for p in config.calibration_patients],
            "evaluation": evaluation, "labels_opened": True, "reason": str(exc),
        }
        _atomic_write_json(result_path, calibration_result)
        return {"status": "CALIBRATION_NOT_EVALUABLE", "result": calibration_result}

    calibration_result = {
        "status": "CALIBRATED",
        "policy_id": EVAL_POLICY_ID,
        "stage": "calibration",
        "patients": [p.patient_id for p in config.calibration_patients],
        "evaluation": evaluation,
        "labels_opened": True,
    }

    _atomic_write_json(result_path, calibration_result)

    lock = {
        "LOCK": "FROZEN_NO_FINAL_LABELS",
        "policy_id": EVAL_POLICY_ID,
        "selected_arm": policy["selected_arm"],
        "ranking": policy["ranking"],
        "selectable_arms": list(SELECTABLE_ARMS),
        "excluded_control_arms": list(EXCLUDED_CONTROL_ARMS),
        "metrics": {a: evaluation["arms"][a] for a in SELECTABLE_ARMS},
        "calibration_manifest_sha256": manifest_hashes,
        "calibration_patients": [p.patient_id for p in config.calibration_patients],
        "contract": contract,
        "calibration_result_sha256": life._sha256(result_path),
    }
    _atomic_write_json(lock_path, lock)
    lock_path.chmod(0o444)
    return {"status": "CALIBRATED", "result": calibration_result, "lock": lock}


# ---------------------------------------------------------------------------
# Final (once-only; requires the calibration lock)
# ---------------------------------------------------------------------------
def run_final(config: EvalConfig) -> dict:
    lock_path = config.output_dir / LOCK_FILENAME
    final_path = config.output_dir / FINAL_RESULT_FILENAME
    result_path = config.output_dir / CALIBRATION_RESULT_FILENAME
    claim_path = config.output_dir / FINAL_CLAIM_FILENAME

    if final_path.exists():                                # once-only: never re-open labels
        try:
            return {"status": "ALREADY_FINALIZED", "result": json.loads(final_path.read_text())}
        except (OSError, ValueError):
            return {"status": "FINAL_RESULT_UNREADABLE"}
    if claim_path.exists():
        return {
            "status": "FINAL_UNSEAL_INCOMPLETE",
            "reason": "a prior attempt crossed the final-label boundary; labels will not be reopened",
        }
    if not lock_path.exists():
        return {"status": "NO_UNIVERSAL_POLICY_LOCK"}
    try:
        lock = json.loads(lock_path.read_text())
    except (OSError, ValueError):
        return {"status": "LOCK_UNREADABLE"}
    problem = _lock_problem(config, lock, result_path)
    if problem is not None:
        return {"status": "LOCK_INVALID", "reason": problem}
    selected_arm = lock["selected_arm"]
    arm_ids = (selected_arm, COMPARATOR_ARM)

    pre = preflight_stage(config, "final")
    if not pre["ok"]:
        return {"status": "PREFLIGHT_FAILED", "preflight": pre}

    manifest_hashes = {pid: detail["manifest_sha256"] for pid, detail in pre["patients"].items()}
    try:
        snapshots = _snapshot_stage_inputs(config, config.final_patients, arm_ids, manifest_hashes)
    except (OSError, KeyError, ValueError) as exc:
        return {"status": "SNAPSHOT_FAILED", "phase": "final", "reason": str(exc)}
    claim = {
        "CLAIM": "FINAL_UNSEAL_STARTED",
        "labels_opened": False,
        "patients": [p.patient_id for p in config.final_patients],
        "final_manifest_sha256": manifest_hashes,
        "universal_policy_lock_sha256": life._sha256(lock_path),
        "contract": _runtime_contract(config),
    }
    if not _exclusive_claim(claim_path, claim):
        return {"status": "FINAL_UNSEAL_INCOMPLETE"}

    try:
        labels = _read_labels_once(config.labels_path)     # <-- the ONLY label read in this call
        _validate_phase_label_support(labels, config.final_patients)
    except (OSError, ValueError) as exc:
        return {"status": "LABELS_INVALID", "phase": "final", "reason": str(exc)}

    evaluation = evaluate_stage(
        config, config.final_patients, labels, arm_ids=arm_ids, snapshots=snapshots
    )
    selected_metrics = evaluation["arms"][selected_arm]
    prime_metrics = evaluation["arms"][COMPARATOR_ARM]
    verdict = _verdict(selected_metrics, prime_metrics)

    final_result = {
        "status": "FINALIZED",
        "verdict": verdict,
        "selected_arm": selected_arm,
        "comparator_arm": COMPARATOR_ARM,
        "final_patients": [p.patient_id for p in config.final_patients],
        "evaluation": evaluation,
        "final_manifest_sha256": manifest_hashes,
        "universal_policy_lock_sha256": life._sha256(lock_path),
    }
    _atomic_write_json(final_path, final_result)
    return {"status": "FINALIZED", "result": final_result}
