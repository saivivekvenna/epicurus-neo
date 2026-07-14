"""Label-isolated Miller multi-patient generalization evaluator: hermetic tripwires against synthetic
fixtures only. No real cohort files, real labels, or real calibration/final CLI runs are touched."""

from __future__ import annotations

import json
import hashlib
from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest

from benchmark import miller_generalization_eval as mg
from benchmark import miller_storage_lifecycle as life
from scripts import miller_generalization_eval as cli
from benchmark.miller_generalization_eval import (
    COMPARATOR_ARM,
    EXCLUDED_CONTROL_ARMS,
    PRODUCT_POLICY_ID,
    REGISTERED_ARMS,
    SELECTABLE_ARMS,
    VERDICTS,
    EvalConfig,
    PatientRef,
    preflight_stage,
    run_calibration,
    run_final,
)

COMMIT = "a" * 40
OTHER_COMMIT = "b" * 40
TRACKED_FILES = ("code/mod.py", "code/frozen_mod.py")

M1, M2, M3 = "1:100:A:T", "1:200:A:T", "1:300:A:T"


# ---------------------------------------------------------------------------
# Fixture builders
# ---------------------------------------------------------------------------
def _blob_store(root: Path) -> tuple[dict, callable]:
    """A historical-blob table: (commit, rel) -> sha, frozen at call time (independent of later disk edits)."""
    store: dict[tuple[str, str], str] = {}

    def snapshot(commit: str, rel: str) -> None:
        store[(commit, rel)] = life._sha256(root / rel)

    def git_blob_sha256(commit: str, rel: str, root_: Path) -> str | None:
        return store.get((commit, rel))

    return store, git_blob_sha256


def _write_tracked_files(root: Path) -> None:
    (root / "code").mkdir(parents=True, exist_ok=True)
    (root / "code/mod.py").write_bytes(b"# pinned code\n")
    (root / "code/frozen_mod.py").write_bytes(b"# frozen scoring module\n")


def _write_patient_freeze(
    freeze_dir: Path,
    patient_id: str,
    *,
    root: Path,
    commit: str,
    arm_selections: dict[str, list[str]],
    non_evaluable_arms: tuple[str, ...] = (),
    extra_mutate=None,
) -> dict:
    freeze_dir.mkdir(parents=True, exist_ok=True)
    all_muts = sorted({m for muts in arm_selections.values() for m in muts})
    variants_df = pd.DataFrame({"key": all_muts, "pass_filters": [True] * len(all_muts)})
    variants_df.to_csv(freeze_dir / "variants.csv", index=False)
    universe_rows = []
    for arm_id in REGISTERED_ARMS:
        muts = [] if arm_id in non_evaluable_arms else arm_selections.get(arm_id, [])
        universe_rows.extend(
            {
                "candidate_id": f"{arm_id}-{i}",
                "mutation_id": mutation_id,
                "prime_rank": (
                    float("nan") if "prime_plain" in non_evaluable_arms else float(i + 1)
                ),
            }
            for i, mutation_id in enumerate(muts)
        )
    pd.DataFrame(universe_rows, columns=["candidate_id", "mutation_id", "prime_rank"]).to_csv(
        freeze_dir / "universe.csv", index=False
    )

    sha256_out = {
        "variants.csv": life._sha256(freeze_dir / "variants.csv"),
        "universe.csv": life._sha256(freeze_dir / "universe.csv"),
    }
    product_hashes: dict[str, str] = {}
    arms_meta: dict[str, dict] = {}
    for arm_id in REGISTERED_ARMS:
        muts = [] if arm_id in non_evaluable_arms else arm_selections.get(arm_id, [])
        rows = [
            {"selection_rank": i + 1, "candidate_id": f"{arm_id}-{i}", "mutation_id": m}
            for i, m in enumerate(muts)
        ]
        df = pd.DataFrame(rows, columns=["selection_rank", "candidate_id", "mutation_id"])
        fn = f"product_select_{arm_id}.csv"
        df.to_csv(freeze_dir / fn, index=False)
        sha256_out[fn] = life._sha256(freeze_dir / fn)
        product_hashes[fn] = sha256_out[fn]
        arms_meta[arm_id] = {
            "selection_file": fn, "n_selected": len(rows),
            "n_unique_mutations": len(set(muts)), "saturated": len(rows) >= 20,
            "ordered_candidate_ids": [r["candidate_id"] for r in rows],
            "ordered_mutation_ids": [r["mutation_id"] for r in rows],
        }

    data_rel = f"data/raw/{patient_id}.vcf"
    (root / "data/raw").mkdir(parents=True, exist_ok=True)
    (root / data_rel).write_bytes(f"vcf-for-{patient_id}".encode())

    input_sha = {rel: life._sha256(root / rel) for rel in TRACKED_FILES}
    input_sha[data_rel] = life._sha256(root / data_rel)
    tracked = {rel: "CLEAN" for rel in TRACKED_FILES}

    man = {
        "patient_id": patient_id,
        "LOCK": "FROZEN_NO_LABELS",
        "labels_opened": False,
        "n_universe_rows": len(universe_rows),
        "git_commit": commit,
        "sha256": sha256_out,
        "input_sha256": input_sha,
        "code_files": list(TRACKED_FILES),
        "git_tracked_clean": tracked,
        "frozen_module_integrity": {"module": TRACKED_FILES[1], "module_sha256": input_sha[TRACKED_FILES[1]]},
        "tool_commits": {
            "PRIME": {
                "dir_head": "7b18d4e11042141e7102f7c69be2b0e03d138dab",
                "adapter_constant": "7b18d4e11042141e7102f7c69be2b0e03d138dab",
                "match": True,
                "tracked_clean": True,
            },
            "MixMHCpred": {
                "dir_head": "0a7f9b9e20d1cf02236f4a0a90d16735be879b38",
                "adapter_constant": "0a7f9b9e20d1cf02236f4a0a90d16735be879b38",
                "match": True,
                "tracked_clean": True,
            },
        },
        "product_portfolios": {
            "policy_id": PRODUCT_POLICY_ID,
            "k": 20,
            "labels_opened": False,
            "prime_provenance": (
                "genuine PRIME percentile rank from frozen universe prime_rank; converted once to "
                "higher-is-better recognition_score"
            ),
            "epicurus_arm_score": (
                "shipped epicurus_lower_evidence_score; legacy frozen_epicurus_score is diagnostic only"
            ),
            "shipped_product_arm": "shipped_epicurus_product",
            "preregistered_arm_ids": [a for a in REGISTERED_ARMS if a != "shipped_epicurus_product"],
            "feature_availability_rows": {
                "translated": len(universe_rows), "presented": 0, "recognized": 0, "coverage": 0,
                "genuine_prime_available": (
                    0 if "prime_plain" in non_evaluable_arms else len(universe_rows)
                ),
                "frozen_epicurus_available": len(universe_rows),
                "shipped_epicurus_score_available": len(universe_rows),
            },
            "arms": arms_meta,
            "sha256": product_hashes,
        },
    }
    if extra_mutate is not None:
        extra_mutate(man, freeze_dir, root)
    (freeze_dir / "FREEZE_MANIFEST.json").write_text(json.dumps(man))
    return man


def _uniform_selections(hits_pattern: dict[str, int], recognized=(M1, M2, M3), pad=20) -> dict[str, list[str]]:
    """Build one patient's per-arm top-20 selection lists: the first ``hits`` slots are recognized keys,
    padded with unique filler mutation ids to a realistic portfolio size."""
    out = {}
    for arm_id, hits in hits_pattern.items():
        chosen = list(recognized[:hits])
        filler = [f"{arm_id}-filler-{i}" for i in range(max(0, pad - len(chosen)))]
        out[arm_id] = chosen + filler
    return out


def _make_stage(
    tmp_path: Path,
    *,
    label: str,
    n_patients: int,
    hits_by_patient: list[dict[str, int]],
    recognized: tuple[str, ...] = (M1, M2, M3),
    commit: str = COMMIT,
    selection_pad: int = 20,
    freeze_kwargs_by_index: dict[int, dict] | None = None,
) -> tuple[tuple[PatientRef, ...], dict, callable, list[str]]:
    root = tmp_path / label
    _write_tracked_files(root)
    store, git_blob_sha256 = _blob_store(root)
    for rel in TRACKED_FILES:
        store_key_snapshot(store, root, commit, rel)

    patient_ids = [f"Hu_{label}{i}" for i in range(n_patients)]
    refs = []
    for i, pid in enumerate(patient_ids):
        freeze_dir = root / "patients" / pid / "freeze"
        kwargs = (freeze_kwargs_by_index or {}).get(i, {})
        arm_selections = _uniform_selections(
            hits_by_patient[i], recognized=recognized, pad=selection_pad
        )
        _write_patient_freeze(freeze_dir, pid, root=root, commit=commit, arm_selections=arm_selections, **kwargs)
        refs.append(PatientRef(patient_id=pid, freeze_dir=freeze_dir))
    return tuple(refs), {"root": root}, git_blob_sha256, patient_ids


def store_key_snapshot(store: dict, root: Path, commit: str, rel: str) -> None:
    store[(commit, rel)] = life._sha256(root / rel)


def _labels_csv(root: Path, patient_recognized: dict[str, tuple[str, ...]]) -> Path:
    rows = []
    for pid, keys in patient_recognized.items():
        for key in keys:
            chrom, pos, ref, alt = key.split(":")
            rows.append({"patient_id": pid, "chrom": chrom, "pos": int(pos), "ref": ref, "alt": alt,
                        "label": "POSITIVE"})
        rows.append({
            "patient_id": pid, "chrom": "9", "pos": 999, "ref": "G", "alt": "C",
            "label": "TESTED_NEGATIVE",
        })
    path = root / "labels.csv"
    pd.DataFrame(rows).to_csv(path, index=False)
    return path


def _config(root: Path, calibration_refs, final_refs, labels_path: Path, git_blob_sha256) -> EvalConfig:
    return EvalConfig(
        root=root,
        output_dir=root / "eval_out",
        labels_path=labels_path,
        calibration_patients=calibration_refs,
        final_patients=final_refs,
        expected_calibration_ids=tuple(p.patient_id for p in calibration_refs),
        expected_final_ids=tuple(p.patient_id for p in final_refs),
        split_path=None,
        expected_split_sha256=None,
        commit_exists=lambda sha, r: sha == COMMIT,
        git_blob_sha256=git_blob_sha256,
    )


def _default_hits_pattern(prime_hits=1, others=1) -> dict[str, int]:
    pattern = {arm: others for arm in SELECTABLE_ARMS}
    pattern[COMPARATOR_ARM] = prime_hits
    pattern["prime_mutation_cap1"] = prime_hits
    return pattern


def _read_labels_spy(monkeypatch):
    calls = {"n": 0}
    real = mg._read_labels_once

    def spy(path):
        calls["n"] += 1
        return real(path)

    monkeypatch.setattr(mg, "_read_labels_once", spy)
    return calls


def test_repeated_assays_collapse_to_mutation_level_any_positive(tmp_path):
    labels = tmp_path / "labels.csv"
    labels.write_text(
        "patient_id,chrom,pos,ref,alt,label\n"
        "Hu_X,1,100,A,T,TESTED_NEGATIVE\n"
        "Hu_X,1,100,A,T,POSITIVE\n"
        "Hu_X,1,200,G,C,TESTED_NEGATIVE\n"
        "Hu_X,1,200,G,C,TESTED_NEGATIVE\n"
    )
    result = mg._read_labels_once(labels).sort_values("pos").reset_index(drop=True)
    assert result[["pos", "label"]].to_dict("records") == [
        {"pos": 100, "label": "POSITIVE"},
        {"pos": 200, "label": "TESTED_NEGATIVE"},
    ]


def _rehash_frozen_output(manifest: dict, freeze_dir: Path, filename: str) -> None:
    digest = life._sha256(freeze_dir / filename)
    manifest["sha256"][filename] = digest
    if filename in manifest["product_portfolios"]["sha256"]:
        manifest["product_portfolios"]["sha256"][filename] = digest


def _rewrite_as_product_v1_raw_universe(freeze_dir: Path) -> None:
    """Mirror the real product-v1 freeze: raw universe has no derived candidate_id."""
    manifest_path = freeze_dir / "FREEZE_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    universe_path = freeze_dir / "universe.csv"
    raw = pd.read_csv(universe_path)
    old_ids = raw["candidate_id"].astype(str).tolist()
    alphabet = "ACDEFGHIKLMNPQRSTVWY"

    def peptide(value: str) -> str:
        digest = hashlib.sha256(value.encode()).digest()
        return "".join(alphabet[byte % len(alphabet)] for byte in digest[:9])

    raw["patient_id"] = manifest["patient_id"]
    raw["mutant_peptide"] = [peptide(value) for value in old_ids]
    raw["hla_allele"] = "HLA-A*02:01"
    adapted = mg._universe_with_candidate_ids(raw.drop(columns="candidate_id"))
    replacements = dict(zip(old_ids, adapted["candidate_id"].astype(str)))
    raw.drop(columns="candidate_id").to_csv(universe_path, index=False)
    _rehash_frozen_output(manifest, freeze_dir, "universe.csv")

    for arm_id, meta in manifest["product_portfolios"]["arms"].items():
        selection_path = freeze_dir / meta["selection_file"]
        selected = pd.read_csv(selection_path)
        selected["candidate_id"] = selected["candidate_id"].astype(str).map(replacements)
        assert selected["candidate_id"].notna().all(), arm_id
        selected.to_csv(selection_path, index=False)
        meta["ordered_candidate_ids"] = selected["candidate_id"].astype(str).tolist()
        _rehash_frozen_output(manifest, freeze_dir, meta["selection_file"])
    manifest_path.write_text(json.dumps(manifest))


# ---------------------------------------------------------------------------
# Preflight: exact ID / hash / arm / order checks
# ---------------------------------------------------------------------------
def test_preflight_stage_passes_for_six_valid_patients(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="cal", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is True
    assert set(out["patients"]) == set(pids)
    for detail in out["patients"].values():
        assert detail["manifest_sha256"]


def test_preflight_accepts_product_v1_raw_universe_and_rederives_candidate_ids(tmp_path):
    refs, meta, blob, _ = _make_stage(
        tmp_path, label="rawv1", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    _rewrite_as_product_v1_raw_universe(refs[0].freeze_dir)
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is True, out["patients"][refs[0].patient_id]


def test_fixture_uses_exact_production_product_arm_schema(tmp_path):
    refs, meta, blob, _ = _make_stage(
        tmp_path, label="schema", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    manifest = json.loads((refs[0].freeze_dir / "FREEZE_MANIFEST.json").read_text())
    arm = manifest["product_portfolios"]["arms"]["shipped_epicurus_product"]
    assert set(arm) == {
        "selection_file", "n_selected", "n_unique_mutations", "saturated",
        "ordered_candidate_ids", "ordered_mutation_ids",
    }
    assert preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                           "calibration")["ok"] is True


def test_preflight_rejects_selection_csv_order_even_when_hashes_are_rewritten(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="csvorder", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    freeze = refs[0].freeze_dir
    man_path = freeze / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    arm = man["product_portfolios"]["arms"]["shipped_epicurus_product"]
    sel_path = freeze / arm["selection_file"]
    selected = pd.read_csv(sel_path)
    selected.loc[0, "candidate_id"] = "tampered-but-rehashed"
    selected.to_csv(sel_path, index=False)
    new_hash = life._sha256(sel_path)
    man["sha256"][arm["selection_file"]] = new_hash
    man["product_portfolios"]["sha256"][arm["selection_file"]] = new_hash
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "ordered IDs disagree" in out["patients"][pids[0]]["reason"]


def test_preflight_rejects_forged_tool_commit_provenance(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="badtool", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    man_path = refs[0].freeze_dir / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    man["tool_commits"]["PRIME"]["tracked_clean"] = False
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "tool commit provenance mismatch for PRIME" in out["patients"][pids[0]]["reason"]


def test_preflight_rejects_two_arms_aliasing_one_selection_file(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="alias", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    man_path = refs[0].freeze_dir / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    arms = man["product_portfolios"]["arms"]
    arms["epicurus_plain"] = dict(arms["shipped_epicurus_product"])
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "distinct selection file" in out["patients"][pids[0]]["reason"]


def test_preflight_rejects_duplicate_candidate_ids_inside_selection(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="dupcandidate", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    freeze = refs[0].freeze_dir
    man_path = freeze / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    arm = man["product_portfolios"]["arms"]["shipped_epicurus_product"]
    selected = pd.read_csv(freeze / arm["selection_file"])
    selected.loc[1, "candidate_id"] = selected.loc[0, "candidate_id"]
    selected.to_csv(freeze / arm["selection_file"], index=False)
    arm["ordered_candidate_ids"] = selected["candidate_id"].astype(str).tolist()
    _rehash_frozen_output(man, freeze, arm["selection_file"])
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "duplicate candidate IDs" in out["patients"][pids[0]]["reason"]


def test_preflight_rejects_null_required_selection_value(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="nullselection", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    freeze = refs[0].freeze_dir
    man_path = freeze / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    arm = man["product_portfolios"]["arms"]["shipped_epicurus_product"]
    selected = pd.read_csv(freeze / arm["selection_file"])
    selected.loc[0, "mutation_id"] = float("nan")
    selected.to_csv(freeze / arm["selection_file"], index=False)
    _rehash_frozen_output(man, freeze, arm["selection_file"])
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "contains null required values" in out["patients"][pids[0]]["reason"]


def test_preflight_rejects_selection_pair_absent_from_universe(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="notinuniverse", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    freeze = refs[0].freeze_dir
    man_path = freeze / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    arm = man["product_portfolios"]["arms"]["shipped_epicurus_product"]
    selected = pd.read_csv(freeze / arm["selection_file"])
    selected.loc[0, "candidate_id"] = "candidate-not-in-universe"
    selected.to_csv(freeze / arm["selection_file"], index=False)
    arm["ordered_candidate_ids"] = selected["candidate_id"].astype(str).tolist()
    _rehash_frozen_output(man, freeze, arm["selection_file"])
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "absent from universe.csv" in out["patients"][pids[0]]["reason"]


def test_preflight_rejects_prime_selection_without_genuine_prime_score(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="fakeprime", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    freeze = refs[0].freeze_dir
    man_path = freeze / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    prime_meta = man["product_portfolios"]["arms"]["prime_plain"]
    target_candidate = prime_meta["ordered_candidate_ids"][0]
    universe = pd.read_csv(freeze / "universe.csv")
    universe.loc[universe["candidate_id"] == target_candidate, "prime_rank"] = float("nan")
    universe.to_csv(freeze / "universe.csv", index=False)
    _rehash_frozen_output(man, freeze, "universe.csv")
    man["product_portfolios"]["feature_availability_rows"]["genuine_prime_available"] -= 1
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "without genuine PRIME scores" in out["patients"][pids[0]]["reason"]


def test_preflight_rejects_availability_count_exceeding_universe(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="badavailability", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    man_path = refs[0].freeze_dir / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    man["product_portfolios"]["feature_availability_rows"]["translated"] = (
        man["n_universe_rows"] + 1
    )
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "exceeds frozen universe rows" in out["patients"][pids[0]]["reason"]


def test_preflight_rejects_universe_without_prime_rank_column(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="noprimecolumn", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    freeze = refs[0].freeze_dir
    man_path = freeze / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    universe = pd.read_csv(freeze / "universe.csv").drop(columns="prime_rank")
    universe.to_csv(freeze / "universe.csv", index=False)
    _rehash_frozen_output(man, freeze, "universe.csv")
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "missing candidate_id/mutation_id/prime_rank" in out["patients"][pids[0]]["reason"]


def test_six_arbitrary_ids_do_not_satisfy_production_locked_split(tmp_path):
    refs, meta, blob, _ = _make_stage(
        tmp_path, label="wrongids", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    config = EvalConfig(
        root=meta["root"], output_dir=meta["root"] / "eval_out",
        labels_path=meta["root"] / "labels.csv", calibration_patients=refs, final_patients=refs,
        split_path=None, expected_split_sha256=None,
        commit_exists=lambda sha, root: True, git_blob_sha256=blob,
    )
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert "locked calibration split" in out["reason"]


def test_preflight_fails_on_missing_patient_manifest(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="miss", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    (refs[3].freeze_dir / "FREEZE_MANIFEST.json").unlink()
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert "no FROZEN_NO_LABELS manifest" in out["patients"][pids[3]]["reason"]


def test_preflight_fails_on_wrong_patient_id(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="wrongid", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    man_path = refs[0].freeze_dir / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    man["patient_id"] = "Hu_someone_else"
    man_path.write_text(json.dumps(man))
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert "patient_id mismatch" in out["patients"][pids[0]]["reason"]


def test_preflight_fails_on_missing_registered_arm(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="armgap", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    man_path = refs[1].freeze_dir / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    del man["product_portfolios"]["arms"]["evidence_lane_portfolio"]
    man_path.write_text(json.dumps(man))
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert "must exactly equal registered arms" in out["patients"][pids[1]]["reason"]


def test_preflight_fails_on_ordered_id_length_mismatch(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="order", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    man_path = refs[2].freeze_dir / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    man["product_portfolios"]["arms"]["prime_plain"]["n_selected"] = 999
    man_path.write_text(json.dumps(man))
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert "ordered id metadata malformed" in out["patients"][pids[2]]["reason"]


def test_preflight_fails_on_tampered_output_hash(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="tamperout", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    (refs[4].freeze_dir / "variants.csv").write_bytes(b"TAMPERED")
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert "output hash mismatch" in out["patients"][pids[4]]["reason"]


def test_preflight_fails_on_non_repo_input_tampered_on_disk(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="tamperin", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    (meta["root"] / f"data/raw/{pids[5]}.vcf").write_bytes(b"TAMPERED")
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert "input hash mismatch" in out["patients"][pids[5]]["reason"]


def test_preflight_requires_exactly_six_distinct_ids(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="count", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    short = refs[:5]
    config = _config(meta["root"], short, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False and "locked calibration split" in out["reason"]


# ---------------------------------------------------------------------------
# Historical provenance vs current disk
# ---------------------------------------------------------------------------
def test_current_code_drift_still_succeeds_when_historical_blob_matches(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="drift", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    # code evolves on disk AFTER the freeze; the historical blob (frozen at fixture build time) is unchanged
    (meta["root"] / "code/mod.py").write_bytes(b"# evolved code, different bytes now\n")
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is True


def test_historical_blob_mismatch_fails_even_if_disk_untouched(tmp_path):
    refs, meta, blob_store_and_fn, pids = _make_stage(
        tmp_path, label="badblob", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )

    def tampered_blob(commit, rel, root):
        if rel == "code/mod.py":
            return "f" * 64
        return blob_store_and_fn(commit, rel, root)

    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", tampered_blob)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert any("historical blob hash mismatch" in d["reason"] for d in out["patients"].values())


def test_missing_historical_blob_fails_closed(tmp_path):
    refs, meta, _, pids = _make_stage(
        tmp_path, label="noblob", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    config = _config(meta["root"], refs, refs, meta["root"] / "labels.csv", lambda c, r, root: None)
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert any("historical blob unreadable" in d["reason"] for d in out["patients"].values())


def test_historical_traversal_key_is_explicitly_rejected(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="historicaltraversal", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    man_path = refs[0].freeze_dir / "FREEZE_MANIFEST.json"
    man = json.loads(man_path.read_text())
    old = "code/mod.py"
    unsafe = "../code/mod.py"
    man["code_files"][man["code_files"].index(old)] = unsafe
    man["input_sha256"][unsafe] = man["input_sha256"].pop(old)
    man["git_tracked_clean"][unsafe] = man["git_tracked_clean"].pop(old)
    man_path.write_text(json.dumps(man))
    out = preflight_stage(_config(meta["root"], refs, refs, meta["root"] / "labels.csv", blob),
                          "calibration")
    assert out["ok"] is False
    assert "unsafe input key" in out["patients"][pids[0]]["reason"]


def test_git_commit_must_resolve_to_existing_commit(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="nocommit", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    config = EvalConfig(
        root=meta["root"], output_dir=meta["root"] / "eval_out", labels_path=meta["root"] / "labels.csv",
        calibration_patients=refs, final_patients=refs,
        expected_calibration_ids=tuple(p.patient_id for p in refs),
        expected_final_ids=tuple(p.patient_id for p in refs),
        split_path=None, expected_split_sha256=None,
        commit_exists=lambda sha, r: False, git_blob_sha256=blob,
    )
    out = preflight_stage(config, "calibration")
    assert out["ok"] is False
    assert any("does not resolve to an existing commit" in d["reason"] for d in out["patients"].values())


# ---------------------------------------------------------------------------
# Label isolation: preflight failure must prevent any label read
# ---------------------------------------------------------------------------
def test_one_missing_patient_prevents_label_read(tmp_path, monkeypatch):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="gate", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    (refs[0].freeze_dir / "FREEZE_MANIFEST.json").unlink()
    labels_path = _labels_csv(meta["root"], {pid: (M1, M2, M3) for pid in pids})
    config = _config(meta["root"], refs, refs, labels_path, blob)

    def forbidden(path):
        raise AssertionError("labels must never be read when preflight fails")

    monkeypatch.setattr(mg, "_read_labels_once", forbidden)
    out = run_calibration(config)
    assert out["status"] == "PREFLIGHT_FAILED"


def test_one_tampered_patient_prevents_label_read(tmp_path, monkeypatch):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="gate2", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    (refs[2].freeze_dir / "variants.csv").write_bytes(b"TAMPERED")
    labels_path = _labels_csv(meta["root"], {pid: (M1, M2, M3) for pid in pids})
    config = _config(meta["root"], refs, refs, labels_path, blob)

    def forbidden(path):
        raise AssertionError("labels must never be read when preflight fails")

    monkeypatch.setattr(mg, "_read_labels_once", forbidden)
    out = run_calibration(config)
    assert out["status"] == "PREFLIGHT_FAILED"


def test_snapshot_failure_after_preflight_still_prevents_label_read(tmp_path, monkeypatch):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="snapshotgate", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
    )
    labels_path = _labels_csv(meta["root"], {pid: (M1, M2, M3) for pid in pids})
    config = _config(meta["root"], refs, refs, labels_path, blob)

    def broken_snapshot(*args, **kwargs):
        raise ValueError("synthetic post-preflight mutation")

    def forbidden(path):
        raise AssertionError("labels must never be read when frozen-byte snapshot fails")

    monkeypatch.setattr(mg, "_snapshot_stage_inputs", broken_snapshot)
    monkeypatch.setattr(mg, "_read_labels_once", forbidden)
    out = run_calibration(config)
    assert out["status"] == "SNAPSHOT_FAILED"


# ---------------------------------------------------------------------------
# Calibration: metrics, lexicographic selection, once-only, exactly-once read
# ---------------------------------------------------------------------------
def _calibration_config_for_selection(tmp_path, label, hits_by_patient, monkeypatch=None):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label=label, n_patients=6, hits_by_patient=hits_by_patient,
    )
    labels_path = _labels_csv(meta["root"], {pid: (M1, M2, M3) for pid in pids})
    config = _config(meta["root"], refs, refs, labels_path, blob)
    return config, pids


def test_calibration_reads_labels_exactly_once_and_locks_atomically(tmp_path, monkeypatch):
    config, pids = _calibration_config_for_selection(
        tmp_path, "onceread", [_default_hits_pattern() for _ in range(6)],
    )
    calls = _read_labels_spy(monkeypatch)
    out = run_calibration(config)
    assert out["status"] == "CALIBRATED"
    assert calls["n"] == 1
    assert (config.output_dir / mg.LOCK_FILENAME).is_file()
    assert (config.output_dir / mg.CALIBRATION_RESULT_FILENAME).is_file()
    lock = json.loads((config.output_dir / mg.LOCK_FILENAME).read_text())
    assert lock["LOCK"] == "FROZEN_NO_FINAL_LABELS"
    assert lock["selected_arm"] in SELECTABLE_ARMS


def test_calibration_is_once_only_and_second_call_never_reads_labels(tmp_path, monkeypatch):
    config, pids = _calibration_config_for_selection(
        tmp_path, "onceonly", [_default_hits_pattern() for _ in range(6)],
    )
    calls = _read_labels_spy(monkeypatch)
    first = run_calibration(config)
    assert first["status"] == "CALIBRATED"
    second = run_calibration(config)
    assert second["status"] == "ALREADY_CALIBRATED"
    assert calls["n"] == 1               # second call did NOT reopen the labels file


def test_calibration_crash_after_label_read_is_fail_closed_on_retry(tmp_path, monkeypatch):
    config, _ = _calibration_config_for_selection(
        tmp_path, "crashsafe", [_default_hits_pattern() for _ in range(6)],
    )
    calls = _read_labels_spy(monkeypatch)

    def crash_after_read(*args, **kwargs):
        raise RuntimeError("synthetic crash after label unseal")

    monkeypatch.setattr(mg, "evaluate_stage", crash_after_read)
    with pytest.raises(RuntimeError, match="synthetic crash"):
        run_calibration(config)
    assert calls["n"] == 1
    second = run_calibration(config)
    assert second["status"] == "CALIBRATION_UNSEAL_INCOMPLETE"
    assert calls["n"] == 1


def test_calibration_crash_between_result_and_lock_never_reopens_labels(tmp_path, monkeypatch):
    config, _ = _calibration_config_for_selection(
        tmp_path, "twowritecrash", [_default_hits_pattern() for _ in range(6)],
    )
    calls = _read_labels_spy(monkeypatch)
    real_write = mg._atomic_write_json

    def fail_lock_write(path, payload):
        if path.name == mg.LOCK_FILENAME:
            raise OSError("synthetic crash before policy lock")
        return real_write(path, payload)

    monkeypatch.setattr(mg, "_atomic_write_json", fail_lock_write)
    with pytest.raises(OSError, match="synthetic crash"):
        run_calibration(config)
    assert (config.output_dir / mg.CALIBRATION_RESULT_FILENAME).is_file()
    assert not (config.output_dir / mg.LOCK_FILENAME).exists()
    second = run_calibration(config)
    assert second["status"] == "CALIBRATION_UNSEAL_INCOMPLETE"
    assert calls["n"] == 1


def test_verified_empty_epicurus_output_is_evaluable_zero_not_missing(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="emptyvalid", n_patients=6,
        hits_by_patient=[{arm: 0 for arm in REGISTERED_ARMS} for _ in range(6)],
        selection_pad=0,
    )
    labels_path = _labels_csv(meta["root"], {pid: (M1,) for pid in pids})
    config = _config(meta["root"], refs, refs, labels_path, blob)
    pre = preflight_stage(config, "calibration")
    assert pre["ok"] is True
    manifest_hashes = {pid: row["manifest_sha256"] for pid, row in pre["patients"].items()}
    snapshots = mg._snapshot_stage_inputs(config, refs, REGISTERED_ARMS, manifest_hashes)
    evaluation = mg.evaluate_stage(
        config, refs, mg._read_labels_once(labels_path), snapshots=snapshots
    )
    per_patient = evaluation["per_arm_per_patient"]["shipped_epicurus_product"]
    assert all(row["evaluable"] is True and row["hits"] == 0 for row in per_patient.values())
    assert evaluation["arms"]["shipped_epicurus_product"]["all_patients_evaluable"] is True
    assert evaluation["arms"]["prime_plain"]["all_patients_evaluable"] is False


def test_calibration_never_selects_a_prime_control_even_when_prime_dominates(tmp_path):
    hits = []
    for _ in range(6):
        pattern = {arm: 1 for arm in SELECTABLE_ARMS}
        pattern[COMPARATOR_ARM] = 3                 # PRIME strictly dominates every selectable arm
        pattern["prime_mutation_cap1"] = 3
        hits.append(pattern)
    config, pids = _calibration_config_for_selection(tmp_path, "noprime", hits)
    out = run_calibration(config)
    assert out["status"] == "CALIBRATED"
    selected = out["lock"]["selected_arm"]
    assert selected in SELECTABLE_ARMS
    assert selected not in EXCLUDED_CONTROL_ARMS


def test_calibration_refuses_partial_genuine_prime_comparator_support(tmp_path):
    refs, meta, blob, pids = _make_stage(
        tmp_path, label="partialprime", n_patients=6,
        hits_by_patient=[_default_hits_pattern() for _ in range(6)],
        freeze_kwargs_by_index={
            2: {"non_evaluable_arms": ("prime_plain", "prime_mutation_cap1")},
        },
    )
    labels_path = _labels_csv(meta["root"], {pid: (M1, M2, M3) for pid in pids})
    out = run_calibration(_config(meta["root"], refs, refs, labels_path, blob))
    assert out["status"] == "CALIBRATION_NOT_EVALUABLE"
    assert "PRIME comparator" in out["result"]["reason"]
    assert not (meta["root"] / "eval_out" / mg.LOCK_FILENAME).exists()


def test_calibration_selects_true_lexicographic_winner_not_simplicity_default(tmp_path):
    hits = []
    for _ in range(6):
        pattern = _default_hits_pattern(prime_hits=1, others=2)
        pattern["evidence_lane_portfolio"] = 3       # strictly best worst/macro hits; listed LAST in simplicity
        hits.append(pattern)
    config, pids = _calibration_config_for_selection(tmp_path, "bestwins", hits)
    out = run_calibration(config)
    assert out["lock"]["selected_arm"] == "evidence_lane_portfolio"
    assert out["lock"]["ranking"][0] == "evidence_lane_portfolio"


def test_calibration_tie_breaks_deterministically_by_simplicity_order(tmp_path):
    hits = []
    for _ in range(6):
        pattern = _default_hits_pattern(prime_hits=1, others=1)
        pattern["epicurus_mutation_cap1"] = 2
        pattern["evidence_lane_portfolio"] = 2        # exact tie with epicurus_mutation_cap1 on every metric
        hits.append(pattern)
    config, pids = _calibration_config_for_selection(tmp_path, "tiebreak", hits)
    out = run_calibration(config)
    assert out["lock"]["selected_arm"] == "epicurus_mutation_cap1"     # earlier in SELECTABLE_ARMS wins


def test_zero_positive_patient_retained_in_hits_but_excluded_from_recall(tmp_path):
    hits = [_default_hits_pattern() for _ in range(6)]
    refs, meta, blob, pids = _make_stage(tmp_path, label="zeropos", n_patients=6, hits_by_patient=hits)
    # patient 0 has NO measured positives at all (recognized set empty for that patient)
    recognized = {pid: (M1, M2, M3) for pid in pids}
    recognized[pids[0]] = ()
    labels_path = _labels_csv(meta["root"], recognized)
    config = _config(meta["root"], refs, refs, labels_path, blob)
    out = run_calibration(config)
    assert out["status"] == "CALIBRATED"
    arm_result = out["result"]["evaluation"]["arms"]["shipped_epicurus_product"]
    assert arm_result["n_evaluable_patients"] == 6          # zero-positive patient still counted (hits=0)
    assert arm_result["n_recall_eligible_patients"] == 5    # but excluded from the recall denominator
    per_patient = out["result"]["evaluation"]["per_arm_per_patient"]["shipped_epicurus_product"]
    assert per_patient[pids[0]]["hits"] == 0
    assert arm_result["mean_duplicate_slot_burden"] == 0
    assert arm_result["worst_patient_duplicate_slot_burden"] == 0
    reachability = out["result"]["evaluation"]["reachability"][pids[1]]
    assert reachability["reachability_generated"] == 1
    assert reachability["reachability_valid"] == reachability["reachability_generated"]
    assert reachability["reachability_eligible"] == reachability["reachability_generated"]
    assert "persisted lossless-universe boundary" in reachability["stage_semantics"]["valid"]


# ---------------------------------------------------------------------------
# Final: requires the lock, once-only, reads labels once, verdict vocabulary
# ---------------------------------------------------------------------------
def _calibrated_config(tmp_path, monkeypatch=None, expected_final_ids=None):
    cal_hits = []
    for _ in range(6):
        pattern = _default_hits_pattern(prime_hits=1, others=1)
        pattern["shipped_epicurus_product"] = 3
        cal_hits.append(pattern)
    cal_config, cal_pids = _calibration_config_for_selection(tmp_path, "cal", cal_hits)
    cal_config = replace(
        cal_config,
        expected_final_ids=tuple(expected_final_ids or [f"Hu_fin{i}" for i in range(6)]),
    )
    out = run_calibration(cal_config)
    assert out["status"] == "CALIBRATED"
    assert out["lock"]["selected_arm"] == "shipped_epicurus_product"
    return cal_config


def _final_refs(tmp_path, root, hits_by_patient, recognized_map, label="fin"):
    patient_ids = [f"Hu_{label}{i}" for i in range(6)]
    refs = []
    for i, pid in enumerate(patient_ids):
        freeze_dir = root / "patients" / pid / "freeze"
        arm_selections = _uniform_selections(hits_by_patient[i])
        _write_patient_freeze(freeze_dir, pid, root=root, commit=COMMIT, arm_selections=arm_selections)
        refs.append(PatientRef(patient_id=pid, freeze_dir=freeze_dir))
    labels_path = _labels_csv(root, {pid: recognized_map.get(pid, (M1, M2, M3)) for pid in patient_ids})
    return tuple(refs), patient_ids, labels_path


def _final_config_from(cal_config: EvalConfig, final_refs, labels_path: Path) -> EvalConfig:
    return EvalConfig(
        root=cal_config.root,
        output_dir=cal_config.output_dir,
        labels_path=labels_path,
        calibration_patients=cal_config.calibration_patients,
        final_patients=tuple(final_refs),
        expected_calibration_ids=cal_config.expected_calibration_ids,
        expected_final_ids=cal_config.expected_final_ids,
        split_path=None,
        expected_split_sha256=None,
        commit_exists=cal_config.commit_exists,
        git_blob_sha256=cal_config.git_blob_sha256,
    )


def test_final_refuses_before_a_lock_exists(tmp_path, monkeypatch):
    refs, meta, blob, pids = _make_stage(tmp_path, label="nolock", n_patients=6,
                                         hits_by_patient=[_default_hits_pattern() for _ in range(6)])
    labels_path = _labels_csv(meta["root"], {pid: (M1, M2, M3) for pid in pids})
    config = _config(meta["root"], refs, refs, labels_path, blob)

    def forbidden(path):
        raise AssertionError("final must never read labels without a valid lock")

    monkeypatch.setattr(mg, "_read_labels_once", forbidden)
    out = run_final(config)
    assert out["status"] == "NO_UNIVERSAL_POLICY_LOCK"


def test_final_reads_labels_exactly_once_and_is_once_only(tmp_path, monkeypatch):
    cal_config = _calibrated_config(tmp_path)
    final_refs, final_pids, labels_path = _final_refs(
        tmp_path, cal_config.root, [_default_hits_pattern(prime_hits=1, others=2) for _ in range(6)],
        recognized_map={},
    )
    final_config = _final_config_from(cal_config, final_refs, labels_path)
    calls = _read_labels_spy(monkeypatch)
    first = run_final(final_config)
    assert first["status"] == "FINALIZED"
    assert calls["n"] == 1
    assert first["result"]["verdict"] in VERDICTS
    second = run_final(final_config)
    assert second["status"] == "ALREADY_FINALIZED"
    assert calls["n"] == 1                # second call did not reopen labels


def test_final_verdict_generalizes_when_selected_beats_prime(tmp_path):
    cal_config = _calibrated_config(tmp_path)
    hits = [_default_hits_pattern(prime_hits=1, others=1) for _ in range(6)]
    for h in hits:
        h["shipped_epicurus_product"] = 3
    final_refs, final_pids, labels_path = _final_refs(tmp_path, cal_config.root, hits, recognized_map={})
    final_config = _final_config_from(cal_config, final_refs, labels_path)
    out = run_final(final_config)
    assert out["status"] == "FINALIZED"
    assert out["result"]["verdict"] == "GENERALIZES"
    paired = out["result"]["evaluation"]["arms"]["shipped_epicurus_product"]["vs_prime_plain"]
    assert paired["macro_mean_delta_hits_at_20"] == 2
    assert set(paired["per_patient_delta_hits_at_20"].values()) == {2}
    for reachability in out["result"]["evaluation"]["reachability"].values():
        assert set(reachability["reachability_selected_by_arm"]) == {
            "shipped_epicurus_product", "prime_plain",
        }


def test_final_verdict_ties_prime_when_equal(tmp_path):
    cal_config = _calibrated_config(tmp_path)
    hits = [_default_hits_pattern(prime_hits=1, others=1) for _ in range(6)]  # shipped == prime everywhere
    final_refs, final_pids, labels_path = _final_refs(tmp_path, cal_config.root, hits, recognized_map={})
    final_config = _final_config_from(cal_config, final_refs, labels_path)
    out = run_final(final_config)
    assert out["result"]["verdict"] == "TIES_PRIME"


def test_final_verdict_does_not_generalize_when_selected_loses(tmp_path):
    cal_config = _calibrated_config(tmp_path)
    hits = [_default_hits_pattern(prime_hits=3, others=1) for _ in range(6)]   # prime now wins on final
    final_refs, final_pids, labels_path = _final_refs(tmp_path, cal_config.root, hits, recognized_map={})
    final_config = _final_config_from(cal_config, final_refs, labels_path)
    out = run_final(final_config)
    assert out["result"]["verdict"] == "DOES_NOT_GENERALIZE"


def test_final_verdict_not_evaluable_when_no_positive_label_support(tmp_path):
    cal_config = _calibrated_config(tmp_path)
    hits = [_default_hits_pattern(prime_hits=1, others=1) for _ in range(6)]
    final_refs, final_pids, labels_path = _final_refs(
        tmp_path, cal_config.root, hits, recognized_map={pid: () for pid in [f"Hu_fin{i}" for i in range(6)]},
    )
    final_config = _final_config_from(cal_config, final_refs, labels_path)
    out = run_final(final_config)
    assert out["result"]["verdict"] == "NOT_EVALUABLE"


def test_final_verdict_not_evaluable_when_genuine_prime_is_unavailable(tmp_path):
    patient_ids = [f"Hu_naE{i}" for i in range(6)]
    cal_config = _calibrated_config(tmp_path, expected_final_ids=patient_ids)
    refs = []
    for i, pid in enumerate(patient_ids):
        freeze_dir = cal_config.root / "patients" / pid / "freeze"
        arm_selections = _uniform_selections(_default_hits_pattern())
        _write_patient_freeze(
            freeze_dir, pid, root=cal_config.root, commit=COMMIT, arm_selections=arm_selections,
            non_evaluable_arms=("prime_plain", "prime_mutation_cap1"),
        )
        refs.append(PatientRef(patient_id=pid, freeze_dir=freeze_dir))
    labels_path = _labels_csv(cal_config.root, {pid: (M1, M2, M3) for pid in patient_ids})
    final_config = _final_config_from(cal_config, refs, labels_path)
    out = run_final(final_config)
    assert out["result"]["verdict"] == "NOT_EVALUABLE"


def test_final_lock_invalid_when_calibration_result_tampered(tmp_path):
    cal_config = _calibrated_config(tmp_path)
    (cal_config.output_dir / mg.CALIBRATION_RESULT_FILENAME).write_text('{"tampered": true}')
    final_refs, final_pids, labels_path = _final_refs(
        tmp_path, cal_config.root, [_default_hits_pattern() for _ in range(6)], recognized_map={},
    )
    final_config = _final_config_from(cal_config, final_refs, labels_path)
    out = run_final(final_config)
    assert out["status"] == "LOCK_INVALID"


def test_final_lock_invalid_when_selected_arm_is_edited(tmp_path):
    cal_config = _calibrated_config(tmp_path)
    lock_path = cal_config.output_dir / mg.LOCK_FILENAME
    lock_path.chmod(0o644)
    lock = json.loads(lock_path.read_text())
    lock["selected_arm"] = "evidence_lane_portfolio"
    lock_path.write_text(json.dumps(lock))
    final_refs, _, labels_path = _final_refs(
        tmp_path, cal_config.root, [_default_hits_pattern() for _ in range(6)], recognized_map={},
    )
    out = run_final(_final_config_from(cal_config, final_refs, labels_path))
    assert out["status"] == "LOCK_INVALID"
    assert "does not match calibration result" in out["reason"]


def test_final_crash_after_label_read_is_fail_closed_on_retry(tmp_path, monkeypatch):
    cal_config = _calibrated_config(tmp_path)
    final_refs, _, labels_path = _final_refs(
        tmp_path, cal_config.root, [_default_hits_pattern() for _ in range(6)], recognized_map={},
    )
    final_config = _final_config_from(cal_config, final_refs, labels_path)
    calls = _read_labels_spy(monkeypatch)

    def crash_after_read(*args, **kwargs):
        raise RuntimeError("synthetic final crash")

    monkeypatch.setattr(mg, "evaluate_stage", crash_after_read)
    with pytest.raises(RuntimeError, match="synthetic final crash"):
        run_final(final_config)
    assert calls["n"] == 1
    second = run_final(final_config)
    assert second["status"] == "FINAL_UNSEAL_INCOMPLETE"
    assert calls["n"] == 1


def test_final_preflight_gates_final_patients_independently_of_calibration(tmp_path, monkeypatch):
    cal_config = _calibrated_config(tmp_path)
    final_refs, final_pids, labels_path = _final_refs(
        tmp_path, cal_config.root, [_default_hits_pattern() for _ in range(6)], recognized_map={},
    )
    (final_refs[1].freeze_dir / "FREEZE_MANIFEST.json").unlink()
    final_config = _final_config_from(cal_config, final_refs, labels_path)

    def forbidden(path):
        raise AssertionError("final must never read labels when a final-patient preflight fails")

    monkeypatch.setattr(mg, "_read_labels_once", forbidden)
    out = run_final(final_config)
    assert out["status"] == "PREFLIGHT_FAILED"


def test_verdict_vocabulary_is_exactly_the_four_values():
    assert set(VERDICTS) == {"GENERALIZES", "TIES_PRIME", "DOES_NOT_GENERALIZE", "NOT_EVALUABLE"}


# ---------------------------------------------------------------------------
# CLI dispatch and exit-code contract
# ---------------------------------------------------------------------------
def test_cli_rejects_unknown_command_before_loading_real_config(monkeypatch, capsys):
    def forbidden_config():
        raise AssertionError("unknown command must not touch split/config")

    monkeypatch.setattr(cli, "default_config", forbidden_config)
    assert cli.main(["not-a-command"]) == 2
    assert "usage:" in capsys.readouterr().out


@pytest.mark.parametrize("status", ["PREFLIGHT_FAILED", "SNAPSHOT_FAILED", "LABELS_INVALID"])
def test_cli_returns_nonzero_for_fail_closed_calibration_status(monkeypatch, status, capsys):
    sentinel = object()
    monkeypatch.setattr(cli, "default_config", lambda: sentinel)
    monkeypatch.setattr(cli, "run_calibration", lambda config: {"status": status})
    assert cli.main(["calibrate"]) == 1
    assert f'"status": "{status}"' in capsys.readouterr().out


def test_cli_returns_zero_for_successful_calibration(monkeypatch, capsys):
    sentinel = object()
    monkeypatch.setattr(cli, "default_config", lambda: sentinel)
    monkeypatch.setattr(cli, "run_calibration", lambda config: {"status": "CALIBRATED"})
    assert cli.main(["calibrate"]) == 0
    assert '"status": "CALIBRATED"' in capsys.readouterr().out
