import json
from pathlib import Path

from epicurus_neo.m6.runner import run


def test_runner_writes_audit_with_standing_verdict(tmp_path):
    audit = run(out_dir=tmp_path, seed=17, bootstrap_n=2000)
    assert audit["corpus_verdict"] == "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"
    assert audit["universal"]["verdict"] in {"ACCEPT", "CONSISTENT_WITH_NO_EFFECT", "REJECT"}
    assert (Path(tmp_path) / "m6a_audit.json").exists()
    assert (Path(tmp_path) / "m6a_audit.md").exists()
    # Universal track always runs on all 4 studies.
    assert set(audit["universal"]["per_fold"]) == {
        "braun_rcc_2025",
        "hu_neovax_2021",
        "mkras_vax_2026",
        "pdac_neovax_2023",
    }


def test_runner_is_deterministic(tmp_path):
    first = run(out_dir=tmp_path / "a", seed=17, bootstrap_n=2000)
    second = run(out_dir=tmp_path / "b", seed=17, bootstrap_n=2000)
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
