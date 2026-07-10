import json
from pathlib import Path

from epicurus_neo.m6.transfer_runner import run_m6b

_VERDICTS = {"ACCEPT_TRANSFER", "REJECT_TRANSFER", "CONSISTENT_WITH_NO_EFFECT_TRANSFER"}
_STUDIES = {"braun_rcc_2025", "hu_neovax_2021", "mkras_vax_2026", "pdac_neovax_2023"}


def test_run_m6b_writes_audit_with_declared_verdict(tmp_path):
    audit = run_m6b(out_dir=tmp_path, seed=17, bootstrap_n=2000)
    assert audit["corpus_verdict"] == "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"
    assert audit["transfer"]["verdict"] in _VERDICTS
    assert set(audit["transfer"]["per_fold"]) == _STUDIES
    assert (Path(tmp_path) / "m6b_audit.json").exists()
    assert (Path(tmp_path) / "m6b_audit.md").exists()
    md = (Path(tmp_path) / "m6b_audit.md").read_text()
    assert audit["transfer"]["verdict"] in md
    assert "auxiliary" in md.lower()


def test_run_m6b_is_deterministic(tmp_path):
    first = run_m6b(out_dir=tmp_path / "a", seed=17, bootstrap_n=2000)
    second = run_m6b(out_dir=tmp_path / "b", seed=17, bootstrap_n=2000)
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
