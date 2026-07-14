from pathlib import Path

from scripts import miller_interim_readout as interim


def test_existing_result_is_returned_without_any_label_or_preflight_access(tmp_path, monkeypatch):
    monkeypatch.setattr(interim, "OUTPUT_ROOT", tmp_path)
    result_dir = tmp_path / "Hu_315"
    result_dir.mkdir()
    (result_dir / "INTERIM_RESULT.json").write_text(
        '{"status":"INTERIM_EVALUATED","patient_id":"Hu_315"}'
    )
    monkeypatch.setattr(
        interim,
        "_read_labels_once",
        lambda path: (_ for _ in ()).throw(AssertionError("outcomes reopened")),
    )
    assert interim.run_interim("Hu_315")["status"] == "INTERIM_EVALUATED"


def test_existing_claim_without_result_fails_closed(tmp_path, monkeypatch):
    monkeypatch.setattr(interim, "OUTPUT_ROOT", tmp_path)
    result_dir = tmp_path / "Hu_315"
    result_dir.mkdir()
    (result_dir / "INTERIM_UNSEAL_STARTED.json").write_text("{}")
    out = interim.run_interim("Hu_315")
    assert out["status"] == "INTERIM_UNSEAL_INCOMPLETE"


def test_recovery_claim_without_result_fails_closed(tmp_path, monkeypatch):
    monkeypatch.setattr(interim, "OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr(interim, "RECOVERY_ADDENDUM", Path(__file__))
    result_dir = tmp_path / "Hu_315"
    result_dir.mkdir()
    (result_dir / "INTERIM_UNSEAL_STARTED.json").write_text("{}")
    (result_dir / "INTERIM_FIRST_ATTEMPT_FAILURE.json").write_text("{}")
    (result_dir / "INTERIM_RECOVERY_STARTED.json").write_text("{}")
    out = interim.run_interim("Hu_315")
    assert out["status"] == "INTERIM_RECOVERY_INCOMPLETE"


def test_final_patient_is_refused_before_outcome_access(tmp_path, monkeypatch):
    monkeypatch.setattr(interim, "OUTPUT_ROOT", tmp_path)
    monkeypatch.setattr(interim, "ADDENDUM", Path(__file__))

    class Config:
        calibration_patients = ()

    monkeypatch.setattr(interim, "default_config", lambda output_dir: Config())
    monkeypatch.setattr(
        interim,
        "_read_labels_once",
        lambda path: (_ for _ in ()).throw(AssertionError("outcomes opened for final patient")),
    )
    out = interim.run_interim("Hu_333")
    assert out["status"] == "REFUSED_NOT_CALIBRATION"
