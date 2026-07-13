from benchmark.miller_patient import load_patient
from scripts.miller_patient_reconstruct import patient_manifest


def test_resolves_hu287_runs_without_labels():
    p = load_patient("Hu_287")
    assert (p.normal_exome_run, p.tumor_exome_run, p.tumor_rna_run) == (
        "SRR24836184", "SRR24836169", "SRR24836183"
    )
    assert p.normal_sample == "Hu_287_N"
    assert p.raw_dir.name == "hu_287"


def test_resolves_new_patient_and_declares_no_label_columns():
    p = load_patient("Hu_315")
    assert len({p.normal_exome_run, p.tumor_exome_run, p.tumor_rna_run}) == 3
    manifest = patient_manifest("Hu_315")
    assert manifest["label_columns_read"] == []
    assert manifest["runs"]["tumor_rna"] == p.tumor_rna_run
    assert "milestone_8_generalization" in manifest["artifact_dir"]
