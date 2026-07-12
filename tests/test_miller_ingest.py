"""Tests for the Miller IPV ingestion scaffold (PRJNA980652 inputs + S1/S2 label contract).

The raw WES/RNA inputs are public (PRJNA980652); the per-peptide ELISpot label table (data files
S1/S2, DOI 10.1126/scitranslmed.abj9905) is paywalled and not yet obtained. So this scaffold covers:
  * parsing the PUBLIC SRA run metadata into a patient input crosswalk + download tranches, and
  * the ingestion CONTRACT the S1/S2 label frame must satisfy once obtained (RUNNABLE BUT BLOCKED ON FILE),
enforcing the north-star invariants: three-state labels, no collapsing of contradictory longitudinal
rows, and a deterministic patient crosswalk to the SRA BioSamples.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from benchmark.miller_ingest import (
    SRA_RUNINFO_FIXTURE,
    build_download_tranches,
    parse_sra_runinfo,
    patient_input_crosswalk,
    validate_recognition_labels,
)

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "miller_ipv_sra_runinfo.csv"


# ---------------------------------------------------------------------------
# Public SRA metadata -> input crosswalk
# ---------------------------------------------------------------------------
def test_parse_sra_runinfo_splits_wes_and_rna():
    runs = parse_sra_runinfo(FIXTURE)
    assert len(runs) == 39
    strat = runs["library_strategy"].value_counts().to_dict()
    assert strat["WXS"] == 26
    assert strat["RNA-Seq"] == 13
    # patient id derived from SampleName Hu_<NNN>_<tumor|normal>
    assert runs["patient_id"].str.match(r"Hu_\d+").all()


def test_patient_crosswalk_is_13_complete_trios():
    runs = parse_sra_runinfo(FIXTURE)
    xwalk = patient_input_crosswalk(runs)
    assert len(xwalk) == 13
    # every patient has a normal exome, a tumor exome, and a tumor RNA run
    assert xwalk["has_normal_exome"].all()
    assert xwalk["has_tumor_exome"].all()
    assert xwalk["has_tumor_rna"].all()
    assert xwalk["complete"].all()


def test_download_tranches_smallest_first_is_one_patient_trio():
    runs = parse_sra_runinfo(FIXTURE)
    tranches = build_download_tranches(runs)
    first = tranches[0]
    # smallest scientifically valid first tranche = ONE patient's full trio
    assert first["n_patients"] == 1
    assert set(first["library_strategies"]) == {"WXS", "RNA-Seq"}
    assert first["n_runs"] == 3  # normal exome + tumor exome + tumor RNA
    # tranche sizes are real, positive, and sum to the full cohort
    assert first["size_gb"] > 0
    total = sum(t["size_gb"] for t in tranches)
    assert abs(total - parse_sra_runinfo(FIXTURE)["size_gb"].sum()) < 1e-6


def test_fixture_path_constant_points_at_committed_metadata():
    assert SRA_RUNINFO_FIXTURE.name == "SRA_RUNINFO.csv"


# ---------------------------------------------------------------------------
# S1/S2 label ingestion contract
# ---------------------------------------------------------------------------
def _labels(rows):
    return pd.DataFrame(rows)


VALID_ROW = {
    "patient_id": "Hu_048", "mutant_peptide": "SLLQHLIGLSNLTHV",
    "hla_allele": "HLA-A*02:01", "assay": "IFNg_ELISpot", "assay_timepoint": "baseline",
    "label": "POSITIVE",
}


def test_valid_label_frame_passes_and_counts_three_states():
    frame = _labels([
        {**VALID_ROW, "label": "POSITIVE"},
        {**VALID_ROW, "mutant_peptide": "AAAAAAAAAAAAAAA", "label": "TESTED_NEGATIVE"},
        {**VALID_ROW, "mutant_peptide": "CCCCCCCCCCCCCCC", "label": "UNTESTED"},
    ])
    rep = validate_recognition_labels(frame, sra_patients={"Hu_048"})
    assert rep["ok"] is True
    assert rep["label_counts"] == {"POSITIVE": 1, "TESTED_NEGATIVE": 1, "UNTESTED": 1}


def test_invalid_label_vocabulary_is_rejected():
    frame = _labels([{**VALID_ROW, "label": "reactive"}])
    rep = validate_recognition_labels(frame, sra_patients={"Hu_048"})
    assert rep["ok"] is False
    assert "REACTIVE" in rep["invalid_labels"]  # normalized to uppercase before vocab check


def test_contradictory_longitudinal_rows_are_preserved_not_collapsed():
    # same peptide, DIFFERENT timepoint, different outcome -> two legitimate rows, NOT a conflict
    frame = _labels([
        {**VALID_ROW, "assay_timepoint": "baseline", "label": "TESTED_NEGATIVE"},
        {**VALID_ROW, "assay_timepoint": "post_treatment", "label": "POSITIVE"},
    ])
    rep = validate_recognition_labels(frame, sra_patients={"Hu_048"})
    assert rep["ok"] is True
    assert rep["n_rows"] == 2  # nothing collapsed
    assert rep["n_conflicting_keys"] == 0


def test_same_key_conflicting_label_is_flagged_not_dropped():
    # identical (patient, peptide, hla, assay, timepoint) with different label = a real conflict to SURFACE
    frame = _labels([
        {**VALID_ROW, "label": "POSITIVE"},
        {**VALID_ROW, "label": "TESTED_NEGATIVE"},
    ])
    rep = validate_recognition_labels(frame, sra_patients={"Hu_048"})
    assert rep["n_conflicting_keys"] == 1
    assert rep["n_rows"] == 2  # both preserved for human adjudication


def test_patient_not_in_sra_crosswalk_is_reported():
    frame = _labels([{**VALID_ROW, "patient_id": "Hu_999"}])
    rep = validate_recognition_labels(frame, sra_patients={"Hu_048"})
    assert "Hu_999" in rep["patients_without_inputs"]


def test_invalid_peptide_is_flagged():
    frame = _labels([{**VALID_ROW, "mutant_peptide": "SLLQ1LIGL"}])  # digit is not an amino acid
    rep = validate_recognition_labels(frame, sra_patients={"Hu_048"})
    assert rep["ok"] is False
    assert rep["n_invalid_peptides"] == 1
