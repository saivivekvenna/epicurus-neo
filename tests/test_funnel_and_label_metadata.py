import json

import pandas as pd
import pytest

from benchmark.funnel import (
    FunnelStage,
    ReachabilityStatus,
    annotate_reachability,
    candidate_reachability_funnel,
    validate_reachability_ledger,
    wilson_ci,
)
from benchmark.label_metadata import Assay, EventType, Timepoint, validate_label_metadata
from benchmark.labels import Label
from epicurus_neo.cli import build_parser, cmd_funnel_report


def _ledger() -> pd.DataFrame:
    stages = [stage.value for stage in FunnelStage]
    rows = []
    patterns = [
        ["reached"] * 8,
        ["lost"] + ["not_assessed"] * 7,
        ["reached", "lost"] + ["not_assessed"] * 6,
        ["reached", "reached", "reached"] + ["not_assessed"] * 5,
    ]
    for index, pattern in enumerate(patterns):
        row = {"positive_id": f"v{index}", "patient_id": f"p{index // 2}"}
        row.update(dict(zip(stages, pattern, strict=True)))
        rows.append(row)
    return pd.DataFrame(rows)


def test_reachability_funnel_accounts_for_loss_and_missing_evidence():
    report = candidate_reachability_funnel(_ledger())
    stages = {row["stage"]: row for row in report["stages"]}
    assert report["total_validated_positives"] == 4
    assert report["patients"] == 2
    assert stages["mutation_called"]["reached"] == 3
    assert stages["mutation_called"]["lost_here"] == 1
    assert stages["transcript_represented"]["reached"] == 2
    assert stages["transcript_represented"]["lost_here"] == 1
    assert stages["survives_gating"]["reached"] == 1
    assert stages["survives_gating"]["cumulative_lost"] == 2
    assert stages["survives_gating"]["not_assessed"] == 1
    assert stages["survives_gating"]["reachability_lower_bound"] == 0.25
    assert stages["survives_gating"]["reachability_upper_bound"] == 0.5
    assert stages["survives_gating"]["ci_lo"] < 0.25 < stages["survives_gating"]["ci_hi"]


def test_reachability_ledger_rejects_downstream_resurrection():
    ledger = _ledger()
    ledger.loc[1, "peptide_generated"] = "reached"
    with pytest.raises(ValueError, match="after an upstream stage was lost"):
        validate_reachability_ledger(ledger)


def test_annotate_reachability_uses_explicit_stage_identity_keys():
    positives = pd.DataFrame(
        {
            "positive_id": ["v1", "v2"],
            "patient_id": ["p1", "p1"],
            "mutation_id": ["m1", "m2"],
            "transcript_id": ["t1", "t2"],
        }
    )
    stage_tables = {
        "mutation_called": pd.DataFrame({"patient_id": ["p1", "p1"], "mutation_id": ["m1", "m2"]}),
        "transcript_represented": pd.DataFrame({"patient_id": ["p1"], "transcript_id": ["t1"]}),
    }
    stage_keys = {
        "mutation_called": ("patient_id", "mutation_id"),
        "transcript_represented": ("patient_id", "transcript_id"),
    }
    ledger = annotate_reachability(
        positives,
        stage_tables,
        stage_keys,
        {"mutation_called": True, "transcript_represented": True},
    )
    assert ledger["mutation_called"].tolist() == [
        ReachabilityStatus.REACHED,
        ReachabilityStatus.REACHED,
    ]
    assert ledger["transcript_represented"].tolist() == [
        ReachabilityStatus.REACHED,
        ReachabilityStatus.LOST,
    ]
    assert ledger.loc[1, "peptide_generated"] is ReachabilityStatus.LOST


def test_incomplete_stage_does_not_convert_absence_to_loss():
    positives = pd.DataFrame(
        {
            "positive_id": ["v1", "v2"],
            "patient_id": ["p1", "p1"],
            "mutation_id": ["m1", "m2"],
        }
    )
    ledger = annotate_reachability(
        positives,
        {"mutation_called": pd.DataFrame({"patient_id": ["p1"], "mutation_id": ["m1"]})},
        {"mutation_called": ("patient_id", "mutation_id")},
        {"mutation_called": False},
    )
    assert ledger["mutation_called"].tolist() == [
        ReachabilityStatus.REACHED,
        ReachabilityStatus.NOT_ASSESSED,
    ]


def test_wilson_ci_is_never_a_bare_point_estimate():
    lo, hi = wilson_ci(14, 14)
    assert lo < 1.0 <= hi
    assert lo == pytest.approx(0.7847, abs=0.001)


def test_label_metadata_keeps_event_assay_label_and_timepoint_distinct():
    frame = pd.DataFrame(
        {
            "record_id": ["r1", "r2", "r3"],
            "dataset_id": ["trial"] * 3,
            "event_type": [
                "pre_existing_reactivity",
                "vaccine_induced_response",
                "unknown",
            ],
            "assay": ["tetramer", "elispot", "unknown"],
            "label": ["positive", "tested_negative", "untested"],
            "timepoint": ["pre_vaccine", "post_boost", "unknown"],
            "provenance": ["doi:table:1", "doi:table:2", "doi:table:3"],
        }
    )
    validated = validate_label_metadata(frame)
    assert validated.loc[0, "event_type"] is EventType.PRE_EXISTING_REACTIVITY
    assert validated.loc[1, "assay"] is Assay.ELISPOT
    assert validated.loc[1, "label"] is Label.TESTED_NEGATIVE
    assert validated.loc[1, "timepoint"] is Timepoint.POST_BOOST


def test_label_metadata_rejects_semantic_contradictions():
    frame = pd.DataFrame(
        {
            "record_id": ["r1"],
            "dataset_id": ["trial"],
            "event_type": ["vaccine_induced_response"],
            "assay": ["elispot"],
            "label": ["positive"],
            "timepoint": ["pre_vaccine"],
            "provenance": ["doi:table:1"],
        }
    )
    with pytest.raises(ValueError, match="cannot be measured pre-vaccine"):
        validate_label_metadata(frame)


def test_funnel_report_cli_writes_confidence_intervals(tmp_path):
    source = tmp_path / "ledger.csv"
    output = tmp_path / "funnel.json"
    _ledger().to_csv(source, index=False)
    args = build_parser().parse_args(["funnel-report", str(source), "--output", str(output)])
    assert cmd_funnel_report(args) == 0
    payload = json.loads(output.read_text())
    assert payload["stages"][0]["ci_lo"] < payload["stages"][0]["recall"]
