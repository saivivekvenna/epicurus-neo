from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from benchmark.ablation import generate_masking_sets
from benchmark.gates import prime_rule, retroactive_prime_checks
from benchmark.headroom import tie_break_canary
from benchmark.improve import regression_values
from benchmark.labels import Label, load_manifest, validate_labels
from benchmark.leakage import assert_no_candidate_overlap, assert_temporal_cutoff
from benchmark.metrics import capture_fraction, hits_at_k, p_at_least_one, precision_at_k
from benchmark.preregistration import lint_preregistration
from benchmark.scorecard import scorecard
from benchmark.stats import bootstrap_ci, mde, n_required, paired_bootstrap


def _toy() -> pd.DataFrame:
    rows = []
    for patient in ("p1", "p2", "p3"):
        for index in range(25):
            rows.append(
                {
                    "patient_id": patient,
                    "mutant_peptide": f"{patient}PEPTIDE{index:02d}",
                    "hla_allele": "HLA-A*02:01",
                    "label": int(patient != "p3" and index in {1, 7}),
                    "candidate": 100 - index if index in {1, 7} else -index,
                    "Prime": -index,
                    "constant": 1.0,
                }
            )
    return pd.DataFrame(rows)


def test_metrics_return_sorted_patient_vectors_and_nan_only_where_registered():
    frame = _toy().sample(frac=1, random_state=4)
    hits = hits_at_k(frame, score_col="candidate", k=5)
    capture = capture_fraction(frame, score_col="candidate", k=5)
    clinical = p_at_least_one(frame, score_col="candidate", k=5)
    precision = precision_at_k(frame, score_col="candidate", k=5)
    assert hits.tolist() == [2.0, 2.0, 0.0]
    assert capture[:2].tolist() == [1.0, 1.0]
    assert np.isnan(capture[2]) and np.isnan(clinical[2])
    assert precision.tolist() == [0.4, 0.4, 0.0]


def test_tie_break_is_independent_of_source_order():
    frame = _toy()
    original = hits_at_k(frame, score_col="constant", k=10)
    shuffled = hits_at_k(frame.sample(frac=1, random_state=99), score_col="constant", k=10)
    assert original.tolist() == shuffled.tolist()
    canary = tie_break_canary(frame, k=10, leaked_value=999)
    assert canary["constant_md5"] != canary["source_order"]
    with pytest.raises(AssertionError, match="Tie-break leakage"):
        tie_break_canary(frame, k=10, leaked_value=canary["constant_md5"])


def test_bootstrap_drops_nan_and_preserves_pairing():
    interval = bootstrap_ci(np.array([1.0, 2.0, np.nan]))
    assert interval.mean == 1.5 and interval.n == 2
    paired = paired_bootstrap(np.array([2.0, 3.0, np.nan]), np.array([1.0, 1.0, 0.0]))
    assert paired.delta == 1.5 and paired.n == 2
    assert paired.lo > 0


def test_mde_is_paired_and_n_required_uses_prospective_power():
    baseline = np.zeros(70)
    candidate = np.linspace(-1.0, 1.0, 70)
    result = mde(candidate, baseline)
    expected = 1.9599639845400536 * np.std(candidate, ddof=1) / np.sqrt(70)
    assert result == pytest.approx(expected)
    assert result.n == 70
    assert n_required(1.013, 0.05) == pytest.approx(3222, abs=2)


def test_scorecard_is_complete_and_computes_verdict():
    report = scorecard(_toy(), "candidate", "Prime", k=5)
    assert report["hits@5"]["n"] == 3
    assert report["capture_fraction"]["n"] == 2
    assert report["unreachable_patients"]["display"] == "1 / 3"
    assert report["candidate_recall"] is None
    assert report["candidate_list"]["rows"] == 75
    assert report["random_baseline"] > 0
    assert report["verdict"] == "CONSISTENT_WITH_NO_EFFECT"
    reachable = _toy().query("patient_id != 'p3'")
    assert prime_rule(reachable, "candidate", k=5)


def test_prime_rule_rejects_a_hand_rule_without_significant_gain():
    frame = _toy()
    frame["hand_rule"] = frame["Prime"]
    assert not prime_rule(frame, "hand_rule", k=5)
    assert retroactive_prime_checks() == {
        "llm_articulated_rule": False,
        "contact_residue_gate": False,
        "anchor_creation_dai": False,
    }


def test_three_state_labels_and_manifest_are_mandatory():
    assert validate_labels([1, 0, -1]).tolist() == [
        Label.POSITIVE,
        Label.TESTED_NEGATIVE,
        Label.UNTESTED,
    ]
    with pytest.raises(ValueError, match="two-valued"):
        validate_labels([1, 0, 1, 0])
    manifest = load_manifest(Path("src/benchmark/manifests/improve.yml"))
    assert set(manifest.label_states) == set(Label)


def test_extraction_and_temporal_leakage_guards():
    left = pd.DataFrame({"mutant_peptide": ["AAAA"], "hla_allele": ["A01"]})
    with pytest.raises(ValueError, match="Extraction leakage"):
        assert_no_candidate_overlap(left, left.copy())
    temporal = pd.DataFrame({"decision_date": ["2025-01-01"], "max_sample_date": ["2025-01-02"]})
    with pytest.raises(ValueError, match="Temporal leakage"):
        assert_temporal_cutoff(temporal)


def test_preregistration_linter_blocks_posthoc_selection():
    config = {
        "primary_metric": "hits@20",
        "kill_condition": "paired lower CI <= 0",
        "minimum_detectable_effect": 0.237,
        "claimed_effect": 0.3,
        "external_sets_touched": [],
        "test_split_identifiers": ["heldout_patient"],
    }
    assert lint_preregistration(config).passed
    result = lint_preregistration(config, changed_source="select_best_score(heldout_patient)")
    assert not result.passed


def test_masking_ablation_sets_are_balanced_blind_and_disjoint():
    rows = []
    for label in (0, 1):
        for index in range(20):
            rows.append(
                {
                    "patient_id": f"p{label}-{index}",
                    "mutant_peptide": f"MUT{label}{index}",
                    "wildtype_peptide": f"WTT{label}{index}",
                    "hla_allele": "HLA-A*02:01",
                    "label": label,
                    "Gene_Symbol": f"G{index}",
                    "Mutation_Consequence": "M",
                    "Protein_position": index,
                }
            )
    generated = generate_masking_sets(pd.DataFrame(rows), seeds=range(2), n_per_condition=10)
    for payload in generated.values():
        a_ids = {row["candidate_id"] for row in payload["condition_a"]}
        b_ids = {row["candidate_id"] for row in payload["condition_b"]}
        assert not a_ids & b_ids
        assert all("label" not in row for row in payload["condition_a"] + payload["condition_b"])
        answers = pd.DataFrame(payload["answer_key"])
        assert answers.groupby("condition")["label"].sum().to_dict() == {
            "A_UNMASKED": 5,
            "B_MASKED": 5,
        }


@pytest.mark.skipif(
    not (
        Path("/tmp/IMPROVE_paper/data.zip").is_file()
        and Path("/tmp/IMPROVE_paper/results.zip").is_file()
    ),
    reason="clone https://github.com/SRHgroup/IMPROVE_paper to /tmp/IMPROVE_paper",
)
def test_official_improve_regressions():
    values = regression_values("/tmp/IMPROVE_paper")
    assert (values["rows"], values["positives"], values["patients"], values["partitions"]) == (
        17_520,
        467,
        70,
        5,
    )
    exact = {
        "oracle_mean_hits_at_20": 6.457142857142857,
        "netmhcpan_rankel_4_1": 0.9142857142857143,
        "prime": 1.2,
        "prioscore": 0.5857142857142857,
        "foreignness": 0.5142857142857142,
        "dai_4_1": 0.6428571428571429,
        "rf": 1.4428571428571428,
        "rf_without_prime": 1.3857142857142857,
        "rf_p_at_least_one": 41 / 61,
        "prime_p_at_least_one": 41 / 61,
        "paired_mde": 0.237374,
    }
    for name, expected in exact.items():
        assert values[name] == pytest.approx(expected, abs=5e-4)
    assert values["random_expectation"] == pytest.approx(0.5818, abs=0.01)
    assert values["positive_count_sd"] == pytest.approx(6.626, abs=0.01)
    assert values["paired_sd_diff"] == pytest.approx(1.013, abs=0.01)
    assert values["unreachable_patients"] == 9

    paired = values["rf_vs_prime"]
    assert paired.delta == pytest.approx(0.242857)
    assert paired.lo == pytest.approx(0.014286, abs=1e-6)
    assert paired.hi == pytest.approx(0.485714, abs=1e-6)
    capture = values["rf_vs_prime_capture_fraction"]
    assert capture.delta == pytest.approx(0.0302, abs=1e-4)
    assert capture.lo < 0 < capture.hi
    clinical = values["rf_vs_prime_p_at_least_one"]
    assert clinical.delta == 0.0

    expected_headroom = [
        (20, 1.2, 1.2, 29),
        (50, 1.2, 2.3, 19),
        (100, 1.2, 3.657142857, 13),
        (200, 1.2, 5.371428571, 10),
        ("all", 1.2, 6.457142857, 9),
    ]
    for row, expected in zip(values["headroom"], expected_headroom, strict=True):
        slate, base, oracle, zero = expected
        assert row["slate"] == slate
        assert row["base"] == pytest.approx(base)
        assert row["oracle"] == pytest.approx(oracle)
        assert row["zero_positive_patients"] == zero

    canary = values["tie_break_canary"]
    assert canary["source_order"] == pytest.approx(2.471428571)
    assert canary["constant_md5"] != pytest.approx(canary["source_order"])
