import json

from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.event_a import load_event_a_frame
from epicurus_neo.m6.transfer import evaluate_transfer_track

_VERDICTS = {"ACCEPT_TRANSFER", "REJECT_TRANSFER", "CONSISTENT_WITH_NO_EFFECT_TRANSFER"}
_STUDIES = {"braun_rcc_2025", "hu_neovax_2021", "mkras_vax_2026", "pdac_neovax_2023"}


def test_transfer_track_reports_declared_gate_and_per_fold_auroc():
    result = evaluate_transfer_track(
        load_label_frame(), load_event_a_frame(), seed=17, bootstrap_n=2000
    )
    assert result["verdict"] in _VERDICTS
    assert set(result["per_fold"]) == _STUDIES
    for entry in result["per_fold"].values():
        assert {"baseline_auroc", "candidate_auroc", "auroc_delta"} <= set(entry)
    assert "macro_auroc_delta" in result
    assert 0 <= result["folds_improved"] <= 4
    # The candidate differs from the baseline by exactly the one Event-A teacher feature.
    assert result["teacher"]["n_event_a"] == 17082
    assert "macro_delta_hits_at_k" in result


def test_transfer_track_is_deterministic():
    first = evaluate_transfer_track(
        load_label_frame(), load_event_a_frame(), seed=17, bootstrap_n=2000
    )
    second = evaluate_transfer_track(
        load_label_frame(), load_event_a_frame(), seed=17, bootstrap_n=2000
    )
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
