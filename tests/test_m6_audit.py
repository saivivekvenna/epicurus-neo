from epicurus_neo.m6.audit import assemble_audit, render_audit_markdown


def _fake_track(track, verdict):
    return {
        "track": track,
        "model": "logistic",
        "baseline": "prevalence",
        "verdict": verdict,
        "macro_hits_at_k": {"delta": 0.1, "delta_ci": [-0.2, 0.4], "per_study": {}},
        "micro_hits_at_k": {"delta_vs_baseline": 0.05, "delta_ci": [-0.1, 0.2]},
        "classification": {
            "auroc": 0.55,
            "brier": 0.2,
            "average_precision": 0.3,
            "calibration": [],
        },
        "ranking_informative_patients": 8,
        "per_fold": {},
    }


def test_audit_stamps_insufficiency_and_carries_both_tracks():
    audit = assemble_audit(
        universal=_fake_track("universal", "CONSISTENT_WITH_NO_EFFECT"),
        presentation=_fake_track("presentation", "REJECT"),
        completeness=[],
        prevalence=[],
        confound={"accuracy": 0.9, "majority_rate": 0.5, "per_study": {}},
        availability=[],
    )
    assert audit["corpus_verdict"] == "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"
    assert audit["universal"]["verdict"] == "CONSISTENT_WITH_NO_EFFECT"
    markdown = render_audit_markdown(audit)
    assert "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA" in markdown
    assert "Universal track" in markdown
    assert "ranking-informative" in markdown.lower()


def test_render_handles_skipped_presentation_track():
    audit = assemble_audit(
        universal=_fake_track("universal", "REJECT"),
        presentation={
            "verdict": "SKIPPED_PRESENTATION_UNAVAILABLE",
            "macro_hits_at_k": {"delta": float("nan"), "delta_ci": [None, None]},
            "classification": {"auroc": float("nan"), "brier": float("nan")},
            "ranking_informative_patients": 0,
        },
        completeness=[],
        prevalence=[{"study_id": "hu_neovax_2021", "positive_rate": 0.24}],
        confound={"accuracy": 0.9, "majority_rate": 0.5, "per_study": {}},
        availability=[],
    )
    # Must not raise on NaN/None CI bounds.
    markdown = render_audit_markdown(audit)
    assert "SKIPPED_PRESENTATION_UNAVAILABLE" in markdown
