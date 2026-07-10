from __future__ import annotations

import json
from pathlib import Path

from epicurus_neo.m6.audit import assemble_audit, render_audit_markdown
from epicurus_neo.m6.confounds import prevalence_by_study, study_only_classifier
from epicurus_neo.m6.dataset import CORPUS_DIR, completeness_report, load_label_frame
from epicurus_neo.m6.evaluate import evaluate_track
from epicurus_neo.m6.presentation import (
    PRESENTATION_STUDIES,
    PresentationUnavailable,
    add_presentation_score,
    presentation_availability,
    resolve_class_i_alleles,
)


MIN_COMPATIBLE_CANDIDATES = 10


def _inert_presentation(verdict: str, reason: str, **extra) -> dict:
    base = {
        "track": "presentation",
        "verdict": verdict,
        "reason": reason,
        "macro_hits_at_k": {"delta": float("nan"), "delta_ci": [None, None]},
        "classification": {"auroc": float("nan"), "brier": float("nan")},
        "ranking_informative_patients": 0,
        "per_fold": {},
    }
    base.update(extra)
    return base


def _presentation_track(frame, *, seed: int, bootstrap_n: int) -> dict:
    """Presentation comparison on the presentation-COMPATIBLE subset only.

    MHCflurry (class-I, <=15mer) cannot score long SLPs, so the baseline is only
    defined for short class-I peptides. Comparing against a baseline blind to most
    candidates is meaningless; we restrict to scored candidates and require >=2
    studies to retain a viable subset, else report the incompatibility honestly.
    """
    try:
        scored = add_presentation_score(frame)
    except PresentationUnavailable as exc:
        return _inert_presentation("SKIPPED_PRESENTATION_UNAVAILABLE", str(exc))

    compatible = scored[
        scored.study_id.isin(PRESENTATION_STUDIES) & scored.presentation_score.notna()
    ].reset_index(drop=True)
    per_study = compatible.groupby("study_id").agg(n=("candidate_id", "size"), pos=("label", "sum"))
    compatibility = {
        study: {
            "n_compatible": int(per_study.n.get(study, 0)),
            "n_positive": int(per_study.pos.get(study, 0)),
        }
        for study in PRESENTATION_STUDIES
    }
    viable = [
        study
        for study in PRESENTATION_STUDIES
        if compatibility[study]["n_positive"] > 0
        and compatibility[study]["n_compatible"] >= MIN_COMPATIBLE_CANDIDATES
    ]
    if len(viable) < 2:
        return _inert_presentation(
            "NOT_VIABLE_PRESENTATION_INCOMPATIBLE",
            "MHCflurry (class-I, <=15mer) scores only short peptides; long SLPs are "
            "incompatible, so fewer than two studies retain a viable presentation-"
            "comparable subset (the comparison would collapse to one study).",
            compatibility=compatibility,
        )
    subset = compatible[compatible.study_id.isin(viable)].reset_index(drop=True)
    result = evaluate_track(
        subset,
        model_name="logistic",
        baseline_name="presentation",
        track="presentation",
        tier="presentation",
        seed=seed,
        bootstrap_n=bootstrap_n,
    )
    result["compatibility"] = compatibility
    return result


def run(
    corpus_dir=CORPUS_DIR,
    out_dir: str | Path = "artifacts/milestone_6",
    *,
    seed: int = 17,
    bootstrap_n: int = 20_000,
) -> dict:
    """Run the M6A two-track LOSO swing and write the audit artifacts."""
    frame = load_label_frame(corpus_dir)
    universal = evaluate_track(
        frame,
        model_name="logistic",
        baseline_name="prevalence",
        track="universal",
        tier="core",
        seed=seed,
        bootstrap_n=bootstrap_n,
    )
    availability = presentation_availability(resolve_class_i_alleles(frame))
    presentation = _presentation_track(frame, seed=seed, bootstrap_n=bootstrap_n)

    audit = assemble_audit(
        universal=universal,
        presentation=presentation,
        completeness=completeness_report(frame),
        prevalence=prevalence_by_study(frame),
        confound=study_only_classifier(frame, seed=seed),
        availability=availability,
    )
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "m6a_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    (out / "m6a_audit.md").write_text(render_audit_markdown(audit))
    return audit
