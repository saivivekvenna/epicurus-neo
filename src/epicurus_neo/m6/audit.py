from __future__ import annotations

import pandas as pd

CORPUS_VERDICT = "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"


def _records(value) -> list[dict]:
    if isinstance(value, pd.DataFrame):
        return value.to_dict("records")
    return list(value)


def assemble_audit(
    *, universal, presentation, completeness, prevalence, confound, availability
) -> dict:
    return {
        "corpus_verdict": CORPUS_VERDICT,
        "note": "Diagnostic swing under the standing insufficiency verdict; not a headline claim.",
        "universal": universal,
        "presentation": presentation,
        "completeness": _records(completeness),
        "prevalence_by_study": _records(prevalence),
        "study_confound": confound,
        "presentation_availability": _records(availability),
    }


def _fmt(value: object) -> str:
    if isinstance(value, bool):
        return str(value)
    if isinstance(value, (int, float)) and value == value:  # not NaN
        return f"{value:.4f}"
    return str(value)


def _fmt_ci(ci: object) -> str:
    if isinstance(ci, (list, tuple)) and len(ci) == 2:
        return f"[{_fmt(ci[0])}, {_fmt(ci[1])}]"
    return str(ci)


def _track_lines(title: str, track: dict, opponent: str) -> list[str]:
    macro = track.get("macro_hits_at_k", {})
    classification = track.get("classification", {})
    lines = [
        f"## {title} (learned vs {opponent})",
        f"- Verdict: **{track.get('verdict', 'n/a')}**",
        f"- Macro-study Δ hits@k_patient: {_fmt(macro.get('delta'))} CI {_fmt_ci(macro.get('delta_ci'))}",
        f"- Macro AUROC (per-study mean): {_fmt(classification.get('auroc'))} | "
        f"Brier: {_fmt(classification.get('brier'))}",
    ]
    per_fold = classification.get("per_fold_auroc")
    if isinstance(per_fold, dict) and per_fold:
        lines.append(
            "- Per-held-out-study AUROC: "
            + ", ".join(f"{study}={_fmt(value)}" for study, value in per_fold.items())
        )
    if "pooled_out_of_fold_auroc" in classification:
        lines.append(
            f"- Pooled OOF AUROC (caveated): {_fmt(classification['pooled_out_of_fold_auroc'])} "
            "— conflates cross-study calibration shift with discrimination"
        )
    lines.append(
        f"- Ranking-informative patients (n_eligible > k): "
        f"{track.get('ranking_informative_patients', 'n/a')}"
    )
    if track.get("reason"):
        lines.append(f"- Reason: {track['reason']}")
    compatibility = track.get("compatibility")
    if isinstance(compatibility, dict) and compatibility:
        lines.append(
            "- Presentation-compatible candidates: "
            + ", ".join(
                f"{study}={info['n_compatible']} ({info['n_positive']} pos)"
                for study, info in compatibility.items()
            )
        )
    return lines


def _completeness_lines(records: list) -> list[str]:
    has_neg = sum(1 for r in records if r.get("denominator_type") == "HAS_TESTED_NEGATIVE")
    no_neg = sum(1 for r in records if r.get("denominator_type") == "NO_TESTED_NEGATIVE")
    return [
        "## Candidate-universe completeness gate",
        f"- Patients with >=1 tested negative (HAS_TESTED_NEGATIVE, rankable): {has_neg}",
        f"- Patients with no tested negative (NO_TESTED_NEGATIVE): {no_neg}",
        "- NO_TESTED_NEGATIVE is a rankability flag (no negative to rank against), not a "
        "denominator-bias claim: e.g. the mKRAS 6/6-responders sit on a complete shared "
        "six-peptide panel. These patients are excluded from primary top-k and kept in "
        "pooled classification.",
    ]


def render_audit_markdown(audit: dict) -> str:
    lines = [
        "# Milestone 6A audit: Event-B-only recognition swing",
        "",
        f"**Corpus verdict (standing):** `{audit['corpus_verdict']}`",
        "",
        audit["note"],
        "",
        *_track_lines("Universal track", audit["universal"], "prevalence (all 4 studies)"),
        "",
        *_track_lines("Presentation track", audit["presentation"], "presentation-only (hu + pdac)"),
        "",
        *_completeness_lines(audit.get("completeness", [])),
        "",
        "## Study confound",
        f"- Study-only classifier accuracy: {_fmt(audit['study_confound'].get('accuracy'))} "
        f"(majority rate {_fmt(audit['study_confound'].get('majority_rate'))})",
    ]
    prevalence = audit.get("prevalence_by_study", [])
    if prevalence:
        lines.append(
            "- Positive rate by study: "
            + ", ".join(f"{row['study_id']}={_fmt(row['positive_rate'])}" for row in prevalence)
        )
    return "\n".join(lines) + "\n"
