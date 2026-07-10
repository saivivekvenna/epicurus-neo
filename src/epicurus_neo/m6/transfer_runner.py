"""M6B runner + audit renderer: Event-A -> Event-B transfer (auxiliary arm)."""

from __future__ import annotations

import json
from pathlib import Path

from epicurus_neo.m6.audit import CORPUS_VERDICT, _fmt, _fmt_ci
from epicurus_neo.m6.dataset import load_label_frame
from epicurus_neo.m6.event_a import load_event_a_frame
from epicurus_neo.m6.transfer import evaluate_transfer_track

_NOTE = (
    "M6B is auxiliary: an Event-A -> Event-B transfer probe that enters no primary success gate. "
    "Diagnostic swing under the standing insufficiency verdict; not a headline claim."
)


def assemble_transfer_audit(transfer: dict) -> dict:
    return {"corpus_verdict": CORPUS_VERDICT, "note": _NOTE, "transfer": transfer}


def render_transfer_audit_markdown(audit: dict) -> str:
    track = audit["transfer"]
    macro_hits = track.get("macro_delta_hits_at_k", {})
    lines = [
        "# Milestone 6B audit: Event-A -> Event-B transfer (auxiliary)",
        "",
        f"**Corpus verdict (standing):** `{audit['corpus_verdict']}`",
        "",
        audit["note"],
        "",
        f"**Question:** {track['question']}",
        "",
        f"## Declared-gate verdict: **{track['verdict']}**",
        f"- Macro per-fold AUROC delta (candidate - baseline): {_fmt(track['macro_auroc_delta'])}",
        f"- Folds improved: {track['folds_improved']} / {track['n_folds_scored']} "
        "(ACCEPT_TRANSFER needs >=3 and a positive macro delta and no harm)",
        f"- Macro Δ hits@k (reported, underpowered): {_fmt(macro_hits.get('delta'))} "
        f"CI {_fmt_ci(macro_hits.get('delta_ci'))}",
        f"- Ranking-informative patients: {track['ranking_informative_patients']}",
        "",
        "## Per-held-out-study AUROC",
        "Baseline = frozen M6A `logistic(core)`; candidate = `logistic(core + event_a_teacher_score)`.",
        "",
        "| Study | baseline | candidate | delta | n_eval |",
        "|---|---|---|---|---|",
    ]
    for study, entry in sorted(track["per_fold"].items()):
        lines.append(
            f"| {study} | {_fmt(entry['baseline_auroc'])} | {_fmt(entry['candidate_auroc'])} | "
            f"{_fmt(entry['auroc_delta'])} | {entry['n_eval']} |"
        )
    teacher = track["teacher"]
    lines += [
        "",
        "## Teacher",
        f"- Frozen Event-A teacher: `{teacher['model']}` on the `{teacher['tier']}` tier, trained on "
        f"{teacher['n_event_a']} IMPROVE Event-A rows ({teacher['n_event_a_positive']} positive). "
        "Labels never merged; the teacher never sees an Event-B row.",
        f"- Sanity: in-distribution 5-fold AUROC on Event-A = "
        f"{_fmt(teacher.get('in_distribution_auroc'))} (a genuine teacher); its score pools to AUROC = "
        f"{_fmt(teacher.get('event_b_pooled_auroc'))} on Event-B. Real Event-A signal that does not "
        "generalize to Event-B - not a weak teacher.",
        "- Event-A is short class-I (8-11mer); most Event-B is long SLP. The teacher score is added as "
        "one feature to the Event-B-only model and is the sole candidate-vs-baseline difference.",
    ]
    return "\n".join(lines) + "\n"


def run_m6b(
    out_dir: str | Path = "artifacts/milestone_6", *, seed: int = 17, bootstrap_n: int = 20_000
) -> dict:
    """Run the M6B transfer arm and write the audit artifacts."""
    transfer = evaluate_transfer_track(
        load_label_frame(), load_event_a_frame(), seed=seed, bootstrap_n=bootstrap_n
    )
    audit = assemble_transfer_audit(transfer)
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    (out / "m6b_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    (out / "m6b_audit.md").write_text(render_transfer_audit_markdown(audit))
    return audit
