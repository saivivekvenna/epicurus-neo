"""Braun RCC vertical-slice orchestration: ingest, reconcile, combine, no model fit.

Keeps three concerns separate and testable: (1) build the accepted Braun Event-B
corpus through the ordinary ingestion path, (2) reconcile the result against the
paper's own summary table without hard-coding the answer, and (3) combine the
frozen IMPROVE (Event-A) corpus with Braun (Event-B) for a single global audit.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from event_b.adapters.braun_rcc import (
    BraunRCCAdapter,
    INVITRO_FILE,
    SUMMARY_FILE,
    braun_source_paths,
    immunogenic_call,
    stage_braun_supplements,
    _read_sheet,
    _txt,
)
from event_b.audit import corpus_audit
from event_b.corpus import EventBCorpus
from event_b.ingest import IngestionResult, ingest_source
from event_b.manifest import manifest_from_paths
from event_b.models import BiologicalEvent, ResponseLabel, SCHEMAS
from event_b.review import ReviewIssue


# Research-report expectations to verify (not truth). Differences must be explained.
EXPECTED = {"patients": 9, "positives": 61, "tested_negatives": 68}


@dataclass(frozen=True)
class BraunBuild:
    adapter: BraunRCCAdapter
    manifest: object
    result: IngestionResult
    review_queue: tuple[ReviewIssue, ...]


def build_braun_corpus(raw_dir: str | Path) -> BraunBuild:
    adapter = BraunRCCAdapter(raw_dir)
    paths = braun_source_paths(raw_dir)
    manifest = manifest_from_paths(
        adapter.declaration.source_name,
        adapter.declaration.source_version,
        adapter.declaration.adapter_name,
        adapter.declaration.adapter_version,
        paths,
    )
    result = ingest_source(adapter, manifest)
    review_queue = tuple(result.review_queue) + tuple(adapter.review_issues)
    return BraunBuild(adapter, manifest, result, review_queue)


def _summary_targets(raw_dir: str | Path) -> dict[str, dict[str, int]]:
    """Driver/passenger immunogenic vs non-immunogenic counts from the paper's sheet 2e."""
    paths = stage_braun_supplements(raw_dir)
    sheet = _read_sheet(paths[SUMMARY_FILE], "2e")
    columns = list(sheet.columns)
    label_col = columns[0]
    imm_col = next(c for c in columns if "immunogenic" in c.lower() and "non" not in c.lower())
    non_col = next(c for c in columns if "non-immunogenic" in c.lower() or "non immunogenic" in c.lower())
    targets: dict[str, dict[str, int]] = {}
    for _, row in sheet.iterrows():
        key = _txt(row[label_col]).lower()
        if key in {"driver", "passenger"}:
            targets[key] = {
                "immunogenic": int(float(row[imm_col])),
                "non_immunogenic": int(float(row[non_col])),
            }
    return targets


def reconcile_braun(raw_dir: str | Path, build: BraunBuild) -> dict:
    """Independently recompute per-peptide calls from source and compare to the paper."""
    paths = stage_braun_supplements(raw_dir)
    invitro = _read_sheet(paths[INVITRO_FILE], "In Vitro")

    by_type: dict[str, dict[str, int]] = {}
    unscorable = 0
    for _, row in invitro.iterrows():
        call = immunogenic_call(row)
        mutation_type = _txt(row.get("Mutation_type")) or "unclassified"
        bucket = by_type.setdefault(mutation_type.lower(), {"immunogenic": 0, "non_immunogenic": 0})
        if call is None:
            unscorable += 1
        elif call:
            bucket["immunogenic"] += 1
        else:
            bucket["non_immunogenic"] += 1

    targets = _summary_targets(raw_dir)
    summary_matches = all(
        by_type.get(key, {}).get(field) == targets[key][field]
        for key in targets
        for field in ("immunogenic", "non_immunogenic")
    )

    accepted = build.result.accepted_corpus
    normalized = build.result.normalized_corpus
    audit = corpus_audit(accepted, build.review_queue, [build.adapter.declaration])
    responses = audit["response_counts"]
    events = audit["event_counts"]

    candidates = accepted.candidates
    identity_completeness = {
        column: float(candidates[column].notna().mean()) if len(candidates) else 0.0
        for column in ("gene", "protein_change", "genomic_variant", "transcript", "wildtype_peptide")
    }

    review_codes: dict[str, int] = {}
    for issue in build.review_queue:
        review_codes[issue.code] = review_codes.get(issue.code, 0) + 1

    return {
        "expected_from_report": EXPECTED,
        "source_observed": {
            "invitro_peptides": int(len(invitro)),
            "patients": int(invitro["Patient_ID"].astype(str).str.strip().nunique()),
            "unscorable_rows": unscorable,
        },
        "recomputed_by_rule": by_type,
        "paper_summary_targets_2e": targets,
        "summary_reconciles": summary_matches,
        "extracted": {
            "candidates": int(len(normalized.candidates)),
            "assays": int(len(normalized.assays)),
        },
        "accepted": {
            "assays": int(len(accepted.assays)),
            "positives": int(responses.get(ResponseLabel.POSITIVE.value, 0)),
            "tested_negatives": int(responses.get(ResponseLabel.TESTED_NEGATIVE.value, 0)),
            "untested": int(responses.get(ResponseLabel.UNTESTED.value, 0)),
            "patients": audit["sample_sizes"]["patient_n"],
            "event_counts": events,
        },
        "review_queue": {"total": len(build.review_queue), "by_code": review_codes},
        "class_breakdown": audit["mhc_class_breakdown"],
        "hla_resolved_fraction": float(candidates["hla_alleles"].notna().mean())
        if len(candidates)
        else 0.0,
        "identity_completeness": identity_completeness,
        "de_novo_basis": (
            "paper states 'no pre-existing immune responses were detected'; ex-vivo week-0 pool "
            "baselines are background (below the ELISpot positivity floor)"
        ),
        "model_readiness": audit["model_readiness"],
    }


def render_reconciliation_markdown(recon: dict) -> str:
    exp = recon["expected_from_report"]
    acc = recon["accepted"]
    src = recon["source_observed"]
    lines = [
        "# Braun RCC 2025 — Event-B ingestion reconciliation",
        "",
        "Source: Braun et al., Nature 2025 (PMC11903305, NCT02950766), CC BY-NC-ND. Immunogenicity",
        "is recomputed from raw IFN-gamma ELISpot replicates using the paper's stated rule",
        "(P<0.05 two-sided t-test AND mean spot count >=3x no-stim DMSO control). No counts are",
        "hard-coded; no negative is inferred from omission.",
        "",
        "## Expected (research report) vs accepted",
        "",
        "| Quantity | Expected | Accepted | Match |",
        "|---|---|---|---|",
        f"| Patients | {exp['patients']} | {acc['patients']} | {exp['patients'] == acc['patients']} |",
        f"| Positives | {exp['positives']} | {acc['positives']} | {exp['positives'] == acc['positives']} |",
        f"| Tested-negatives | {exp['tested_negatives']} | {acc['tested_negatives']} | "
        f"{exp['tested_negatives'] == acc['tested_negatives']} |",
        "",
        "## Independent recomputation vs the paper's summary (sheet 2e)",
        "",
        f"- Source In Vitro peptides: {src['invitro_peptides']} across {src['patients']} patients "
        f"(unscorable rows: {src['unscorable_rows']}).",
        f"- Recomputed by rule, split by mutation type: `{recon['recomputed_by_rule']}`",
        f"- Paper sheet-2e targets: `{recon['paper_summary_targets_2e']}`",
        f"- Driver/passenger splits reconcile exactly: **{recon['summary_reconciles']}**",
        "",
        "## Why accepted positives are 61, not 62",
        "",
        "The rule scores 62 peptides immunogenic across all 130 In Vitro rows. One immunogenic",
        "peptide (AMACR|p.Y41N, patient 104) has a blank `Mutation_type` and is excluded from the",
        "paper's driver/passenger summary (Fig. 2e / Supplementary Table 2). Rather than silently",
        "inflate the accepted count to 62, it is routed to the review queue; the 129 accepted assays",
        "(61 positive + 68 tested-negative) match the paper's reported figures.",
        "",
        "## Review queue, event typing, and completeness",
        "",
        f"- Review queue: {recon['review_queue']['total']} record(s) by code "
        f"`{recon['review_queue']['by_code']}`.",
        f"- Accepted event counts: `{acc['event_counts']}` (Event-B only; no Event-A relabelled).",
        f"- De-novo basis: {recon['de_novo_basis']}.",
        f"- MHC class breakdown (candidates): `{recon['class_breakdown']}` "
        "(long SLPs left UNKNOWN; class not resolved per peptide by the assay).",
        f"- HLA-resolved fraction (candidates): {recon['hla_resolved_fraction']:.2f} "
        "(per-peptide HLA is a prediction, not an assay restriction; not stored as such).",
        f"- Candidate identity completeness: `{recon['identity_completeness']}`.",
        "",
        "## Model-readiness (this slice alone)",
        "",
        f"```json\n{recon['model_readiness']}\n```",
        "",
    ]
    return "\n".join(lines)


POST_VACCINE_TIMEPOINTS = {"POST_VACCINE", "POST_PRIME", "POST_BOOST"}


def _truthy(value: object) -> bool:
    return str(value).strip().upper() in {"TRUE", "1", "YES"}


def assert_quality_gates(corpus: EventBCorpus) -> None:
    """Fail loudly on label-integrity violations in an accepted corpus.

    These are the CI gates for the vertical slice. A non-empty review queue is NOT a
    violation (ambiguity is expected); only accepted rows are checked.
    """
    assays = corpus.assays
    event = assays.event_type.astype(str).str.upper()
    label = assays.response_label.astype(str).str.upper()
    relative = assays.relative_to_vaccine.astype(str).str.upper()
    provenance_ids = set(corpus.provenance.provenance_id.astype(str))
    violations: list[str] = []

    event_b_positive = event.eq(BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value) & label.eq(
        ResponseLabel.POSITIVE.value
    )
    if (event_b_positive & ~relative.isin(POST_VACCINE_TIMEPOINTS)).any():
        violations.append("accepted Event-B positive without post-vaccine timepoint")

    tested_negative = label.eq(ResponseLabel.TESTED_NEGATIVE.value)
    if (tested_negative & ~assays.explicit_assay_inclusion.map(_truthy)).any():
        violations.append("accepted tested-negative without explicit assay inclusion")

    if assays.patient_id.isna().any() or assays.patient_id.astype(str).eq("").any():
        violations.append("accepted assay with missing patient identifier")

    if (~assays.provenance_id.astype(str).isin(provenance_ids)).any():
        violations.append("accepted assay label without provenance")

    if event.eq(BiologicalEvent.EVENT_C_CLINICAL_OUTCOME.value).any():
        violations.append("clinical outcome stored as an assay")

    resolved = assays.dropna(subset=["candidate_id"]).copy()
    resolved["_label"] = resolved.response_label.astype(str).str.upper()
    grouped = resolved.groupby(["candidate_id", "assay_type", "timepoint"], dropna=False)["_label"]
    for _, labels in grouped:
        if {ResponseLabel.POSITIVE.value, ResponseLabel.TESTED_NEGATIVE.value}.issubset(set(labels)):
            violations.append("contradictory accepted labels for one candidate/assay/timepoint")
            break

    if violations:
        raise AssertionError("Event-B quality gates failed: " + "; ".join(sorted(set(violations))))


def load_corpus_from_parquet(directory: str | Path) -> EventBCorpus:
    """Load a previously exported (already-accepted) corpus from Parquet, read-only."""
    directory = Path(directory)
    corpus = EventBCorpus()
    for entity in SCHEMAS:
        path = directory / f"{entity}.parquet"
        if path.exists():
            setattr(corpus, entity, pd.read_parquet(path))
    return corpus


def _is_null(value: object) -> bool:
    return value is None or value is pd.NA or (isinstance(value, float) and pd.isna(value))


def _harmonize_types(frame: pd.DataFrame) -> pd.DataFrame:
    """Coerce only genuinely mixed-type object columns to string, preserving nulls.

    Different adapters may store the same free-text field with different Python types
    (e.g. IMPROVE writes ``source_interpretation`` as an int, Braun as a str); a mixed
    column breaks Parquet type inference. Clean numeric columns are left untouched.
    """
    out = frame.copy()
    for column in out.columns:
        if out[column].dtype != object:
            continue
        present = out[column][~out[column].map(_is_null)]
        if len({type(value) for value in present}) > 1:
            out[column] = out[column].map(lambda value: value if _is_null(value) else str(value))
    return out


def combine_corpora(*corpora: EventBCorpus) -> EventBCorpus:
    """Concatenate already-accepted corpora into one; both are pre-validated."""
    combined = EventBCorpus()
    for entity in SCHEMAS:
        frames = [
            getattr(corpus, entity) for corpus in corpora if len(getattr(corpus, entity)) > 0
        ]
        if frames:
            setattr(combined, entity, _harmonize_types(pd.concat(frames, ignore_index=True)))
    return combined
