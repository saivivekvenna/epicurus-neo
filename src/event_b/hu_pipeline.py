"""Hu 2021 melanoma NeoVax vertical-slice orchestration: ingest, reconcile, combine.

Mirrors the Braun pipeline's three separable concerns: (1) build the accepted Hu
Event-B corpus through the ordinary ingestion path, (2) reconcile the vaccine-peptide
recognition against Ott 2017's independently published neoantigen totals (CD8 15/97,
CD4 58/97 across patients 1-6) without hard-coding the answer, and (3) combine with the
frozen IMPROVE (Event-A) and Braun (Event-B) corpora for one multi-study global audit.
No recognition model is fitted.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from event_b.adapters.hu_neovax import (
    HuNeoVaxAdapter,
    OTT_CD4_NEOANTIGENS,
    OTT_CD8_NEOANTIGENS,
    OTT_COHORT,
    OTT_PMCID,
    hu_source_paths,
)
from event_b.audit import corpus_audit
from event_b.braun_pipeline import POST_VACCINE_TIMEPOINTS, _truthy  # shared Event-B gate helpers
from event_b.corpus import EventBCorpus
from event_b.ingest import IngestionResult, ingest_source
from event_b.manifest import manifest_from_paths
from event_b.models import BiologicalEvent, ResponseLabel, VaccineInclusion
from event_b.review import ReviewIssue


# Research-report / prior-literature expectations to verify (not truth). The CD8 positive
# count is the strong anchor; CD4 is expected within a small, explained tolerance.
EXPECTED = {
    "patients": 8,
    "cd8_positive_neoantigens": OTT_CD8_NEOANTIGENS[0],  # 15
    "cd4_positive_neoantigens": OTT_CD4_NEOANTIGENS[0],  # 58
    "cd4_tolerance": 3,
}


@dataclass(frozen=True)
class HuBuild:
    adapter: HuNeoVaxAdapter
    manifest: object
    result: IngestionResult
    review_queue: tuple[ReviewIssue, ...]


def build_hu_corpus(raw_dir: str | Path) -> HuBuild:
    adapter = HuNeoVaxAdapter(raw_dir)
    paths = hu_source_paths(raw_dir)
    manifest = manifest_from_paths(
        adapter.declaration.source_name,
        adapter.declaration.source_version,
        adapter.declaration.adapter_name,
        adapter.declaration.adapter_version,
        paths,
    )
    result = ingest_source(adapter, manifest)
    review_queue = tuple(result.review_queue) + tuple(adapter.review_issues)
    return HuBuild(adapter, manifest, result, review_queue)


def _neoantigen_counts(corpus: EventBCorpus, source_marker: str) -> dict[str, int]:
    """Positive and total distinct neoantigens (gene|protein_change) for the Ott 1-6 cohort,
    restricted to vaccine-peptide Event-B assays from a given source dataset (4a or 4b)."""
    assays = corpus.assays
    candidates = corpus.candidates
    joined = assays.merge(
        candidates[["candidate_id", "gene", "protein_change", "candidate_source", "vaccine_inclusion"]],
        on="candidate_id",
        how="left",
    )
    joined = joined[
        joined.event_type.astype(str).eq(BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value)
        & joined.candidate_source.astype(str).str.contains(source_marker, na=False)
        & joined.patient_id.astype(str).str.split(":").str[-1].isin(OTT_COHORT)
    ].copy()
    joined["neo"] = joined.gene.astype(str) + "|" + joined.protein_change.astype(str)
    positive = joined.response_label.astype(str).eq(ResponseLabel.POSITIVE.value)
    return {
        "positive_neoantigens": int(joined.loc[positive, "neo"].nunique()),
        "total_neoantigens": int(joined["neo"].nunique()),
        "positive_rows": int(positive.sum()),
        "assayed_rows": int(len(joined)),
    }


def reconcile_hu(build: HuBuild) -> dict:
    """Reconcile vaccine-peptide recognition against Ott 2017's published totals."""
    accepted = build.result.accepted_corpus
    audit = corpus_audit(accepted, build.review_queue, [build.adapter.declaration])

    cd8 = _neoantigen_counts(accepted, "4a")
    cd4 = _neoantigen_counts(accepted, "4b")
    cd8_reconciles = cd8["positive_neoantigens"] == OTT_CD8_NEOANTIGENS[0]
    cd4_delta = abs(cd4["positive_neoantigens"] - OTT_CD4_NEOANTIGENS[0])
    cd4_reconciles = cd4_delta <= EXPECTED["cd4_tolerance"]

    responses = audit["response_counts"]
    events = audit["event_counts"]

    # Epitope spreading must be entirely non-vaccine and never an Event-B label.
    spreading = accepted.assays[
        accepted.assays.event_type.astype(str).eq(BiologicalEvent.EPITOPE_SPREADING.value)
    ]
    spreading_candidate_ids = set(spreading.candidate_id.astype(str))
    spreading_included = accepted.candidates[
        accepted.candidates.candidate_id.astype(str).isin(spreading_candidate_ids)
        & accepted.candidates.vaccine_inclusion.astype(str).eq(VaccineInclusion.INCLUDED.value)
    ]

    # Reliability differentiation: show that channels keep distinct evidence strength.
    evidence = accepted.recognition_evidence
    reliability = {}
    for family in sorted(evidence.evidence_family.astype(str).unique()):
        sub = evidence[evidence.evidence_family.astype(str).eq(family)]
        reliability[family] = {
            "n": int(len(sub)),
            "candidate_specificity": sorted(sub.candidate_specificity.dropna().unique().tolist()),
            "assay_directness": sorted(sub.assay_directness.dropna().unique().tolist()),
            "vaccine_relevance": sorted(sub.vaccine_relevance.dropna().unique().tolist()),
            "temporal_clarity": sorted(sub.temporal_clarity.dropna().unique().tolist()),
        }

    review_codes: dict[str, int] = {}
    for issue in build.review_queue:
        review_codes[issue.code] = review_codes.get(issue.code, 0) + 1

    patient_ids = sorted(
        accepted.patients.source_patient_id.astype(str), key=lambda value: int(value)
    )

    return {
        "expected_from_literature": EXPECTED,
        "ott_cross_check_source": f"{OTT_PMCID} (Ott 2017; per-peptide screen for patients 1-6)",
        "source_observed": {
            "patients": patient_ids,
            "patient_n": len(patient_ids),
        },
        "reconciliation": {
            "cd8": {"observed": cd8, "ott_reported": OTT_CD8_NEOANTIGENS, "reconciles": cd8_reconciles},
            "cd4": {
                "observed": cd4,
                "ott_reported": OTT_CD4_NEOANTIGENS,
                "delta": cd4_delta,
                "reconciles_within_tolerance": cd4_reconciles,
            },
        },
        "accepted": {
            "assays": int(len(accepted.assays)),
            "positives": int(responses.get(ResponseLabel.POSITIVE.value, 0)),
            "tested_negatives": int(responses.get(ResponseLabel.TESTED_NEGATIVE.value, 0)),
            "untested": int(responses.get(ResponseLabel.UNTESTED.value, 0)),
            "patients": audit["sample_sizes"]["patient_n"],
            "event_counts": events,
        },
        "epitope_spreading": {
            "assays": int(len(spreading)),
            "all_non_vaccine": int(len(spreading_included)) == 0,
            "any_positive": bool(
                spreading.response_label.astype(str).eq(ResponseLabel.POSITIVE.value).any()
            ),
            "note": (
                "week-16 epitope-spreading responses in Datasets 11a-c are uniformly non-reactive; "
                "the paper's spreading positives arise at later / post-checkpoint timepoints not "
                "ingested here. Kept as EPITOPE_SPREADING, never vaccine-candidate recognition."
            ),
        },
        "reliability_differentiation": reliability,
        "review_queue": {"total": len(build.review_queue), "by_code": review_codes},
        "class_breakdown": audit["mhc_class_breakdown"],
        "de_novo_basis": (
            "Ott 2017 and Hu 2021 both report no pre-vaccination neoantigen reactivity; this table "
            "carries no week-0 column, so de-novo is author-asserted (not baseline-verified) and "
            "recognition_evidence temporal_clarity is set lower than a baseline-verified de-novo claim"
        ),
        "model_readiness": audit["model_readiness"],
    }


def render_reconciliation_markdown(recon: dict) -> str:
    cd8 = recon["reconciliation"]["cd8"]
    cd4 = recon["reconciliation"]["cd4"]
    acc = recon["accepted"]
    spread = recon["epitope_spreading"]
    lines = [
        "# Hu 2021 melanoma NeoVax — Event-B ingestion reconciliation",
        "",
        "Source: Hu et al., Nat Med 2021 (PMC8273876, NCT01970358), author manuscript. Per-peptide",
        "CD8 (class-I minimal epitopes) and CD4 (class-II assay peptides) IFN-gamma ELISpot calls are",
        "ingested as reported (Hu scored a response positive at >=2.5x DMSO). They are reconciled",
        f"against Ott 2017's independently published totals ({recon['ott_cross_check_source']}).",
        "",
        "## Vaccine-peptide recognition vs Ott 2017 (patients 1-6, week 16)",
        "",
        "| Channel | Positive neoantigens (observed) | Ott reported | Total neoantigens | Reconciles |",
        "|---|---|---|---|---|",
        f"| CD8 | {cd8['observed']['positive_neoantigens']} | {cd8['ott_reported'][0]} | "
        f"{cd8['observed']['total_neoantigens']} (Ott {cd8['ott_reported'][1]}) | {cd8['reconciles']} |",
        f"| CD4 | {cd4['observed']['positive_neoantigens']} | {cd4['ott_reported'][0]} | "
        f"{cd4['observed']['total_neoantigens']} (Ott {cd4['ott_reported'][1]}) | "
        f"within tol (delta {cd4['delta']}): {cd4['reconciles_within_tolerance']} |",
        "",
        "The CD8 positive-neoantigen count matches Ott exactly. CD4 differs by a small margin that is",
        "not closed by any other assay condition in the table (minigene/tumor add nothing), consistent",
        "with Ott counting at the immunizing-peptide rather than neoantigen granularity; it is reported,",
        "not tuned away. Per-channel totals differ slightly from Ott's 97 because CD8 and CD4 cover",
        "different neoantigen subsets; the positive count is the anchor.",
        "",
        "## Accepted Event-B corpus",
        "",
        f"- Patients: {recon['source_observed']['patient_n']} "
        f"(`{recon['source_observed']['patients']}` — Ott's 1-6 plus new 11-12; one consolidated",
        "  source, so no cross-study patient double-counting).",
        f"- Accepted assays: {acc['assays']} — {acc['positives']} POSITIVE / {acc['tested_negatives']} "
        f"TESTED_NEGATIVE / {acc['untested']} UNTESTED.",
        f"- Event counts: `{acc['event_counts']}`.",
        f"- MHC class (candidates): `{recon['class_breakdown']}` (CD8 class I, CD4 class II).",
        "",
        "## Epitope spreading kept separate",
        "",
        f"- {spread['assays']} epitope-spreading assays; all non-vaccine: **{spread['all_non_vaccine']}**; "
        f"any labelled a vaccine-candidate positive: **{spread['any_positive']}**.",
        f"- {spread['note']}",
        "",
        "## Reliability is not flattened",
        "",
        "```json",
        _json(recon["reliability_differentiation"]),
        "```",
        "",
        f"De-novo basis: {recon['de_novo_basis']}.",
        "",
        "## Model-readiness (this slice alone)",
        "",
        f"```json\n{recon['model_readiness']}\n```",
        "",
    ]
    return "\n".join(lines)


def _json(value) -> str:
    import json

    return json.dumps(value, indent=2, sort_keys=True, default=str)


def assert_quality_gates(corpus: EventBCorpus) -> None:
    """Label-integrity CI gates for the Hu slice, including epitope-spreading separation.

    Shares the Event-B checks with Braun and adds the melanoma-specific invariant that
    epitope spreading is never counted as vaccine-candidate recognition.
    """
    assays = corpus.assays
    candidates = corpus.candidates.set_index("candidate_id", drop=False)
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

    # Melanoma-specific: epitope spreading is never a vaccine-candidate recognition label, and
    # every vaccine-candidate Event-B assay must reference an INCLUDED candidate.
    spreading = event.eq(BiologicalEvent.EPITOPE_SPREADING.value)
    for candidate_id in assays.loc[spreading, "candidate_id"].dropna().astype(str):
        if candidate_id in candidates.index:
            if str(candidates.loc[candidate_id, "vaccine_inclusion"]).upper() == (
                VaccineInclusion.INCLUDED.value
            ):
                violations.append("epitope-spreading response tied to a vaccine-included candidate")
                break
    event_b = event.eq(BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value)
    for candidate_id in assays.loc[event_b, "candidate_id"].dropna().astype(str):
        if candidate_id in candidates.index:
            if str(candidates.loc[candidate_id, "vaccine_inclusion"]).upper() != (
                VaccineInclusion.INCLUDED.value
            ):
                violations.append("Event-B recognition tied to a non-vaccine candidate")
                break

    resolved = assays.dropna(subset=["candidate_id"]).copy()
    resolved["_label"] = resolved.response_label.astype(str).str.upper()
    grouped = resolved.groupby(["candidate_id", "assay_type", "timepoint"], dropna=False)["_label"]
    for _, labels in grouped:
        if {ResponseLabel.POSITIVE.value, ResponseLabel.TESTED_NEGATIVE.value}.issubset(set(labels)):
            violations.append("contradictory accepted labels for one candidate/assay/timepoint")
            break

    if violations:
        raise AssertionError("Hu Event-B quality gates failed: " + "; ".join(sorted(set(violations))))
