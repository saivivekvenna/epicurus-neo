"""Link Event-B observations to the existing candidate-reachability ledger."""

from __future__ import annotations

from hashlib import sha256

import pandas as pd

from benchmark.funnel import STAGES, ReachabilityStatus, validate_reachability_ledger
from event_b.corpus import EventBCorpus
from event_b.models import SCHEMAS


def _status_value(value) -> str:
    return value.value if isinstance(value, ReachabilityStatus) else str(value)


def link_event_b_to_funnel(corpus: EventBCorpus, ledger: pd.DataFrame) -> pd.DataFrame:
    """Attach canonical candidates without converting missing upstream data to loss."""
    candidate_ids = set(corpus.candidates.candidate_id.astype(str))
    missing = set(ledger.candidate_id.astype(str)).difference(candidate_ids)
    if missing:
        raise ValueError(
            f"Funnel ledger contains unknown Event-B candidates: {sorted(missing)[:5]}"
        )
    work = ledger.rename(columns={"candidate_id": "positive_id"}).copy()
    canonical = validate_reachability_ledger(
        work, positive_id_col="positive_id", patient_col="patient_id"
    ).rename(columns={"positive_id": "candidate_id"})
    assays_by_candidate = set(corpus.assays.candidate_id.dropna().astype(str))
    included = set(
        corpus.candidates.loc[
            corpus.candidates.vaccine_inclusion.astype(str).str.upper().eq("INCLUDED"),
            "candidate_id",
        ].astype(str)
    )
    rows = []
    for _, row in canonical.iterrows():
        candidate_id = str(row.candidate_id)
        output = {
            "funnel_link_id": "funnel:" + sha256(candidate_id.encode()).hexdigest()[:20],
            "candidate_id": candidate_id,
            "patient_id": row.patient_id,
            "study_id": row.get("study_id", pd.NA),
            **{stage: _status_value(row[stage]) for stage in STAGES},
            "recognition_scored": row.get(
                "recognition_scored", ReachabilityStatus.NOT_ASSESSED.value
            ),
            "vaccine_inclusion": (
                ReachabilityStatus.REACHED.value
                if candidate_id in included
                else ReachabilityStatus.NOT_ASSESSED.value
            ),
            "functional_assay": (
                ReachabilityStatus.REACHED.value
                if candidate_id in assays_by_candidate
                else ReachabilityStatus.NOT_ASSESSED.value
            ),
            "provenance_id": row.get("provenance_id", pd.NA),
        }
        rows.append(output)
    return SCHEMAS["candidate_funnel_links"].normalize(pd.DataFrame(rows))
