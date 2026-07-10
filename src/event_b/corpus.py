"""In-memory canonical Event-B corpus."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd

from event_b.models import SCHEMAS


def _empty(entity: str) -> pd.DataFrame:
    return pd.DataFrame(columns=SCHEMAS[entity].columns)


@dataclass
class EventBCorpus:
    studies: pd.DataFrame = field(default_factory=lambda: _empty("studies"))
    patients: pd.DataFrame = field(default_factory=lambda: _empty("patients"))
    vaccines: pd.DataFrame = field(default_factory=lambda: _empty("vaccines"))
    candidates: pd.DataFrame = field(default_factory=lambda: _empty("candidates"))
    assays: pd.DataFrame = field(default_factory=lambda: _empty("assays"))
    clinical_outcomes: pd.DataFrame = field(default_factory=lambda: _empty("clinical_outcomes"))
    recognition_evidence: pd.DataFrame = field(
        default_factory=lambda: _empty("recognition_evidence")
    )
    candidate_funnel_links: pd.DataFrame = field(
        default_factory=lambda: _empty("candidate_funnel_links")
    )
    provenance: pd.DataFrame = field(default_factory=lambda: _empty("provenance"))

    def normalized(self) -> "EventBCorpus":
        return EventBCorpus(
            **{entity: SCHEMAS[entity].normalize(getattr(self, entity)) for entity in SCHEMAS}
        )

    def tables(self) -> dict[str, pd.DataFrame]:
        return {entity: getattr(self, entity) for entity in SCHEMAS}
