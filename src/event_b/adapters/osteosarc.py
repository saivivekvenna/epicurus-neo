"""Declaration-only adapter for the deep n=1 Osteosarc acceptance case."""

from __future__ import annotations

import pandas as pd

from event_b.adapters.base import AdapterDeclaration
from event_b.corpus import EventBCorpus
from event_b.manifest import SourceManifest
from event_b.models import BiologicalEvent, SCHEMAS


class OsteosarcCaseStudyAdapter:
    declaration = AdapterDeclaration(
        "Osteosarc longitudinal case",
        "source-supplied",
        "osteosarc_case_study",
        "1.0.0",
        (
            "studies",
            "patients",
            "vaccines",
            "candidates",
            "assays",
            "clinical_outcomes",
            "provenance",
        ),
        (
            BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value,
            BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
            BiologicalEvent.EVENT_C_CLINICAL_OUTCOME.value,
        ),
        (
            "One patient; never a population-scale training cohort",
            "Candidate testing is selection-conditioned",
            "Most generated variants have no candidate-resolved assay label",
        ),
        ("Only explicitly reported peptide-specific assays receive response labels",),
        ("Population-level effect estimates",),
    )

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        raise RuntimeError("Osteosarc source files were not supplied; no records were fabricated")

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        del manifest
        corpus = EventBCorpus()
        for entity, records in extracted.items():
            if entity not in SCHEMAS:
                raise ValueError(f"Unknown Osteosarc entity: {entity}")
            setattr(corpus, entity, SCHEMAS[entity].normalize(pd.DataFrame(records)))
        return corpus
