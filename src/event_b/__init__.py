"""Canonical Event-B vaccine-response corpus and recognition-evidence substrate."""

from event_b.corpus import EventBCorpus
from event_b.models import BiologicalEvent, ResponseLabel, stable_candidate_id
from event_b.validation import validate_corpus

__all__ = [
    "BiologicalEvent",
    "EventBCorpus",
    "ResponseLabel",
    "stable_candidate_id",
    "validate_corpus",
]
