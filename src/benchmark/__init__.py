"""Leakage-safe, patient-level benchmark instrumentation.

This package is deliberately independent of the historical
``epicurus_neo.benchmark`` training module.  Milestone 1 trains no models.
"""

from benchmark.gates import prime_rule
from benchmark.labels import Label, validate_labels
from benchmark.metrics import (
    capture_fraction,
    hits_at_k,
    mrr,
    p_at_least_one,
    precision_at_k,
)
from benchmark.scorecard import scorecard
from benchmark.stats import bootstrap_ci, mde, n_required, paired_bootstrap

__all__ = [
    "Label",
    "bootstrap_ci",
    "capture_fraction",
    "hits_at_k",
    "mde",
    "mrr",
    "n_required",
    "p_at_least_one",
    "paired_bootstrap",
    "precision_at_k",
    "prime_rule",
    "scorecard",
    "validate_labels",
]
