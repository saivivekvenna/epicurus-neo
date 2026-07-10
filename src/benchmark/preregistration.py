"""Deterministic portion of the plan §8.4 pre-registration linter."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any


ALLOWED_METRICS = {"hits@20", "capture_fraction", "p_at_least_one", "precision@20", "mrr"}


@dataclass(frozen=True)
class LintResult:
    passed: bool
    errors: tuple[str, ...]


def lint_preregistration(config: dict[str, Any], *, changed_source: str = "") -> LintResult:
    errors: list[str] = []
    metric = config.get("primary_metric")
    if metric not in ALLOWED_METRICS:
        errors.append(f"primary_metric must be one of {sorted(ALLOWED_METRICS)}")
    if not str(config.get("kill_condition", "")).strip():
        errors.append("kill_condition must be present and falsifiable")
    try:
        mde_value = float(config["minimum_detectable_effect"])
        claimed = float(config["claimed_effect"])
        if claimed < mde_value:
            errors.append("claimed_effect is smaller than minimum_detectable_effect")
    except (KeyError, TypeError, ValueError):
        errors.append("minimum_detectable_effect and claimed_effect must be numeric")
    external = config.get("external_sets_touched", [])
    if external and not bool(config.get("milestone_run")):
        errors.append("external_sets_touched must be empty outside a milestone run")
    lowered_source = changed_source.lower()
    forbidden = tuple(str(value).lower() for value in config.get("test_split_identifiers", []))
    if any(identifier and identifier in lowered_source for identifier in forbidden):
        errors.append("test-split identifier appears in changed feature/selection source")
    selector_tokens = ("score_selection", "gridsearch", "selector.fit", "select_best_score")
    if any(token in lowered_source for token in selector_tokens):
        errors.append("score-column selection/search is forbidden in Milestone 1")
    return LintResult(not errors, tuple(errors))
