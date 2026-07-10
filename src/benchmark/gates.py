"""Automated gates for hand-reasoned and LLM-derived scores."""

from __future__ import annotations

import pandas as pd

from benchmark.scorecard import scorecard


def retroactive_prime_checks() -> dict[str, bool]:
    """Apply the registered gates to the three historical hand-reasoned rules.

    The first two are paired hits@20 comparisons; the third was registered as
    a balanced-subsample AUROC experiment and therefore uses its own 0.5 gate.
    """
    checks = {
        "llm_articulated_rule": {"delta": -0.357, "ci": (-0.614, -0.100)},
        "contact_residue_gate": {"delta": -0.071, "ci": (float("nan"), float("nan"))},
    }
    result = {
        name: bool(values["delta"] > 0 and values["ci"][0] > 0) for name, values in checks.items()
    }
    result["anchor_creation_dai"] = bool(0.403 > 0.5)
    return result


def prime_rule(
    df: pd.DataFrame,
    candidate_score: str,
    group_col: str = "patient_id",
    *,
    prime_col: str = "Prime",
    label_col: str = "label",
    k: int = 20,
) -> bool:
    """Return True only when candidate hits@k beats PRIME with paired CI > 0."""
    if prime_col not in df.columns:
        aliases = ("prime_source_score", "PRIME", "prime")
        prime_col = next((column for column in aliases if column in df.columns), prime_col)
    if prime_col not in df.columns:
        raise ValueError("PRIME baseline column is required by the merge gate")
    report = scorecard(
        df,
        score_col=candidate_score,
        baseline_col=prime_col,
        group_col=group_col,
        k=k,
        label_col=label_col,
    )
    return report["verdict"] == "ACCEPT"
