from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class FailureExample:
    candidate_id: str
    group_id: str
    rank: int
    score: float
    label: str
    peptide: str
    hla: str
    gene: str


@dataclass(frozen=True)
class FailureReport:
    k: int
    score_col: str
    group_col: str
    groups: int
    positives_total: int
    positives_missed_at_k: int
    negatives_in_top_k: int
    false_negatives: tuple[FailureExample, ...]
    false_positives: tuple[FailureExample, ...]
    numeric_feature_means: dict[str, dict[str, float]]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _ranked(frame: pd.DataFrame, *, group_col: str, score_col: str) -> pd.DataFrame:
    ranked = frame.copy()
    ranked["_rank"] = (
        ranked.sort_values(score_col, ascending=False, kind="mergesort")
        .groupby(group_col)
        .cumcount()
        + 1
    )
    return ranked


def _example(row: pd.Series, *, group_col: str, score_col: str) -> FailureExample:
    return FailureExample(
        candidate_id=str(row.get("candidate_id", "")),
        group_id=str(row.get(group_col, "")),
        rank=int(row["_rank"]),
        score=float(row[score_col]),
        label=str(row.get("label", "")),
        peptide=str(row.get("mutant_peptide", "")),
        hla=str(row.get("hla_allele", "")),
        gene=str(row.get("gene_symbol", "")),
    )


def _numeric_feature_means(frame: pd.DataFrame, masks: dict[str, pd.Series]) -> dict[str, dict[str, float]]:
    numeric_cols = [
        col
        for col in frame.columns
        if pd.api.types.is_numeric_dtype(frame[col])
        and col not in {"label_weight", "_rank"}
        and not col.startswith("epicurus_")
    ]
    means: dict[str, dict[str, float]] = {}
    for name, mask in masks.items():
        subset = frame.loc[mask, numeric_cols]
        if subset.empty:
            means[name] = {}
            continue
        means[name] = {
            col: float(value)
            for col, value in subset.mean(numeric_only=True).replace({np.nan: None}).dropna().items()
        }
    return means


def build_failure_report(
    scored: pd.DataFrame,
    *,
    group_col: str = "patient_id",
    score_col: str = "epicurus_score",
    k: int = 20,
    max_examples: int = 20,
) -> FailureReport:
    required = {group_col, score_col, "label"}
    missing = required - set(scored.columns)
    if missing:
        raise ValueError(f"Missing required columns for failure report: {sorted(missing)}")

    ranked = _ranked(scored[scored["label"].isin(["positive", "negative"])], group_col=group_col, score_col=score_col)
    in_top_k = ranked["_rank"] <= k
    is_positive = ranked["label"] == "positive"
    is_negative = ranked["label"] == "negative"
    false_negative_mask = is_positive & ~in_top_k
    false_positive_mask = is_negative & in_top_k

    false_negatives = tuple(
        _example(row, group_col=group_col, score_col=score_col)
        for _, row in ranked.loc[false_negative_mask]
        .sort_values(["_rank", score_col], ascending=[True, False])
        .head(max_examples)
        .iterrows()
    )
    false_positives = tuple(
        _example(row, group_col=group_col, score_col=score_col)
        for _, row in ranked.loc[false_positive_mask]
        .sort_values(score_col, ascending=False)
        .head(max_examples)
        .iterrows()
    )

    masks = {
        "true_positive_top_k": is_positive & in_top_k,
        "false_negative_missed": false_negative_mask,
        "false_positive_top_k": false_positive_mask,
        "true_negative_below_k": is_negative & ~in_top_k,
    }

    return FailureReport(
        k=k,
        score_col=score_col,
        group_col=group_col,
        groups=int(ranked[group_col].nunique()),
        positives_total=int(is_positive.sum()),
        positives_missed_at_k=int(false_negative_mask.sum()),
        negatives_in_top_k=int(false_positive_mask.sum()),
        false_negatives=false_negatives,
        false_positives=false_positives,
        numeric_feature_means=_numeric_feature_means(ranked, masks),
    )


def make_hypothesis_prompt(report: FailureReport) -> str:
    payload = json.dumps(report.to_dict(), indent=2)
    return f"""You are helping improve Epicurus Neo for the exact neoantigen ranking task.

Goal:
- Improve held-out hits@{report.k}, precision@{report.k}, recall@{report.k}, and nDCG@{report.k}.
- Do not optimize pooled AUC at the expense of top-k ranking.
- Do not propose using locked test labels.
- Do not propose LLM direct peptide scoring.

Given the failure report below, propose up to five concrete modeling hypotheses.
Each hypothesis must be testable by code and benchmarked with leakage-aware grouped CV.

Return YAML with this schema:

hypotheses:
  - id: short_snake_case
    rationale: one sentence
    expected_metric_effect: hits@k / precision@k / recall@k / calibration
    implementation: exact feature/model/scoring change
    risk: leakage, overfit, missing data, biology caveat
    validation: grouped split and acceptance criterion

Failure report:

```json
{payload}
```
"""


def write_research_artifacts(
    report: FailureReport,
    *,
    output_dir: str | Path,
) -> tuple[Path, Path]:
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    report_path = root / "failure_report.json"
    prompt_path = root / "hypothesis_prompt.md"
    report_path.write_text(json.dumps(report.to_dict(), indent=2))
    prompt_path.write_text(make_hypothesis_prompt(report))
    return report_path, prompt_path

