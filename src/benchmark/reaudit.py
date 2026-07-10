"""Re-audit frozen out-of-fold score files without fitting or selecting."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from benchmark.improve import load_improve_rf
from benchmark.scorecard import scorecard


@dataclass(frozen=True)
class AuditRow:
    iteration: str
    artifact: str
    claimed_delta: str
    measured_delta: float
    ci_lo: float
    ci_hi: float
    verdict: str


def audit_frame(
    frame: pd.DataFrame,
    *,
    iteration: str,
    artifact: str,
    claimed_delta: str,
    score_col: str,
    baseline_col: str,
    group_col: str = "patient_id",
    ascending: bool = False,
    baseline_ascending: bool = False,
) -> AuditRow:
    report = scorecard(
        frame,
        score_col,
        baseline_col,
        group_col=group_col,
        ascending=ascending,
        baseline_ascending=baseline_ascending,
    )
    primary = report["hits@20"]
    return AuditRow(
        iteration=iteration,
        artifact=artifact,
        claimed_delta=claimed_delta,
        measured_delta=primary["delta_vs_baseline"],
        ci_lo=primary["delta_ci"][0],
        ci_hi=primary["delta_ci"][1],
        verdict=report["verdict"],
    )


def render_markdown(rows: list[AuditRow]) -> str:
    lines = [
        "| iteration | artifact | claimed Δ | measured Δ | 95% CI | verdict |",
        "|---:|---|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(
            f"| {row.iteration} | `{row.artifact}` | {row.claimed_delta} | "
            f"{row.measured_delta:+.4f} | [{row.ci_lo:+.4f}, {row.ci_hi:+.4f}] | "
            f"{row.verdict} |"
        )
    return "\n".join(lines) + "\n"


def write_markdown(rows: list[AuditRow], path: str | Path) -> None:
    Path(path).write_text(render_markdown(rows))


def reaudit_repository(repo_root: str | Path, improve_repo: str | Path) -> list[AuditRow]:
    """Re-audit every retained patient-level OOF score artifact from iterations 001–029."""
    root = Path(repo_root)
    rows: list[AuditRow] = []
    gartner_artifact = root / "outputs/benchmarks/gartner_patient_oof_scored.csv"
    gartner = pd.read_csv(gartner_artifact)
    rows.append(
        audit_frame(
            gartner,
            iteration="002",
            artifact=gartner_artifact.name,
            claimed_delta="+0.0769 vs prior learned ranker",
            score_col="epicurus_score",
            baseline_col="baseline_gartner_nmer_score",
        )
    )

    keys = ["patient_id", "mutant_peptide", "hla_allele"]
    epicurus_artifact = root / "outputs/improve_cv.oof_scored.csv"
    epicurus = pd.read_csv(epicurus_artifact)
    rows.append(
        audit_frame(
            epicurus,
            iteration="026",
            artifact=epicurus_artifact.name,
            claimed_delta="+0.2714 vs PRIME",
            score_col="epicurus_score",
            baseline_col="prime_source_score",
        )
    )

    rf = load_improve_rf(improve_repo)
    rows.append(
        audit_frame(
            rf,
            iteration="028",
            artifact="results.zip:pred_df_TME_excluded.txt",
            claimed_delta="+0.2429 vs PRIME",
            score_col="prediction_rf",
            baseline_col="Prime",
        )
    )

    xgb_artifact = root / "outputs/benchmarks/improve_xgb_slot_oof.scored.csv"
    xgb = pd.read_csv(xgb_artifact).rename(columns={"Patient": "patient_id"})
    xgb["patient_id"] = "improve:" + xgb["patient_id"].astype(str)
    xgb = xgb.merge(epicurus[keys + ["epicurus_score"]], on=keys, validate="one_to_one")
    rows.append(
        audit_frame(
            xgb,
            iteration="028",
            artifact=xgb_artifact.name,
            claimed_delta="+0.0143 vs Epicurus",
            score_col="xgb_slot_score",
            baseline_col="epicurus_score",
        )
    )

    prefixed_rf = rf.copy()
    prefixed_rf["patient_id"] = "improve:" + prefixed_rf["patient_id"].astype(str)
    blend = epicurus[keys + ["epicurus_score", "label"]].merge(
        prefixed_rf[keys + ["prediction_rf"]], on=keys, validate="one_to_one"
    )
    blend["rf_rank"] = blend.groupby("patient_id")["prediction_rf"].rank(method="average", pct=True)
    blend["epicurus_rank"] = blend.groupby("patient_id")["epicurus_score"].rank(
        method="average", pct=True
    )
    blend["rf_epicurus_blend"] = (blend["rf_rank"] + blend["epicurus_rank"]) / 2
    rows.append(
        audit_frame(
            blend,
            iteration="028",
            artifact="reconstructed frozen 50/50 RF-Epicurus blend",
            claimed_delta="+0.0714 vs RF",
            score_col="rf_epicurus_blend",
            baseline_col="prediction_rf",
        )
    )

    base = pd.read_csv(
        root / "data/processed/improve_cv.normalized.csv",
        usecols=keys + ["prime_source_score"],
    )
    esm_cases = (
        ("improve_none_esm_oof.csv", "+0.2143 vs PRIME"),
        ("improve_delta_esm_oof.csv", "-0.0143 vs PRIME"),
        ("improve_paired_esm_oof.csv", "+0.2143 vs PRIME"),
    )
    for artifact, claimed in esm_cases:
        frame = pd.read_csv(root / "outputs/benchmarks/paired_esm" / artifact).merge(
            base, on=keys, validate="one_to_one"
        )
        rows.append(
            audit_frame(
                frame,
                iteration="028",
                artifact=artifact,
                claimed_delta=claimed,
                score_col="score",
                baseline_col="prime_source_score",
            )
        )

    neoguider_artifact = root / "outputs/benchmarks/improve_neoguider_official_cv.scored.csv"
    neoguider = pd.read_csv(neoguider_artifact)
    neoguider = neoguider.loc[neoguider["feature_set"] == "tme_excluded"].rename(
        columns={
            "Patient": "patient_id",
            "Mut_peptide": "mutant_peptide",
            "HLA_allele": "hla_allele",
        }
    )
    neoguider["patient_id"] = "improve:" + neoguider["patient_id"].astype(str)
    neoguider = neoguider.merge(base, on=keys, validate="one_to_one")
    rows.append(
        audit_frame(
            neoguider,
            iteration="028",
            artifact=f"{neoguider_artifact.name}:tme_excluded",
            claimed_delta="+0.1429 vs PRIME",
            score_col="neoguider_score",
            baseline_col="prime_source_score",
        )
    )

    neoprecis_artifact = root / "outputs/benchmarks/improve_neoprecis_approx.scored.csv"
    neoprecis = pd.read_csv(neoprecis_artifact).merge(base, on=keys, validate="one_to_one")
    rows.append(
        audit_frame(
            neoprecis,
            iteration="029",
            artifact=neoprecis_artifact.name,
            claimed_delta="-0.6286 vs PRIME",
            score_col="neoprecis_immuno_score",
            baseline_col="prime_source_score",
        )
    )
    return rows
