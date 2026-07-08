from __future__ import annotations

import math

import numpy as np
import pandas as pd


NON_FEATURE_COLUMNS = {
    "candidate_id",
    "source_dataset",
    "study_id",
    "patient_id",
    "hla_allele",
    "hla_allele_norm",
    "mutant_peptide",
    "mutant_peptide_norm",
    "wildtype_peptide",
    "wildtype_peptide_norm",
    "mutant_hla_key",
    "wildtype_hla_key",
    "label",
    "label_weight",
    "assay_type",
    "split",
}


def infer_numeric_feature_columns(frame: pd.DataFrame) -> list[str]:
    """Infer usable numeric feature columns from a canonical candidate table."""
    columns: list[str] = []
    for column in frame.columns:
        if column in NON_FEATURE_COLUMNS:
            continue
        if pd.api.types.is_numeric_dtype(frame[column]):
            columns.append(column)
    return sorted(columns)


def safe_log_inverse(value: object) -> float:
    """Score a positive quantity where smaller is better, such as binding nM."""
    if value is None or pd.isna(value):
        return float("nan")
    value_float = float(value)
    if value_float <= 0:
        return float("nan")
    return -math.log10(value_float)


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    std = values.std(skipna=True)
    if pd.isna(std) or std == 0:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (values - values.mean(skipna=True)) / std


def add_baseline_scores(frame: pd.DataFrame) -> pd.DataFrame:
    """Add deterministic baseline rank scores when source feature columns exist."""
    out = frame.copy()

    if "binding_affinity_nm" in out.columns:
        out["baseline_binding_score"] = out["binding_affinity_nm"].map(safe_log_inverse)
    elif "binding_score" in out.columns:
        out["baseline_binding_score"] = pd.to_numeric(out["binding_score"], errors="coerce")

    presentation_candidates = [
        "presentation_score",
        "bigmhc_el_score",
        "mhcflurry_presentation_score",
    ]
    for column in presentation_candidates:
        if column in out.columns:
            out["baseline_presentation_score"] = pd.to_numeric(out[column], errors="coerce")
            break

    components: list[pd.Series] = []
    if "baseline_binding_score" in out.columns:
        components.append(zscore(out["baseline_binding_score"]))
    if "baseline_presentation_score" in out.columns:
        components.append(zscore(out["baseline_presentation_score"]))
    if "expression_tpm" in out.columns:
        components.append(zscore(np.log1p(pd.to_numeric(out["expression_tpm"], errors="coerce"))))
    if "mutant_wildtype_binding_delta" in out.columns:
        components.append(zscore(out["mutant_wildtype_binding_delta"]))
    if "foreignness_score" in out.columns:
        components.append(zscore(out["foreignness_score"]))

    if components:
        out["baseline_pvac_style_score"] = sum(components) / len(components)

    return out

