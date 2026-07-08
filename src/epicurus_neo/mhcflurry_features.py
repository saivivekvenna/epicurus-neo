from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from epicurus_neo.schema import add_normalized_columns


def _load_default_predictor() -> Any:
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError as exc:
        raise RuntimeError(
            "MHCflurry is not installed. Install the optional dependency with "
            "`pip install -e '.[mhc]'`."
        ) from exc
    return Class1PresentationPredictor.load()


def add_mhcflurry_predictions(
    frame: pd.DataFrame,
    *,
    predictor: Any | None = None,
    peptide_col: str = "mutant_peptide",
    allele_col: str = "hla_allele",
) -> pd.DataFrame:
    """Add MHCflurry class-I affinity/presentation features to a canonical table."""
    if peptide_col not in frame.columns or allele_col not in frame.columns:
        raise ValueError(f"Missing required columns: {peptide_col}, {allele_col}")

    out = add_normalized_columns(frame)
    valid = (out["mutant_peptide_norm"] != "") & (out["hla_allele_norm"] != "HLA-")
    if predictor is None:
        predictor = _load_default_predictor()

    predictions = pd.DataFrame(index=out.index)
    predictions[
        [
            "mhcflurry_affinity",
            "mhcflurry_processing_score",
            "mhcflurry_presentation_score",
            "mhcflurry_presentation_percentile",
        ]
    ] = pd.NA

    if valid.any():
        for allele, group in out.loc[valid].groupby("hla_allele_norm", sort=False):
            predicted = predictor.predict(
                peptides=group["mutant_peptide_norm"].tolist(),
                alleles=[allele],
                verbose=0,
                throw=False,
            )
            predicted.index = group.index
            predictions.loc[group.index, "mhcflurry_affinity"] = predicted["affinity"]
            predictions.loc[group.index, "mhcflurry_processing_score"] = predicted["processing_score"]
            predictions.loc[group.index, "mhcflurry_presentation_score"] = predicted[
                "presentation_score"
            ]
            predictions.loc[group.index, "mhcflurry_presentation_percentile"] = predicted[
                "presentation_percentile"
            ]

    for column in predictions.columns:
        out[column] = pd.to_numeric(predictions[column], errors="coerce")
    return out


def add_mhcflurry_predictions_file(input_path: str | Path, output_path: str | Path) -> Path:
    frame = pd.read_csv(input_path)
    out = add_mhcflurry_predictions(frame)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output
