"""Genuine MHCflurry class-I presentation features for peptide x HLA candidate frames.

Why this exists: lossless-recovered candidates are generated de novo (raw allele -> Ensembl -> windows)
and were only scored by genuine PRIME (whose backbone is MixMHCpred). To attribute the four-arm scorer
stage FAIRLY, recovered candidates need a genuine presentation feature — not a 0.5 neutral impute and
not the MixMHCpred proxy that double-counts PRIME. MHCflurry 2.x is an independent, learned class-I
presentation predictor (affinity + antigen-processing + presentation), so it supplies real evidence.

NetMHCpan (the frozen Epicurus ``el`` feature) is not locally runnable (licensed binary); MHCflurry is
the genuine substitute. Callers that recompute ``el`` with MHCflurry MUST disclose the substitution and
its agreement with NetMHCpan-EL — do not silently relabel MHCflurry output as NetMHCpan-EL.

Deterministic: MHCflurry inference is a fixed forward pass. Invalid peptides (empty / non-standard AA /
wrong length) return NaN rather than raising, so a mixed candidate frame scores in one pass.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

# Output columns (genuine MHCflurry, oriented as the predictor emits them):
#   mhcflurry_el_percentile : presentation %rank, LOWER = stronger presenter (EL-analog)
#   mhcflurry_affinity_nm   : predicted binding affinity in nM, LOWER = tighter binder (BA-analog)
#   mhcflurry_processing    : antigen-processing score in [0, 1], HIGHER = better processed
PRESENTATION_COLUMNS = (
    "mhcflurry_el_percentile",
    "mhcflurry_affinity_nm",
    "mhcflurry_processing",
)

_STD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def load_presentation_predictor() -> Any:
    """Load the default MHCflurry Class1PresentationPredictor (fail closed if the package is absent)."""
    try:
        from mhcflurry import Class1PresentationPredictor
    except ImportError as exc:  # pragma: no cover - environment guard
        raise RuntimeError(
            "MHCflurry is not installed. Install the optional dependency: pip install -e '.[mhc]' "
            "and run `mhcflurry-downloads fetch models_class1_presentation`."
        ) from exc
    return Class1PresentationPredictor.load()


def _normalize_allele(value: object) -> str:
    text = str(value).strip()
    if not text:
        return ""
    return text if text.upper().startswith("HLA-") else f"HLA-{text}"


def _valid_peptide(value: object) -> bool:
    pep = str(value).strip().upper()
    return bool(pep) and 8 <= len(pep) <= 15 and not (set(pep) - _STD_AA)


def add_presentation_features(
    frame: pd.DataFrame,
    *,
    peptide_col: str = "mutant_peptide",
    allele_col: str = "hla_allele",
    predictor: Any | None = None,
) -> pd.DataFrame:
    """Attach genuine MHCflurry presentation features to ``frame`` (row order preserved).

    Invalid peptide/allele rows get NaN for every presentation column. Grouped by allele so MHCflurry
    runs one batched forward pass per HLA.
    """
    if peptide_col not in frame.columns or allele_col not in frame.columns:
        raise ValueError(f"missing required columns: {peptide_col!r}, {allele_col!r}")
    if predictor is None:
        predictor = load_presentation_predictor()

    out = frame.copy()
    peptide = out[peptide_col].map(lambda p: str(p).strip().upper())
    allele = out[allele_col].map(_normalize_allele)
    valid = peptide.map(_valid_peptide) & (allele != "")

    for col in PRESENTATION_COLUMNS:
        out[col] = pd.NA

    for al, idx in out.index[valid].to_series().groupby(allele[valid], sort=False).groups.items():
        peps = peptide.loc[idx].tolist()
        predicted = predictor.predict(peptides=peps, alleles=[al], verbose=0, throw=False)
        predicted.index = idx
        out.loc[idx, "mhcflurry_el_percentile"] = predicted["presentation_percentile"].to_numpy()
        out.loc[idx, "mhcflurry_affinity_nm"] = predicted["affinity"].to_numpy()
        out.loc[idx, "mhcflurry_processing"] = predicted["processing_score"].to_numpy()

    for col in PRESENTATION_COLUMNS:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    return out
