from __future__ import annotations

from pathlib import Path

import pandas as pd

from epicurus_neo.m6.dataset import CORPUS_DIR, parse_alleles

PRESENTATION_STUDIES = ("hu_neovax_2021", "pdac_neovax_2023")


class PresentationUnavailable(RuntimeError):
    """Raised when the MHCflurry presentation baseline cannot be computed."""


def _is_class_i_allele(allele: str) -> bool:
    """True for HLA class-I loci (A/B/C); MHCflurry rejects class-II (DP/DQ/DR)."""
    locus = allele.upper().replace("HLA-", "").lstrip()
    return bool(locus) and locus[0] in {"A", "B", "C"} and "/" not in allele


def _predicted_antigen_alleles(corpus_dir: Path) -> pd.DataFrame:
    antigens = pd.read_parquet(corpus_dir / "antigens.parquet")
    predicted = antigens[antigens.hla_evidence_type == "PREDICTED_BEST_BINDER"].copy()
    predicted["antigen_alleles"] = predicted.hla_alleles.map(parse_alleles)
    keyed = predicted.groupby(["study_id", "gene", "protein_change"])["antigen_alleles"].agg(
        lambda lists: sorted({allele for sublist in lists for allele in sublist})
    )
    return keyed.rename("antigen_alleles").reset_index()


def resolve_class_i_alleles(frame: pd.DataFrame, corpus_dir: str | Path = CORPUS_DIR) -> pd.DataFrame:
    """Attach a per-candidate ``class_i_alleles`` list (hu direct; pdac by antigen join)."""
    candidate_alleles = frame.hla_alleles.map(parse_alleles).tolist()
    predicted = _predicted_antigen_alleles(Path(corpus_dir))
    out = frame.merge(
        predicted, on=["study_id", "gene", "protein_change"], how="left", validate="many_to_one"
    )

    def _resolve(direct: list[str], antigen: object) -> list[str]:
        source = direct if direct else (list(antigen) if isinstance(antigen, list) else [])
        return [allele for allele in source if _is_class_i_allele(allele)]

    out["class_i_alleles"] = [
        _resolve(direct, antigen)
        for direct, antigen in zip(candidate_alleles, out["antigen_alleles"])
    ]
    return out.drop(columns=["antigen_alleles"])


def presentation_availability(frame: pd.DataFrame) -> pd.DataFrame:
    if "class_i_alleles" not in frame.columns:
        frame = resolve_class_i_alleles(frame)
    has = frame.class_i_alleles.map(lambda alleles: len(alleles) > 0)
    table = (
        frame.assign(_has=has)
        .groupby("study_id")["_has"]
        .agg(resolved="sum", total="size")
        .reset_index()
    )
    table["resolved"] = table.resolved.astype(int)
    return table


def add_presentation_score(frame: pd.DataFrame, *, predictor: object | None = None) -> pd.DataFrame:
    """Add ``presentation_score`` = best MHCflurry presentation across class-I alleles."""
    try:
        from epicurus_neo.mhcflurry_features import add_mhcflurry_predictions
    except ImportError as exc:  # pragma: no cover - exercised only without the extra
        raise PresentationUnavailable(str(exc)) from exc

    resolved = resolve_class_i_alleles(frame).reset_index(drop=True)
    resolved["presentation_score"] = float("nan")
    exploded = resolved[["class_i_alleles"]].explode("class_i_alleles")
    exploded = exploded[exploded.class_i_alleles.map(lambda a: isinstance(a, str) and bool(a))]
    if exploded.empty:
        return resolved

    scoreable = resolved.loc[exploded.index].copy()
    scoreable["hla_allele"] = exploded.class_i_alleles.to_numpy()
    try:
        scored = add_mhcflurry_predictions(scoreable, predictor=predictor, allele_col="hla_allele")
    except Exception as exc:  # model-load / prediction failure: degrade, never crash the swing
        raise PresentationUnavailable(str(exc)) from exc
    scored = scored.assign(_row=exploded.index.to_numpy())
    best = scored.groupby("_row")["mhcflurry_presentation_score"].max()
    resolved.loc[best.index, "presentation_score"] = best.to_numpy()
    return resolved
