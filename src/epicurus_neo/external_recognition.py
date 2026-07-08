from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from epicurus_neo.plm_retrieval import load_plm_embedding_cache
from epicurus_neo.schema import normalize_hla, normalize_peptide


AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")


def normalize_vdjdb(
    frame: pd.DataFrame,
    *,
    min_score: int = 0,
) -> pd.DataFrame:
    required = {
        "species",
        "antigen.epitope",
        "antigen.species",
        "mhc.a",
        "mhc.class",
        "vdjdb.score",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"VDJdb input is missing columns: {sorted(missing)}")

    work = frame[
        (frame["species"] == "HomoSapiens")
        & (frame["mhc.class"] == "MHCI")
        & (pd.to_numeric(frame["vdjdb.score"], errors="coerce").fillna(0) >= min_score)
    ].copy()
    work["mutant_peptide"] = work["antigen.epitope"].map(normalize_peptide)
    work["hla_allele"] = work["mhc.a"].map(normalize_hla)
    work["vdjdb_score"] = pd.to_numeric(work["vdjdb.score"], errors="coerce").fillna(0)
    work = work[
        work["mutant_peptide"].map(
            lambda peptide: 8 <= len(peptide) <= 14 and set(peptide).issubset(AMINO_ACIDS)
        )
        & work["hla_allele"].str.startswith("HLA-")
    ]
    work["recognition_origin"] = np.where(
        work["antigen.species"].fillna("").astype(str).str.lower() == "homosapiens",
        "human",
        "pathogen",
    )

    grouped = (
        work.groupby(["mutant_peptide", "hla_allele", "recognition_origin"], as_index=False)
        .agg(
            recognition_support=("mutant_peptide", "size"),
            recognition_max_evidence=("vdjdb_score", "max"),
        )
        .sort_values(["hla_allele", "mutant_peptide", "recognition_origin"])
        .reset_index(drop=True)
    )
    grouped.insert(
        0,
        "candidate_id",
        [f"vdjdb:{idx}" for idx in range(len(grouped))],
    )
    grouped["source_dataset"] = "vdjdb"
    grouped["label"] = "positive"
    return grouped


def normalize_vdjdb_file(
    input_path: str | Path,
    output_path: str | Path,
    *,
    min_score: int = 0,
) -> Path:
    frame = pd.read_csv(input_path, sep="\t", low_memory=False)
    out = normalize_vdjdb(frame, min_score=min_score)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output


def _similarity_summary(
    similarities: np.ndarray,
    weights: np.ndarray,
    *,
    prefix: str,
    top_k: int,
) -> dict[str, float]:
    descending = np.sort(similarities)[::-1]
    weighted_descending = np.sort(similarities * weights)[::-1]
    return {
        f"{prefix}_max_similarity": float(descending[0]) if len(descending) else float("nan"),
        f"{prefix}_topk_similarity_mean": (
            float(np.mean(descending[:top_k])) if len(descending) else float("nan")
        ),
        f"{prefix}_max_weighted_similarity": (
            float(weighted_descending[0]) if len(weighted_descending) else float("nan")
        ),
        f"{prefix}_reference_count": float(len(similarities)),
    }


def add_external_recognition_features(
    frame: pd.DataFrame,
    reference: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    *,
    top_k: int = 5,
) -> pd.DataFrame:
    out = frame.copy()
    ref = reference.copy()
    ref["mutant_peptide"] = ref["mutant_peptide"].map(normalize_peptide)
    ref["hla_allele"] = ref["hla_allele"].map(normalize_hla)
    ref = ref.drop_duplicates(["mutant_peptide", "hla_allele", "recognition_origin"])
    embedded = ref["mutant_peptide"].map(embeddings.get)
    ref = ref.loc[embedded.notna()].reset_index(drop=True)
    reference_matrix = np.stack(embedded.loc[embedded.notna()].to_list()).astype(np.float32)
    embedding_center = reference_matrix.mean(axis=0)
    centered_reference_matrix = reference_matrix - embedding_center
    centered_norms = np.linalg.norm(centered_reference_matrix, axis=1, keepdims=True)
    centered_reference_matrix = centered_reference_matrix / np.maximum(centered_norms, 1e-12)
    weights = (
        np.log1p(pd.to_numeric(ref["recognition_support"], errors="coerce").fillna(0).to_numpy())
        * (
            1.0
            + 0.25
            * pd.to_numeric(ref["recognition_max_evidence"], errors="coerce")
            .fillna(0)
            .to_numpy()
        )
    )
    hla_values = ref["hla_allele"].to_numpy()
    origin_values = ref["recognition_origin"].to_numpy()

    rows: list[dict[str, float]] = []
    for row in out.itertuples():
        peptide = normalize_peptide(row.mutant_peptide)
        hla = normalize_hla(row.hla_allele)
        query_embedding = embeddings.get(peptide)
        similarities = reference_matrix @ query_embedding if query_embedding is not None else None
        centered_similarities = None
        if query_embedding is not None:
            centered_query = query_embedding - embedding_center
            centered_query = centered_query / max(float(np.linalg.norm(centered_query)), 1e-12)
            centered_similarities = centered_reference_matrix @ centered_query
        hla_mask = hla_values == hla
        pathogen_mask = origin_values == "pathogen"
        human_mask = origin_values == "human"
        empty = np.asarray([], dtype=np.float32)
        summary = {}
        summary.update(
            _similarity_summary(
                similarities[hla_mask] if similarities is not None else empty,
                weights[hla_mask] if similarities is not None else empty,
                prefix="recognition_hla",
                top_k=top_k,
            )
        )
        summary.update(
            _similarity_summary(
                centered_similarities[hla_mask] if centered_similarities is not None else empty,
                weights[hla_mask] if centered_similarities is not None else empty,
                prefix="recognition_hla_centered",
                top_k=top_k,
            )
        )
        summary.update(
            _similarity_summary(
                similarities[pathogen_mask] if similarities is not None else empty,
                weights[pathogen_mask] if similarities is not None else empty,
                prefix="recognition_pathogen",
                top_k=top_k,
            )
        )
        summary.update(
            _similarity_summary(
                centered_similarities[pathogen_mask]
                if centered_similarities is not None
                else empty,
                weights[pathogen_mask] if centered_similarities is not None else empty,
                prefix="recognition_pathogen_centered",
                top_k=top_k,
            )
        )
        summary.update(
            _similarity_summary(
                similarities[human_mask] if similarities is not None else empty,
                weights[human_mask] if similarities is not None else empty,
                prefix="recognition_human",
                top_k=top_k,
            )
        )
        summary.update(
            _similarity_summary(
                centered_similarities[human_mask] if centered_similarities is not None else empty,
                weights[human_mask] if centered_similarities is not None else empty,
                prefix="recognition_human_centered",
                top_k=top_k,
            )
        )
        summary["recognition_pathogen_minus_human_similarity"] = (
            summary["recognition_pathogen_max_similarity"]
            - summary["recognition_human_max_similarity"]
        )
        summary["recognition_centered_pathogen_minus_human_similarity"] = (
            summary["recognition_pathogen_centered_max_similarity"]
            - summary["recognition_human_centered_max_similarity"]
        )
        rows.append(summary)

    feature_frame = pd.DataFrame(rows, index=out.index)
    for column in feature_frame:
        out[column] = feature_frame[column]
    return out


def add_external_recognition_features_file(
    input_path: str | Path,
    reference_path: str | Path,
    embedding_cache_path: str | Path,
    output_path: str | Path,
    *,
    top_k: int = 5,
) -> Path:
    frame = pd.read_csv(input_path)
    reference = pd.read_csv(reference_path)
    embeddings, _ = load_plm_embedding_cache(embedding_cache_path)
    out = add_external_recognition_features(frame, reference, embeddings, top_k=top_k)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output
