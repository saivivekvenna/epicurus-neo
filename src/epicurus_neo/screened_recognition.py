from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from epicurus_neo.plm_retrieval import load_plm_embedding_cache
from epicurus_neo.schema import normalize_hla, normalize_peptide


def _hla_family(value: object) -> str:
    allele = normalize_hla(value)
    if "*" not in allele:
        return allele
    locus, fields = allele.split("*", maxsplit=1)
    return f"{locus}*{fields.split(':', maxsplit=1)[0]}"


def aggregate_screened_reference(reference: pd.DataFrame) -> pd.DataFrame:
    required = {"mutant_peptide", "hla_allele", "label"}
    missing = required.difference(reference.columns)
    if missing:
        raise ValueError(f"Screened reference is missing columns: {sorted(missing)}")

    work = reference[reference["label"].isin(["positive", "negative"])].copy()
    work["mutant_peptide"] = work["mutant_peptide"].map(normalize_peptide)
    work["hla_allele"] = work["hla_allele"].map(normalize_hla)
    work = work[work["mutant_peptide"].ne("") & work["hla_allele"].ne("")]
    work["positive_count"] = (work["label"] == "positive").astype(int)
    work["negative_count"] = (work["label"] == "negative").astype(int)
    grouped = (
        work.groupby(["mutant_peptide", "hla_allele"], as_index=False)
        .agg(
            positive_count=("positive_count", "sum"),
            negative_count=("negative_count", "sum"),
        )
        .sort_values(["hla_allele", "mutant_peptide"])
        .reset_index(drop=True)
    )
    grouped["screen_count"] = grouped["positive_count"] + grouped["negative_count"]
    grouped["response_rate"] = grouped["positive_count"] / grouped["screen_count"]
    grouped["hla_family"] = grouped["hla_allele"].map(_hla_family)
    return grouped


def _empty_summary(prefix: str, top_ks: tuple[int, ...]) -> dict[str, float]:
    summary = {
        f"{prefix}_max_positive_similarity": float("nan"),
        f"{prefix}_max_negative_similarity": float("nan"),
        f"{prefix}_positive_minus_negative_similarity": float("nan"),
        f"{prefix}_reference_count": 0.0,
        f"{prefix}_positive_reference_count": 0.0,
    }
    for top_k in top_ks:
        summary[f"{prefix}_top{top_k}_response_rate"] = float("nan")
        summary[f"{prefix}_top{top_k}_weighted_response_rate"] = float("nan")
        summary[f"{prefix}_top{top_k}_mean_similarity"] = float("nan")
    return summary


def _neighborhood_summary(
    similarities: np.ndarray,
    positive_counts: np.ndarray,
    negative_counts: np.ndarray,
    *,
    prefix: str,
    top_ks: tuple[int, ...],
) -> dict[str, float]:
    if len(similarities) == 0:
        return _empty_summary(prefix, top_ks)

    response_rates = positive_counts / np.maximum(positive_counts + negative_counts, 1.0)
    order = np.argsort(similarities)[::-1]
    positive_mask = positive_counts > 0
    negative_mask = negative_counts > 0
    max_positive = (
        float(np.max(similarities[positive_mask])) if positive_mask.any() else float("nan")
    )
    max_negative = (
        float(np.max(similarities[negative_mask])) if negative_mask.any() else float("nan")
    )
    summary = {
        f"{prefix}_max_positive_similarity": max_positive,
        f"{prefix}_max_negative_similarity": max_negative,
        f"{prefix}_positive_minus_negative_similarity": max_positive - max_negative,
        f"{prefix}_reference_count": float(len(similarities)),
        f"{prefix}_positive_reference_count": float(positive_mask.sum()),
    }
    for top_k in top_ks:
        selected = order[:top_k]
        selected_similarity = similarities[selected]
        selected_positive = positive_counts[selected]
        selected_negative = negative_counts[selected]
        selected_support = selected_positive + selected_negative
        similarity_weights = np.exp(
            8.0 * (selected_similarity - float(np.max(selected_similarity)))
        )
        evidence_weights = similarity_weights * np.log1p(selected_support)
        summary[f"{prefix}_top{top_k}_response_rate"] = float(
            np.mean(response_rates[selected])
        )
        summary[f"{prefix}_top{top_k}_weighted_response_rate"] = float(
            np.average(response_rates[selected], weights=evidence_weights)
        )
        summary[f"{prefix}_top{top_k}_mean_similarity"] = float(
            np.mean(selected_similarity)
        )
    return summary


def add_screened_recognition_features(
    frame: pd.DataFrame,
    reference: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    *,
    top_ks: tuple[int, ...] = (1, 3, 5, 10, 20),
    prefix: str = "screened_recognition",
) -> pd.DataFrame:
    if not top_ks or any(top_k < 1 for top_k in top_ks):
        raise ValueError("top_ks must contain positive integers")

    out = frame.copy()
    ref = aggregate_screened_reference(reference)
    embedded = ref["mutant_peptide"].map(embeddings.get)
    ref = ref.loc[embedded.notna()].reset_index(drop=True)
    if ref.empty:
        raise ValueError("No screened reference peptides were found in the embedding cache")

    matrix = np.stack(embedded.loc[embedded.notna()].to_list()).astype(np.float32)
    center = matrix.mean(axis=0)
    centered_matrix = matrix - center
    centered_matrix /= np.maximum(
        np.linalg.norm(centered_matrix, axis=1, keepdims=True),
        1e-12,
    )
    alleles = ref["hla_allele"].to_numpy()
    families = ref["hla_family"].to_numpy()
    positive_counts = ref["positive_count"].to_numpy(dtype=np.float32)
    negative_counts = ref["negative_count"].to_numpy(dtype=np.float32)

    rows: list[dict[str, float]] = []
    for row in out.itertuples(index=False):
        peptide = normalize_peptide(row.mutant_peptide)
        allele = normalize_hla(row.hla_allele)
        family = _hla_family(allele)
        embedding = embeddings.get(peptide)
        if embedding is None:
            summary: dict[str, float] = {}
            for scope in ("hla", "family", "global"):
                summary.update(_empty_summary(f"{prefix}_{scope}", top_ks))
                summary.update(_empty_summary(f"{prefix}_{scope}_centered", top_ks))
            rows.append(summary)
            continue

        similarities = matrix @ embedding
        centered_embedding = embedding - center
        centered_embedding /= max(float(np.linalg.norm(centered_embedding)), 1e-12)
        centered_similarities = centered_matrix @ centered_embedding
        masks = {
            "hla": alleles == allele,
            "family": families == family,
            "global": np.ones(len(ref), dtype=bool),
        }
        summary = {}
        for scope, mask in masks.items():
            summary.update(
                _neighborhood_summary(
                    similarities[mask],
                    positive_counts[mask],
                    negative_counts[mask],
                    prefix=f"{prefix}_{scope}",
                    top_ks=top_ks,
                )
            )
            summary.update(
                _neighborhood_summary(
                    centered_similarities[mask],
                    positive_counts[mask],
                    negative_counts[mask],
                    prefix=f"{prefix}_{scope}_centered",
                    top_ks=top_ks,
                )
            )
        rows.append(summary)

    return pd.concat([out, pd.DataFrame(rows, index=out.index)], axis=1)


def add_screened_recognition_features_file(
    input_path: str | Path,
    reference_path: str | Path,
    embedding_cache_path: str | Path,
    output_path: str | Path,
    *,
    top_ks: tuple[int, ...] = (1, 3, 5, 10, 20),
) -> Path:
    frame = pd.read_csv(input_path)
    reference = pd.read_csv(reference_path)
    embeddings, _ = load_plm_embedding_cache(embedding_cache_path)
    out = add_screened_recognition_features(
        frame,
        reference,
        embeddings,
        top_ks=top_ks,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output
