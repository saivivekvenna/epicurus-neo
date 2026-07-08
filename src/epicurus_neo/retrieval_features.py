from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from epicurus_neo.schema import add_normalized_columns


def peptide_similarity(left: str, right: str) -> float:
    """Return a simple length-aware peptide similarity in [0, 1]."""
    if not left or not right:
        return float("nan")
    if len(left) == len(right):
        return sum(a == b for a, b in zip(left, right, strict=True)) / len(left)

    shorter, longer = (left, right) if len(left) < len(right) else (right, left)
    best = 0
    for offset in range(len(longer) - len(shorter) + 1):
        matches = sum(
            a == b
            for a, b in zip(shorter, longer[offset : offset + len(shorter)], strict=True)
        )
        best = max(best, matches)
    return best / len(longer)


def _topk_mean(values: list[float], k: int) -> float:
    clean = [value for value in values if not pd.isna(value)]
    if not clean:
        return float("nan")
    return float(np.mean(sorted(clean, reverse=True)[:k]))


def _reference_subset(reference: pd.DataFrame, hla: str) -> pd.DataFrame:
    same_hla = reference[reference["hla_allele_norm"] == hla]
    return same_hla if not same_hla.empty else reference


def _ensure_canonical_minimum(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for column in ["hla_allele", "mutant_peptide", "wildtype_peptide"]:
        if column not in out.columns:
            out[column] = ""
    return out


def add_retrieval_features(
    frame: pd.DataFrame,
    reference: pd.DataFrame,
    *,
    top_k: int = 5,
    exclude_self: bool = True,
) -> pd.DataFrame:
    """Add nearest-neighbor peptide-label features from a labeled reference table."""
    out = add_normalized_columns(_ensure_canonical_minimum(frame))
    ref = add_normalized_columns(_ensure_canonical_minimum(reference))
    ref = ref[ref["label"].isin(["positive", "negative"])].copy()
    if ref.empty:
        return out

    rows: list[dict[str, float]] = []
    for _, row in out.iterrows():
        subset = _reference_subset(ref, str(row["hla_allele_norm"]))
        if exclude_self and "candidate_id" in row and "candidate_id" in subset.columns:
            subset = subset[subset["candidate_id"].astype(str) != str(row["candidate_id"])]

        peptide = str(row["mutant_peptide_norm"])
        pos_sims: list[float] = []
        neg_sims: list[float] = []
        all_pairs: list[tuple[float, str]] = []
        for _, ref_row in subset.iterrows():
            similarity = peptide_similarity(peptide, str(ref_row["mutant_peptide_norm"]))
            if pd.isna(similarity):
                continue
            label = str(ref_row["label"])
            all_pairs.append((similarity, label))
            if label == "positive":
                pos_sims.append(similarity)
            elif label == "negative":
                neg_sims.append(similarity)

        nearest = sorted(all_pairs, reverse=True)[:top_k]
        nearest_pos_fraction = (
            sum(1 for _, label in nearest if label == "positive") / len(nearest)
            if nearest
            else float("nan")
        )
        max_pos = max(pos_sims) if pos_sims else float("nan")
        max_neg = max(neg_sims) if neg_sims else float("nan")
        rows.append(
            {
                "retrieval_max_positive_similarity": max_pos,
                "retrieval_max_negative_similarity": max_neg,
                "retrieval_positive_minus_negative_similarity": max_pos - max_neg,
                "retrieval_topk_positive_similarity_mean": _topk_mean(pos_sims, top_k),
                "retrieval_topk_negative_similarity_mean": _topk_mean(neg_sims, top_k),
                "retrieval_topk_positive_fraction": nearest_pos_fraction,
                "retrieval_reference_count": float(len(subset)),
            }
        )

    feature_frame = pd.DataFrame(rows, index=out.index)
    for column in feature_frame.columns:
        out[column] = feature_frame[column]
    return out


def add_retrieval_features_file(
    input_path: str | Path,
    reference_path: str | Path,
    output_path: str | Path,
    *,
    top_k: int = 5,
) -> Path:
    frame = pd.read_csv(input_path)
    reference = pd.read_csv(reference_path)
    out = add_retrieval_features(frame, reference, top_k=top_k)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output
