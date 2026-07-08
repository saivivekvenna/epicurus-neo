from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from epicurus_neo.features import CHARGE, HYDROPHOBICITY
from epicurus_neo.schema import add_normalized_columns


AA_GROUPS = [
    set("AVLIM"),
    set("FYW"),
    set("STNQ"),
    set("KRH"),
    set("DE"),
    set("CGP"),
]


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


def residue_biochemical_similarity(left: str, right: str) -> float:
    if left == right:
        return 1.0
    hydrophobicity_distance = abs(HYDROPHOBICITY.get(left, 0.0) - HYDROPHOBICITY.get(right, 0.0))
    charge_distance = abs(CHARGE.get(left, 0.0) - CHARGE.get(right, 0.0))
    same_group = any(left in group and right in group for group in AA_GROUPS)
    score = 1.0 - 0.10 * hydrophobicity_distance - 0.25 * charge_distance
    if same_group:
        score += 0.15
    return float(min(1.0, max(0.0, score)))


def peptide_biochemical_similarity(left: str, right: str) -> float:
    """Return length-aware peptide similarity with conservative substitutions."""
    if not left or not right:
        return float("nan")

    def aligned_score(shorter: str, longer: str, offset: int) -> float:
        total = sum(
            residue_biochemical_similarity(a, b)
            for a, b in zip(shorter, longer[offset : offset + len(shorter)], strict=True)
        )
        return total / len(longer)

    if len(left) == len(right):
        return aligned_score(left, right, 0)

    shorter, longer = (left, right) if len(left) < len(right) else (right, left)
    return max(aligned_score(shorter, longer, offset) for offset in range(len(longer) - len(shorter) + 1))


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


def _retrieval_summary(prefix: str, pos_sims: list[float], neg_sims: list[float], top_k: int) -> dict[str, float]:
    all_pairs = [(similarity, "positive") for similarity in pos_sims]
    all_pairs.extend((similarity, "negative") for similarity in neg_sims)
    nearest = sorted(all_pairs, reverse=True)[:top_k]
    nearest_pos_fraction = (
        sum(1 for _, label in nearest if label == "positive") / len(nearest)
        if nearest
        else float("nan")
    )
    max_pos = max(pos_sims) if pos_sims else float("nan")
    max_neg = max(neg_sims) if neg_sims else float("nan")
    return {
        f"{prefix}_max_positive_similarity": max_pos,
        f"{prefix}_max_negative_similarity": max_neg,
        f"{prefix}_positive_minus_negative_similarity": max_pos - max_neg,
        f"{prefix}_topk_positive_similarity_mean": _topk_mean(pos_sims, top_k),
        f"{prefix}_topk_negative_similarity_mean": _topk_mean(neg_sims, top_k),
        f"{prefix}_topk_positive_fraction": nearest_pos_fraction,
    }


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
        exact_pos_sims: list[float] = []
        exact_neg_sims: list[float] = []
        biochemical_pos_sims: list[float] = []
        biochemical_neg_sims: list[float] = []
        for _, ref_row in subset.iterrows():
            ref_peptide = str(ref_row["mutant_peptide_norm"])
            exact_similarity = peptide_similarity(peptide, ref_peptide)
            biochemical_similarity = peptide_biochemical_similarity(peptide, ref_peptide)
            if pd.isna(exact_similarity):
                continue
            label = str(ref_row["label"])
            if label == "positive":
                exact_pos_sims.append(exact_similarity)
                biochemical_pos_sims.append(biochemical_similarity)
            elif label == "negative":
                exact_neg_sims.append(exact_similarity)
                biochemical_neg_sims.append(biochemical_similarity)

        summary = _retrieval_summary("retrieval", exact_pos_sims, exact_neg_sims, top_k)
        summary.update(
            _retrieval_summary(
                "retrieval_biochemical",
                biochemical_pos_sims,
                biochemical_neg_sims,
                top_k,
            )
        )
        summary["retrieval_reference_count"] = float(len(subset))
        rows.append(summary)

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
