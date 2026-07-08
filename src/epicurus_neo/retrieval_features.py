from __future__ import annotations

import hashlib
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
AA_ALPHABET = tuple("ACDEFGHIKLMNPQRSTVWY")
AA_TO_INDEX = {aa: idx for idx, aa in enumerate(AA_ALPHABET)}


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


def peptide_motif_embedding(peptide: str) -> np.ndarray:
    """Return a compact deterministic peptide embedding for retrieval prototypes.

    The embedding intentionally avoids supervised labels. It combines amino-acid
    composition, terminal residue identity, T-cell-face composition, length, and
    simple biochemical summaries. This gives the retrieval layer a smoother
    neighborhood than exact sequence matching while staying reproducible and
    lightweight.
    """
    normalized = "".join(aa for aa in str(peptide).upper() if aa in AA_TO_INDEX)
    if not normalized:
        return np.zeros(87, dtype=float)

    composition = np.zeros(len(AA_ALPHABET), dtype=float)
    first = np.zeros(len(AA_ALPHABET), dtype=float)
    last = np.zeros(len(AA_ALPHABET), dtype=float)
    tcr_face = np.zeros(len(AA_ALPHABET), dtype=float)

    for aa in normalized:
        composition[AA_TO_INDEX[aa]] += 1.0
    composition /= len(normalized)
    first[AA_TO_INDEX[normalized[0]]] = 1.0
    last[AA_TO_INDEX[normalized[-1]]] = 1.0

    face_residues = normalized[2:-1] if len(normalized) >= 4 else normalized
    for aa in face_residues:
        tcr_face[AA_TO_INDEX[aa]] += 1.0
    if face_residues:
        tcr_face /= len(face_residues)

    hydrophobicity = np.array([HYDROPHOBICITY.get(aa, 0.0) for aa in normalized])
    charge = np.array([CHARGE.get(aa, 0.0) for aa in normalized])
    biochemical = np.array(
        [
            min(len(normalized), 30) / 30.0,
            float(hydrophobicity.mean()),
            float(hydrophobicity.std()),
            float(charge.mean()),
            float(charge.sum() / max(len(normalized), 1)),
            float(sum(aa in "FYW" for aa in normalized) / len(normalized)),
            float(sum(aa in "STNQ" for aa in normalized) / len(normalized)),
        ],
        dtype=float,
    )
    return np.concatenate([composition, first, last, tcr_face, biochemical])


def embedding_cosine_similarity(left: np.ndarray, right: np.ndarray) -> float:
    denom = float(np.linalg.norm(left) * np.linalg.norm(right))
    if denom == 0.0:
        return float("nan")
    return float(np.dot(left, right) / denom)


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


def _fold_key(row: pd.Series) -> str:
    parts = [
        str(row.get("candidate_id", "")),
        str(row.get("hla_allele", "")),
        str(row.get("mutant_peptide", "")),
    ]
    return "|".join(parts)


def _stable_fold(row: pd.Series, n_folds: int) -> int:
    digest = hashlib.sha256(_fold_key(row).encode("utf-8")).hexdigest()
    return int(digest[:12], 16) % n_folds


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


def _prototype_similarity_summary(
    prefix: str,
    query_embedding: np.ndarray,
    pos_embeddings: list[np.ndarray],
    neg_embeddings: list[np.ndarray],
) -> dict[str, float]:
    pos_prototype = np.mean(pos_embeddings, axis=0) if pos_embeddings else None
    neg_prototype = np.mean(neg_embeddings, axis=0) if neg_embeddings else None
    pos_similarity = (
        embedding_cosine_similarity(query_embedding, pos_prototype)
        if pos_prototype is not None
        else float("nan")
    )
    neg_similarity = (
        embedding_cosine_similarity(query_embedding, neg_prototype)
        if neg_prototype is not None
        else float("nan")
    )
    return {
        f"{prefix}_positive_prototype_similarity": pos_similarity,
        f"{prefix}_negative_prototype_similarity": neg_similarity,
        f"{prefix}_positive_minus_negative_prototype_similarity": pos_similarity - neg_similarity,
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
        embedding_pos_sims: list[float] = []
        embedding_neg_sims: list[float] = []
        positive_embeddings: list[np.ndarray] = []
        negative_embeddings: list[np.ndarray] = []
        query_embedding = peptide_motif_embedding(peptide)
        for _, ref_row in subset.iterrows():
            ref_peptide = str(ref_row["mutant_peptide_norm"])
            exact_similarity = peptide_similarity(peptide, ref_peptide)
            biochemical_similarity = peptide_biochemical_similarity(peptide, ref_peptide)
            ref_embedding = peptide_motif_embedding(ref_peptide)
            embedding_similarity = embedding_cosine_similarity(query_embedding, ref_embedding)
            if pd.isna(exact_similarity):
                continue
            label = str(ref_row["label"])
            if label == "positive":
                exact_pos_sims.append(exact_similarity)
                biochemical_pos_sims.append(biochemical_similarity)
                embedding_pos_sims.append(embedding_similarity)
                positive_embeddings.append(ref_embedding)
            elif label == "negative":
                exact_neg_sims.append(exact_similarity)
                biochemical_neg_sims.append(biochemical_similarity)
                embedding_neg_sims.append(embedding_similarity)
                negative_embeddings.append(ref_embedding)

        summary = _retrieval_summary("retrieval", exact_pos_sims, exact_neg_sims, top_k)
        summary.update(
            _retrieval_summary(
                "retrieval_biochemical",
                biochemical_pos_sims,
                biochemical_neg_sims,
                top_k,
            )
        )
        summary.update(
            _retrieval_summary(
                "retrieval_motif",
                embedding_pos_sims,
                embedding_neg_sims,
                top_k,
            )
        )
        summary.update(
            _prototype_similarity_summary(
                "retrieval_motif",
                query_embedding,
                positive_embeddings,
                negative_embeddings,
            )
        )
        summary["retrieval_reference_count"] = float(len(subset))
        rows.append(summary)

    feature_frame = pd.DataFrame(rows, index=out.index)
    for column in feature_frame.columns:
        out[column] = feature_frame[column]
    return out


def add_crossfit_retrieval_features(
    frame: pd.DataFrame,
    *,
    top_k: int = 5,
    n_folds: int = 5,
    fold_col: str = "retrieval_fold",
) -> pd.DataFrame:
    """Add retrieval features using out-of-fold labeled references.

    This is for training rows. Each row retrieves only from rows assigned to
    other deterministic folds, avoiding direct label memorization while keeping
    the same feature schema as `add_retrieval_features`.
    """
    if n_folds < 2:
        raise ValueError("n_folds must be at least 2 for crossfit retrieval features")

    base = add_normalized_columns(_ensure_canonical_minimum(frame))
    if fold_col in base.columns:
        folds = pd.to_numeric(base[fold_col], errors="coerce").fillna(0).astype(int) % n_folds
    else:
        folds = base.apply(lambda row: _stable_fold(row, n_folds), axis=1)

    pieces: list[pd.DataFrame] = []
    for fold_id in range(n_folds):
        query = base[folds == fold_id]
        if query.empty:
            continue
        reference = base[folds != fold_id]
        scored = add_retrieval_features(
            query,
            reference,
            top_k=top_k,
            exclude_self=True,
        )
        scored[fold_col] = fold_id
        pieces.append(scored)

    if not pieces:
        out = base.copy()
        out[fold_col] = folds
        return out
    out = pd.concat(pieces, axis=0).sort_index()
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


def add_crossfit_retrieval_features_file(
    input_path: str | Path,
    output_path: str | Path,
    *,
    top_k: int = 5,
    n_folds: int = 5,
    fold_col: str = "retrieval_fold",
) -> Path:
    frame = pd.read_csv(input_path)
    out = add_crossfit_retrieval_features(
        frame,
        top_k=top_k,
        n_folds=n_folds,
        fold_col=fold_col,
    )
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output
