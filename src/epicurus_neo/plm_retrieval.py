from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from epicurus_neo.retrieval_features import _retrieval_summary, embedding_cosine_similarity
from epicurus_neo.schema import add_normalized_columns


def _mean_pool_embeddings(hidden_states: object, attention_mask: object) -> np.ndarray:
    mask = attention_mask.unsqueeze(-1).to(hidden_states.dtype)
    masked = hidden_states * mask
    denom = mask.sum(dim=1).clamp(min=1)
    pooled = masked.sum(dim=1) / denom
    return pooled.detach().cpu().numpy()


def _ensure_canonical_minimum(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    for column in ["hla_allele", "mutant_peptide", "wildtype_peptide"]:
        if column not in out.columns:
            out[column] = ""
    return out


def compute_plm_embeddings(
    peptides: list[str],
    *,
    model_name: str = "facebook/esm2_t6_8M_UR50D",
    batch_size: int = 64,
    device: str | None = None,
) -> dict[str, np.ndarray]:
    try:
        import torch
        from transformers import AutoModel, AutoTokenizer
    except ImportError as exc:
        raise ImportError(
            "PLM retrieval requires optional dependencies: transformers and torch"
        ) from exc

    unique = sorted({str(peptide).upper() for peptide in peptides if str(peptide)})
    if not unique:
        return {}
    if device is None:
        device = "mps" if torch.backends.mps.is_available() else "cpu"

    tokenizer = AutoTokenizer.from_pretrained(model_name)
    model = AutoModel.from_pretrained(model_name).to(device).eval()
    embeddings: dict[str, np.ndarray] = {}

    for start in range(0, len(unique), batch_size):
        batch = unique[start : start + batch_size]
        inputs = tokenizer(batch, return_tensors="pt", padding=True)
        inputs = {key: value.to(device) for key, value in inputs.items()}
        with torch.no_grad():
            outputs = model(**inputs)
        pooled = _mean_pool_embeddings(outputs.last_hidden_state, inputs["attention_mask"])
        for peptide, vector in zip(batch, pooled, strict=True):
            norm = np.linalg.norm(vector)
            embeddings[peptide] = vector / norm if norm else vector
    return embeddings


def save_plm_embedding_cache(
    embeddings: dict[str, np.ndarray],
    output_path: str | Path,
    *,
    model_name: str,
) -> Path:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    peptides = sorted(embeddings)
    matrix = (
        np.stack([embeddings[peptide] for peptide in peptides])
        if peptides
        else np.empty((0, 0), dtype=np.float32)
    )
    np.savez_compressed(
        output,
        peptides=np.asarray(peptides),
        embeddings=matrix.astype(np.float32),
        model_name=np.asarray(model_name),
    )
    return output


def load_plm_embedding_cache(path: str | Path) -> tuple[dict[str, np.ndarray], str]:
    with np.load(path, allow_pickle=False) as payload:
        peptides = payload["peptides"].astype(str).tolist()
        matrix = payload["embeddings"]
        model_name = str(payload["model_name"].item())
    embeddings = {
        peptide: vector.astype(np.float32, copy=False)
        for peptide, vector in zip(peptides, matrix, strict=True)
    }
    return embeddings, model_name


def build_plm_embedding_cache(
    input_paths: list[str | Path],
    output_path: str | Path,
    *,
    model_name: str = "facebook/esm2_t6_8M_UR50D",
    batch_size: int = 64,
    device: str | None = None,
) -> Path:
    peptides: list[str] = []
    for input_path in input_paths:
        frame = pd.read_csv(input_path, usecols=["mutant_peptide"])
        peptides.extend(frame["mutant_peptide"].astype(str).tolist())
    embeddings = compute_plm_embeddings(
        peptides,
        model_name=model_name,
        batch_size=batch_size,
        device=device,
    )
    return save_plm_embedding_cache(embeddings, output_path, model_name=model_name)


def add_embedding_retrieval_features(
    frame: pd.DataFrame,
    reference: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    *,
    top_k: int = 5,
    prefix: str = "retrieval_plm",
    exclude_self: bool = True,
) -> pd.DataFrame:
    out = add_normalized_columns(_ensure_canonical_minimum(frame))
    ref = add_normalized_columns(_ensure_canonical_minimum(reference))
    ref = ref[ref["label"].isin(["positive", "negative"])].copy()
    if ref.empty:
        return out

    rows: list[dict[str, float]] = []
    for _, row in out.iterrows():
        hla = str(row.get("hla_allele_norm", ""))
        subset = ref[ref["hla_allele_norm"] == hla]
        if subset.empty:
            subset = ref
        if exclude_self and "candidate_id" in row and "candidate_id" in subset.columns:
            subset = subset[subset["candidate_id"].astype(str) != str(row["candidate_id"])]

        peptide = str(row["mutant_peptide_norm"])
        query_embedding = embeddings.get(peptide)
        pos_sims: list[float] = []
        neg_sims: list[float] = []
        pos_embeddings: list[np.ndarray] = []
        neg_embeddings: list[np.ndarray] = []
        if query_embedding is not None:
            for _, ref_row in subset.iterrows():
                ref_embedding = embeddings.get(str(ref_row["mutant_peptide_norm"]))
                if ref_embedding is None:
                    continue
                similarity = embedding_cosine_similarity(query_embedding, ref_embedding)
                if str(ref_row["label"]) == "positive":
                    pos_sims.append(similarity)
                    pos_embeddings.append(ref_embedding)
                elif str(ref_row["label"]) == "negative":
                    neg_sims.append(similarity)
                    neg_embeddings.append(ref_embedding)

        summary = _retrieval_summary(prefix, pos_sims, neg_sims, top_k)
        pos_prototype = np.mean(pos_embeddings, axis=0) if pos_embeddings else None
        neg_prototype = np.mean(neg_embeddings, axis=0) if neg_embeddings else None
        pos_proto_sim = (
            embedding_cosine_similarity(query_embedding, pos_prototype)
            if query_embedding is not None and pos_prototype is not None
            else float("nan")
        )
        neg_proto_sim = (
            embedding_cosine_similarity(query_embedding, neg_prototype)
            if query_embedding is not None and neg_prototype is not None
            else float("nan")
        )
        summary[f"{prefix}_positive_prototype_similarity"] = pos_proto_sim
        summary[f"{prefix}_negative_prototype_similarity"] = neg_proto_sim
        summary[f"{prefix}_positive_minus_negative_prototype_similarity"] = (
            pos_proto_sim - neg_proto_sim
        )
        summary[f"{prefix}_reference_count"] = float(len(subset))
        rows.append(summary)

    feature_frame = pd.DataFrame(rows, index=out.index)
    for column in feature_frame.columns:
        out[column] = feature_frame[column]
    return out


def add_plm_retrieval_features_file(
    input_path: str | Path,
    reference_path: str | Path,
    output_path: str | Path,
    *,
    model_name: str = "facebook/esm2_t6_8M_UR50D",
    batch_size: int = 64,
    device: str | None = None,
    top_k: int = 5,
) -> Path:
    frame = pd.read_csv(input_path)
    reference = pd.read_csv(reference_path)
    peptides = pd.concat(
        [frame["mutant_peptide"].astype(str), reference["mutant_peptide"].astype(str)],
        ignore_index=True,
    ).tolist()
    embeddings = compute_plm_embeddings(
        peptides,
        model_name=model_name,
        batch_size=batch_size,
        device=device,
    )
    out = add_embedding_retrieval_features(frame, reference, embeddings, top_k=top_k)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(output, index=False)
    return output
