from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from epicurus_neo.metrics import group_metrics, summarize_group_metrics
from epicurus_neo.plm_retrieval import load_plm_embedding_cache
from epicurus_neo.schema import normalize_hla, normalize_peptide


DEFAULT_BASE_FEATURES = (
    "mhcflurry_presentation_score",
    "mhcflurry_processing_score",
    "mutation_count",
    "mutation_anchor_count",
    "mutation_tcr_face_count",
    "mutation_hydrophobicity_delta",
    "mutation_charge_delta",
    "retrieval_max_positive_similarity",
    "retrieval_max_negative_similarity",
    "retrieval_positive_minus_negative_similarity",
    "retrieval_topk_positive_similarity_mean",
    "retrieval_topk_negative_similarity_mean",
    "retrieval_topk_positive_fraction",
    "retrieval_biochemical_max_positive_similarity",
    "retrieval_biochemical_max_negative_similarity",
    "retrieval_biochemical_positive_minus_negative_similarity",
    "retrieval_biochemical_topk_positive_similarity_mean",
    "retrieval_biochemical_topk_negative_similarity_mean",
    "retrieval_biochemical_topk_positive_fraction",
    "recognition_hla_max_similarity",
    "recognition_hla_topk_similarity_mean",
    "recognition_hla_max_weighted_similarity",
    "recognition_hla_centered_max_similarity",
    "recognition_hla_centered_topk_similarity_mean",
    "recognition_hla_centered_max_weighted_similarity",
    "recognition_pathogen_max_similarity",
    "recognition_pathogen_topk_similarity_mean",
    "recognition_pathogen_max_weighted_similarity",
    "recognition_pathogen_centered_max_similarity",
    "recognition_pathogen_centered_topk_similarity_mean",
    "recognition_pathogen_centered_max_weighted_similarity",
    "recognition_human_max_similarity",
    "recognition_human_topk_similarity_mean",
    "recognition_human_max_weighted_similarity",
    "recognition_human_centered_max_similarity",
    "recognition_human_centered_topk_similarity_mean",
    "recognition_human_centered_max_weighted_similarity",
    "recognition_pathogen_minus_human_similarity",
    "recognition_centered_pathogen_minus_human_similarity",
    "screened_recognition_hla_centered_positive_minus_negative_similarity",
    "screened_recognition_hla_centered_top10_response_rate",
    "screened_recognition_hla_centered_top10_weighted_response_rate",
    "screened_recognition_family_centered_positive_minus_negative_similarity",
    "screened_recognition_family_centered_top10_response_rate",
    "screened_recognition_family_centered_top10_weighted_response_rate",
    "screened_recognition_global_centered_positive_minus_negative_similarity",
    "screened_recognition_global_centered_top10_response_rate",
    "screened_recognition_global_centered_top10_weighted_response_rate",
)


@dataclass(frozen=True)
class XGBRankerConfig:
    objective: str = "rank:ndcg"
    n_estimators: int = 200
    learning_rate: float = 0.05
    max_depth: int = 3
    min_child_weight: float = 1.0
    subsample: float = 0.9
    colsample_bytree: float = 0.8
    reg_lambda: float = 1.0
    random_state: int = 17


@dataclass(frozen=True)
class FrozenPLMSelection:
    config: XGBRankerConfig
    feature_columns: tuple[str, ...]
    allele_vocabulary: tuple[str, ...]
    validation_summary: dict[str, float]
    candidate_summaries: tuple[dict[str, object], ...]
    model_name: str


def default_ranker_configs() -> tuple[XGBRankerConfig, ...]:
    configs: list[XGBRankerConfig] = []
    for objective in ("rank:ndcg", "rank:pairwise"):
        for max_depth in (2, 3, 4):
            for learning_rate, n_estimators in ((0.03, 300), (0.06, 180)):
                configs.append(
                    XGBRankerConfig(
                        objective=objective,
                        n_estimators=n_estimators,
                        learning_rate=learning_rate,
                        max_depth=max_depth,
                        min_child_weight=2.0 if max_depth >= 4 else 1.0,
                        subsample=0.9,
                        colsample_bytree=0.7,
                    )
                )
    return tuple(configs)


def _embedding_dimension(embeddings: dict[str, np.ndarray]) -> int:
    if not embeddings:
        raise ValueError("PLM embedding cache is empty")
    dimensions = {int(np.asarray(vector).shape[0]) for vector in embeddings.values()}
    if len(dimensions) != 1:
        raise ValueError("PLM embedding vectors have inconsistent dimensions")
    return dimensions.pop()


def _allele_vocabulary(frame: pd.DataFrame) -> tuple[str, ...]:
    return tuple(sorted({normalize_hla(value) for value in frame["hla_allele"] if normalize_hla(value)}))


def build_frozen_plm_features(
    frame: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    *,
    allele_vocabulary: tuple[str, ...],
    base_features: tuple[str, ...] = DEFAULT_BASE_FEATURES,
) -> pd.DataFrame:
    dimension = _embedding_dimension(embeddings)
    rows: list[np.ndarray] = []
    for peptide in frame["mutant_peptide"]:
        vector = embeddings.get(normalize_peptide(peptide))
        if vector is None:
            rows.append(np.zeros(dimension, dtype=np.float32))
        else:
            rows.append(np.asarray(vector, dtype=np.float32))
    matrix = np.stack(rows)
    features = pd.DataFrame(
        matrix,
        index=frame.index,
        columns=[f"plm_embedding_{idx:03d}" for idx in range(dimension)],
    )

    extra_columns: dict[str, pd.Series] = {}
    for column in base_features:
        if column in frame.columns:
            extra_columns[column] = pd.to_numeric(frame[column], errors="coerce")

    if "mhcflurry_affinity" in frame.columns:
        affinity = pd.to_numeric(frame["mhcflurry_affinity"], errors="coerce")
        extra_columns["mhcflurry_affinity_inverse_score"] = -np.log10(
            affinity.where(affinity > 0)
        )
    if "mhcflurry_presentation_percentile" in frame.columns:
        percentile = pd.to_numeric(frame["mhcflurry_presentation_percentile"], errors="coerce")
        extra_columns["mhcflurry_presentation_percentile_inverse_score"] = -percentile

    normalized_alleles = frame["hla_allele"].map(normalize_hla)
    allele_columns: dict[str, pd.Series] = {}
    for allele in allele_vocabulary:
        safe_name = allele.replace("*", "_").replace(":", "_").replace("-", "_")
        allele_columns[f"hla_{safe_name}"] = (normalized_alleles == allele).astype(np.float32)
    return pd.concat(
        [
            features,
            pd.DataFrame(extra_columns, index=frame.index),
            pd.DataFrame(allele_columns, index=frame.index),
        ],
        axis=1,
    )


def _sorted_rank_data(
    frame: pd.DataFrame,
    features: pd.DataFrame,
    *,
    group_col: str,
) -> tuple[pd.DataFrame, np.ndarray, list[int], pd.Index]:
    labeled = frame["label"].isin(["positive", "negative"])
    work = frame.loc[labeled].copy()
    work["__row_index"] = work.index
    work = work.sort_values(group_col, kind="mergesort")
    index = pd.Index(work["__row_index"])
    x = features.loc[index]
    y = (work["label"] == "positive").astype(int).to_numpy()
    groups = work.groupby(group_col, sort=False).size().astype(int).tolist()
    return x, y, groups, index


def _make_ranker(config: XGBRankerConfig):
    try:
        from xgboost import XGBRanker
    except ImportError as exc:
        raise ImportError("Frozen PLM ranking requires xgboost") from exc

    return XGBRanker(
        objective=config.objective,
        n_estimators=config.n_estimators,
        learning_rate=config.learning_rate,
        max_depth=config.max_depth,
        min_child_weight=config.min_child_weight,
        subsample=config.subsample,
        colsample_bytree=config.colsample_bytree,
        reg_lambda=config.reg_lambda,
        tree_method="hist",
        random_state=config.random_state,
        n_jobs=-1,
    )


def _fit_ranker(
    frame: pd.DataFrame,
    features: pd.DataFrame,
    *,
    group_col: str,
    config: XGBRankerConfig,
):
    x, y, groups, _ = _sorted_rank_data(frame, features, group_col=group_col)
    model = _make_ranker(config)
    model.fit(x, y, group=groups, verbose=False)
    return model


def _score_frame(
    frame: pd.DataFrame,
    features: pd.DataFrame,
    model: object,
    *,
    output_col: str = "epicurus_frozen_plm_score",
) -> pd.DataFrame:
    out = frame.copy()
    out[output_col] = model.predict(features)
    return out


def _selection_key(summary: dict[str, float]) -> tuple[float, float, float, float, float]:
    return (
        summary["mean_hits_at_k"],
        summary["mean_precision_at_k"],
        summary["mean_recall_at_k"],
        summary["mean_ndcg_at_k"],
        summary["mean_mrr"],
    )


def select_frozen_plm_ranker(
    train: pd.DataFrame,
    validation: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    *,
    model_name: str,
    group_col: str = "hla_allele",
    k: int = 20,
    configs: tuple[XGBRankerConfig, ...] | None = None,
) -> FrozenPLMSelection:
    configs = configs or default_ranker_configs()
    if not configs:
        raise ValueError("At least one XGBRanker configuration is required")

    allele_vocabulary = _allele_vocabulary(train)
    train_features = build_frozen_plm_features(
        train,
        embeddings,
        allele_vocabulary=allele_vocabulary,
    )
    validation_features = build_frozen_plm_features(
        validation,
        embeddings,
        allele_vocabulary=allele_vocabulary,
    )

    candidates: list[dict[str, object]] = []
    best: tuple[tuple[float, float, float, float, float], XGBRankerConfig, dict[str, float]] | None = None
    for config in configs:
        model = _fit_ranker(train, train_features, group_col=group_col, config=config)
        scored = _score_frame(validation, validation_features, model)
        summary = summarize_group_metrics(
            group_metrics(
                scored,
                group_col=group_col,
                score_col="epicurus_frozen_plm_score",
                k=k,
            )
        )
        candidates.append({"config": asdict(config), "summary": summary})
        key = _selection_key(summary)
        if best is None or key > best[0]:
            best = (key, config, summary)

    assert best is not None
    _, selected_config, selected_summary = best
    return FrozenPLMSelection(
        config=selected_config,
        feature_columns=tuple(train_features.columns),
        allele_vocabulary=allele_vocabulary,
        validation_summary=selected_summary,
        candidate_summaries=tuple(candidates),
        model_name=model_name,
    )


def refit_and_score_frozen_plm_ranker(
    train_validation: pd.DataFrame,
    target: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    selection: FrozenPLMSelection,
    *,
    group_col: str = "hla_allele",
) -> pd.DataFrame:
    train_features = build_frozen_plm_features(
        train_validation,
        embeddings,
        allele_vocabulary=selection.allele_vocabulary,
    ).reindex(columns=list(selection.feature_columns))
    target_features = build_frozen_plm_features(
        target,
        embeddings,
        allele_vocabulary=selection.allele_vocabulary,
    ).reindex(columns=list(selection.feature_columns))
    model = _fit_ranker(
        train_validation,
        train_features,
        group_col=group_col,
        config=selection.config,
    )
    return _score_frame(target, target_features, model)


def run_frozen_plm_ranker_files(
    train_path: str | Path,
    validation_path: str | Path,
    train_validation_path: str | Path,
    target_path: str | Path,
    embedding_cache_path: str | Path,
    scored_output_path: str | Path,
    selection_output_path: str | Path,
    *,
    group_col: str = "hla_allele",
    k: int = 20,
) -> tuple[Path, Path, FrozenPLMSelection]:
    train = pd.read_csv(train_path)
    validation = pd.read_csv(validation_path)
    train_validation = pd.read_csv(train_validation_path)
    target = pd.read_csv(target_path)
    embeddings, model_name = load_plm_embedding_cache(embedding_cache_path)

    selection = select_frozen_plm_ranker(
        train,
        validation,
        embeddings,
        model_name=model_name,
        group_col=group_col,
        k=k,
    )
    scored = refit_and_score_frozen_plm_ranker(
        train_validation,
        target,
        embeddings,
        selection,
        group_col=group_col,
    )

    scored_output = Path(scored_output_path)
    scored_output.parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(scored_output, index=False)

    selection_output = Path(selection_output_path)
    selection_output.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "model_name": selection.model_name,
        "config": asdict(selection.config),
        "feature_columns": list(selection.feature_columns),
        "allele_vocabulary": list(selection.allele_vocabulary),
        "validation_summary": selection.validation_summary,
        "candidate_summaries": list(selection.candidate_summaries),
    }
    selection_output.write_text(json.dumps(payload, indent=2) + "\n")
    return scored_output, selection_output, selection
