from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from epicurus_neo.frozen_plm_ranker import build_frozen_plm_features
from epicurus_neo.leakage import detect_exact_leakage
from epicurus_neo.metrics import group_metrics, summarize_group_metrics
from epicurus_neo.plm_retrieval import load_plm_embedding_cache
from epicurus_neo.schema import normalize_hla


COMMON_TRANSFER_FEATURES = (
    "mhcflurry_affinity",
    "mhcflurry_processing_score",
    "mhcflurry_presentation_score",
    "mhcflurry_presentation_percentile",
    "mutation_count",
    "mutation_anchor_count",
    "mutation_tcr_face_count",
    "mutation_hydrophobicity_delta",
    "mutation_charge_delta",
)


@dataclass(frozen=True)
class TransferRankerConfig:
    strategy: str
    objective: str = "rank:ndcg"
    pretrain_estimators: int = 180
    target_estimators: int = 180
    learning_rate: float = 0.05
    max_depth: int = 2
    random_state: int = 17


@dataclass(frozen=True)
class TransferRankerSelection:
    config: TransferRankerConfig
    feature_columns: tuple[str, ...]
    allele_vocabulary: tuple[str, ...]
    validation_summary: dict[str, float]
    candidate_summaries: tuple[dict[str, object], ...]
    model_name: str


def default_transfer_configs() -> tuple[TransferRankerConfig, ...]:
    configs: list[TransferRankerConfig] = []
    for strategy in (
        "external_only",
        "target_only",
        "pooled",
        "pretrain_finetune",
        "teacher",
    ):
        for max_depth in (2, 3):
            if strategy == "pretrain_finetune":
                configs.append(
                    TransferRankerConfig(
                        strategy=strategy,
                        pretrain_estimators=180,
                        target_estimators=60,
                        max_depth=max_depth,
                    )
                )
            elif strategy == "teacher":
                configs.append(
                    TransferRankerConfig(
                        strategy=strategy,
                        pretrain_estimators=180,
                        target_estimators=180,
                        max_depth=max_depth,
                    )
                )
            else:
                configs.append(
                    TransferRankerConfig(
                        strategy=strategy,
                        pretrain_estimators=180,
                        target_estimators=180,
                        max_depth=max_depth,
                    )
                )
    return tuple(configs)


def _allele_vocabulary(*frames: pd.DataFrame) -> tuple[str, ...]:
    return tuple(
        sorted(
            {
                allele
                for frame in frames
                for allele in frame["hla_allele"].map(normalize_hla)
                if allele
            }
        )
    )


def _features(
    frame: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    allele_vocabulary: tuple[str, ...],
) -> pd.DataFrame:
    return build_frozen_plm_features(
        frame,
        embeddings,
        allele_vocabulary=allele_vocabulary,
        base_features=COMMON_TRANSFER_FEATURES,
    )


def _sorted_rank_data(
    frame: pd.DataFrame,
    features: pd.DataFrame,
    *,
    group_col: str,
) -> tuple[pd.DataFrame, np.ndarray, list[int]]:
    labeled = frame["label"].isin(["positive", "negative"])
    work = frame.loc[labeled].copy()
    work["__row_index"] = work.index
    work = work.sort_values(group_col, kind="mergesort")
    index = pd.Index(work["__row_index"])
    x = features.loc[index]
    y = (work["label"] == "positive").astype(int).to_numpy()
    groups = work.groupby(group_col, sort=False).size().astype(int).tolist()
    return x, y, groups


def _make_ranker(config: TransferRankerConfig, *, n_estimators: int):
    try:
        from xgboost import XGBRanker
    except ImportError as exc:
        raise ImportError("Transfer ranking requires the optional xgboost dependency") from exc

    return XGBRanker(
        objective=config.objective,
        n_estimators=n_estimators,
        learning_rate=config.learning_rate,
        max_depth=config.max_depth,
        min_child_weight=1.0,
        subsample=0.9,
        colsample_bytree=0.7,
        reg_lambda=1.0,
        tree_method="hist",
        random_state=config.random_state,
        n_jobs=-1,
    )


def _fit(
    frame: pd.DataFrame,
    features: pd.DataFrame,
    *,
    group_col: str,
    config: TransferRankerConfig,
    n_estimators: int,
    initial_model: object | None = None,
):
    x, y, groups = _sorted_rank_data(frame, features, group_col=group_col)
    model = _make_ranker(config, n_estimators=n_estimators)
    kwargs = {}
    if initial_model is not None:
        kwargs["xgb_model"] = initial_model.get_booster()
    model.fit(x, y, group=groups, verbose=False, **kwargs)
    return model


def _pooled_training_data(
    external: pd.DataFrame,
    target: pd.DataFrame,
    external_features: pd.DataFrame,
    target_features: pd.DataFrame,
    *,
    external_group_col: str,
    target_group_col: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    external_work = external.copy()
    target_work = target.copy()
    external_work["__transfer_group"] = (
        "external:" + external_work[external_group_col].astype(str)
    )
    target_work["__transfer_group"] = "target:" + target_work[target_group_col].astype(str)
    pooled_frame = pd.concat([external_work, target_work], ignore_index=True)
    pooled_features = pd.concat(
        [
            external_features.reset_index(drop=True),
            target_features.reset_index(drop=True),
        ],
        ignore_index=True,
    )
    return pooled_frame, pooled_features


def _fit_strategy(
    external: pd.DataFrame,
    target_train: pd.DataFrame,
    external_features: pd.DataFrame,
    target_features: pd.DataFrame,
    *,
    external_group_col: str,
    target_group_col: str,
    config: TransferRankerConfig,
):
    if config.strategy == "external_only":
        return _fit(
            external,
            external_features,
            group_col=external_group_col,
            config=config,
            n_estimators=config.pretrain_estimators,
        ), None
    if config.strategy == "target_only":
        return _fit(
            target_train,
            target_features,
            group_col=target_group_col,
            config=config,
            n_estimators=config.target_estimators,
        ), None
    if config.strategy == "pooled":
        pooled_frame, pooled_features = _pooled_training_data(
            external,
            target_train,
            external_features,
            target_features,
            external_group_col=external_group_col,
            target_group_col=target_group_col,
        )
        return _fit(
            pooled_frame,
            pooled_features,
            group_col="__transfer_group",
            config=config,
            n_estimators=config.target_estimators,
        ), None

    teacher = _fit(
        external,
        external_features,
        group_col=external_group_col,
        config=config,
        n_estimators=config.pretrain_estimators,
    )
    if config.strategy == "pretrain_finetune":
        return _fit(
            target_train,
            target_features,
            group_col=target_group_col,
            config=config,
            n_estimators=config.target_estimators,
            initial_model=teacher,
        ), None
    if config.strategy == "teacher":
        teacher_train = target_features.copy()
        teacher_train["external_teacher_score"] = teacher.predict(target_features)
        student = _fit(
            target_train,
            teacher_train,
            group_col=target_group_col,
            config=config,
            n_estimators=config.target_estimators,
        )
        return student, teacher
    raise ValueError(f"Unsupported transfer strategy: {config.strategy}")


def _score(
    frame: pd.DataFrame,
    features: pd.DataFrame,
    model: object,
    teacher: object | None,
    *,
    output_col: str = "epicurus_transfer_ranker_score",
) -> pd.DataFrame:
    prediction_features = features.copy()
    if teacher is not None:
        prediction_features["external_teacher_score"] = teacher.predict(features)
    out = frame.copy()
    out[output_col] = model.predict(prediction_features)
    return out


def _selection_key(summary: dict[str, float]) -> tuple[float, float, float, float, float]:
    return (
        summary["mean_hits_at_k"],
        summary["mean_precision_at_k"],
        summary["mean_recall_at_k"],
        summary["mean_ndcg_at_k"],
        summary["mean_mrr"],
    )


def select_transfer_ranker(
    external: pd.DataFrame,
    target_train: pd.DataFrame,
    validation: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    *,
    model_name: str,
    external_group_col: str = "patient_id",
    target_group_col: str = "hla_allele",
    k: int = 20,
    configs: tuple[TransferRankerConfig, ...] | None = None,
) -> tuple[TransferRankerSelection, pd.DataFrame]:
    overlap = detect_exact_leakage(external, validation)
    if overlap.shared_mutant_hla or overlap.shared_wildtype_hla:
        raise ValueError(f"External/validation peptide leakage detected: {overlap}")

    configs = configs or default_transfer_configs()
    vocabulary = _allele_vocabulary(external, target_train)
    external_features = _features(external, embeddings, vocabulary)
    target_features = _features(target_train, embeddings, vocabulary)
    validation_features = _features(validation, embeddings, vocabulary)

    candidates: list[dict[str, object]] = []
    best: tuple[
        tuple[float, float, float, float, float],
        TransferRankerConfig,
        dict[str, float],
        pd.DataFrame,
        tuple[str, ...],
    ] | None = None
    for config in configs:
        model, teacher = _fit_strategy(
            external,
            target_train,
            external_features,
            target_features,
            external_group_col=external_group_col,
            target_group_col=target_group_col,
            config=config,
        )
        scored = _score(validation, validation_features, model, teacher)
        summary = summarize_group_metrics(
            group_metrics(
                scored,
                group_col=target_group_col,
                score_col="epicurus_transfer_ranker_score",
                k=k,
            )
        )
        candidates.append({"config": asdict(config), "summary": summary})
        key = _selection_key(summary)
        prediction_features = validation_features.copy()
        if teacher is not None:
            prediction_features["external_teacher_score"] = teacher.predict(validation_features)
        feature_columns = tuple(prediction_features.columns)
        if best is None or key > best[0]:
            best = (key, config, summary, scored, feature_columns)

    assert best is not None
    _, config, summary, scored, feature_columns = best
    selection = TransferRankerSelection(
        config=config,
        feature_columns=feature_columns,
        allele_vocabulary=vocabulary,
        validation_summary=summary,
        candidate_summaries=tuple(candidates),
        model_name=model_name,
    )
    return selection, scored


def refit_and_score_transfer_ranker(
    external: pd.DataFrame,
    target_train_validation: pd.DataFrame,
    target: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    selection: TransferRankerSelection,
    *,
    external_group_col: str = "patient_id",
    target_group_col: str = "hla_allele",
) -> pd.DataFrame:
    overlap = detect_exact_leakage(external, target)
    if overlap.shared_mutant_hla or overlap.shared_wildtype_hla:
        raise ValueError(f"External/target peptide leakage detected: {overlap}")

    external_features = _features(external, embeddings, selection.allele_vocabulary)
    train_features = _features(
        target_train_validation,
        embeddings,
        selection.allele_vocabulary,
    )
    target_features = _features(target, embeddings, selection.allele_vocabulary)
    model, teacher = _fit_strategy(
        external,
        target_train_validation,
        external_features,
        train_features,
        external_group_col=external_group_col,
        target_group_col=target_group_col,
        config=selection.config,
    )
    return _score(target, target_features, model, teacher)


def run_transfer_ranker_files(
    external_path: str | Path,
    train_path: str | Path,
    validation_path: str | Path,
    embedding_cache_path: str | Path,
    validation_output_path: str | Path,
    selection_output_path: str | Path,
    *,
    external_group_col: str = "patient_id",
    target_group_col: str = "hla_allele",
    k: int = 20,
) -> tuple[Path, Path, TransferRankerSelection]:
    external = pd.read_csv(external_path)
    train = pd.read_csv(train_path)
    validation = pd.read_csv(validation_path)
    embeddings, model_name = load_plm_embedding_cache(embedding_cache_path)
    selection, scored = select_transfer_ranker(
        external,
        train,
        validation,
        embeddings,
        model_name=model_name,
        external_group_col=external_group_col,
        target_group_col=target_group_col,
        k=k,
    )

    validation_output = Path(validation_output_path)
    validation_output.parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(validation_output, index=False)
    selection_output = Path(selection_output_path)
    selection_output.parent.mkdir(parents=True, exist_ok=True)
    selection_output.write_text(
        json.dumps(
            {
                "model_name": selection.model_name,
                "config": asdict(selection.config),
                "feature_columns": list(selection.feature_columns),
                "allele_vocabulary": list(selection.allele_vocabulary),
                "validation_summary": selection.validation_summary,
                "candidate_summaries": list(selection.candidate_summaries),
            },
            indent=2,
        )
        + "\n"
    )
    return validation_output, selection_output, selection
