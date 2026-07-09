from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from epicurus_neo.metrics import group_metrics, summarize_group_metrics
from epicurus_neo.schema import normalize_hla, normalize_peptide


NUMERIC_FEATURES = (
    "mhcflurry_affinity_inverse_score",
    "mhcflurry_processing_score",
    "mhcflurry_presentation_score",
    "mhcflurry_presentation_percentile_inverse_score",
)


@dataclass(frozen=True)
class PLMFineTuneConfig:
    strategy: str
    pretrain_epochs: int = 3
    target_epochs: int = 4
    encoder_learning_rate: float = 1e-5
    head_learning_rate: float = 1e-3
    weight_decay: float = 1e-4
    negatives_per_positive: int = 2
    hard_negative_pool: int = 100
    batch_pairs: int = 32
    random_state: int = 17


@dataclass(frozen=True)
class PLMFineTuneSelection:
    config: PLMFineTuneConfig
    selected_target_epoch: int
    validation_summary: dict[str, float]
    candidate_summaries: tuple[dict[str, object], ...]
    model_name: str


def default_finetune_configs() -> tuple[PLMFineTuneConfig, ...]:
    return (
        PLMFineTuneConfig(strategy="target_only"),
        PLMFineTuneConfig(strategy="external_pretrain"),
    )


def _numeric_frame(frame: pd.DataFrame) -> pd.DataFrame:
    out = pd.DataFrame(index=frame.index)
    missing = pd.Series(np.nan, index=frame.index)
    affinity = pd.to_numeric(
        frame.get("mhcflurry_affinity", missing),
        errors="coerce",
    )
    out["mhcflurry_affinity_inverse_score"] = -np.log10(affinity.where(affinity > 0))
    out["mhcflurry_processing_score"] = pd.to_numeric(
        frame.get("mhcflurry_processing_score", missing), errors="coerce"
    )
    out["mhcflurry_presentation_score"] = pd.to_numeric(
        frame.get("mhcflurry_presentation_score", missing), errors="coerce"
    )
    percentile = pd.to_numeric(
        frame.get("mhcflurry_presentation_percentile", missing), errors="coerce"
    )
    out["mhcflurry_presentation_percentile_inverse_score"] = -percentile
    return out


def fit_numeric_stats(*frames: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    matrix = pd.concat([_numeric_frame(frame) for frame in frames], ignore_index=True)
    means = matrix.mean(skipna=True).fillna(0.0).to_numpy(dtype=np.float32)
    stds = matrix.std(skipna=True).replace(0.0, 1.0).fillna(1.0).to_numpy(dtype=np.float32)
    return means, stds


def numeric_matrix(
    frame: pd.DataFrame,
    means: np.ndarray,
    stds: np.ndarray,
) -> np.ndarray:
    values = _numeric_frame(frame).to_numpy(dtype=np.float32)
    values = np.where(np.isnan(values), means, values)
    return ((values - means) / stds).astype(np.float32)


def make_hard_pairs(
    frame: pd.DataFrame,
    *,
    group_col: str,
    negatives_per_positive: int,
    hard_negative_pool: int,
    random_state: int,
    hard_score_col: str = "mhcflurry_presentation_score",
) -> tuple[np.ndarray, np.ndarray]:
    if negatives_per_positive < 1:
        raise ValueError("negatives_per_positive must be at least 1")
    if hard_negative_pool < 1:
        raise ValueError("hard_negative_pool must be at least 1")

    rng = np.random.default_rng(random_state)
    positive_indices: list[int] = []
    negative_indices: list[int] = []
    for _, group in frame.groupby(group_col, sort=False):
        positives = group.index[group["label"] == "positive"].to_numpy()
        negatives = group[group["label"] == "negative"].copy()
        if len(positives) == 0 or negatives.empty:
            continue
        if hard_score_col in negatives:
            negatives["__hard_score"] = pd.to_numeric(
                negatives[hard_score_col], errors="coerce"
            ).fillna(-np.inf)
            negatives = negatives.sort_values("__hard_score", ascending=False)
        pool = negatives.index.to_numpy()[:hard_negative_pool]
        for positive_index in positives:
            sampled = rng.choice(
                pool,
                size=negatives_per_positive,
                replace=len(pool) < negatives_per_positive,
            )
            positive_indices.extend([int(positive_index)] * negatives_per_positive)
            negative_indices.extend(int(index) for index in sampled)
    if not positive_indices:
        raise ValueError("No within-group positive/negative pairs are available")
    order = rng.permutation(len(positive_indices))
    return (
        np.asarray(positive_indices, dtype=int)[order],
        np.asarray(negative_indices, dtype=int)[order],
    )


def _device_name(requested: str | None) -> str:
    import torch

    if requested:
        return requested
    return "mps" if torch.backends.mps.is_available() else "cpu"


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


def _allele_indices(frame: pd.DataFrame, vocabulary: tuple[str, ...]) -> np.ndarray:
    mapping = {allele: index + 1 for index, allele in enumerate(vocabulary)}
    return frame["hla_allele"].map(normalize_hla).map(mapping).fillna(0).to_numpy(dtype=np.int64)


def _make_model(model_name: str, allele_count: int, device: str):
    import torch
    from transformers import AutoModel

    class PeptideHLARanker(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.encoder = AutoModel.from_pretrained(model_name)
            hidden_size = int(self.encoder.config.hidden_size)
            allele_dim = 16
            self.allele_embedding = torch.nn.Embedding(allele_count + 1, allele_dim)
            self.head = torch.nn.Sequential(
                torch.nn.LayerNorm(hidden_size + allele_dim + len(NUMERIC_FEATURES)),
                torch.nn.Linear(hidden_size + allele_dim + len(NUMERIC_FEATURES), 128),
                torch.nn.GELU(),
                torch.nn.Dropout(0.15),
                torch.nn.Linear(128, 1),
            )

        def forward(self, inputs, allele_index, numeric):
            special_tokens_mask = inputs.pop("special_tokens_mask")
            outputs = self.encoder(**inputs)
            residue_mask = inputs["attention_mask"] * (1 - special_tokens_mask)
            mask = residue_mask.unsqueeze(-1).to(outputs.last_hidden_state.dtype)
            pooled = (outputs.last_hidden_state * mask).sum(dim=1) / mask.sum(dim=1).clamp(
                min=1
            )
            allele = self.allele_embedding(allele_index)
            return self.head(torch.cat([pooled, allele, numeric], dim=1)).squeeze(1)

    return PeptideHLARanker().to(device)


def _batch_scores(
    model: Any,
    tokenizer: Any,
    peptides: list[str],
    allele_indices: np.ndarray,
    numeric: np.ndarray,
    *,
    device: str,
):
    import torch

    inputs = tokenizer(
        peptides,
        return_tensors="pt",
        padding=True,
        return_special_tokens_mask=True,
    )
    inputs = {key: value.to(device) for key, value in inputs.items()}
    allele_tensor = torch.as_tensor(allele_indices, dtype=torch.long, device=device)
    numeric_tensor = torch.as_tensor(numeric, dtype=torch.float32, device=device)
    return model(inputs, allele_tensor, numeric_tensor)


def _train_epoch(
    model: Any,
    tokenizer: Any,
    frame: pd.DataFrame,
    positive_indices: np.ndarray,
    negative_indices: np.ndarray,
    allele_indices: np.ndarray,
    numeric: np.ndarray,
    optimizer: Any,
    *,
    batch_pairs: int,
    device: str,
) -> float:
    import torch

    model.train()
    losses: list[float] = []
    peptides = frame["mutant_peptide"].map(normalize_peptide)
    for start in range(0, len(positive_indices), batch_pairs):
        positive = positive_indices[start : start + batch_pairs]
        negative = negative_indices[start : start + batch_pairs]
        indices = np.concatenate([positive, negative])
        scores = _batch_scores(
            model,
            tokenizer,
            peptides.loc[indices].tolist(),
            allele_indices[indices],
            numeric[indices],
            device=device,
        )
        pair_count = len(positive)
        positive_scores = scores[:pair_count]
        negative_scores = scores[pair_count:]
        ranking_loss = torch.nn.functional.softplus(
            negative_scores - positive_scores
        ).mean()
        calibration_loss = 0.5 * (
            torch.nn.functional.binary_cross_entropy_with_logits(
                positive_scores, torch.ones_like(positive_scores)
            )
            + torch.nn.functional.binary_cross_entropy_with_logits(
                negative_scores, torch.zeros_like(negative_scores)
            )
        )
        loss = ranking_loss + 0.1 * calibration_loss
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        optimizer.step()
        losses.append(float(loss.detach().cpu()))
    return float(np.mean(losses))


def _predict(
    model: Any,
    tokenizer: Any,
    frame: pd.DataFrame,
    allele_indices: np.ndarray,
    numeric: np.ndarray,
    *,
    batch_size: int,
    device: str,
) -> np.ndarray:
    import torch

    model.eval()
    peptides = frame["mutant_peptide"].map(normalize_peptide).tolist()
    predictions: list[np.ndarray] = []
    with torch.no_grad():
        for start in range(0, len(frame), batch_size):
            stop = start + batch_size
            scores = _batch_scores(
                model,
                tokenizer,
                peptides[start:stop],
                allele_indices[start:stop],
                numeric[start:stop],
                device=device,
            )
            predictions.append(scores.detach().cpu().numpy())
    return np.concatenate(predictions)


def _optimizer(model: Any, config: PLMFineTuneConfig):
    import torch

    encoder_parameters = list(model.encoder.parameters())
    head_parameters = [
        parameter
        for name, parameter in model.named_parameters()
        if not name.startswith("encoder.")
    ]
    return torch.optim.AdamW(
        [
            {"params": encoder_parameters, "lr": config.encoder_learning_rate},
            {"params": head_parameters, "lr": config.head_learning_rate},
        ],
        weight_decay=config.weight_decay,
    )


def _summary(
    validation: pd.DataFrame,
    predictions: np.ndarray,
    *,
    group_col: str,
    k: int,
) -> dict[str, float]:
    scored = validation.copy()
    scored["epicurus_finetuned_plm_score"] = predictions
    return summarize_group_metrics(
        group_metrics(
            scored,
            group_col=group_col,
            score_col="epicurus_finetuned_plm_score",
            k=k,
        )
    )


def _selection_key(summary: dict[str, float]) -> tuple[float, float, float, float, float]:
    return (
        summary["mean_hits_at_k"],
        summary["mean_precision_at_k"],
        summary["mean_recall_at_k"],
        summary["mean_ndcg_at_k"],
        summary["mean_mrr"],
    )


def select_finetuned_plm_ranker(
    external: pd.DataFrame,
    target_train: pd.DataFrame,
    validation: pd.DataFrame,
    *,
    model_name: str = "facebook/esm2_t6_8M_UR50D",
    external_group_col: str = "hla_allele",
    target_group_col: str = "hla_allele",
    k: int = 20,
    configs: tuple[PLMFineTuneConfig, ...] | None = None,
    device: str | None = None,
) -> tuple[PLMFineTuneSelection, pd.DataFrame]:
    import torch
    from transformers import AutoTokenizer

    configs = configs or default_finetune_configs()
    device = _device_name(device)
    vocabulary = _allele_vocabulary(external, target_train)
    means, stds = fit_numeric_stats(external, target_train)
    external_numeric = numeric_matrix(external, means, stds)
    target_numeric = numeric_matrix(target_train, means, stds)
    validation_numeric = numeric_matrix(validation, means, stds)
    external_alleles = _allele_indices(external, vocabulary)
    target_alleles = _allele_indices(target_train, vocabulary)
    validation_alleles = _allele_indices(validation, vocabulary)
    tokenizer = AutoTokenizer.from_pretrained(model_name)

    candidates: list[dict[str, object]] = []
    best: tuple[
        tuple[float, float, float, float, float],
        PLMFineTuneConfig,
        int,
        dict[str, float],
        np.ndarray,
    ] | None = None
    for config in configs:
        torch.manual_seed(config.random_state)
        np.random.seed(config.random_state)
        model = _make_model(model_name, len(vocabulary), device)
        optimizer = _optimizer(model, config)
        if config.strategy == "external_pretrain":
            for epoch in range(config.pretrain_epochs):
                positive, negative = make_hard_pairs(
                    external,
                    group_col=external_group_col,
                    negatives_per_positive=config.negatives_per_positive,
                    hard_negative_pool=config.hard_negative_pool,
                    random_state=config.random_state + epoch,
                )
                _train_epoch(
                    model,
                    tokenizer,
                    external,
                    positive,
                    negative,
                    external_alleles,
                    external_numeric,
                    optimizer,
                    batch_pairs=config.batch_pairs,
                    device=device,
                )
        elif config.strategy != "target_only":
            raise ValueError(f"Unsupported fine-tuning strategy: {config.strategy}")

        optimizer = _optimizer(model, config)
        for epoch in range(1, config.target_epochs + 1):
            positive, negative = make_hard_pairs(
                target_train,
                group_col=target_group_col,
                negatives_per_positive=config.negatives_per_positive,
                hard_negative_pool=config.hard_negative_pool,
                random_state=config.random_state + 100 + epoch,
            )
            loss = _train_epoch(
                model,
                tokenizer,
                target_train,
                positive,
                negative,
                target_alleles,
                target_numeric,
                optimizer,
                batch_pairs=config.batch_pairs,
                device=device,
            )
            predictions = _predict(
                model,
                tokenizer,
                validation,
                validation_alleles,
                validation_numeric,
                batch_size=config.batch_pairs * 2,
                device=device,
            )
            summary = _summary(
                validation,
                predictions,
                group_col=target_group_col,
                k=k,
            )
            candidates.append(
                {
                    "config": asdict(config),
                    "target_epoch": epoch,
                    "training_loss": loss,
                    "summary": summary,
                }
            )
            key = _selection_key(summary)
            if best is None or key > best[0]:
                best = (key, config, epoch, summary, predictions.copy())

        del model
        if device == "mps":
            torch.mps.empty_cache()

    assert best is not None
    _, config, epoch, summary, predictions = best
    selection = PLMFineTuneSelection(
        config=config,
        selected_target_epoch=epoch,
        validation_summary=summary,
        candidate_summaries=tuple(candidates),
        model_name=model_name,
    )
    scored = validation.copy()
    scored["epicurus_finetuned_plm_score"] = predictions
    return selection, scored


def run_finetuned_plm_ranker_files(
    external_path: str | Path,
    train_path: str | Path,
    validation_path: str | Path,
    validation_output_path: str | Path,
    selection_output_path: str | Path,
    *,
    model_name: str = "facebook/esm2_t6_8M_UR50D",
    external_group_col: str = "hla_allele",
    target_group_col: str = "hla_allele",
    k: int = 20,
    device: str | None = None,
) -> tuple[Path, Path, PLMFineTuneSelection]:
    external = pd.read_csv(external_path)
    train = pd.read_csv(train_path)
    validation = pd.read_csv(validation_path)
    selection, scored = select_finetuned_plm_ranker(
        external,
        train,
        validation,
        model_name=model_name,
        external_group_col=external_group_col,
        target_group_col=target_group_col,
        k=k,
        device=device,
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
                "selected_target_epoch": selection.selected_target_epoch,
                "validation_summary": selection.validation_summary,
                "candidate_summaries": list(selection.candidate_summaries),
            },
            indent=2,
        )
        + "\n"
    )
    return validation_output, selection_output, selection
