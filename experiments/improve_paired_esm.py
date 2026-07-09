"""Evaluate paired mutant/wild-type ESM features on official IMPROVE folds."""

from __future__ import annotations

import argparse
import json
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import OneHotEncoder
from xgboost import XGBClassifier

from epicurus_neo.metrics import group_metrics, summarize_group_metrics
from epicurus_neo.plm_retrieval import load_plm_embedding_cache

DEFAULT_MEMBER = (
    "data/03_data_for_CV/IMPROVE/"
    "03_3_final_peptide_features_Partition_for_CV.txt"
)
CATEGORICAL_COLUMNS = ("HLA_allele", "Mutation_Consequence", "IB_CB_cat", "cohort")
EXCLUDED_NUMERIC_COLUMNS = {"response", "validation", "Partition"}


def read_improve(path: Path, member: str) -> pd.DataFrame:
    with zipfile.ZipFile(path) as archive, archive.open(member) as source:
        return pd.read_csv(source, sep="\t")


def safe_numeric_columns(frame: pd.DataFrame) -> list[str]:
    return [
        column
        for column in frame.select_dtypes(include=[np.number]).columns
        if column not in EXCLUDED_NUMERIC_COLUMNS
    ]


def canonical_frame(raw: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": [f"improve:{index}" for index in range(len(raw))],
            "patient_id": "improve:" + raw["Patient"].astype(str),
            "mutant_peptide": raw["Mut_peptide"].astype(str),
            "wildtype_peptide": raw["Norm_peptide"].astype(str),
            "hla_allele": raw["HLA_allele"].astype(str),
            "label": raw["response"].map({1: "positive", 0: "negative"}),
            "official_partition": raw["Partition"].astype(str),
        }
    )


def slot_weights(frame: pd.DataFrame) -> np.ndarray:
    weights = np.ones(len(frame), dtype=np.float32)
    for patient in frame["Patient"].unique():
        patient_mask = frame["Patient"].eq(patient).to_numpy()
        negative_mask = patient_mask & frame["response"].eq(0).to_numpy()
        negative_count = int(negative_mask.sum())
        if negative_count:
            weights[negative_mask] = 20.0 / negative_count
    return weights


def paired_embedding_matrix(
    frame: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    feature_set: str,
) -> np.ndarray:
    mutant = np.stack([embeddings[str(value).upper()] for value in frame["Mut_peptide"]])
    wildtype = np.stack([embeddings[str(value).upper()] for value in frame["Norm_peptide"]])
    delta = mutant - wildtype
    absolute_delta = np.abs(delta)
    cosine = np.sum(mutant * wildtype, axis=1, keepdims=True)

    if feature_set == "none":
        return np.empty((len(frame), 0), dtype=np.float32)
    if feature_set == "delta":
        return np.concatenate([delta, absolute_delta, cosine], axis=1)
    if feature_set == "paired":
        return np.concatenate([mutant, wildtype, delta, absolute_delta, cosine], axis=1)
    raise ValueError(f"Unknown feature set: {feature_set}")


def fit_fold_features(
    train: pd.DataFrame,
    test: pd.DataFrame,
    numeric_columns: list[str],
    embeddings: dict[str, np.ndarray],
    feature_set: str,
) -> tuple[np.ndarray, np.ndarray]:
    numeric_imputer = SimpleImputer(strategy="median")
    train_numeric = numeric_imputer.fit_transform(train[numeric_columns])
    test_numeric = numeric_imputer.transform(test[numeric_columns])

    categorical_encoder = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    train_categorical = categorical_encoder.fit_transform(
        train[list(CATEGORICAL_COLUMNS)].fillna("missing").astype(str)
    )
    test_categorical = categorical_encoder.transform(
        test[list(CATEGORICAL_COLUMNS)].fillna("missing").astype(str)
    )

    train_paired = paired_embedding_matrix(train, embeddings, feature_set)
    test_paired = paired_embedding_matrix(test, embeddings, feature_set)
    return (
        np.concatenate([train_numeric, train_categorical, train_paired], axis=1),
        np.concatenate([test_numeric, test_categorical, test_paired], axis=1),
    )


def evaluate_feature_set(
    raw: pd.DataFrame,
    embeddings: dict[str, np.ndarray],
    feature_set: str,
) -> tuple[pd.DataFrame, dict[str, float]]:
    numeric_columns = safe_numeric_columns(raw)
    scored_parts: list[pd.DataFrame] = []

    for partition in sorted(raw["Partition"].astype(str).unique()):
        test_mask = raw["Partition"].astype(str).eq(partition)
        train = raw.loc[~test_mask].copy()
        test = raw.loc[test_mask].copy()

        test_peptides = set(test["Mut_peptide"].astype(str))
        train = train.loc[~train["Mut_peptide"].astype(str).isin(test_peptides)].copy()
        train_features, test_features = fit_fold_features(
            train,
            test,
            numeric_columns,
            embeddings,
            feature_set,
        )

        model = XGBClassifier(
            n_estimators=400,
            learning_rate=0.03,
            max_depth=2,
            min_child_weight=3,
            subsample=0.9,
            colsample_bytree=0.8,
            reg_lambda=3.0,
            objective="binary:logistic",
            eval_metric="logloss",
            tree_method="hist",
            random_state=17,
            n_jobs=-1,
        )
        model.fit(
            train_features,
            train["response"].astype(int),
            sample_weight=slot_weights(train),
        )
        part = canonical_frame(test)
        part["score"] = model.predict_proba(test_features)[:, 1]
        scored_parts.append(part)

    scored = pd.concat(scored_parts, ignore_index=True)
    summary = summarize_group_metrics(
        group_metrics(scored, group_col="patient_id", score_col="score", k=20)
    )
    return scored, summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--embedding-cache", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--member", default=DEFAULT_MEMBER)
    parser.add_argument(
        "--feature-set",
        action="append",
        choices=("none", "delta", "paired"),
        dest="feature_sets",
    )
    args = parser.parse_args()

    raw = read_improve(args.data, args.member)
    embeddings, model_name = load_plm_embedding_cache(args.embedding_cache)
    feature_sets = args.feature_sets or ["none", "delta", "paired"]
    args.output_dir.mkdir(parents=True, exist_ok=True)
    results: list[dict[str, object]] = []

    for feature_set in feature_sets:
        scored, summary = evaluate_feature_set(raw, embeddings, feature_set)
        scored.to_csv(args.output_dir / f"improve_{feature_set}_esm_oof.csv", index=False)
        results.append(
            {
                "feature_set": feature_set,
                "embedding_model": model_name if feature_set != "none" else None,
                "summary": summary,
            }
        )
        print(feature_set, json.dumps(summary, sort_keys=True))

    (args.output_dir / "improve_paired_esm.json").write_text(
        json.dumps(results, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
