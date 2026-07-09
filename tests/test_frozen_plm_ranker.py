import numpy as np
import pandas as pd

from epicurus_neo.cli import build_parser
from epicurus_neo.frozen_plm_ranker import (
    XGBRankerConfig,
    build_frozen_plm_features,
    refit_and_score_frozen_plm_ranker,
    select_frozen_plm_ranker,
)


def _frame(prefix: str) -> pd.DataFrame:
    rows = []
    for allele in ("HLA-A*01:01", "HLA-A*02:01"):
        for idx in range(6):
            positive = idx < 2
            rows.append(
                {
                    "candidate_id": f"{prefix}-{allele}-{idx}",
                    "hla_allele": allele,
                    "mutant_peptide": f"{prefix}{idx}PEPTIDE",
                    "label": "positive" if positive else "negative",
                    "mhcflurry_presentation_score": 0.9 if positive else 0.1,
                    "mutation_tcr_face_count": 3 if positive else 0,
                }
            )
    return pd.DataFrame(rows)


def _embeddings(*frames: pd.DataFrame) -> dict[str, np.ndarray]:
    result = {}
    for frame in frames:
        for row in frame.itertuples():
            positive = row.label == "positive"
            result[row.mutant_peptide] = np.array(
                [1.0, 0.0, 0.5] if positive else [0.0, 1.0, -0.5],
                dtype=np.float32,
            )
    return result


def _tiny_config() -> tuple[XGBRankerConfig, ...]:
    return (
        XGBRankerConfig(
            n_estimators=5,
            learning_rate=0.2,
            max_depth=2,
            subsample=1.0,
            colsample_bytree=1.0,
        ),
    )


def test_build_frozen_plm_features_combines_embeddings_and_alleles():
    frame = _frame("train")
    embeddings = _embeddings(frame)

    features = build_frozen_plm_features(
        frame,
        embeddings,
        allele_vocabulary=("HLA-A*01:01", "HLA-A*02:01"),
    )

    assert features.shape[0] == len(frame)
    assert "plm_embedding_002" in features
    assert "mhcflurry_presentation_score" in features
    assert features.filter(like="hla_").shape[1] == 2
    assert (features.filter(like="hla_").sum(axis=1) == 1).all()


def test_build_frozen_plm_features_includes_screened_recognition_signal():
    frame = _frame("train")
    frame["screened_recognition_hla_centered_top10_response_rate"] = [
        0.8,
        0.1,
        0.6,
        0.2,
        0.7,
        0.3,
    ] * 2
    embeddings = _embeddings(frame)

    features = build_frozen_plm_features(
        frame,
        embeddings,
        allele_vocabulary=("HLA-A*01:01", "HLA-A*02:01"),
    )

    assert (
        features["screened_recognition_hla_centered_top10_response_rate"].tolist()
        == frame["screened_recognition_hla_centered_top10_response_rate"].tolist()
    )


def test_frozen_plm_ranker_selects_on_validation_and_refits():
    train = _frame("train")
    validation = _frame("validation")
    target = _frame("target")
    embeddings = _embeddings(train, validation, target)

    selection = select_frozen_plm_ranker(
        train,
        validation,
        embeddings,
        model_name="test/model",
        configs=_tiny_config(),
        k=3,
    )
    scored = refit_and_score_frozen_plm_ranker(
        pd.concat([train, validation], ignore_index=True),
        target,
        embeddings,
        selection,
    )

    assert selection.model_name == "test/model"
    assert selection.validation_summary["mean_hits_at_k"] == 2.0
    assert len(selection.candidate_summaries) == 1
    assert scored["epicurus_frozen_plm_score"].notna().all()


def test_frozen_plm_rank_cli_parses_paths():
    parser = build_parser()
    args = parser.parse_args(
        [
            "frozen-plm-rank",
            "--train",
            "train.csv",
            "--validation",
            "validation.csv",
            "--train-validation",
            "train-validation.csv",
            "--target",
            "target.csv",
            "--embedding-cache",
            "embeddings.npz",
            "--output",
            "scored.csv",
            "--selection-output",
            "selection.json",
        ]
    )

    assert args.command == "frozen-plm-rank"
    assert args.group_col == "hla_allele"
