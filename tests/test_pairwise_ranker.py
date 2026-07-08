import pandas as pd

from epicurus_neo.pairwise_ranker import fit_pairwise_ranker


def test_pairwise_ranker_scores_positive_like_rows_higher():
    train = pd.DataFrame(
        {
            "patient_id": ["p1"] * 6,
            "label": ["positive", "positive", "negative", "negative", "negative", "negative"],
            "label_weight": [1.0] * 6,
            "feature": [0.9, 0.8, 0.2, 0.1, 0.3, 0.4],
        }
    )
    target = pd.DataFrame(
        {
            "patient_id": ["p2", "p2"],
            "label": ["unknown", "unknown"],
            "feature": [0.95, 0.05],
        }
    )

    ranker = fit_pairwise_ranker(train, group_col="patient_id")
    scored = ranker.predict_scores(target)

    assert "epicurus_pairwise_score" in scored.columns
    assert scored.loc[0, "epicurus_pairwise_score"] > scored.loc[1, "epicurus_pairwise_score"]
