import pandas as pd

from epicurus_neo.score_selection import apply_score_selection, select_score_columns_by_group


def test_select_score_columns_by_group_uses_validation_winner():
    validation = pd.DataFrame(
        {
            "hla_allele": ["A", "A", "A", "B", "B", "B"],
            "label": ["positive", "negative", "negative", "positive", "negative", "negative"],
            "score_a": [0.9, 0.2, 0.1, 0.1, 0.2, 0.9],
            "score_b": [0.1, 0.2, 0.9, 0.9, 0.2, 0.1],
        }
    )

    selection = select_score_columns_by_group(
        validation,
        group_col="hla_allele",
        score_columns=["score_a", "score_b"],
        k=1,
    )

    assert selection.group_score_cols["A"] == "score_a"
    assert selection.group_score_cols["B"] == "score_b"


def test_apply_score_selection_tracks_score_source():
    validation = pd.DataFrame(
        {
            "hla_allele": ["A", "A", "B", "B"],
            "label": ["positive", "negative", "positive", "negative"],
            "score_a": [0.9, 0.1, 0.1, 0.9],
            "score_b": [0.1, 0.9, 0.9, 0.1],
        }
    )
    target = pd.DataFrame(
        {
            "hla_allele": ["A", "B", "C"],
            "score_a": [0.2, 0.3, 0.4],
            "score_b": [0.8, 0.7, 0.6],
        }
    )
    selection = select_score_columns_by_group(
        validation,
        group_col="hla_allele",
        score_columns=["score_a", "score_b"],
        k=1,
    )

    out = apply_score_selection(target, selection, group_col="hla_allele")

    assert out.loc[0, "epicurus_selected_score"] == 0.2
    assert out.loc[1, "epicurus_selected_score"] == 0.7
    assert out.loc[2, "epicurus_selected_score_source"] == selection.default_score_col


def test_select_score_columns_uses_default_when_group_has_too_few_positives():
    validation = pd.DataFrame(
        {
            "hla_allele": ["A", "A", "B", "B"],
            "label": ["positive", "negative", "negative", "negative"],
            "score_a": [0.9, 0.1, 0.1, 0.9],
            "score_b": [0.1, 0.9, 0.9, 0.1],
        }
    )

    selection = select_score_columns_by_group(
        validation,
        group_col="hla_allele",
        score_columns=["score_a", "score_b"],
        k=1,
        min_positive=1,
    )

    assert selection.group_score_cols["B"] == selection.default_score_col
