import pandas as pd

from epicurus_neo.cli import build_parser, cmd_apply_score_selector
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


def test_select_score_columns_can_prioritize_mrr():
    validation = pd.DataFrame(
        {
            "hla_allele": ["A", "A", "A", "A", "A"],
            "label": ["negative", "positive", "positive", "positive", "negative"],
            "score_hits": [0.9, 0.8, 0.7, 0.1, 0.0],
            "score_mrr": [0.8, 0.9, 0.1, 0.0, 0.7],
        }
    )

    hits_selection = select_score_columns_by_group(
        validation,
        group_col="hla_allele",
        score_columns=["score_hits", "score_mrr"],
        k=3,
        objective="hits",
    )
    mrr_selection = select_score_columns_by_group(
        validation,
        group_col="hla_allele",
        score_columns=["score_hits", "score_mrr"],
        k=3,
        objective="mrr",
    )

    assert hits_selection.group_score_cols["A"] == "score_hits"
    assert mrr_selection.group_score_cols["A"] == "score_mrr"


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


def test_apply_score_selector_cli_writes_outputs(tmp_path):
    validation = tmp_path / "validation.csv"
    target = tmp_path / "target.csv"
    output = tmp_path / "scored.csv"
    selection_output = tmp_path / "selection.json"
    pd.DataFrame(
        {
            "hla_allele": ["A", "A"],
            "label": ["positive", "negative"],
            "score_a": [0.9, 0.1],
            "score_b": [0.1, 0.9],
        }
    ).to_csv(validation, index=False)
    pd.DataFrame(
        {
            "hla_allele": ["A"],
            "score_a": [0.4],
            "score_b": [0.6],
        }
    ).to_csv(target, index=False)

    parser = build_parser()
    args = parser.parse_args(
        [
            "apply-score-selector",
            "--validation",
            str(validation),
            "--target",
            str(target),
            "--output",
            str(output),
            "--selection-output",
            str(selection_output),
            "--group-col",
            "hla_allele",
            "--score-col",
            "score_a",
            "--score-col",
            "score_b",
            "-k",
            "1",
            "--objective",
            "hits",
        ]
    )

    assert cmd_apply_score_selector(args) == 0
    assert "epicurus_selected_score" in pd.read_csv(output).columns
    assert "score_a" in selection_output.read_text()
