import pandas as pd

from epicurus_neo.cli import build_parser, cmd_apply_blend_selector
from epicurus_neo.score_blending import (
    add_rank_blend_score,
    apply_blend_selection,
    select_blends_by_group,
)


def test_add_rank_blend_score_uses_group_local_percentile_ranks():
    frame = pd.DataFrame(
        {
            "hla_allele": ["A", "A", "B", "B"],
            "score_a": [10.0, 1.0, 2.0, 1.0],
            "score_b": [1.0, 10.0, 1.0, 2.0],
        }
    )

    out = add_rank_blend_score(
        frame,
        group_col="hla_allele",
        weights={"score_a": 0.5, "score_b": 0.5},
        output_col="blend",
    )

    assert out.loc[0, "blend"] == 0.75
    assert out.loc[1, "blend"] == 0.75
    assert out.loc[2, "blend"] == 0.75
    assert out.loc[3, "blend"] == 0.75


def test_select_blends_by_group_can_choose_pairwise_blend():
    validation = pd.DataFrame(
        {
            "hla_allele": ["A"] * 5,
            "label": ["positive", "negative", "positive", "negative", "negative"],
            "score_a": [0.9, 0.8, 0.1, 0.2, 0.0],
            "score_b": [0.2, 0.0, 0.9, 0.8, 0.1],
        }
    )

    selection = select_blends_by_group(
        validation,
        group_col="hla_allele",
        score_columns=["score_a", "score_b"],
        k=2,
        objective="hits",
    )

    assert selection.group_weights["A"] == {"score_a": 0.5, "score_b": 0.5}


def test_apply_blend_selection_tracks_blend_source():
    validation = pd.DataFrame(
        {
            "hla_allele": ["A"] * 5,
            "label": ["positive", "negative", "positive", "negative", "negative"],
            "score_a": [0.9, 0.8, 0.1, 0.2, 0.0],
            "score_b": [0.2, 0.0, 0.9, 0.8, 0.1],
        }
    )
    target = pd.DataFrame(
        {
            "hla_allele": ["A", "A", "B"],
            "score_a": [0.9, 0.1, 0.4],
            "score_b": [0.2, 0.8, 0.6],
        }
    )
    selection = select_blends_by_group(
        validation,
        group_col="hla_allele",
        score_columns=["score_a", "score_b"],
        k=2,
        objective="hits",
    )

    out = apply_blend_selection(target, selection, group_col="hla_allele")

    assert "epicurus_blend_score" in out.columns
    assert out.loc[0, "epicurus_blend_score_source"] == "0.50*score_a+0.50*score_b"
    assert out.loc[2, "epicurus_blend_score_source"] == "0.50*score_a+0.50*score_b"


def test_apply_blend_selector_cli_writes_outputs(tmp_path):
    validation = tmp_path / "validation.csv"
    target = tmp_path / "target.csv"
    output = tmp_path / "scored.csv"
    selection_output = tmp_path / "selection.json"
    pd.DataFrame(
        {
            "hla_allele": ["A"] * 5,
            "label": ["positive", "negative", "positive", "negative", "negative"],
            "score_a": [0.9, 0.8, 0.1, 0.2, 0.0],
            "score_b": [0.2, 0.0, 0.9, 0.8, 0.1],
        }
    ).to_csv(validation, index=False)
    pd.DataFrame(
        {
            "hla_allele": ["A", "A"],
            "score_a": [0.4, 0.6],
            "score_b": [0.6, 0.4],
        }
    ).to_csv(target, index=False)

    parser = build_parser()
    args = parser.parse_args(
        [
            "apply-blend-selector",
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
            "2",
            "--objective",
            "hits",
        ]
    )

    assert cmd_apply_blend_selector(args) == 0
    assert "epicurus_blend_score" in pd.read_csv(output).columns
    assert "default_blend_name" in selection_output.read_text()
