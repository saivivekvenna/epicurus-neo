import pandas as pd

from epicurus_neo.cli import build_parser, cmd_apply_blend_selector, cmd_apply_guarded_blend_selector
from epicurus_neo.score_blending import (
    add_rank_blend_score,
    apply_blend_selection,
    guarded_best_blend,
    select_guarded_blends_by_group,
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


def test_select_blends_by_group_accepts_custom_pair_weights():
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
        pair_weights=(0.4,),
    )

    assert selection.group_weights["A"] == {"score_a": 0.4, "score_b": 0.6}


def test_guarded_best_blend_blocks_candidates_that_lose_hits():
    validation = pd.DataFrame(
        {
            "hla_allele": ["A"] * 4,
            "label": ["positive", "positive", "negative", "negative"],
            "score_hits": [0.8, 0.7, 0.9, 0.1],
            "score_mrr": [0.95, 0.2, 0.9, 0.8],
        }
    )

    weights, summary = guarded_best_blend(
        validation,
        group_col="hla_allele",
        score_columns=["score_hits", "score_mrr"],
        k=3,
        baseline_objective="hits",
        objective="mrr",
        guard_metric="mean_hits_at_k",
        pair_weights=(0.9,),
    )

    assert weights == {"score_hits": 1.0}
    assert summary["baseline"]["mean_hits_at_k"] == 2.0
    assert summary["selected"]["mean_hits_at_k"] == 2.0


def test_select_guarded_blends_by_group_tracks_baseline_and_selected_summaries():
    validation = pd.DataFrame(
        {
            "hla_allele": ["A"] * 4,
            "label": ["positive", "positive", "negative", "negative"],
            "score_hits": [0.8, 0.7, 0.9, 0.1],
            "score_mrr": [0.95, 0.2, 0.9, 0.8],
        }
    )

    selection = select_guarded_blends_by_group(
        validation,
        group_col="hla_allele",
        score_columns=["score_hits", "score_mrr"],
        k=3,
        baseline_objective="hits",
        objective="mrr",
        pair_weights=(0.9,),
    )

    assert selection.group_weights["A"] == {"score_hits": 1.0}
    assert "baseline" in selection.validation_summary["A"]
    assert "selected" in selection.validation_summary["A"]


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
            "--pair-weight",
            "0.4",
            "-k",
            "2",
            "--objective",
            "hits",
        ]
    )

    assert cmd_apply_blend_selector(args) == 0
    assert "epicurus_blend_score" in pd.read_csv(output).columns
    selection_text = selection_output.read_text()
    assert "default_blend_name" in selection_text
    assert "0.4" in selection_text


def test_apply_guarded_blend_selector_cli_writes_outputs(tmp_path):
    validation = tmp_path / "validation.csv"
    target = tmp_path / "target.csv"
    output = tmp_path / "scored.csv"
    selection_output = tmp_path / "selection.json"
    pd.DataFrame(
        {
            "hla_allele": ["A"] * 4,
            "label": ["positive", "positive", "negative", "negative"],
            "score_hits": [0.8, 0.7, 0.9, 0.1],
            "score_mrr": [0.95, 0.2, 0.9, 0.8],
        }
    ).to_csv(validation, index=False)
    pd.DataFrame(
        {
            "hla_allele": ["A", "A"],
            "score_hits": [0.4, 0.6],
            "score_mrr": [0.6, 0.4],
        }
    ).to_csv(target, index=False)

    parser = build_parser()
    args = parser.parse_args(
        [
            "apply-guarded-blend-selector",
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
            "score_hits",
            "--score-col",
            "score_mrr",
            "-k",
            "3",
            "--baseline-objective",
            "hits",
            "--objective",
            "mrr",
            "--guard-metric",
            "mean_hits_at_k",
        ]
    )

    assert cmd_apply_guarded_blend_selector(args) == 0
    assert "epicurus_blend_score" in pd.read_csv(output).columns
    selection_text = selection_output.read_text()
    assert "baseline_objective" in selection_text
    assert "guard_metric" in selection_text
