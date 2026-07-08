import pandas as pd
import pytest

from epicurus_neo.benchmark import (
    add_groupwise_ensemble_scores,
    add_retrieval_score,
    add_transferable_presentation_score,
    add_weighted_groupwise_score,
    train_and_evaluate,
)
from epicurus_neo.cli import build_parser, cmd_score_report, cmd_train_eval
from epicurus_neo.features import add_baseline_scores, infer_numeric_feature_columns
from epicurus_neo.model import fit_ranker


def _toy_rows(patient: str, study: str, offset: int) -> list[dict]:
    rows = []
    for idx in range(12):
        positive = idx in {2, 7}
        rows.append(
            {
                "candidate_id": f"{patient}-{idx}",
                "source_dataset": "toy",
                "study_id": study,
                "patient_id": patient,
                "hla_allele": "HLA-A*02:01",
                "mutant_peptide": f"PEPTIDE{offset}{idx}",
                "wildtype_peptide": f"SELTIDE{offset}{idx}",
                "label": "positive" if positive else "negative",
                "label_weight": 1.0,
                "assay_type": "synthetic",
                "binding_affinity_nm": 20.0 + idx,
                "presentation_score": 0.95 if positive else 0.2 + idx * 0.01,
                "expression_tpm": 50.0 if positive else 5.0,
                "foreignness_score": 0.8 if positive else 0.1,
                "mutation_tcr_face": 1 if positive else 0,
                "Screening Status": 1 if positive else 0,
                "Nmer score": 0.9 if positive else 0.1,
                "Top netMHCpan4.0 EL ranked minimal": 0.1 if positive else 2.0,
            }
        )
    return rows


def test_add_baseline_scores_and_feature_inference():
    frame = pd.DataFrame(_toy_rows("p1", "s1", 1))
    frame["all_null_numeric"] = pd.NA
    scored = add_baseline_scores(frame)
    assert "baseline_binding_score" in scored.columns
    assert "baseline_pvac_style_score" in scored.columns
    assert "baseline_gartner_nmer_score" in scored.columns
    assert "baseline_netmhcpan_el_score" in scored.columns
    features = infer_numeric_feature_columns(scored)
    assert "presentation_score" in features
    assert "all_null_numeric" not in features
    assert "label_weight" not in features
    assert "Screening Status" not in features
    assert "rand" not in features
    assert "retrieval_fold" not in features


def test_fit_ranker_scores_candidates():
    train = pd.DataFrame(_toy_rows("p1", "s1", 1) + _toy_rows("p2", "s2", 2))
    test = pd.DataFrame(_toy_rows("p3", "s3", 3))
    model = fit_ranker(add_baseline_scores(train))
    scored = model.predict_scores(add_baseline_scores(test))
    assert "epicurus_score" in scored.columns
    assert scored["epicurus_score"].between(0, 1).all()


def test_fit_ranker_can_emit_uncertainty_scores():
    train = pd.DataFrame(_toy_rows("p1", "s1", 1) + _toy_rows("p2", "s2", 2))
    test = pd.DataFrame(_toy_rows("p3", "s3", 3))
    model = fit_ranker(
        add_baseline_scores(train),
        uncertainty_ensemble_size=3,
        uncertainty_penalty=1.0,
    )
    scored = model.predict_scores(add_baseline_scores(test))

    assert "epicurus_score_std" in scored.columns
    assert "epicurus_lower_confidence_score" in scored.columns
    assert scored["epicurus_lower_confidence_score"].between(0, 1).all()
    assert (scored["epicurus_lower_confidence_score"] <= scored["epicurus_score"]).all()


def test_ranker_predict_scores_tolerates_missing_test_features():
    train = pd.DataFrame(_toy_rows("p1", "s1", 1) + _toy_rows("p2", "s2", 2))
    test = pd.DataFrame(_toy_rows("p3", "s3", 3)).drop(columns=["Nmer score"])
    model = fit_ranker(add_baseline_scores(train))
    scored = model.predict_scores(add_baseline_scores(test))
    assert scored["epicurus_score"].between(0, 1).all()


def test_train_and_evaluate_blocks_exact_leakage():
    train = pd.DataFrame(_toy_rows("p1", "s1", 1))
    test = pd.DataFrame(_toy_rows("p1", "s1", 1))
    with pytest.raises(ValueError, match="leakage"):
        train_and_evaluate(train, test)


def test_train_and_evaluate_reports_ranker_and_baselines():
    train = pd.DataFrame(_toy_rows("p1", "s1", 1) + _toy_rows("p2", "s2", 2))
    test = pd.DataFrame(_toy_rows("p3", "s3", 3) + _toy_rows("p4", "s4", 4))
    result = train_and_evaluate(train, test, k=5, uncertainty_ensemble_size=3)
    score_cols = {item.score_col for item in result.benchmark_results}
    assert "epicurus_hits20_score" in score_cols
    assert "epicurus_blend_score" in score_cols
    assert "epicurus_pairwise_score" in score_cols
    assert "epicurus_lower_confidence_score" in score_cols
    assert "epicurus_score" in score_cols
    assert "baseline_gartner_nmer_score" in score_cols
    assert "baseline_pvac_style_score" in score_cols
    assert result.feature_columns


def test_add_groupwise_ensemble_scores_rank_normalizes_components():
    frame = pd.DataFrame(
        {
            "patient_id": ["p1", "p1", "p1"],
            "a": [0.1, 0.2, 0.3],
            "b": [3.0, 2.0, 1.0],
        }
    )
    out = add_groupwise_ensemble_scores(
        frame,
        group_col="patient_id",
        component_cols=["a", "b"],
    )
    assert "epicurus_blend_score" in out.columns
    assert out["epicurus_blend_score"].between(0, 1).all()


def test_add_weighted_groupwise_score_uses_available_components():
    frame = pd.DataFrame(
        {
            "patient_id": ["p1", "p1", "p1"],
            "baseline_gartner_nmer_score": [0.1, 0.2, 0.3],
            "baseline_netmhcpan_el_score": [3.0, 2.0, 1.0],
        }
    )
    out = add_weighted_groupwise_score(
        frame,
        group_col="patient_id",
        weights={"baseline_gartner_nmer_score": 0.9, "baseline_netmhcpan_el_score": 0.1},
        output_col="epicurus_hits20_score",
    )
    assert "epicurus_hits20_score" in out.columns
    assert out["epicurus_hits20_score"].between(0, 1).all()


def test_add_transferable_presentation_score():
    frame = pd.DataFrame(
        {
            "patient_id": ["p1", "p1", "p1"],
            "mhcflurry_presentation_score": [0.9, 0.8, 0.1],
            "seq_hydrophobicity_mean": [-1.0, 2.0, 0.0],
            "seq_cysteine_fraction": [0.0, 0.1, 0.0],
            "seq_aromatic_fraction": [0.1, 0.2, 0.0],
        }
    )
    out = add_transferable_presentation_score(frame, group_col="patient_id")
    assert "epicurus_transfer_score" in out.columns
    assert out["epicurus_transfer_score"].between(0, 1).all()
    assert out.loc[0, "epicurus_transfer_score"] > out.loc[2, "epicurus_transfer_score"]


def test_add_retrieval_score_uses_positive_similarity():
    frame = pd.DataFrame({"retrieval_max_positive_similarity": [0.2, 0.9]})
    out = add_retrieval_score(frame)
    assert out["epicurus_retrieval_score"].tolist() == [0.2, 0.9]


def test_train_eval_cli_can_ignore_shared_study_and_purge(tmp_path):
    train = pd.DataFrame(_toy_rows("p1", "s1", 1) + _toy_rows("p2", "s1", 2))
    test = pd.DataFrame(_toy_rows("p3", "s1", 2))
    train_path = tmp_path / "train.csv"
    test_path = tmp_path / "test.csv"
    scored_path = tmp_path / "scored.csv"
    train.to_csv(train_path, index=False)
    test.to_csv(test_path, index=False)

    parser = build_parser()
    args = parser.parse_args(
        [
            "train-eval",
            "--train",
            str(train_path),
            "--test",
            str(test_path),
            "--ignore-shared-study",
            "--purge-exact-overlaps",
            "--uncertainty-ensemble-size",
            "2",
            "--write-scored",
            str(scored_path),
        ]
    )
    assert cmd_train_eval(args) == 0
    assert scored_path.exists()


def test_score_report_cli_writes_multiple_scores(tmp_path):
    table = pd.DataFrame(
        {
            "patient_id": ["p1", "p1", "p1"],
            "label": ["positive", "negative", "positive"],
            "score_a": [0.9, 0.1, 0.8],
            "score_b": [0.1, 0.9, 0.8],
        }
    )
    table_path = tmp_path / "scores.csv"
    output_path = tmp_path / "report.json"
    table.to_csv(table_path, index=False)

    parser = build_parser()
    args = parser.parse_args(
        [
            "score-report",
            str(table_path),
            "--score-col",
            "score_a",
            "--score-col",
            "score_b",
            "--output",
            str(output_path),
        ]
    )
    assert cmd_score_report(args) == 0
    assert output_path.exists()
    assert "score_a" in output_path.read_text()
