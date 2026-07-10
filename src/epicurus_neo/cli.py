from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from benchmark.scorecard import scorecard as build_scorecard

from epicurus_neo.benchmark import train_and_evaluate
from epicurus_neo.auto_research import build_failure_report, write_research_artifacts
from epicurus_neo.data_manifest import load_dataset_manifest
from epicurus_neo.download import dataset_file_plans, download_file
from epicurus_neo.experiment import grouped_cross_validate, summarize_cross_validation
from epicurus_neo.leakage import detect_exact_leakage, purge_train_overlaps
from epicurus_neo.normalize import (
    normalize_bigmhc_table,
    normalize_candidate_table,
    normalize_cd8_multimer_2025,
    normalize_gartner_table,
    normalize_improve_cv,
    normalize_neoranking_neopep,
    write_normalized,
)
from epicurus_neo.portfolio import PortfolioConstraints, select_portfolio
from epicurus_neo.schema import validate_schema
from epicurus_neo.splits import assign_holdout_split


def _load_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    if path.suffix.lower() in {".tsv", ".txt"}:
        return pd.read_csv(path, sep="\t")
    raise ValueError(f"Unsupported table format: {path}")


def cmd_validate(args: argparse.Namespace) -> int:
    frame = _load_table(Path(args.table))
    report = validate_schema(frame)
    print(json.dumps(report.__dict__, indent=2))
    return 0 if report.ok else 1


def cmd_metrics(args: argparse.Namespace) -> int:
    frame = _load_table(Path(args.table))
    baseline_col = args.baseline_col or args.score_col
    report = build_scorecard(
        frame,
        args.score_col,
        baseline_col,
        group_col=args.group_col,
        k=args.k,
    )
    print(json.dumps(report, indent=2))
    return 0


def cmd_score_report(args: argparse.Namespace) -> int:
    frame = _load_table(Path(args.table))
    baseline_col = args.baseline_col or args.score_col[0]
    results = [
        {
            "score_col": score_col,
            "scorecard": build_scorecard(
                frame,
                score_col,
                baseline_col,
                group_col=args.group_col,
                k=args.k,
            ),
        }
        for score_col in args.score_col
    ]
    payload = {
        "table": args.table,
        "group_col": args.group_col,
        "baseline_col": baseline_col,
        "k": args.k,
        "benchmarks": results,
    }
    text = json.dumps(payload, indent=2)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(text + "\n")
    print(text)
    return 0


def cmd_compare_metrics(args: argparse.Namespace) -> int:
    from epicurus_neo.metric_compare import compare_metric_reports, compare_metric_reports_file

    if args.output:
        output = compare_metric_reports_file(args.report, args.output, sort_by=args.sort_by)
        print(output)
    else:
        frame = compare_metric_reports(args.report, sort_by=args.sort_by)
        print(frame.to_csv(index=False))
    return 0


def cmd_precision_filter(args: argparse.Namespace) -> int:
    from epicurus_neo.precision_filter import (
        GroupedPrecisionThreshold,
        apply_precision_threshold_files,
    )

    output, threshold, target_summary = apply_precision_threshold_files(
        args.validation,
        args.target,
        args.output,
        score_col=args.score_col,
        target_precision=args.target_precision,
        min_selected=args.min_selected,
        group_col=args.group_col,
        min_group_positives=args.min_group_positives,
        min_group_selected=args.min_group_selected,
    )
    if isinstance(threshold, GroupedPrecisionThreshold):
        threshold_payload = {
            "score_col": threshold.score_col,
            "group_col": threshold.group_col,
            "default_threshold": threshold.default_threshold.__dict__,
            "group_thresholds": {
                group: group_threshold.__dict__
                for group, group_threshold in threshold.group_thresholds.items()
            },
            "target_precision": threshold.target_precision,
            "min_group_positives": threshold.min_group_positives,
            "min_group_selected": threshold.min_group_selected,
        }
    else:
        threshold_payload = threshold.__dict__
    payload = {
        "output": str(output),
        "threshold": threshold_payload,
        "target_summary": target_summary,
    }
    if args.report_output:
        Path(args.report_output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.report_output).write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


def cmd_leakage(args: argparse.Namespace) -> int:
    train = _load_table(Path(args.train))
    test = _load_table(Path(args.test))
    report = detect_exact_leakage(train, test)
    print(json.dumps(report.__dict__, indent=2))
    return 1 if report.has_leakage else 0


def cmd_train_eval(args: argparse.Namespace) -> int:
    train = _load_table(Path(args.train))
    test = _load_table(Path(args.test))
    if args.purge_exact_overlaps:
        train = purge_train_overlaps(train, test)
    result = train_and_evaluate(
        train,
        test,
        group_col=args.group_col,
        k=args.k,
        allow_exact_leakage=args.allow_exact_leakage,
        include_shared_studies_as_leakage=not args.ignore_shared_study,
        uncertainty_ensemble_size=args.uncertainty_ensemble_size,
        uncertainty_penalty=args.uncertainty_penalty,
    )
    payload = {
        "feature_columns": result.feature_columns,
        "leakage": result.leakage.__dict__,
        "benchmarks": [
            {"score_col": item.score_col, "summary": item.summary}
            for item in result.benchmark_results
        ],
    }
    print(json.dumps(payload, indent=2))
    if args.write_scored:
        result.scored_test.to_csv(args.write_scored, index=False)
    return 0


def cmd_list_datasets(args: argparse.Namespace) -> int:
    sources = load_dataset_manifest(args.manifest)
    payload = [source.__dict__ for source in sources]
    print(json.dumps(payload, indent=2))
    return 0


def cmd_download_plan(args: argparse.Namespace) -> int:
    plans = dataset_file_plans(
        args.manifest,
        output_dir=args.output_dir,
        dataset_key=args.dataset,
    )
    print(
        json.dumps(
            [plan.__dict__ | {"output_path": str(plan.output_path)} for plan in plans], indent=2
        )
    )
    return 0


def cmd_download_file(args: argparse.Namespace) -> int:
    path = download_file(args.url, args.output, overwrite=args.overwrite)
    print(path)
    return 0


def cmd_normalize(args: argparse.Namespace) -> int:
    if args.kind == "neoranking-neopep":
        normalized = normalize_neoranking_neopep(args.input)
    elif args.kind == "gartner":
        normalized = normalize_gartner_table(args.input)
    elif args.kind == "tesla":
        from epicurus_neo.normalize import normalize_tesla_table

        normalized = normalize_tesla_table(args.input)
    elif args.kind == "bigmhc":
        normalized = normalize_bigmhc_table(args.input, zip_member=args.zip_member)
    elif args.kind == "cd8-multimer-2025":
        normalized = normalize_cd8_multimer_2025(args.input)
    elif args.kind == "improve-cv":
        normalized = normalize_improve_cv(args.input, zip_member=args.zip_member)
    else:
        normalized = normalize_candidate_table(
            _load_table(Path(args.input)),
            source_dataset=args.source_dataset,
            study_default=args.study_default,
        )
    write_normalized(normalized, args.output)
    print(json.dumps({"rows": len(normalized), "output": args.output}, indent=2))
    return 0


def cmd_make_holdout_split(args: argparse.Namespace) -> int:
    frame = _load_table(Path(args.table))
    assigned = assign_holdout_split(
        frame,
        group_col=args.group_col,
        holdout_values=args.holdout,
    )
    assigned.to_csv(args.output, index=False)
    return 0


def cmd_select_portfolio(args: argparse.Namespace) -> int:
    frame = _load_table(Path(args.table))
    selected = select_portfolio(
        frame,
        score_col=args.score_col,
        constraints=PortfolioConstraints(
            k=args.k,
            max_per_hla=args.max_per_hla,
            max_per_gene=args.max_per_gene,
            min_score=args.min_score,
        ),
    )
    selected.to_csv(args.output, index=False)
    return 0


def cmd_group_cv(args: argparse.Namespace) -> int:
    frame = _load_table(Path(args.table))
    folds = grouped_cross_validate(
        frame,
        group_col=args.group_col,
        metric_group_col=args.metric_group_col,
        k=args.k,
        max_splits=args.max_splits,
        purge_exact_overlaps=not args.no_purge_exact_overlaps,
    )
    payload = {
        "summary": summarize_cross_validation(folds),
        "folds": [
            {
                "name": fold.name,
                "status": fold.status,
                "test_groups": fold.test_groups,
                "feature_columns": fold.feature_columns,
                "reason": fold.reason,
                "benchmarks": [
                    {"score_col": item.score_col, "summary": item.summary}
                    for item in fold.benchmark_results
                ],
            }
            for fold in folds
        ],
    }
    print(json.dumps(payload, indent=2))
    return 0 if all(fold.status != "leakage_blocked" for fold in folds) else 1


def cmd_research_report(args: argparse.Namespace) -> int:
    scored = _load_table(Path(args.scored))
    report = build_failure_report(
        scored,
        group_col=args.group_col,
        score_col=args.score_col,
        k=args.k,
        max_examples=args.max_examples,
    )
    report_path, prompt_path = write_research_artifacts(report, output_dir=args.output_dir)
    print(json.dumps({"report": str(report_path), "prompt": str(prompt_path)}, indent=2))
    return 0


def cmd_mhcflurry_features(args: argparse.Namespace) -> int:
    from epicurus_neo.mhcflurry_features import add_mhcflurry_predictions_file

    output = add_mhcflurry_predictions_file(args.input, args.output)
    print(output)
    return 0


def cmd_retrieval_features(args: argparse.Namespace) -> int:
    from epicurus_neo.retrieval_features import add_retrieval_features_file

    output = add_retrieval_features_file(
        args.input,
        args.reference,
        args.output,
        top_k=args.top_k,
    )
    print(output)
    return 0


def cmd_multik_retrieval_features(args: argparse.Namespace) -> int:
    from epicurus_neo.retrieval_features import add_multik_retrieval_features_file

    output = add_multik_retrieval_features_file(
        args.input,
        args.reference,
        args.output,
        top_ks=tuple(args.top_k or [1, 3, 5, 10, 20]),
    )
    print(output)
    return 0


def cmd_crossfit_retrieval_features(args: argparse.Namespace) -> int:
    from epicurus_neo.retrieval_features import add_crossfit_retrieval_features_file

    output = add_crossfit_retrieval_features_file(
        args.input,
        args.output,
        top_k=args.top_k,
        n_folds=args.n_folds,
        fold_col=args.fold_col,
    )
    print(output)
    return 0


def cmd_plm_retrieval_features(args: argparse.Namespace) -> int:
    from epicurus_neo.plm_retrieval import add_plm_retrieval_features_file

    output = add_plm_retrieval_features_file(
        args.input,
        args.reference,
        args.output,
        model_name=args.model_name,
        batch_size=args.batch_size,
        device=args.device,
        top_k=args.top_k,
    )
    print(output)
    return 0


def cmd_build_plm_embedding_cache(args: argparse.Namespace) -> int:
    from epicurus_neo.plm_retrieval import build_plm_embedding_cache

    output = build_plm_embedding_cache(
        args.input,
        args.output,
        model_name=args.model_name,
        batch_size=args.batch_size,
        device=args.device,
    )
    print(output)
    return 0


def cmd_frozen_plm_rank(args: argparse.Namespace) -> int:
    from epicurus_neo.frozen_plm_ranker import run_frozen_plm_ranker_files

    scored_output, selection_output, selection = run_frozen_plm_ranker_files(
        args.train,
        args.validation,
        args.train_validation,
        args.target,
        args.embedding_cache,
        args.output,
        args.selection_output,
        group_col=args.group_col,
        k=args.k,
    )
    payload = {
        "output": str(scored_output),
        "selection_output": str(selection_output),
        "model_name": selection.model_name,
        "config": selection.config.__dict__,
        "validation_summary": selection.validation_summary,
    }
    print(json.dumps(payload, indent=2))
    return 0


def cmd_normalize_vdjdb(args: argparse.Namespace) -> int:
    from epicurus_neo.external_recognition import normalize_vdjdb_file

    output = normalize_vdjdb_file(args.input, args.output, min_score=args.min_score)
    print(output)
    return 0


def cmd_external_recognition_features(args: argparse.Namespace) -> int:
    from epicurus_neo.external_recognition import add_external_recognition_features_file

    output = add_external_recognition_features_file(
        args.input,
        args.reference,
        args.embedding_cache,
        args.output,
        top_k=args.top_k,
    )
    print(output)
    return 0


def cmd_screened_recognition_features(args: argparse.Namespace) -> int:
    from epicurus_neo.screened_recognition import add_screened_recognition_features_file

    output = add_screened_recognition_features_file(
        args.input,
        args.reference,
        args.embedding_cache,
        args.output,
        top_ks=tuple(args.top_k or [1, 3, 5, 10, 20]),
    )
    print(output)
    return 0


def cmd_transfer_rank(args: argparse.Namespace) -> int:
    from epicurus_neo.transfer_ranker import run_transfer_ranker_files

    validation_output, selection_output, selection = run_transfer_ranker_files(
        args.external,
        args.train,
        args.validation,
        args.embedding_cache,
        args.output,
        args.selection_output,
        external_group_col=args.external_group_col,
        target_group_col=args.group_col,
        k=args.k,
    )
    print(
        json.dumps(
            {
                "output": str(validation_output),
                "selection_output": str(selection_output),
                "config": selection.config.__dict__,
                "validation_summary": selection.validation_summary,
            },
            indent=2,
        )
    )
    return 0


def cmd_finetune_plm_rank(args: argparse.Namespace) -> int:
    from epicurus_neo.plm_finetune import run_finetuned_plm_ranker_files

    validation_output, selection_output, selection = run_finetuned_plm_ranker_files(
        args.external,
        args.train,
        args.validation,
        args.output,
        args.selection_output,
        model_name=args.model_name,
        external_group_col=args.external_group_col,
        target_group_col=args.group_col,
        k=args.k,
        device=args.device,
    )
    print(
        json.dumps(
            {
                "output": str(validation_output),
                "selection_output": str(selection_output),
                "config": selection.config.__dict__,
                "selected_target_epoch": selection.selected_target_epoch,
                "validation_summary": selection.validation_summary,
            },
            indent=2,
        )
    )
    return 0


def cmd_apply_score_selector(args: argparse.Namespace) -> int:
    from epicurus_neo.score_selection import apply_score_selection_files

    output, selection = apply_score_selection_files(
        args.validation,
        args.target,
        args.output,
        group_col=args.group_col,
        score_columns=args.score_col,
        k=args.k,
        min_positive=args.min_positive,
        objective=args.objective,
    )
    payload = {
        "output": str(output),
        "default_score_col": selection.default_score_col,
        "group_score_cols": selection.group_score_cols,
        "validation_summary": selection.validation_summary,
    }
    if args.selection_output:
        Path(args.selection_output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.selection_output).write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


def cmd_apply_blend_selector(args: argparse.Namespace) -> int:
    from epicurus_neo.score_blending import apply_blend_selection_files, blend_name

    pair_weights = _parse_pair_weights(args.pair_weight)
    output, selection = apply_blend_selection_files(
        args.validation,
        args.target,
        args.output,
        group_col=args.group_col,
        score_columns=args.score_col,
        k=args.k,
        min_positive=args.min_positive,
        objective=args.objective,
        pair_weights=pair_weights,
    )
    payload = {
        "output": str(output),
        "pair_weights": list(pair_weights),
        "default_blend": selection.default_weights,
        "default_blend_name": blend_name(selection.default_weights),
        "group_blends": selection.group_weights,
        "validation_summary": selection.validation_summary,
    }
    if args.selection_output:
        Path(args.selection_output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.selection_output).write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


def _parse_pair_weights(raw_weights: list[float] | None) -> tuple[float, ...]:
    pair_weights = tuple(raw_weights) if raw_weights else (0.25, 0.5, 0.75)
    for weight in pair_weights:
        if weight <= 0.0 or weight >= 1.0:
            raise ValueError("--pair-weight values must be between 0 and 1")
    return pair_weights


def cmd_apply_guarded_blend_selector(args: argparse.Namespace) -> int:
    from epicurus_neo.score_blending import apply_guarded_blend_selection_files, blend_name

    pair_weights = _parse_pair_weights(args.pair_weight)

    output, selection = apply_guarded_blend_selection_files(
        args.validation,
        args.target,
        args.output,
        group_col=args.group_col,
        score_columns=args.score_col,
        k=args.k,
        min_positive=args.min_positive,
        baseline_objective=args.baseline_objective,
        objective=args.objective,
        guard_metric=args.guard_metric,
        min_guard_delta=args.min_guard_delta,
        pair_weights=pair_weights,
    )
    payload = {
        "output": str(output),
        "pair_weights": list(pair_weights),
        "baseline_objective": args.baseline_objective,
        "objective": args.objective,
        "guard_metric": args.guard_metric,
        "min_guard_delta": args.min_guard_delta,
        "default_blend": selection.default_weights,
        "default_blend_name": blend_name(selection.default_weights),
        "group_blends": selection.group_weights,
        "validation_summary": selection.validation_summary,
    }
    if args.selection_output:
        Path(args.selection_output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.selection_output).write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="epicurus")
    sub = parser.add_subparsers(dest="command", required=True)

    validate = sub.add_parser("validate-schema")
    validate.add_argument("table")
    validate.set_defaults(func=cmd_validate)

    metrics = sub.add_parser("metrics")
    metrics.add_argument("table")
    metrics.add_argument("--group-col", default="patient_id")
    metrics.add_argument("--score-col", required=True)
    metrics.add_argument("--baseline-col")
    metrics.add_argument("-k", type=int, default=20)
    metrics.set_defaults(func=cmd_metrics)

    score_report = sub.add_parser("score-report")
    score_report.add_argument("table")
    score_report.add_argument("--group-col", default="patient_id")
    score_report.add_argument("--score-col", action="append", required=True)
    score_report.add_argument("--baseline-col")
    score_report.add_argument("-k", type=int, default=20)
    score_report.add_argument("--output")
    score_report.set_defaults(func=cmd_score_report)

    compare_metrics = sub.add_parser("compare-metrics")
    compare_metrics.add_argument("report", nargs="+")
    compare_metrics.add_argument("--sort-by", default="mean_hits_at_k")
    compare_metrics.add_argument("--output")
    compare_metrics.set_defaults(func=cmd_compare_metrics)

    precision_filter = sub.add_parser("precision-filter")
    precision_filter.add_argument("--validation", required=True)
    precision_filter.add_argument("--target", required=True)
    precision_filter.add_argument("--output", required=True)
    precision_filter.add_argument("--report-output")
    precision_filter.add_argument("--score-col", required=True)
    precision_filter.add_argument("--target-precision", type=float, default=0.5)
    precision_filter.add_argument("--min-selected", type=int, default=1)
    precision_filter.add_argument("--group-col")
    precision_filter.add_argument("--min-group-positives", type=int, default=2)
    precision_filter.add_argument("--min-group-selected", type=int, default=1)
    precision_filter.set_defaults(func=cmd_precision_filter)

    leakage = sub.add_parser("detect-leakage")
    leakage.add_argument("--train", required=True)
    leakage.add_argument("--test", required=True)
    leakage.set_defaults(func=cmd_leakage)

    train_eval = sub.add_parser("train-eval")
    train_eval.add_argument("--train", required=True)
    train_eval.add_argument("--test", required=True)
    train_eval.add_argument("--group-col", default="patient_id")
    train_eval.add_argument("-k", type=int, default=20)
    train_eval.add_argument("--allow-exact-leakage", action="store_true")
    train_eval.add_argument("--ignore-shared-study", action="store_true")
    train_eval.add_argument("--purge-exact-overlaps", action="store_true")
    train_eval.add_argument("--uncertainty-ensemble-size", type=int, default=1)
    train_eval.add_argument("--uncertainty-penalty", type=float, default=1.0)
    train_eval.add_argument("--write-scored")
    train_eval.set_defaults(func=cmd_train_eval)

    list_datasets = sub.add_parser("list-datasets")
    list_datasets.add_argument("--manifest", default="configs/datasets.yml")
    list_datasets.set_defaults(func=cmd_list_datasets)

    download_plan = sub.add_parser("download-plan")
    download_plan.add_argument("--manifest", default="configs/datasets.yml")
    download_plan.add_argument("--dataset")
    download_plan.add_argument("--output-dir", default="data/raw")
    download_plan.set_defaults(func=cmd_download_plan)

    download = sub.add_parser("download-file")
    download.add_argument("--url", required=True)
    download.add_argument("--output", required=True)
    download.add_argument("--overwrite", action="store_true")
    download.set_defaults(func=cmd_download_file)

    normalize = sub.add_parser("normalize")
    normalize.add_argument(
        "--kind",
        choices=[
            "generic",
            "neoranking-neopep",
            "gartner",
            "tesla",
            "bigmhc",
            "cd8-multimer-2025",
            "improve-cv",
        ],
        default="generic",
    )
    normalize.add_argument("--input", required=True)
    normalize.add_argument("--output", required=True)
    normalize.add_argument("--zip-member")
    normalize.add_argument("--source-dataset", default="external")
    normalize.add_argument("--study-default")
    normalize.set_defaults(func=cmd_normalize)

    holdout = sub.add_parser("make-holdout-split")
    holdout.add_argument("table")
    holdout.add_argument("--group-col", required=True)
    holdout.add_argument("--holdout", action="append", required=True)
    holdout.add_argument("--output", required=True)
    holdout.set_defaults(func=cmd_make_holdout_split)

    portfolio = sub.add_parser("select-portfolio")
    portfolio.add_argument("table")
    portfolio.add_argument("--score-col", default="epicurus_score")
    portfolio.add_argument("-k", type=int, default=20)
    portfolio.add_argument("--max-per-hla", type=int)
    portfolio.add_argument("--max-per-gene", type=int)
    portfolio.add_argument("--min-score", type=float)
    portfolio.add_argument("--output", required=True)
    portfolio.set_defaults(func=cmd_select_portfolio)

    group_cv = sub.add_parser("group-cv")
    group_cv.add_argument("table")
    group_cv.add_argument("--group-col", required=True)
    group_cv.add_argument("--metric-group-col", default="patient_id")
    group_cv.add_argument("-k", type=int, default=20)
    group_cv.add_argument("--max-splits", type=int)
    group_cv.add_argument("--no-purge-exact-overlaps", action="store_true")
    group_cv.set_defaults(func=cmd_group_cv)

    research = sub.add_parser("research-report")
    research.add_argument("--scored", required=True)
    research.add_argument("--group-col", default="patient_id")
    research.add_argument("--score-col", default="epicurus_score")
    research.add_argument("-k", type=int, default=20)
    research.add_argument("--max-examples", type=int, default=20)
    research.add_argument("--output-dir", required=True)
    research.set_defaults(func=cmd_research_report)

    mhcflurry = sub.add_parser("add-mhcflurry-features")
    mhcflurry.add_argument("--input", required=True)
    mhcflurry.add_argument("--output", required=True)
    mhcflurry.set_defaults(func=cmd_mhcflurry_features)

    retrieval = sub.add_parser("add-retrieval-features")
    retrieval.add_argument("--input", required=True)
    retrieval.add_argument("--reference", required=True)
    retrieval.add_argument("--output", required=True)
    retrieval.add_argument("--top-k", type=int, default=5)
    retrieval.set_defaults(func=cmd_retrieval_features)

    multik_retrieval = sub.add_parser("add-multik-retrieval-features")
    multik_retrieval.add_argument("--input", required=True)
    multik_retrieval.add_argument("--reference", required=True)
    multik_retrieval.add_argument("--output", required=True)
    multik_retrieval.add_argument("--top-k", action="append", type=int)
    multik_retrieval.set_defaults(func=cmd_multik_retrieval_features)

    crossfit_retrieval = sub.add_parser("add-crossfit-retrieval-features")
    crossfit_retrieval.add_argument("--input", required=True)
    crossfit_retrieval.add_argument("--output", required=True)
    crossfit_retrieval.add_argument("--top-k", type=int, default=5)
    crossfit_retrieval.add_argument("--n-folds", type=int, default=5)
    crossfit_retrieval.add_argument("--fold-col", default="retrieval_fold")
    crossfit_retrieval.set_defaults(func=cmd_crossfit_retrieval_features)

    plm_retrieval = sub.add_parser("add-plm-retrieval-features")
    plm_retrieval.add_argument("--input", required=True)
    plm_retrieval.add_argument("--reference", required=True)
    plm_retrieval.add_argument("--output", required=True)
    plm_retrieval.add_argument("--model-name", default="facebook/esm2_t6_8M_UR50D")
    plm_retrieval.add_argument("--batch-size", type=int, default=64)
    plm_retrieval.add_argument("--device")
    plm_retrieval.add_argument("--top-k", type=int, default=5)
    plm_retrieval.set_defaults(func=cmd_plm_retrieval_features)

    plm_cache = sub.add_parser("build-plm-embedding-cache")
    plm_cache.add_argument("--input", action="append", required=True)
    plm_cache.add_argument("--output", required=True)
    plm_cache.add_argument("--model-name", default="facebook/esm2_t6_8M_UR50D")
    plm_cache.add_argument("--batch-size", type=int, default=64)
    plm_cache.add_argument("--device")
    plm_cache.set_defaults(func=cmd_build_plm_embedding_cache)

    frozen_plm_rank = sub.add_parser("frozen-plm-rank")
    frozen_plm_rank.add_argument("--train", required=True)
    frozen_plm_rank.add_argument("--validation", required=True)
    frozen_plm_rank.add_argument("--train-validation", required=True)
    frozen_plm_rank.add_argument("--target", required=True)
    frozen_plm_rank.add_argument("--embedding-cache", required=True)
    frozen_plm_rank.add_argument("--output", required=True)
    frozen_plm_rank.add_argument("--selection-output", required=True)
    frozen_plm_rank.add_argument("--group-col", default="hla_allele")
    frozen_plm_rank.add_argument("-k", type=int, default=20)
    frozen_plm_rank.set_defaults(func=cmd_frozen_plm_rank)

    normalize_vdjdb = sub.add_parser("normalize-vdjdb")
    normalize_vdjdb.add_argument("--input", required=True)
    normalize_vdjdb.add_argument("--output", required=True)
    normalize_vdjdb.add_argument("--min-score", type=int, default=0)
    normalize_vdjdb.set_defaults(func=cmd_normalize_vdjdb)

    external_recognition = sub.add_parser("add-external-recognition-features")
    external_recognition.add_argument("--input", required=True)
    external_recognition.add_argument("--reference", required=True)
    external_recognition.add_argument("--embedding-cache", required=True)
    external_recognition.add_argument("--output", required=True)
    external_recognition.add_argument("--top-k", type=int, default=5)
    external_recognition.set_defaults(func=cmd_external_recognition_features)

    screened_recognition = sub.add_parser("add-screened-recognition-features")
    screened_recognition.add_argument("--input", required=True)
    screened_recognition.add_argument("--reference", required=True)
    screened_recognition.add_argument("--embedding-cache", required=True)
    screened_recognition.add_argument("--output", required=True)
    screened_recognition.add_argument("--top-k", action="append", type=int)
    screened_recognition.set_defaults(func=cmd_screened_recognition_features)

    transfer_rank = sub.add_parser("transfer-rank")
    transfer_rank.add_argument("--external", required=True)
    transfer_rank.add_argument("--train", required=True)
    transfer_rank.add_argument("--validation", required=True)
    transfer_rank.add_argument("--embedding-cache", required=True)
    transfer_rank.add_argument("--output", required=True)
    transfer_rank.add_argument("--selection-output", required=True)
    transfer_rank.add_argument("--external-group-col", default="patient_id")
    transfer_rank.add_argument("--group-col", default="hla_allele")
    transfer_rank.add_argument("-k", type=int, default=20)
    transfer_rank.set_defaults(func=cmd_transfer_rank)

    finetune_plm_rank = sub.add_parser("finetune-plm-rank")
    finetune_plm_rank.add_argument("--external", required=True)
    finetune_plm_rank.add_argument("--train", required=True)
    finetune_plm_rank.add_argument("--validation", required=True)
    finetune_plm_rank.add_argument("--output", required=True)
    finetune_plm_rank.add_argument("--selection-output", required=True)
    finetune_plm_rank.add_argument("--model-name", default="facebook/esm2_t6_8M_UR50D")
    finetune_plm_rank.add_argument("--external-group-col", default="hla_allele")
    finetune_plm_rank.add_argument("--group-col", default="hla_allele")
    finetune_plm_rank.add_argument("--device")
    finetune_plm_rank.add_argument("-k", type=int, default=20)
    finetune_plm_rank.set_defaults(func=cmd_finetune_plm_rank)

    selector = sub.add_parser("apply-score-selector")
    selector.add_argument("--validation", required=True)
    selector.add_argument("--target", required=True)
    selector.add_argument("--output", required=True)
    selector.add_argument("--selection-output")
    selector.add_argument("--group-col", default="patient_id")
    selector.add_argument("--score-col", action="append", required=True)
    selector.add_argument("-k", type=int, default=20)
    selector.add_argument("--min-positive", type=int, default=1)
    selector.add_argument(
        "--objective",
        choices=["hits", "recall", "ndcg", "mrr", "balanced"],
        default="hits",
    )
    selector.set_defaults(func=cmd_apply_score_selector)

    blend_selector = sub.add_parser("apply-blend-selector")
    blend_selector.add_argument("--validation", required=True)
    blend_selector.add_argument("--target", required=True)
    blend_selector.add_argument("--output", required=True)
    blend_selector.add_argument("--selection-output")
    blend_selector.add_argument("--group-col", default="patient_id")
    blend_selector.add_argument("--score-col", action="append", required=True)
    blend_selector.add_argument("--pair-weight", action="append", type=float)
    blend_selector.add_argument("-k", type=int, default=20)
    blend_selector.add_argument("--min-positive", type=int, default=1)
    blend_selector.add_argument(
        "--objective",
        choices=["hits", "recall", "ndcg", "mrr", "balanced"],
        default="hits",
    )
    blend_selector.set_defaults(func=cmd_apply_blend_selector)

    guarded_blend = sub.add_parser("apply-guarded-blend-selector")
    guarded_blend.add_argument("--validation", required=True)
    guarded_blend.add_argument("--target", required=True)
    guarded_blend.add_argument("--output", required=True)
    guarded_blend.add_argument("--selection-output")
    guarded_blend.add_argument("--group-col", default="patient_id")
    guarded_blend.add_argument("--score-col", action="append", required=True)
    guarded_blend.add_argument("--pair-weight", action="append", type=float)
    guarded_blend.add_argument("-k", type=int, default=20)
    guarded_blend.add_argument("--min-positive", type=int, default=1)
    guarded_blend.add_argument(
        "--baseline-objective",
        choices=["hits", "recall", "ndcg", "mrr", "balanced"],
        default="ndcg",
    )
    guarded_blend.add_argument(
        "--objective",
        choices=["hits", "recall", "ndcg", "mrr", "balanced"],
        default="mrr",
    )
    guarded_blend.add_argument("--guard-metric", default="mean_hits_at_k")
    guarded_blend.add_argument("--min-guard-delta", type=float, default=0.0)
    guarded_blend.set_defaults(func=cmd_apply_guarded_blend_selector)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
