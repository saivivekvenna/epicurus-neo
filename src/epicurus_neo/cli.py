from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from epicurus_neo.benchmark import evaluate_score_columns, train_and_evaluate
from epicurus_neo.auto_research import build_failure_report, write_research_artifacts
from epicurus_neo.data_manifest import load_dataset_manifest
from epicurus_neo.download import dataset_file_plans, download_file
from epicurus_neo.experiment import grouped_cross_validate, summarize_cross_validation
from epicurus_neo.leakage import detect_exact_leakage, purge_train_overlaps
from epicurus_neo.metrics import group_metrics, summarize_group_metrics
from epicurus_neo.normalize import (
    normalize_bigmhc_table,
    normalize_candidate_table,
    normalize_gartner_table,
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
    per_group = group_metrics(frame, group_col=args.group_col, score_col=args.score_col, k=args.k)
    print(json.dumps(summarize_group_metrics(per_group), indent=2))
    return 0


def cmd_score_report(args: argparse.Namespace) -> int:
    frame = _load_table(Path(args.table))
    results = [
        {"score_col": item.score_col, "summary": item.summary}
        for item in evaluate_score_columns(
            frame,
            group_col=args.group_col,
            score_columns=args.score_col,
            k=args.k,
        )
    ]
    payload = {
        "table": args.table,
        "group_col": args.group_col,
        "k": args.k,
        "benchmarks": results,
    }
    text = json.dumps(payload, indent=2)
    if args.output:
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        Path(args.output).write_text(text + "\n")
    print(text)
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
    print(json.dumps([plan.__dict__ | {"output_path": str(plan.output_path)} for plan in plans], indent=2))
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
    metrics.add_argument("-k", type=int, default=20)
    metrics.set_defaults(func=cmd_metrics)

    score_report = sub.add_parser("score-report")
    score_report.add_argument("table")
    score_report.add_argument("--group-col", default="patient_id")
    score_report.add_argument("--score-col", action="append", required=True)
    score_report.add_argument("-k", type=int, default=20)
    score_report.add_argument("--output")
    score_report.set_defaults(func=cmd_score_report)

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
        choices=["generic", "neoranking-neopep", "gartner", "tesla", "bigmhc"],
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

    selector = sub.add_parser("apply-score-selector")
    selector.add_argument("--validation", required=True)
    selector.add_argument("--target", required=True)
    selector.add_argument("--output", required=True)
    selector.add_argument("--selection-output")
    selector.add_argument("--group-col", default="patient_id")
    selector.add_argument("--score-col", action="append", required=True)
    selector.add_argument("-k", type=int, default=20)
    selector.add_argument("--min-positive", type=int, default=1)
    selector.set_defaults(func=cmd_apply_score_selector)

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
