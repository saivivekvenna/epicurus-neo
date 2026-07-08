from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from epicurus_neo.benchmark import train_and_evaluate
from epicurus_neo.data_manifest import load_dataset_manifest
from epicurus_neo.leakage import detect_exact_leakage
from epicurus_neo.metrics import group_metrics, summarize_group_metrics
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


def cmd_leakage(args: argparse.Namespace) -> int:
    train = _load_table(Path(args.train))
    test = _load_table(Path(args.test))
    report = detect_exact_leakage(train, test)
    print(json.dumps(report.__dict__, indent=2))
    return 1 if report.has_leakage else 0


def cmd_train_eval(args: argparse.Namespace) -> int:
    train = _load_table(Path(args.train))
    test = _load_table(Path(args.test))
    result = train_and_evaluate(
        train,
        test,
        group_col=args.group_col,
        k=args.k,
        allow_exact_leakage=args.allow_exact_leakage,
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
    train_eval.add_argument("--write-scored")
    train_eval.set_defaults(func=cmd_train_eval)

    list_datasets = sub.add_parser("list-datasets")
    list_datasets.add_argument("--manifest", default="configs/datasets.yml")
    list_datasets.set_defaults(func=cmd_list_datasets)

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

    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    return int(args.func(args))


if __name__ == "__main__":
    raise SystemExit(main())
