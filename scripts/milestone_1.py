#!/usr/bin/env python3
"""Reproduce IMPROVE regressions or generate the ten masking question sets."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path

from benchmark.ablation import generate_masking_sets, write_masking_sets
from benchmark.improve import load_improve_data, regression_values
from benchmark.reaudit import reaudit_repository, write_markdown


def _jsonable(value):
    if is_dataclass(value):
        return {key: _jsonable(item) for key, item in asdict(value).items()}
    if isinstance(value, dict):
        return {str(key): _jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_jsonable(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("improve_repo", type=Path, help="clone of SRHgroup/IMPROVE_paper")
    subparsers = parser.add_subparsers(dest="command", required=True)
    subparsers.add_parser("verify")
    generate = subparsers.add_parser("generate-ablation")
    generate.add_argument("output_dir", type=Path)
    reaudit = subparsers.add_parser("reaudit")
    reaudit.add_argument("output", type=Path)
    args = parser.parse_args()

    if args.command == "verify":
        print(json.dumps(_jsonable(regression_values(args.improve_repo)), indent=2, sort_keys=True))
        return 0
    if args.command == "generate-ablation":
        data = load_improve_data(args.improve_repo)
        sets = generate_masking_sets(data)
        write_masking_sets(sets, args.output_dir)
        print(f"wrote {len(sets)} seeds to {args.output_dir}")
        return 0
    rows = reaudit_repository(Path.cwd(), args.improve_repo)
    write_markdown(rows, args.output)
    print(f"wrote {len(rows)} frozen-score audits to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
