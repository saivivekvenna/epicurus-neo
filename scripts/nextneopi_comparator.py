#!/usr/bin/env python3
"""CLI for the locked nextNEOpi full-pipeline comparator boundary."""

from __future__ import annotations

import argparse
import json

from benchmark.nextneopi_comparator import freeze_native_portfolio, prepare_batch


def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare")
    prepare.add_argument("patient_id")
    for name in ("tumor_r1", "tumor_r2", "normal_r1", "normal_r2", "rna_r1", "rna_r2"):
        prepare.add_argument(f"--{name.replace('_', '-')}", required=True)
    prepare.add_argument("--conversion-provenance", required=True)
    prepare.add_argument("--output-dir", required=True)
    prepare.add_argument("--sex", choices=("female", "male", "NA"), default="NA")
    freeze = sub.add_parser("freeze-portfolio")
    freeze.add_argument("aggregate_tsv")
    freeze.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    if args.command == "prepare":
        result = prepare_batch(
            patient_id=args.patient_id,
            tumor_r1=args.tumor_r1,
            tumor_r2=args.tumor_r2,
            normal_r1=args.normal_r1,
            normal_r2=args.normal_r2,
            rna_r1=args.rna_r1,
            rna_r2=args.rna_r2,
            conversion_provenance=args.conversion_provenance,
            output_dir=args.output_dir,
            sex=args.sex,
        )
    else:
        result = freeze_native_portfolio(args.aggregate_tsv, args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
