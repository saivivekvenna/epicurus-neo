#!/usr/bin/env python3
"""CLI for the frozen pVACtools/Vaxrank controlled-input comparator track."""

from __future__ import annotations

import argparse
import json

from benchmark.controlled_comparators import freeze_vaxrank_portfolio, prepare_common_bundle
from benchmark.nextneopi_comparator import freeze_native_portfolio


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    prepare = sub.add_parser("prepare")
    prepare.add_argument("patient_id")
    for name in (
        "pass_vcf", "pass_vcf_index", "rna_bam", "rna_bam_index", "hla_panel",
        "pvac_ready_vcf", "pvac_ready_vcf_index",
    ):
        prepare.add_argument(f"--{name.replace('_', '-')}", required=True)
    prepare.add_argument("--output-dir", required=True)
    vaxrank = sub.add_parser("freeze-vaxrank")
    vaxrank.add_argument("ranked_csv")
    vaxrank.add_argument("--output-dir", required=True)
    pvac = sub.add_parser("freeze-pvacseq")
    pvac.add_argument("aggregate_tsv")
    pvac.add_argument("--output-dir", required=True)
    args = parser.parse_args(argv)

    if args.command == "prepare":
        result = prepare_common_bundle(
            patient_id=args.patient_id,
            pass_vcf=args.pass_vcf,
            pass_vcf_index=args.pass_vcf_index,
            rna_bam=args.rna_bam,
            rna_bam_index=args.rna_bam_index,
            hla_panel=args.hla_panel,
            pvac_ready_vcf=args.pvac_ready_vcf,
            pvac_ready_vcf_index=args.pvac_ready_vcf_index,
            output_dir=args.output_dir,
        )
    elif args.command == "freeze-vaxrank":
        result = freeze_vaxrank_portfolio(args.ranked_csv, args.output_dir)
    else:
        result = freeze_native_portfolio(args.aggregate_tsv, args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
