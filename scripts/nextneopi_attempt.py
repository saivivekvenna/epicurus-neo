#!/usr/bin/env python3
"""Finalize a checksum-bound nextNEOpi Track-A execution attempt."""

from __future__ import annotations

import argparse
import json

from benchmark.nextneopi_attempt import finalize_non_success, finalize_success


def _common(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("patient_id")
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--source", required=True, help="exact pinned nextNEOpi.nf")
    parser.add_argument("--config", required=True)
    parser.add_argument("--instrumentation", required=True)
    parser.add_argument("--execution-log", required=True)
    parser.add_argument("--execution-trace", required=True)
    parser.add_argument("--execution-command", required=True)
    parser.add_argument("--started-at", required=True, help="ISO-8601 timestamp with timezone")
    parser.add_argument("--finished-at", required=True, help="ISO-8601 timestamp with timezone")
    parser.add_argument("--runtime", required=True, help="container/host runtime identity")
    parser.add_argument("--output-dir", required=True)


def main() -> int:
    parser = argparse.ArgumentParser()
    commands = parser.add_subparsers(dest="command", required=True)
    success = commands.add_parser("success")
    _common(success)
    success.add_argument("--aggregate", required=True)
    success.add_argument("--pvac-input-vcf", required=True)
    success.add_argument("--pvac-input-vcf-index", required=True)
    success.add_argument("--reference", required=True)
    failed = commands.add_parser("non-success")
    _common(failed)
    failed.add_argument("--status", choices=("FAILED", "ABSTAINED"), required=True)
    failed.add_argument("--exit-code", type=int)
    failed.add_argument("--reason", required=True)
    args = parser.parse_args()
    common = dict(
        patient_id=args.patient_id, input_manifest=args.input_manifest, source=args.source,
        config=args.config, instrumentation=args.instrumentation,
        execution_log=args.execution_log, execution_trace=args.execution_trace,
        execution_command=args.execution_command, started_at=args.started_at,
        finished_at=args.finished_at, runtime=args.runtime, output_dir=args.output_dir,
    )
    if args.command == "success":
        result = finalize_success(
            **common, aggregate=args.aggregate, pvac_input_vcf=args.pvac_input_vcf,
            pvac_input_vcf_index=args.pvac_input_vcf_index, reference=args.reference,
        )
    else:
        result = finalize_non_success(
            **common, execution_status=args.status, exit_code=args.exit_code,
            reason=args.reason,
        )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
