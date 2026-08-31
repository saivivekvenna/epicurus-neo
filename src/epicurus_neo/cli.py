from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from epicurus_neo.portfolio import PortfolioConstraints, select_portfolio
from epicurus_neo.schema import validate_schema


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


def cmd_validate_product_input(args: argparse.Namespace) -> int:
    from epicurus_neo.contracts import validate_candidate_contract
    from epicurus_neo.product import load_product_candidates

    try:
        candidates = load_product_candidates(
            args.input,
            patient_id=args.patient_id,
            rna_evidence_path=args.rna_evidence,
        )
    except (ValueError, TypeError) as exc:
        print(json.dumps({"ok": False, "error": str(exc)}, indent=2))
        return 1
    report = validate_candidate_contract(candidates)
    print(json.dumps(report.__dict__ | {"ok": report.ok}, indent=2))
    return 0 if report.ok else 1


def cmd_run_patient(args: argparse.Namespace) -> int:
    from epicurus_neo.product import InferenceConfig, run_product_inference

    outputs = run_product_inference(
        args.input,
        args.output_dir,
        patient_id=args.patient_id,
        rna_evidence_path=args.rna_evidence,
        config=InferenceConfig(
            k=args.k,
            max_per_mutation=args.max_per_mutation,
            max_per_gene=args.max_per_gene,
            max_per_hla=args.max_per_hla,
            core_threshold=args.core_threshold,
            supporting_threshold=args.supporting_threshold,
            apply_validity_gate=not args.disable_validity_gate,
        ),
    )
    print(json.dumps(outputs, indent=2))
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


def cmd_run_pipeline(args: argparse.Namespace) -> int:
    from epicurus_neo.pipeline import load_pipeline_config, run_pipeline

    config = load_pipeline_config(args.config)
    result = run_pipeline(
        config,
        output_dir=args.output_dir,
        start=args.start,
        stop=args.stop,
        force=args.force,
    )
    print(json.dumps(result.to_summary(), indent=2))
    return 0 if result.ok else 1


def cmd_doctor(args: argparse.Namespace) -> int:
    from epicurus_neo.pipeline import readiness_report

    report = readiness_report(bundle_dir=args.bundle_dir)
    print(json.dumps(report, indent=2))
    return 0 if report["ready"] else 1


def cmd_references(args: argparse.Namespace) -> int:
    from epicurus_neo.pipeline import scaffold_references

    result = scaffold_references(args.dest)
    print(json.dumps(result, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="epicurus-neo",
        description="Prioritize personalized cancer-vaccine neoantigens from WES/RNA-seq.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    pipeline = sub.add_parser(
        "run-pipeline",
        help="Run the full A-to-Z pipeline: raw WES/RNA -> ranked <=20 portfolio.",
    )
    pipeline.add_argument("--config", required=True, help="Patient pipeline config (YAML).")
    pipeline.add_argument("--output-dir", required=True)
    pipeline.add_argument("--start", help="Optional first stage to run (resume).")
    pipeline.add_argument("--stop", help="Optional last stage to run.")
    pipeline.add_argument(
        "--force",
        action="store_true",
        help="Recompute stages even if a valid cached artifact exists.",
    )
    pipeline.set_defaults(func=cmd_run_pipeline)

    doctor = sub.add_parser(
        "doctor",
        help="Check external tools and reference data are installed and ready.",
    )
    doctor.add_argument("--bundle-dir", help="Reference bundle directory to check.")
    doctor.set_defaults(func=cmd_doctor)

    references = sub.add_parser(
        "references",
        help="Scaffold the reference-bundle directory and write install instructions.",
    )
    references.add_argument("--dest", required=True, help="Directory to create the bundle in.")
    references.set_defaults(func=cmd_references)

    validate = sub.add_parser("validate-schema")
    validate.add_argument("table")
    validate.set_defaults(func=cmd_validate)

    validate_product = sub.add_parser("validate-patient-input")
    validate_product.add_argument("--input", required=True)
    validate_product.add_argument("--patient-id")
    validate_product.add_argument("--rna-evidence")
    validate_product.set_defaults(func=cmd_validate_product_input)

    run_patient = sub.add_parser(
        "run-patient",
        help="Prioritize an existing pVACseq candidate table (start at the ranking stage).",
    )
    run_patient.add_argument("--input", required=True)
    run_patient.add_argument("--patient-id")
    run_patient.add_argument("--rna-evidence")
    run_patient.add_argument("--output-dir", required=True)
    run_patient.add_argument("-k", type=int, default=20)
    run_patient.add_argument(
        "--max-per-mutation",
        type=int,
        default=1,
        help="Maximum selected peptide-HLA routes per mutation (default: 1).",
    )
    run_patient.add_argument("--max-per-gene", type=int, default=4)
    run_patient.add_argument("--max-per-hla", type=int)
    run_patient.add_argument("--core-threshold", type=float, default=0.55)
    run_patient.add_argument("--supporting-threshold", type=float, default=0.35)
    run_patient.add_argument(
        "--disable-validity-gate",
        action="store_true",
        help="Disable the default deterministic biological-validity gate for an audit comparison.",
    )
    run_patient.set_defaults(func=cmd_run_patient)

    portfolio = sub.add_parser("select-portfolio")
    portfolio.add_argument("table")
    portfolio.add_argument("--score-col", default="epicurus_neo_score")
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
