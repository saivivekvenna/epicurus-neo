"""Reusable, label-blind Miller patient download/FASTQ/expression driver.

This is the parameterized replacement for the Hu_287-only plumbing. It deliberately
does not import or open the Miller recognition-label table.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from benchmark.miller_download import download_patient
from benchmark.miller_ingest import SRA_RUNINFO_FIXTURE
from benchmark.miller_patient import load_patient
from scripts.miller_reconstruct import convert_sra, quant_rna


ROOT = Path(__file__).resolve().parents[1]
MAX_LOCAL_THREADS = 12

# Reconstruction stages delegate to the (now patient-parameterized) Hu_287 shell scripts. They read the
# per-patient runs/sample names ONLY from the environment; the defaults inside each script keep Hu_287
# backward-compatible when invoked with no env.
_STAGE_SCRIPTS = {
    "wes": ROOT / "scripts/miller_hu287_somatic.sh",
    "hla": ROOT / "scripts/miller_hu287_hla.sh",
    "rna": ROOT / "scripts/miller_hu287_rna.sh",
}


def script_for(command: str) -> Path:
    """Absolute path of the shell script that runs a reconstruction stage (wes/hla/rna)."""
    return _STAGE_SCRIPTS[command]


def script_env(patient) -> dict:
    """The label-blind env contract passed to every reconstruction shell script for ``patient``.

    These four variables fully parameterize the frozen Hu_287 scripts: PATIENT_ID drives sample names,
    slugged raw dirs, and the milestone-7-vs-8 provenance branch; the three run accessions drive alignment
    inputs. NO recognition-label column is ever read to build this env.
    """
    return {
        "PATIENT_ID": patient.patient_id,
        "NORMAL_EXOME_RUN": patient.normal_exome_run,
        "TUMOR_EXOME_RUN": patient.tumor_exome_run,
        "TUMOR_RNA_RUN": patient.tumor_rna_run,
    }


def patient_manifest(patient_id: str) -> dict:
    p = load_patient(patient_id)
    return {
        "patient_id": p.patient_id,
        "label_columns_read": [],
        "runs": {
            "normal_exome": p.normal_exome_run,
            "tumor_exome": p.tumor_exome_run,
            "tumor_rna": p.tumor_rna_run,
        },
        "raw_dir": str(p.raw_dir),
        "artifact_dir": str(p.artifact_dir),
        "sample_names": {"normal": p.normal_sample, "tumor": p.tumor_sample},
    }


def default_threads() -> int:
    """Use the high-core local profile while retaining headroom for macOS and piped sort stages."""
    override = os.environ.get("EPICURUS_THREADS")
    if override is not None:
        value = int(override)
        if value < 1:
            raise ValueError("EPICURUS_THREADS must be positive")
        return value
    return max(1, min(MAX_LOCAL_THREADS, os.cpu_count() or 4))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("patient_id")
    parser.add_argument("command", choices=("metadata", "download", "convert", "quant", "wes", "hla", "rna"))
    parser.add_argument("--threads", type=int, default=None)
    args = parser.parse_args(argv)
    threads = default_threads() if args.threads is None else args.threads
    if threads < 1:
        parser.error("--threads must be positive")
    p = load_patient(args.patient_id)
    p.artifact_dir.mkdir(parents=True, exist_ok=True)

    if args.command == "metadata":
        result = patient_manifest(args.patient_id)
        (p.artifact_dir / "PATIENT_INPUTS.json").write_text(json.dumps(result, indent=2) + "\n")
    elif args.command == "download":
        result = download_patient(p.patient_id, p.raw_dir, SRA_RUNINFO_FIXTURE)
    elif args.command == "convert":
        fqdir = p.raw_dir / "fastq"
        expected = (p.normal_exome_run, p.tumor_exome_run, p.tumor_rna_run)
        result = [convert_sra(p.raw_dir / f"{run}.sra", fqdir, threads=threads) for run in expected]
        (p.artifact_dir / "CONVERT_PROVENANCE.json").write_text(json.dumps(result, indent=2) + "\n")
    elif args.command == "quant":
        result = quant_rna(
            p.raw_dir / "fastq",
            p.raw_dir / "salmon_quant",
            threads=threads,
            rna_run=p.tumor_rna_run,
        )
        (p.artifact_dir / "EXPRESSION_QUANT.json").write_text(json.dumps(result, indent=2) + "\n")

    else:
        script = script_for(args.command)
        env = os.environ.copy()
        env.update(script_env(p))
        env["THREADS"] = str(threads)
        subprocess.run([str(script)], cwd=ROOT, env=env, check=True)
        result = {"patient_id": p.patient_id, "stage": args.command, "status": "OK"}

    print(json.dumps(result, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
