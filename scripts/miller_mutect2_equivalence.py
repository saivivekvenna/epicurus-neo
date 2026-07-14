"""Strict, outcome-blind serial-vs-scattered Mutect2 equivalence audit.

The audit filters the scattered raw call with the preregistered FilterMutectCalls
path, derives its PASS VCF, and then compares serial and scattered raw, filtered,
and PASS records byte-for-byte after VCF header removal.  Header timestamps and
command lines are intentionally outside the comparison; every variant/sample
field remains inside it.  Any difference refuses acceleration deployment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(payload, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def require_vcf(path: Path, bcftools: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing/empty VCF: {path}")
    probe = subprocess.run(
        [bcftools, "view", "-h", str(path)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
    )
    if probe.returncode != 0:
        raise ValueError(f"unreadable VCF {path}: {probe.stderr.strip()}")


def record_lines(path: Path, bcftools: str) -> list[str]:
    completed = subprocess.run(
        [bcftools, "view", "-H", str(path)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=True,
    )
    return completed.stdout.splitlines()


def variant_key(record: str) -> str:
    fields = record.split("\t")
    if len(fields) < 8:
        raise ValueError("malformed VCF record")
    return ":".join((fields[0], fields[1], fields[3], fields[4]))


def compare_records(serial: list[str], scattered: list[str]) -> dict:
    serial_keys = [variant_key(record) for record in serial]
    scattered_keys = [variant_key(record) for record in scattered]
    serial_set, scattered_set = set(serial_keys), set(scattered_keys)
    mismatched = [
        {"index": index, "serial": left, "scattered": right}
        for index, (left, right) in enumerate(zip(serial, scattered))
        if left != right
    ]
    if len(serial) != len(scattered):
        mismatched.append(
            {"length_mismatch": {"serial": len(serial), "scattered": len(scattered)}}
        )
    serial_digest = hashlib.sha256(("\n".join(serial) + ("\n" if serial else "")).encode()).hexdigest()
    scattered_digest = hashlib.sha256(
        ("\n".join(scattered) + ("\n" if scattered else "")).encode()
    ).hexdigest()
    return {
        "serial_records": len(serial),
        "scattered_records": len(scattered),
        "site_keys_equal_in_order": serial_keys == scattered_keys,
        "missing_from_scattered": sorted(serial_set - scattered_set)[:20],
        "extra_in_scattered": sorted(scattered_set - serial_set)[:20],
        "record_stream_sha256": {"serial": serial_digest, "scattered": scattered_digest},
        "records_exact": serial == scattered,
        "first_record_mismatches": mismatched[:20],
    }


def filter_scattered(
    *,
    gatk: Path,
    reference: Path,
    scattered_raw: Path,
    scattered_filtered: Path,
    scattered_pass: Path,
    bcftools: str,
) -> None:
    stats = Path(f"{scattered_raw}.stats")
    if not stats.is_file() or stats.stat().st_size == 0:
        raise ValueError(f"missing scattered Mutect2 stats: {stats}")
    if not scattered_filtered.exists():
        subprocess.run(
            [
                str(gatk),
                "FilterMutectCalls",
                "-R",
                str(reference),
                "-V",
                str(scattered_raw),
                "-O",
                str(scattered_filtered),
            ],
            check=True,
        )
    require_vcf(scattered_filtered, bcftools)
    if not scattered_pass.exists():
        subprocess.run(
            [bcftools, "view", "-f", "PASS", str(scattered_filtered), "-Oz", "-o", str(scattered_pass)],
            check=True,
        )
        subprocess.run([bcftools, "index", "-t", "-f", str(scattered_pass)], check=True)
    require_vcf(scattered_pass, bcftools)


def run(args: argparse.Namespace) -> dict:
    for path in (
        args.gatk,
        args.reference,
        args.serial_raw,
        args.serial_filtered,
        args.serial_pass,
        args.scattered_raw,
    ):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing/empty input: {path}")
    for path in (args.serial_raw, args.serial_filtered, args.serial_pass, args.scattered_raw):
        require_vcf(path, args.bcftools)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    scattered_filtered = args.output_dir / "scattered.filtered.vcf.gz"
    scattered_pass = args.output_dir / "scattered.pass.vcf.gz"
    filter_scattered(
        gatk=args.gatk,
        reference=args.reference,
        scattered_raw=args.scattered_raw,
        scattered_filtered=scattered_filtered,
        scattered_pass=scattered_pass,
        bcftools=args.bcftools,
    )

    pairs = {
        "raw": (args.serial_raw, args.scattered_raw),
        "filtered": (args.serial_filtered, scattered_filtered),
        "pass": (args.serial_pass, scattered_pass),
    }
    comparisons = {
        name: compare_records(record_lines(left, args.bcftools), record_lines(right, args.bcftools))
        for name, (left, right) in pairs.items()
    }
    scattered_pass_records = record_lines(scattered_pass, args.bcftools)
    nonpass = [record for record in scattered_pass_records if record.split("\t")[6] != "PASS"]
    equivalent = all(result["records_exact"] for result in comparisons.values()) and not nonpass
    report = {
        "status": "EQUIVALENT" if equivalent else "NOT_EQUIVALENT",
        "deploy_scatter": equivalent,
        "isolation": "LOCKED_TEST: no recognition/outcome data read",
        "comparison_contract": "all non-header VCF fields exact, in coordinate order",
        "inputs": {
            name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
            for name, path in {
                "serial_raw": args.serial_raw,
                "serial_filtered": args.serial_filtered,
                "serial_pass": args.serial_pass,
                "scattered_raw": args.scattered_raw,
                "scattered_raw_stats": Path(f"{args.scattered_raw}.stats"),
            }.items()
        },
        "comparisons": comparisons,
        "scattered_pass_nonpass_records": len(nonpass),
        "derived": {
            "scattered_filtered": str(scattered_filtered.resolve()),
            "scattered_pass": str(scattered_pass.resolve()),
        },
    }
    atomic_json(args.output_dir / "MUTECT2_SCATTER_EQUIVALENCE.json", report)
    return report


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gatk", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--serial-raw", type=Path, required=True)
    p.add_argument("--serial-filtered", type=Path, required=True)
    p.add_argument("--serial-pass", type=Path, required=True)
    p.add_argument("--scattered-raw", type=Path, required=True)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--bcftools", default="bcftools")
    return p


def main(argv: list[str] | None = None) -> int:
    report = run(parser().parse_args(argv))
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["deploy_scatter"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
