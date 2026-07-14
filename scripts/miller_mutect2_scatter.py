"""Label-blind, resumable Mutect2 scatter/gather for Miller WES reconstruction.

This module changes execution only, not the biological calling policy: every shard
receives the same tumor, matched normal, reference, CDS intervals and zero interval
padding as the serial preregistered call.  Intervals are never subdivided, which
keeps an original CDS interval wholly within one Mutect2 assembly traversal.

The caller is deliberately not enabled by default yet.  It must first reproduce the
already-frozen Hu_315 serial callset before it may be selected for future patients.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def split_command(gatk: Path, reference: Path, intervals: Path, count: int, output: Path) -> list[str]:
    return [
        str(gatk),
        "SplitIntervals",
        "-R",
        str(reference),
        "-L",
        str(intervals),
        "--interval-padding",
        "0",
        "--scatter-count",
        str(count),
        "--subdivision-mode",
        "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
        "-O",
        str(output),
    ]


def mutect_command(
    gatk: Path,
    reference: Path,
    tumor_bam: Path,
    normal_bam: Path,
    normal_sample: str,
    interval_list: Path,
    output: Path,
    java_heap_gb: int,
) -> list[str]:
    return [
        str(gatk),
        "--java-options",
        f"-Xmx{java_heap_gb}g",
        "Mutect2",
        "-R",
        str(reference),
        "-I",
        str(tumor_bam),
        "-I",
        str(normal_bam),
        "--normal",
        normal_sample,
        "-L",
        str(interval_list),
        "--interval-padding",
        "0",
        "-O",
        str(output),
    ]


def gather_command(gatk: Path, shards: list[Path], output: Path) -> list[str]:
    command = [str(gatk), "GatherVcfs"]
    for shard in shards:
        command.extend(["-I", str(shard)])
    command.extend(["-O", str(output)])
    return command


def merge_stats_command(gatk: Path, shards: list[Path], output: Path) -> list[str]:
    command = [str(gatk), "MergeMutectStats"]
    for shard in shards:
        command.extend(["--stats", str(Path(f"{shard}.stats"))])
    command.extend(["-O", str(output)])
    return command


def run_logged(command: list[str], log: Path) -> None:
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("ab") as handle:
        subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)


def vcf_ready(vcf: Path, bcftools: str) -> bool:
    if not vcf.is_file() or vcf.stat().st_size == 0:
        return False
    if not Path(f"{vcf}.stats").is_file() or Path(f"{vcf}.stats").stat().st_size == 0:
        return False
    return subprocess.run(
        [bcftools, "view", "-h", str(vcf)],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    ).returncode == 0


def require_file(path: Path, description: str) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing/empty {description}: {path}")


def run(args: argparse.Namespace) -> dict:
    for path, description in (
        (args.gatk, "GATK launcher"),
        (args.reference, "reference"),
        (args.intervals, "interval source"),
        (args.tumor_bam, "tumor BAM"),
        (args.normal_bam, "normal BAM"),
    ):
        require_file(path, description)
    if args.scatter_count < 2:
        raise ValueError("scatter-count must be >=2; use the preregistered serial caller for one shard")
    if args.java_heap_gb < 1:
        raise ValueError("java-heap-gb must be positive")

    output = args.output.resolve()
    final_stats = Path(f"{output}.stats")
    if output.exists():
        if vcf_ready(output, args.bcftools):
            return {"status": "ALREADY_COMPLETE", "output": str(output)}
        raise RuntimeError(f"refusing to replace incomplete/corrupt authoritative output: {output}")
    # The VCF is the publication sentinel and is always moved last below.  Sidecars
    # without it can only be private debris from an interrupted publication.
    Path(f"{output}.tbi").unlink(missing_ok=True)
    final_stats.unlink(missing_ok=True)

    work = output.parent / f".{output.name}.scatter-work"
    interval_dir = work / "intervals"
    shard_dir = work / "shards"
    log_dir = work / "logs"
    work.mkdir(parents=True, exist_ok=True)
    shard_dir.mkdir(exist_ok=True)
    log_dir.mkdir(exist_ok=True)

    interval_lists = sorted(interval_dir.glob("*.interval_list")) if interval_dir.exists() else []
    if not interval_lists:
        if interval_dir.exists():
            shutil.rmtree(interval_dir)
        interval_dir.mkdir()
        run_logged(
            split_command(args.gatk, args.reference, args.intervals, args.scatter_count, interval_dir),
            log_dir / "split.log",
        )
        interval_lists = sorted(interval_dir.glob("*.interval_list"))
    if len(interval_lists) != args.scatter_count:
        raise RuntimeError(
            f"SplitIntervals produced {len(interval_lists)} shards; expected {args.scatter_count}"
        )

    shards = [shard_dir / f"shard_{index:04d}.vcf.gz" for index in range(len(interval_lists))]

    def call_one(item: tuple[int, Path, Path]) -> None:
        index, interval_list, shard = item
        if vcf_ready(shard, args.bcftools):
            return
        for stale in (shard, Path(f"{shard}.tbi"), Path(f"{shard}.stats")):
            stale.unlink(missing_ok=True)
        run_logged(
            mutect_command(
                args.gatk,
                args.reference,
                args.tumor_bam,
                args.normal_bam,
                args.normal_sample,
                interval_list,
                shard,
                args.java_heap_gb,
            ),
            log_dir / f"mutect_{index:04d}.log",
        )
        if not vcf_ready(shard, args.bcftools):
            raise RuntimeError(f"Mutect2 shard failed validation: {shard}")

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.scatter_count) as pool:
        futures = [pool.submit(call_one, item) for item in zip(range(len(shards)), interval_lists, shards)]
        for future in futures:
            future.result()

    partial_dir = Path(tempfile.mkdtemp(prefix="gather-", dir=output.parent))
    try:
        partial_vcf = partial_dir / output.name
        partial_stats = Path(f"{partial_vcf}.stats")
        run_logged(gather_command(args.gatk, shards, partial_vcf), log_dir / "gather.log")
        run_logged(merge_stats_command(args.gatk, shards, partial_stats), log_dir / "merge_stats.log")
        if not vcf_ready(partial_vcf, args.bcftools):
            raise RuntimeError("gathered Mutect2 VCF/stats failed validation")
        subprocess.run([args.bcftools, "index", "-t", "-f", str(partial_vcf)], check=True)
        # Publish sidecars first and the authoritative VCF last.  A crash can never
        # leave a final-named VCF that a resume might mistake for a complete call.
        os.replace(partial_stats, final_stats)
        os.replace(Path(f"{partial_vcf}.tbi"), Path(f"{output}.tbi"))
        os.replace(partial_vcf, output)
    finally:
        shutil.rmtree(partial_dir, ignore_errors=True)

    manifest = {
        "status": "COMPLETE_UNVALIDATED_ACCELERATION",
        "isolation": "LOCKED_TEST: no recognition/outcome data read",
        "semantic_contract": {
            "caller": "GATK Mutect2 matched-normal",
            "interval_source": str(args.intervals.resolve()),
            "interval_padding": 0,
            "subdivision_mode": "BALANCING_WITHOUT_INTERVAL_SUBDIVISION",
            "scatter_count": args.scatter_count,
        },
        "inputs": {
            "reference": str(args.reference.resolve()),
            "tumor_bam": str(args.tumor_bam.resolve()),
            "normal_bam": str(args.normal_bam.resolve()),
            "normal_sample": args.normal_sample,
        },
        "interval_sha256": {p.name: sha256_file(p) for p in interval_lists},
        "output": str(output),
        "output_sha256": sha256_file(output),
        "validation_required": "Exact Hu_315 serial-vs-scattered callset equivalence before deployment",
    }
    manifest_path = output.parent / f"{output.name}.scatter-manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    if not args.keep_shards:
        shutil.rmtree(work)
    return manifest


def parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser()
    p.add_argument("--gatk", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--intervals", type=Path, required=True)
    p.add_argument("--tumor-bam", type=Path, required=True)
    p.add_argument("--normal-bam", type=Path, required=True)
    p.add_argument("--normal-sample", required=True)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--scatter-count", type=int, default=4)
    p.add_argument("--java-heap-gb", type=int, default=4)
    p.add_argument("--bcftools", default="bcftools")
    p.add_argument("--keep-shards", action="store_true")
    return p


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    print(json.dumps(run(args), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
