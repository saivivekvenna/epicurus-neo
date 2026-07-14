from pathlib import Path

import pytest

from scripts import miller_mutect2_scatter as scatter


def test_split_never_subdivides_original_intervals():
    command = scatter.split_command(Path("gatk"), Path("ref.fa"), Path("cds.bed"), 4, Path("out"))
    assert command[command.index("--subdivision-mode") + 1] == "BALANCING_WITHOUT_INTERVAL_SUBDIVISION"
    assert command[command.index("--interval-padding") + 1] == "0"
    assert command[command.index("--scatter-count") + 1] == "4"


def test_mutect_shards_preserve_matched_normal_and_zero_padding():
    command = scatter.mutect_command(
        Path("gatk"),
        Path("ref.fa"),
        Path("tumor.bam"),
        Path("normal.bam"),
        "Hu_X_N",
        Path("one.interval_list"),
        Path("one.vcf.gz"),
        3,
    )
    assert command[:3] == ["gatk", "--java-options", "-Xmx3g"]
    assert command[command.index("--normal") + 1] == "Hu_X_N"
    assert command[command.index("--interval-padding") + 1] == "0"
    assert command.count("-I") == 2


def test_gather_and_stats_preserve_shard_order():
    shards = [Path("00.vcf.gz"), Path("01.vcf.gz")]
    gather = scatter.gather_command(Path("gatk"), shards, Path("all.vcf.gz"))
    stats = scatter.merge_stats_command(Path("gatk"), shards, Path("all.vcf.gz.stats"))
    assert gather == ["gatk", "GatherVcfs", "-I", "00.vcf.gz", "-I", "01.vcf.gz", "-O", "all.vcf.gz"]
    assert stats == [
        "gatk",
        "MergeMutectStats",
        "--stats",
        "00.vcf.gz.stats",
        "--stats",
        "01.vcf.gz.stats",
        "-O",
        "all.vcf.gz.stats",
    ]


def test_run_fails_closed_before_subprocess_on_bad_inputs(tmp_path):
    args = scatter.parser().parse_args(
        [
            "--gatk", str(tmp_path / "gatk"),
            "--reference", str(tmp_path / "ref.fa"),
            "--intervals", str(tmp_path / "cds.bed"),
            "--tumor-bam", str(tmp_path / "t.bam"),
            "--normal-bam", str(tmp_path / "n.bam"),
            "--normal-sample", "Hu_X_N",
            "--output", str(tmp_path / "out.vcf.gz"),
        ]
    )
    with pytest.raises(ValueError, match="GATK launcher"):
        scatter.run(args)
