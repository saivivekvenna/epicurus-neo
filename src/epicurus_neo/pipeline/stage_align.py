"""Align stage: BWA-MEM2 alignment + duplicate marking -> tumor/normal BAM."""

from __future__ import annotations

from pathlib import Path

from epicurus_neo.pipeline.stages import PipelineContext, Stage
from epicurus_neo.pipeline.tools import BWA_MEM2


class AlignStage(Stage):
    name = "align"
    tool = BWA_MEM2  # primary tool checked for readiness (also uses samtools + gatk)

    def _genome(self, ctx: PipelineContext) -> Path:
        return Path(ctx.config.references_bundle_dir).expanduser() / "genome.fa"

    def declared_outputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {
            "tumor_bam": ctx.paths["tumor_bam"],
            "normal_bam": ctx.paths["normal_bam"],
        }

    def required_inputs(self, ctx: PipelineContext) -> dict[str, Path]:
        reads = ctx.config.inputs
        return {
            "genome": self._genome(ctx),
            **{f"tumor_wes_{i}": Path(p) for i, p in enumerate(reads.tumor_wes)},
            **{f"normal_wes_{i}": Path(p) for i, p in enumerate(reads.normal_wes)},
        }

    def _sample_commands(
        self, ctx: PipelineContext, sample: str, fastqs: list[str], out_bam: Path
    ) -> list[list[str]]:
        genome = str(self._genome(ctx))
        threads = str(ctx.config.threads)
        stage = self.stage_dir(ctx)
        sam = str(stage / f"{sample}.sam")
        sorted_bam = str(stage / f"{sample}.sorted.bam")
        metrics = str(stage / f"{sample}.dup_metrics.txt")
        return [
            ["bwa-mem2", "mem", "-t", threads, "-o", sam, genome, *fastqs],
            ["samtools", "sort", "-@", threads, "-o", sorted_bam, sam],
            ["gatk", "MarkDuplicates", "-I", sorted_bam, "-O", str(out_bam), "-M", metrics],
            ["samtools", "index", str(out_bam)],
        ]

    def build_commands(self, ctx: PipelineContext) -> list[list[str]]:
        reads = ctx.config.inputs
        return [
            *self._sample_commands(ctx, "tumor", reads.tumor_wes, ctx.paths["tumor_bam"]),
            *self._sample_commands(ctx, "normal", reads.normal_wes, ctx.paths["normal_bam"]),
        ]
