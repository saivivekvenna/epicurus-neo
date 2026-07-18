"""Annotate stage: Ensembl VEP annotation of the somatic VCF."""

from __future__ import annotations

from pathlib import Path

from epicurus_neo.pipeline.stages import PipelineContext, Stage
from epicurus_neo.pipeline.tools import VEP


class AnnotateStage(Stage):
    name = "annotate"
    tool = VEP

    def _bundle(self, ctx: PipelineContext) -> Path:
        return Path(ctx.config.references_bundle_dir).expanduser()

    def declared_outputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {"annotated_vcf": ctx.paths["annotated_vcf"]}

    def required_inputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {"somatic_vcf": ctx.paths["somatic_vcf"]}

    def build_commands(self, ctx: PipelineContext) -> list[list[str]]:
        bundle = self._bundle(ctx)
        return [
            [
                "vep",
                "--input_file", str(ctx.paths["somatic_vcf"]),
                "--output_file", str(ctx.paths["annotated_vcf"]),
                "--vcf",
                "--compress_output", "bgzip",
                "--offline",
                "--cache",
                "--dir_cache", str(bundle / "vep"),
                "--fasta", str(bundle / "genome.fa"),
                "--force_overwrite",
                "--symbol",
                "--transcript_version",
            ]
        ]
