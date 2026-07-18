"""Call stage: GATK Mutect2 + FilterMutectCalls -> filtered somatic VCF.

The exact resource-bundle arguments are finalized against a real GRCh38 bundle on
the Linux acceptance run; the command shape here is the frozen contract.
"""

from __future__ import annotations

from pathlib import Path

from epicurus_neo.pipeline.stages import PipelineContext, Stage
from epicurus_neo.pipeline.tools import GATK


class CallStage(Stage):
    name = "call"
    tool = GATK

    def _genome(self, ctx: PipelineContext) -> Path:
        return Path(ctx.config.references_bundle_dir).expanduser() / "genome.fa"

    def declared_outputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {"somatic_vcf": ctx.paths["somatic_vcf"]}

    def required_inputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {
            "genome": self._genome(ctx),
            "tumor_bam": ctx.paths["tumor_bam"],
            "normal_bam": ctx.paths["normal_bam"],
        }

    def build_commands(self, ctx: PipelineContext) -> list[list[str]]:
        genome = str(self._genome(ctx))
        unfiltered = str(self.stage_dir(ctx) / "unfiltered.vcf.gz")
        return [
            [
                "gatk", "Mutect2",
                "-R", genome,
                "-I", str(ctx.paths["tumor_bam"]),
                "-I", str(ctx.paths["normal_bam"]),
                "-normal", "normal",
                "-O", unfiltered,
            ],
            [
                "gatk", "FilterMutectCalls",
                "-R", genome,
                "-V", unfiltered,
                "-O", str(ctx.paths["somatic_vcf"]),
            ],
        ]
