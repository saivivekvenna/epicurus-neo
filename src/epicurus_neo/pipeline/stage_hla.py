"""HLA stage: class-I typing with OptiType, or a config-provided allele override."""

from __future__ import annotations

from pathlib import Path

from epicurus_neo.pipeline.stages import PipelineContext, Stage
from epicurus_neo.pipeline.tools import OPTITYPE


class HlaStage(Stage):
    name = "hla"
    tool = OPTITYPE

    def declared_outputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {"hla_alleles": ctx.paths["hla_alleles"]}

    def required_inputs(self, ctx: PipelineContext) -> dict[str, Path]:
        if ctx.config.inputs.hla_alleles:
            return {}
        return {f"normal_wes_{i}": Path(p) for i, p in enumerate(ctx.config.inputs.normal_wes)}

    def build_commands(self, ctx: PipelineContext) -> list[list[str]]:
        # Only used when no allele override is supplied.
        return [
            [
                "OptiTypePipeline.py",
                "--input", *ctx.config.inputs.normal_wes,
                "--dna",
                "--outdir", str(self.stage_dir(ctx)),
            ]
        ]

    def execute(self, ctx: PipelineContext) -> None:
        override = ctx.config.inputs.hla_alleles
        self.stage_dir(ctx).mkdir(parents=True, exist_ok=True)
        if override:
            # No typing needed: record the caller-provided alleles directly.
            ctx.paths["hla_alleles"].write_text("\n".join(override) + "\n")
            return
        # Otherwise run OptiType, then the caller-supplied post-parse writes
        # hla_alleles.txt from the OptiType result on the Linux acceptance run.
        super().execute(ctx)
