"""Express stage: Salmon RNA quantification -> transcript TPM (quant.sf)."""

from __future__ import annotations

from pathlib import Path

from epicurus_neo.pipeline.stages import PipelineContext, Stage
from epicurus_neo.pipeline.tools import SALMON


class ExpressStage(Stage):
    name = "express"
    tool = SALMON

    def _index(self, ctx: PipelineContext) -> Path:
        return Path(ctx.config.references_bundle_dir).expanduser() / "salmon_index"

    def declared_outputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {"rna_quant": ctx.paths["rna_quant"]}

    def required_inputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {
            "salmon_index": self._index(ctx),
            **{f"tumor_rna_{i}": Path(p) for i, p in enumerate(ctx.config.inputs.tumor_rna)},
        }

    def build_commands(self, ctx: PipelineContext) -> list[list[str]]:
        rna = ctx.config.inputs.tumor_rna
        base = [
            "salmon", "quant",
            "-i", str(self._index(ctx)),
            "-l", "A",
            "-p", str(ctx.config.threads),
            "--validateMappings",
            "-o", str(self.stage_dir(ctx)),
        ]
        if len(rna) >= 2:
            reads = ["-1", rna[0], "-2", rna[1]]
        else:
            reads = ["-r", rna[0]]
        return [[*base, *reads]]
