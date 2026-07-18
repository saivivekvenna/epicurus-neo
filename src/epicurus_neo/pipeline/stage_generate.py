"""Generate stage: pVACseq peptide generation + MHC presentation -> candidate table."""

from __future__ import annotations

import shutil
from pathlib import Path

from epicurus_neo.pipeline.stages import PipelineContext, Stage
from epicurus_neo.pipeline.tools import PVACSEQ


class GenerateStage(Stage):
    name = "generate"
    tool = PVACSEQ

    def declared_outputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {"candidates": ctx.paths["candidates"]}

    def required_inputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {
            "annotated_vcf": ctx.paths["annotated_vcf"],
            "hla_alleles": ctx.paths["hla_alleles"],
        }

    def _alleles(self, ctx: PipelineContext) -> str:
        if ctx.config.inputs.hla_alleles:
            alleles = ctx.config.inputs.hla_alleles
        else:
            hla_file = ctx.paths["hla_alleles"]
            alleles = (
                [line.strip() for line in hla_file.read_text().splitlines() if line.strip()]
                if hla_file.exists()
                else []
            )
        return ",".join(alleles)

    def build_commands(self, ctx: PipelineContext) -> list[list[str]]:
        lengths = ",".join(str(x) for x in ctx.config.epitope_lengths)
        command = [
            "pvacseq", "run",
            str(ctx.paths["annotated_vcf"]),
            ctx.config.patient_id,
            self._alleles(ctx),
            *ctx.config.predictors,
            str(self.stage_dir(ctx)),
            "-e1", lengths,
            "-t", str(ctx.config.threads),
        ]
        return [command]

    def execute(self, ctx: PipelineContext) -> None:
        super().execute(ctx)
        # pVACseq writes MHC_Class_I/<sample>.all_epitopes.tsv; expose it as the
        # canonical candidate artifact the prioritize stage consumes.
        produced = (
            self.stage_dir(ctx)
            / "MHC_Class_I"
            / f"{ctx.config.patient_id}.all_epitopes.tsv"
        )
        if produced.exists():
            shutil.copyfile(produced, ctx.paths["candidates"])
