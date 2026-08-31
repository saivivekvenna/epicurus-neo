"""Prioritize stage: the Epicurus Neo differentiator (gate -> rank -> portfolio)."""

from __future__ import annotations

from pathlib import Path

from epicurus_neo.pipeline.stages import PipelineContext, Stage
from epicurus_neo.product import InferenceConfig, run_product_inference


class PrioritizeStage(Stage):
    """Run the deterministic validity gate, calibrated ranking, and portfolio selection.

    This is the only stage Epicurus Neo implements itself; it consumes the generated
    candidate table and writes the ranked candidate list plus a report.
    """

    name = "prioritize"
    tool = None  # pure in-repo Python; no external tool

    def declared_outputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {"ranked_candidates": ctx.paths["ranked_candidates"]}

    def required_inputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {"candidates": ctx.paths["candidates"]}

    def _config(self, ctx: PipelineContext) -> InferenceConfig:
        p = ctx.config.prioritize
        return InferenceConfig(
            k=p.k,
            max_per_mutation=p.max_per_mutation,
            max_per_gene=p.max_per_gene,
            max_per_hla=p.max_per_hla,
        )

    def execute(self, ctx: PipelineContext) -> None:
        candidates = ctx.paths["candidates"]
        if not candidates.exists():
            raise FileNotFoundError(
                f"prioritize stage needs generated candidates at {candidates}; "
                "run the 'generate' stage first"
            )
        self.stage_dir(ctx).mkdir(parents=True, exist_ok=True)
        run_product_inference(
            candidates,
            self.stage_dir(ctx),
            patient_id=ctx.config.patient_id,
            config=self._config(ctx),
        )
