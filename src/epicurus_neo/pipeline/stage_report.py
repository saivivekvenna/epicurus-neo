"""Report stage: emit the final top-20 portfolio and copy the patient report."""

from __future__ import annotations

import shutil
from pathlib import Path

import pandas as pd

from epicurus_neo.pipeline.stages import PipelineContext, Stage


class ReportStage(Stage):
    """Extract the selected portfolio from the ranked candidates and finalize outputs."""

    name = "report"
    tool = None

    def declared_outputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {
            "portfolio": ctx.paths["portfolio"],
            "report_json": ctx.paths["report_json"],
            "report_md": ctx.paths["report_md"],
        }

    def required_inputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {"ranked_candidates": ctx.paths["ranked_candidates"]}

    def execute(self, ctx: PipelineContext) -> None:
        ranked_path = ctx.paths["ranked_candidates"]
        if not ranked_path.exists():
            raise FileNotFoundError(
                f"report stage needs ranked candidates at {ranked_path}; "
                "run the 'prioritize' stage first"
            )
        self.stage_dir(ctx).mkdir(parents=True, exist_ok=True)

        ranked = pd.read_csv(ranked_path)
        selected = ranked
        if "selected" in ranked.columns:
            selected = ranked[ranked["selected"].astype(bool)]
        if "rank" in selected.columns:
            selected = selected.sort_values("rank")
        selected.to_csv(ctx.paths["portfolio"], index=False)

        # Carry the human/machine report emitted by the prioritize stage forward.
        prioritize_dir = ctx.stage_dir("prioritize")
        for src_name, dest in (
            ("report.json", ctx.paths["report_json"]),
            ("report.md", ctx.paths["report_md"]),
        ):
            src = prioritize_dir / src_name
            if src.exists():
                shutil.copyfile(src, dest)
