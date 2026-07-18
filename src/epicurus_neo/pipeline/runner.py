"""Pipeline orchestrator: run stages A -> Z with resume, provenance, and fail-stop."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

from epicurus_neo.pipeline.config import PipelineConfig
from epicurus_neo.pipeline.provenance import stage_provenance
from epicurus_neo.pipeline.stage_align import AlignStage
from epicurus_neo.pipeline.stage_annotate import AnnotateStage
from epicurus_neo.pipeline.stage_call import CallStage
from epicurus_neo.pipeline.stage_express import ExpressStage
from epicurus_neo.pipeline.stage_generate import GenerateStage
from epicurus_neo.pipeline.stage_hla import HlaStage
from epicurus_neo.pipeline.stage_prioritize import PrioritizeStage
from epicurus_neo.pipeline.stage_report import ReportStage
from epicurus_neo.pipeline.stages import (
    CACHED,
    COMPLETED,
    FAILED,
    STAGE_ORDER,
    PipelineContext,
    Stage,
    StageResult,
    ToolUnavailableError,
)


def build_stages() -> list[Stage]:
    """Instantiate the eight pipeline stages in fixed execution order."""
    registry: dict[str, Stage] = {
        "align": AlignStage(),
        "call": CallStage(),
        "annotate": AnnotateStage(),
        "express": ExpressStage(),
        "hla": HlaStage(),
        "generate": GenerateStage(),
        "prioritize": PrioritizeStage(),
        "report": ReportStage(),
    }
    return [registry[name] for name in STAGE_ORDER]


def _selected_stages(stages: list[Stage], start: str | None, stop: str | None) -> list[Stage]:
    names = [s.name for s in stages]
    for boundary in (start, stop):
        if boundary is not None and boundary not in names:
            raise ValueError(f"unknown stage '{boundary}'; valid stages: {', '.join(names)}")
    lo = names.index(start) if start else 0
    hi = names.index(stop) if stop else len(names) - 1
    if lo > hi:
        raise ValueError(f"start stage '{start}' comes after stop stage '{stop}'")
    return stages[lo : hi + 1]


@dataclass
class PipelineResult:
    patient_id: str
    output_dir: Path
    results: list[StageResult] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return all(result.ok for result in self.results)

    @property
    def portfolio_path(self) -> Path:
        return self.output_dir / "report" / "portfolio_top20.csv"

    def to_summary(self) -> dict:
        return {
            "patient_id": self.patient_id,
            "ok": self.ok,
            "output_dir": str(self.output_dir),
            "portfolio": str(self.portfolio_path) if self.portfolio_path.exists() else None,
            "stages": [
                {
                    "stage": result.name,
                    "status": result.status,
                    "outputs": result.outputs,
                    "message": result.message,
                }
                for result in self.results
            ],
        }


def run_pipeline(
    config: PipelineConfig,
    *,
    output_dir: str | Path,
    start: str | None = None,
    stop: str | None = None,
    force: bool = False,
) -> PipelineResult:
    """Run the pipeline. Stops at the first failed stage (fail-closed)."""
    ctx = PipelineContext(config=config, output_dir=Path(output_dir))
    ctx.output_dir.mkdir(parents=True, exist_ok=True)
    ctx.provenance_dir.mkdir(parents=True, exist_ok=True)

    stages = _selected_stages(build_stages(), start, stop)
    result = PipelineResult(patient_id=config.patient_id, output_dir=ctx.output_dir)

    for stage in stages:
        declared = {name: str(path) for name, path in stage.declared_outputs(ctx).items()}

        if not force and stage.is_satisfied(ctx):
            result.results.append(StageResult(stage.name, CACHED, outputs=declared))
            continue

        try:
            stage.execute(ctx)
        except ToolUnavailableError as exc:
            result.results.append(StageResult(stage.name, FAILED, message=str(exc)))
            break
        except Exception as exc:  # noqa: BLE001 - any stage failure stops the run, recorded
            result.results.append(
                StageResult(stage.name, FAILED, message=f"{type(exc).__name__}: {exc}")
            )
            break

        provenance: dict = {}
        try:
            commands = stage.build_commands(ctx) if stage.tool is not None else []
            provenance = stage_provenance(
                stage=stage.name,
                tool=stage.tool,
                command=[" ".join(cmd) for cmd in commands],
                inputs=dict(stage.required_inputs(ctx)),
                outputs=stage.declared_outputs(ctx),
            )
            (ctx.provenance_dir / f"{stage.name}.json").write_text(
                json.dumps(provenance, indent=2) + "\n"
            )
        except Exception:  # noqa: BLE001 - provenance must never fail a completed stage
            provenance = {}

        result.results.append(
            StageResult(stage.name, COMPLETED, outputs=declared, provenance=provenance)
        )

    return result
