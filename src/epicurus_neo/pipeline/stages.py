"""Stage abstraction, execution context, and result types for the pipeline."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from epicurus_neo.pipeline.config import PipelineConfig
from epicurus_neo.pipeline.tools import ToolSpec

# Fixed left-to-right execution order (A -> Z).
STAGE_ORDER: tuple[str, ...] = (
    "align",
    "call",
    "annotate",
    "express",
    "hla",
    "generate",
    "prioritize",
    "report",
)

COMPLETED = "completed"
CACHED = "cached"
SKIPPED = "skipped"
FAILED = "failed"


class ToolUnavailableError(RuntimeError):
    """Raised when a stage's required external tool is not installed."""


@dataclass
class PipelineContext:
    """Everything a stage needs: config + where artifacts live."""

    config: PipelineConfig
    output_dir: Path

    def stage_dir(self, name: str) -> Path:
        return self.output_dir / name

    @property
    def provenance_dir(self) -> Path:
        return self.output_dir / "provenance"

    def artifact(self, stage: str, filename: str) -> Path:
        return self.stage_dir(stage) / filename

    # Canonical artifact locations, shared so downstream stages agree with upstream.
    @property
    def paths(self) -> dict[str, Path]:
        return {
            "tumor_bam": self.artifact("align", "tumor.md.bam"),
            "normal_bam": self.artifact("align", "normal.md.bam"),
            "somatic_vcf": self.artifact("call", "somatic.filtered.vcf.gz"),
            "annotated_vcf": self.artifact("annotate", "annotated.vcf.gz"),
            "rna_quant": self.artifact("express", "quant.sf"),
            "hla_alleles": self.artifact("hla", "hla_alleles.txt"),
            "candidates": self.artifact("generate", "candidates.tsv"),
            "ranked_candidates": self.artifact("prioritize", "ranked_candidates.csv"),
            "portfolio": self.artifact("report", "portfolio_top20.csv"),
            "report_json": self.artifact("report", "report.json"),
            "report_md": self.artifact("report", "report.md"),
        }


@dataclass(frozen=True)
class StageResult:
    name: str
    status: str
    outputs: dict[str, str] = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)
    message: str = ""

    @property
    def ok(self) -> bool:
        return self.status in {COMPLETED, CACHED, SKIPPED}


class Stage:
    """Base class for a pipeline stage.

    External-tool stages implement :meth:`build_command` and inherit the default
    :meth:`execute`, which verifies the tool then shells out. Pure in-repo stages
    (prioritize, report) set ``tool = None`` and override :meth:`execute`.
    """

    name: str = ""
    tool: ToolSpec | None = None

    def declared_outputs(self, ctx: PipelineContext) -> dict[str, Path]:
        raise NotImplementedError

    def required_inputs(self, ctx: PipelineContext) -> dict[str, Path]:
        return {}

    def build_commands(self, ctx: PipelineContext) -> list[list[str]]:
        """The ordered external commands (argv lists) this stage runs.

        External-tool stages implement this; pure in-repo stages leave it empty.
        """
        raise NotImplementedError(f"stage '{self.name}' has no external command")

    def is_satisfied(self, ctx: PipelineContext) -> bool:
        """True if every declared output already exists (resume support)."""
        outputs = self.declared_outputs(ctx)
        return bool(outputs) and all(path.exists() for path in outputs.values())

    def check_tool(self) -> None:
        if self.tool is not None and not self.tool.is_available():
            raise ToolUnavailableError(
                f"stage '{self.name}' requires '{self.tool.binary}' "
                f"({self.tool.name}), which is not installed or not on PATH"
            )

    def execute(self, ctx: PipelineContext) -> None:
        """Default execution for external-tool stages: check tool, run each command."""
        self.check_tool()
        self.stage_dir(ctx).mkdir(parents=True, exist_ok=True)
        for command in self.build_commands(ctx):
            subprocess.run(command, check=True)  # noqa: S603 - fixed argv, no shell

    def stage_dir(self, ctx: PipelineContext) -> Path:
        return ctx.stage_dir(self.name)
