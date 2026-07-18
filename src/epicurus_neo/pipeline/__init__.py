"""Epicurus full pipeline: raw WES/RNA -> ranked <=20 neoantigen portfolio.

The pipeline orchestrates established external tools (BWA-MEM2, GATK, VEP, Salmon,
OptiType, pVACtools, MHCflurry) for read alignment through candidate generation,
and applies the in-repo Epicurus prioritize stage (validity gate -> calibrated
ranking -> portfolio) as the final step.
"""

from __future__ import annotations

from epicurus_neo.pipeline.config import (
    ConfigError,
    PipelineConfig,
    load_pipeline_config,
    parse_pipeline_config,
)
from epicurus_neo.pipeline.doctor import readiness_report
from epicurus_neo.pipeline.references import references_manifest, scaffold_references
from epicurus_neo.pipeline.runner import PipelineResult, build_stages, run_pipeline
from epicurus_neo.pipeline.stages import STAGE_ORDER, PipelineContext, StageResult

__all__ = [
    "ConfigError",
    "PipelineConfig",
    "PipelineContext",
    "PipelineResult",
    "STAGE_ORDER",
    "StageResult",
    "build_stages",
    "load_pipeline_config",
    "parse_pipeline_config",
    "readiness_report",
    "references_manifest",
    "run_pipeline",
    "scaffold_references",
]
