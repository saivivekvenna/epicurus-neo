"""Patient pipeline configuration: parse and validate ``patient.yaml``."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when a patient config is missing or malformed."""


@dataclass(frozen=True)
class PipelineInputs:
    tumor_wes: list[str]
    normal_wes: list[str]
    tumor_rna: list[str]
    hla_alleles: list[str] | None = None


@dataclass(frozen=True)
class PrioritizeConfig:
    k: int = 20
    max_per_mutation: int = 1
    max_per_gene: int = 4
    max_per_hla: int | None = None


@dataclass(frozen=True)
class PipelineConfig:
    patient_id: str
    inputs: PipelineInputs
    references_bundle_dir: str
    predictors: list[str] = field(default_factory=lambda: ["MHCflurry"])
    epitope_lengths: list[int] = field(default_factory=lambda: [8, 9, 10, 11])
    prioritize: PrioritizeConfig = field(default_factory=PrioritizeConfig)
    threads: int = 8


def _require(mapping: dict[str, Any], key: str, context: str) -> Any:
    if key not in mapping or mapping[key] is None:
        raise ConfigError(f"config is missing required field '{context}{key}'")
    return mapping[key]


def _as_str_list(value: Any, context: str) -> list[str]:
    if isinstance(value, str):
        return [value]
    if isinstance(value, (list, tuple)) and all(isinstance(v, str) for v in value):
        return [str(v) for v in value]
    raise ConfigError(f"config field '{context}' must be a string or list of strings")


def parse_pipeline_config(data: dict[str, Any]) -> PipelineConfig:
    """Build a :class:`PipelineConfig` from a parsed mapping, validating required fields."""
    if not isinstance(data, dict):
        raise ConfigError("config root must be a mapping")

    patient_id = str(_require(data, "patient_id", ""))

    raw_inputs = _require(data, "inputs", "")
    if not isinstance(raw_inputs, dict):
        raise ConfigError("config field 'inputs' must be a mapping")
    inputs = PipelineInputs(
        tumor_wes=_as_str_list(_require(raw_inputs, "tumor_wes", "inputs."), "inputs.tumor_wes"),
        normal_wes=_as_str_list(_require(raw_inputs, "normal_wes", "inputs."), "inputs.normal_wes"),
        tumor_rna=_as_str_list(_require(raw_inputs, "tumor_rna", "inputs."), "inputs.tumor_rna"),
        hla_alleles=(
            _as_str_list(raw_inputs["hla_alleles"], "inputs.hla_alleles")
            if raw_inputs.get("hla_alleles")
            else None
        ),
    )

    references = _require(data, "references", "")
    if not isinstance(references, dict):
        raise ConfigError("config field 'references' must be a mapping")
    bundle_dir = str(_require(references, "bundle_dir", "references."))

    generate = data.get("generate") or {}
    predictors = (
        _as_str_list(generate["predictors"], "generate.predictors")
        if generate.get("predictors")
        else ["MHCflurry"]
    )
    epitope_lengths = generate.get("epitope_lengths") or [8, 9, 10, 11]
    if not all(isinstance(x, int) for x in epitope_lengths):
        raise ConfigError("config field 'generate.epitope_lengths' must be a list of integers")

    raw_prioritize = data.get("prioritize") or {}
    prioritize = PrioritizeConfig(
        k=int(raw_prioritize.get("k", 20)),
        max_per_mutation=int(raw_prioritize.get("max_per_mutation", 1)),
        max_per_gene=int(raw_prioritize.get("max_per_gene", 4)),
        max_per_hla=(
            int(raw_prioritize["max_per_hla"])
            if raw_prioritize.get("max_per_hla") is not None
            else None
        ),
    )

    resources = data.get("resources") or {}
    threads = int(resources.get("threads", 8))

    return PipelineConfig(
        patient_id=patient_id,
        inputs=inputs,
        references_bundle_dir=bundle_dir,
        predictors=predictors,
        epitope_lengths=list(epitope_lengths),
        prioritize=prioritize,
        threads=threads,
    )


def load_pipeline_config(path: str | Path) -> PipelineConfig:
    """Load and validate a patient pipeline config from a YAML file."""
    config_path = Path(path)
    if not config_path.exists():
        raise ConfigError(f"config file not found: {config_path}")
    try:
        data = yaml.safe_load(config_path.read_text())
    except yaml.YAMLError as exc:  # pragma: no cover - passthrough of parser detail
        raise ConfigError(f"config file is not valid YAML: {exc}") from exc
    return parse_pipeline_config(data)
