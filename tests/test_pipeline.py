from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

from epicurus_neo.pipeline import (
    ConfigError,
    build_stages,
    load_pipeline_config,
    parse_pipeline_config,
    readiness_report,
    run_pipeline,
)
from epicurus_neo.pipeline.provenance import hash_file
from epicurus_neo.pipeline.stages import STAGE_ORDER, PipelineContext
from epicurus_neo.pipeline import tools as tools_mod

DEMO_CANDIDATES = Path("examples/demo_patient/pvacseq_all_epitopes.tsv")


def _config_dict() -> dict:
    return {
        "patient_id": "DEMO-001",
        "inputs": {
            "tumor_wes": ["tumor_R1.fastq.gz", "tumor_R2.fastq.gz"],
            "normal_wes": ["normal_R1.fastq.gz", "normal_R2.fastq.gz"],
            "tumor_rna": ["rna_R1.fastq.gz", "rna_R2.fastq.gz"],
        },
        "references": {"bundle_dir": "/refs/GRCh38"},
        "prioritize": {"k": 20, "max_per_mutation": 1, "max_per_gene": 4},
    }


# ---------------------------------------------------------------- config

def test_parse_config_reads_all_sections():
    config = parse_pipeline_config(_config_dict())
    assert config.patient_id == "DEMO-001"
    assert config.inputs.tumor_wes == ["tumor_R1.fastq.gz", "tumor_R2.fastq.gz"]
    assert config.references_bundle_dir == "/refs/GRCh38"
    assert config.prioritize.k == 20
    assert config.prioritize.max_per_mutation == 1
    assert config.predictors == ["MHCflurry"]


def test_parse_config_missing_required_field_raises():
    data = _config_dict()
    del data["inputs"]["tumor_wes"]
    with pytest.raises(ConfigError, match="tumor_wes"):
        parse_pipeline_config(data)


def test_parse_config_hla_override_parsed():
    data = _config_dict()
    data["inputs"]["hla_alleles"] = ["HLA-A*02:01", "HLA-B*07:02"]
    config = parse_pipeline_config(data)
    assert config.inputs.hla_alleles == ["HLA-A*02:01", "HLA-B*07:02"]


def test_load_config_missing_file_raises(tmp_path: Path):
    with pytest.raises(ConfigError, match="not found"):
        load_pipeline_config(tmp_path / "nope.yaml")


def test_example_config_template_loads():
    config = load_pipeline_config("examples/demo_patient/patient.yaml")
    assert config.patient_id == "DEMO-001"
    assert config.prioritize.max_per_mutation == 1


# ---------------------------------------------------------------- stage graph

def test_stage_order_is_complete_and_fixed():
    assert [s.name for s in build_stages()] == list(STAGE_ORDER)
    assert STAGE_ORDER[-2:] == ("prioritize", "report")


def _ctx(tmp_path: Path, **overrides) -> PipelineContext:
    data = _config_dict()
    data.update(overrides)
    return PipelineContext(config=parse_pipeline_config(data), output_dir=tmp_path)


def test_align_command_uses_bwa_and_reads(tmp_path: Path):
    ctx = _ctx(tmp_path)
    align = build_stages()[0]
    commands = align.build_commands(ctx)
    flat = [tok for cmd in commands for tok in cmd]
    assert commands[0][0] == "bwa-mem2"
    assert "tumor_R1.fastq.gz" in flat and "normal_R2.fastq.gz" in flat
    # produces both tumor and normal marked-dup BAMs
    assert any(str(ctx.paths["tumor_bam"]) in cmd for cmd in commands)
    assert any(str(ctx.paths["normal_bam"]) in cmd for cmd in commands)


def test_call_command_is_mutect2_then_filter(tmp_path: Path):
    ctx = _ctx(tmp_path)
    call = build_stages()[1]
    commands = call.build_commands(ctx)
    assert commands[0][:2] == ["gatk", "Mutect2"]
    assert commands[1][:2] == ["gatk", "FilterMutectCalls"]
    assert any(str(ctx.paths["somatic_vcf"]) in cmd for cmd in commands)


def test_generate_command_includes_predictors_and_lengths(tmp_path: Path):
    ctx = _ctx(tmp_path, inputs={**_config_dict()["inputs"], "hla_alleles": ["HLA-A*02:01"]})
    generate = {s.name: s for s in build_stages()}["generate"]
    command = generate.build_commands(ctx)[0]
    assert command[:2] == ["pvacseq", "run"]
    assert "MHCflurry" in command
    assert "HLA-A*02:01" in command
    assert "8,9,10,11" in command


# ---------------------------------------------------------------- tool availability

def test_missing_tool_fails_stage_and_stops_pipeline(tmp_path: Path, monkeypatch):
    # No external tools available -> align fails immediately and the run stops.
    monkeypatch.setattr(tools_mod.shutil, "which", lambda _binary: None)
    config = parse_pipeline_config(_config_dict())
    result = run_pipeline(config, output_dir=tmp_path)
    assert not result.ok
    assert result.results[0].name == "align"
    assert result.results[0].status == "failed"
    assert "bwa-mem2" in result.results[0].message
    # pipeline stopped at the first failure; no later stages ran
    assert len(result.results) == 1


# ---------------------------------------------------------------- hla override (real execute)

def test_hla_override_writes_alleles_without_tool(tmp_path: Path):
    result = run_pipeline(
        parse_pipeline_config(
            {**_config_dict(), "inputs": {**_config_dict()["inputs"], "hla_alleles": ["HLA-A*02:01", "HLA-B*07:02"]}}
        ),
        output_dir=tmp_path,
        start="hla",
        stop="hla",
    )
    assert result.results[0].status == "completed"
    hla_file = tmp_path / "hla" / "hla_alleles.txt"
    assert hla_file.read_text().splitlines() == ["HLA-A*02:01", "HLA-B*07:02"]


# ---------------------------------------------------------------- prioritize + report (real end to end)

def test_prioritize_and_report_produce_portfolio_from_candidates(tmp_path: Path):
    # Seed the generate-stage output, then run the in-repo tail of the pipeline.
    generate_dir = tmp_path / "generate"
    generate_dir.mkdir(parents=True)
    shutil.copyfile(DEMO_CANDIDATES, generate_dir / "candidates.tsv")

    result = run_pipeline(
        parse_pipeline_config(_config_dict()),
        output_dir=tmp_path,
        start="prioritize",
    )
    assert result.ok, result.to_summary()
    assert [r.name for r in result.results] == ["prioritize", "report"]

    portfolio = tmp_path / "report" / "portfolio_top20.csv"
    assert portfolio.exists()
    import pandas as pd

    frame = pd.read_csv(portfolio)
    assert len(frame) >= 1
    assert len(frame) <= 20
    # one route per mutation by default
    assert frame["mutation_id"].is_unique
    # the unexpressed gene (CTNNB1, TPM 0) is not in the portfolio
    assert "CTNNB1" not in set(frame["gene_symbol"])
    # report artifacts carried forward
    assert (tmp_path / "report" / "report.md").exists()
    assert (tmp_path / "report" / "report.json").exists()


def test_resume_skips_completed_stage(tmp_path: Path):
    generate_dir = tmp_path / "generate"
    generate_dir.mkdir(parents=True)
    shutil.copyfile(DEMO_CANDIDATES, generate_dir / "candidates.tsv")

    first = run_pipeline(parse_pipeline_config(_config_dict()), output_dir=tmp_path, start="prioritize")
    assert all(r.status == "completed" for r in first.results)

    # second run without --force: outputs already exist -> cached
    second = run_pipeline(parse_pipeline_config(_config_dict()), output_dir=tmp_path, start="prioritize")
    assert all(r.status == "cached" for r in second.results)

    # with force: recompute
    forced = run_pipeline(
        parse_pipeline_config(_config_dict()), output_dir=tmp_path, start="prioritize", force=True
    )
    assert all(r.status == "completed" for r in forced.results)


def test_provenance_written_for_completed_inrepo_stage(tmp_path: Path):
    generate_dir = tmp_path / "generate"
    generate_dir.mkdir(parents=True)
    shutil.copyfile(DEMO_CANDIDATES, generate_dir / "candidates.tsv")
    run_pipeline(parse_pipeline_config(_config_dict()), output_dir=tmp_path, start="prioritize")
    prov = json.loads((tmp_path / "provenance" / "prioritize.json").read_text())
    assert prov["stage"] == "prioritize"
    assert prov["tool"] is None  # in-repo stage
    assert "candidates" in prov["inputs"]


# ---------------------------------------------------------------- provenance + doctor

def test_hash_file_absent_and_present(tmp_path: Path):
    assert hash_file(tmp_path / "missing") == "absent"
    f = tmp_path / "a.txt"
    f.write_text("hello")
    assert hash_file(f) == hash_file(f)  # deterministic
    assert len(hash_file(f)) == 64


def test_readiness_report_structure_and_reference_check(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(tools_mod.shutil, "which", lambda _binary: None)
    report = readiness_report(bundle_dir=tmp_path)
    assert report["ready"] is False
    assert report["tools_ready"] is False
    assert len(report["tools"]) == 8
    assert {r["item"] for r in report["references"]} == {"genome.fa", "vep", "salmon_index"}
    # all references absent in an empty bundle dir
    assert all(not r["present"] for r in report["references"])


def test_readiness_report_ready_when_tools_and_refs_present(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(tools_mod.shutil, "which", lambda binary: f"/usr/bin/{binary}")
    (tmp_path / "genome.fa").write_text("")
    (tmp_path / "vep").mkdir()
    (tmp_path / "salmon_index").mkdir()
    report = readiness_report(bundle_dir=tmp_path)
    assert report["tools_ready"] is True
    assert report["references_ready"] is True
    assert report["ready"] is True
