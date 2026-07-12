"""Tests for the Miller Hu_287 reconstruction dependency manifest (pure; injected probes, no tools/network)."""

from __future__ import annotations

import importlib

mr = importlib.import_module("scripts.miller_reconstruct")


def test_dependency_manifest_marks_runnable_and_blocked_with_hints():
    present = {"fasterq-dump": "/x/fasterq-dump", "bwa": "/x/bwa", "samtools": "/x/samtools",
               "bcftools": "/x/bcftools", "salmon": "/x/salmon"}
    man = mr.dependency_manifest(which=lambda t: present.get(t), exists=lambda p: False)
    by = {s["stage"]: s for s in man["reconstruction_stages"]}
    assert by["sra_to_fastq"]["status"] == "RUNNABLE"
    assert by["rna_quant"]["status"] == "RUNNABLE"
    # blocked stages carry a concrete install hint (machine-actionable), not a vague blocker
    assert by["hla_typing_classI"]["status"] == "NOT_EVALUABLE"
    assert "conda install" in by["hla_typing_classI"]["install_hint"]
    assert "gatk" in by["somatic_calling"]["install_hint"]
    assert "pvactools" in by["mutanome_enumeration"]["install_hint"]
    assert man["isolation"].startswith("LOCKED_TEST")


def test_reference_status_present_flag_reflects_disk():
    refs = mr.reference_status(exists=lambda p: p.endswith("gencode.v44.transcripts.fa"))
    by = {r["name"]: r for r in refs}
    assert by["GENCODE_v44_transcripts"]["present"] is True
    assert by["GRCh38_primary_assembly"]["present"] is False
    assert all("acquire" in r and r["acquire"] for r in refs)   # every reference has an acquire command


def test_all_references_have_required_fields():
    for r in mr.REFERENCES:
        assert {"name", "purpose", "dest", "acquire"} <= set(r)
