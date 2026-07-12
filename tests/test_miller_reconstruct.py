"""Tests for the Miller Hu_287 reconstruction dependency manifest (pure; injected probes, no tools/network)."""

from __future__ import annotations

import importlib

mr = importlib.import_module("scripts.miller_reconstruct")


def test_dependency_manifest_gates_on_refs_and_carries_install_hints():
    present = {"fasterq-dump": "/x/fasterq-dump", "bwa": "/x/bwa", "samtools": "/x/samtools",
               "bcftools": "/x/bcftools", "salmon": "/x/salmon"}
    # no references present at all -> even tool-satisfied stages are NOT_EVALUABLE (missing reference)
    man = mr.dependency_manifest(which=lambda t: present.get(t), exists=lambda p: False)
    by = {s["stage"]: s for s in man["reconstruction_stages"]}
    assert by["sra_to_fastq"]["status"] == "RUNNABLE"                       # needs no reference
    assert by["rna_quant"]["status"] == "NOT_EVALUABLE"                     # salmon present but index absent
    assert by["rna_quant"]["method_status"][0]["missing_refs"]              # method-level missing reference
    # with the salmon index present, rna_quant becomes RUNNABLE
    man2 = mr.dependency_manifest(which=lambda t: present.get(t), exists=lambda p: "salmon_index" in p)
    by2 = {s["stage"]: s for s in man2["reconstruction_stages"]}
    assert by2["rna_quant"]["status"] == "RUNNABLE"
    # blocked stages carry a concrete install hint (machine-actionable), not a vague blocker
    assert "conda install" in by["hla_typing_classI"]["install_hint"]
    assert "gatk" in by["somatic_calling"]["install_hint"]
    assert "pvactools" in by["mutanome_enumeration"]["install_hint"]
    assert man["isolation"].startswith("LOCKED_TEST")


def test_tool_versions_honors_injected_which():
    # Fix 3: the injected resolver must be used (previously _probe_version called shutil.which directly)
    tv = mr.tool_versions(which=lambda t: None)
    assert tv["fasterq-dump"] is None and tv["samtools"] is None and tv["salmon"] is None


def test_reference_status_present_flag_reflects_disk():
    refs = mr.reference_status(exists=lambda p: p.endswith("gencode.v44.transcripts.fa"))
    by = {r["name"]: r for r in refs}
    assert by["GENCODE_v44_transcripts"]["present"] is True
    assert by["GRCh38_primary_assembly"]["present"] is False
    assert all("acquire" in r and r["acquire"] for r in refs)   # every reference has an acquire command


def test_all_references_have_required_fields():
    for r in mr.REFERENCES:
        assert {"name", "purpose", "dest", "acquire"} <= set(r)


def test_summarize_quant_aggregates_transcript_tpm_to_gene(tmp_path):
    # two transcripts of gene G1, one of G2 -> gene TPM sums; GENCODE header field 5 is the gene symbol
    sf = tmp_path / "quant.sf"
    sf.write_text(
        "Name\tLength\tEffectiveLength\tTPM\tNumReads\n"
        "ENST1|ENSG1|-|-|G1-201|G1|100|x\t100\t80\t60.0\t10\n"
        "ENST2|ENSG1|-|-|G1-202|G1|100|x\t100\t80\t40.0\t8\n"
        "ENST3|ENSG2|-|-|G2-201|G2|100|x\t100\t80\t100.0\t20\n")
    s = mr.summarize_quant(sf)
    assert s["status"] == "OK" and s["n_transcripts"] == 3 and s["n_genes"] == 2
    assert s["top10_genes_tpm"]["G1"] == 100.0 and s["top10_genes_tpm"]["G2"] == 100.0
    assert s["sum_tpm"] == 200.0
