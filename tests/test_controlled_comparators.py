from __future__ import annotations

import csv
from pathlib import Path

import pytest

from benchmark.controlled_comparators import freeze_vaxrank_portfolio, prepare_common_bundle


def _touch(path: Path, data: bytes = b"x") -> Path:
    path.write_bytes(data)
    return path


def test_common_bundle_hashes_inputs_and_keeps_hla_tool_native(tmp_path, monkeypatch):
    monkeypatch.chdir(Path(__file__).parents[1])
    hla = _touch(tmp_path / "hla.txt", b"HLA-A*01:01\nHLA-B*07:02\n")
    files = [_touch(tmp_path / name) for name in ("pass.vcf.gz", "pass.vcf.gz.tbi", "rna.bam",
                                                   "rna.bam.bai", "pvac.vcf.gz", "pvac.vcf.gz.tbi")]
    manifest = prepare_common_bundle(
        patient_id="Hu_test", pass_vcf=files[0], pass_vcf_index=files[1], rna_bam=files[2],
        rna_bam_index=files[3], hla_panel=hla, pvac_ready_vcf=files[4],
        pvac_ready_vcf_index=files[5], output_dir=tmp_path / "out",
    )
    assert manifest["labels_opened"] is False
    assert manifest["allow_dna_only_fallback"] is False
    assert "--mhc-predictor mhcflurry" in manifest["commands"]["vaxrank"]
    assert "--normal-sample-name Hu_test_N" in manifest["commands"]["pvacseq"]
    assert "--pVACseq_filter_set" not in manifest["commands"]["pvacseq"]
    assert set(manifest["inputs"]) == {
        "pass_vcf", "pass_vcf_index", "rna_bam", "rna_bam_index", "hla_panel",
        "pvac_ready_vcf", "pvac_ready_vcf_index",
    }


def test_vaxrank_freeze_uses_first_peptide_per_native_variant_rank(tmp_path):
    source = tmp_path / "ranked.csv"
    fields = ("chr", "pos", "ref", "alt", "variant_rank", "peptide_rank", "amino_acids",
              "combined_score", "expression_score", "target_epitope_score")
    with source.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for rank in range(1, 24):
            for peptide_rank in (2, 1):
                writer.writerow({"chr": "1", "pos": rank, "ref": "A", "alt": "T",
                                 "variant_rank": rank, "peptide_rank": peptide_rank,
                                 "amino_acids": f"PEPTIDE{peptide_rank}", "combined_score": 1,
                                 "expression_score": 1, "target_epitope_score": 1})
    result = freeze_vaxrank_portfolio(source, tmp_path / "out")
    assert result["portfolio"]["size"] == 20
    with (tmp_path / "out" / "vaxrank_native_top20.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert [int(row["rank"]) for row in rows] == list(range(1, 21))
    assert all(row["peptide"] == "PEPTIDE1" for row in rows)


def test_vaxrank_freeze_rejects_noncontiguous_native_ranks(tmp_path):
    source = tmp_path / "ranked.csv"
    source.write_text(
        "chr,pos,ref,alt,variant_rank,peptide_rank,amino_acids\n"
        "1,1,A,T,1,1,AAAAAAAA\n1,2,A,T,3,1,BBBBBBBB\n"
    )
    with pytest.raises(ValueError, match="contiguous"):
        freeze_vaxrank_portfolio(source, tmp_path / "out")
