"""Focused schema/integrity validation for the Miller Hu_287 HLA provenance JSON.

Validates the provenance chain recorded by scripts/miller_hu287_hla.sh (LOCKED_TEST; no label read):
per-file relative path/sha256/size for the normal MD BAM+index, extracted FASTQs, and OptiType result,
plus environment identity and the OptiType objective/read-count/region. Skips cleanly when the artifact
is absent (e.g. a fresh checkout without the local reconstruction data)."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
PROV = (ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv"
        / "hu_287_reconstruction/HLA_PROVENANCE.json")

_ALLELE = re.compile(r"^HLA-[ABC]\*\d{2,3}:\d{2,3}$")
_HEX64 = re.compile(r"^[0-9a-f]{64}$")
_REQUIRED_FILES = ("normal_md_bam", "normal_md_bam_index", "hla_fastq_1", "hla_fastq_2", "optitype_result_tsv")


def _doc() -> dict:
    if not PROV.is_file():
        pytest.skip("HLA_PROVENANCE.json not present (no local reconstruction artifact)")
    return json.loads(PROV.read_text())


def test_hla_provenance_core_fields():
    d = _doc()
    assert d["patient_id"] == "Hu_287"
    assert "no label read" in d["isolation"]                  # LOCKED_TEST flag intact
    assert d["extraction_region_grch38"] == "6:29800000-33600000"
    assert isinstance(d["objective"], float) and d["objective"] > 0
    assert isinstance(d["reads_used"], float) and d["reads_used"] > 0
    assert isinstance(d["result_path"], str) and d["result_path"].endswith("_result.tsv")


def test_hla_provenance_six_class_i_alleles_two_per_locus():
    alleles = _doc()["class_i_alleles"]
    assert len(alleles) == 6 and alleles == sorted(alleles)   # deterministic sorted set
    assert all(_ALLELE.match(a) for a in alleles), alleles
    for locus in ("A", "B", "C"):
        assert sum(a.startswith(f"HLA-{locus}*") for a in alleles) == 2, locus


def test_hla_provenance_environment_identity():
    env = _doc()["environment"]
    for k in ("micromamba_version", "optitype_pkg", "glpk_pkg", "razers3_pkg", "glpk_solver",
              "mamba_env", "mamba_root_prefix"):
        assert env.get(k), k                                  # non-empty identity fields
    assert "GLPK" in env["glpk_solver"]                       # solver version string recorded


def test_hla_provenance_file_records_are_wellformed():
    files = _doc()["provenance_files"]
    assert set(files) == set(_REQUIRED_FILES)
    for name, rec in files.items():
        assert not Path(rec["path"]).is_absolute(), name      # relative repo path only
        assert _HEX64.match(rec["sha256"]), name              # 64-hex sha256
        assert isinstance(rec["size_bytes"], int) and rec["size_bytes"] > 0, name


def test_hla_provenance_result_tsv_hash_and_size_match_on_disk():
    # spot-check the small result TSV: recorded sha256 + size must match the actual file bytes
    rec = _doc()["provenance_files"]["optitype_result_tsv"]
    p = ROOT / rec["path"]
    if not p.is_file():
        pytest.skip("result TSV not present locally")
    assert p.stat().st_size == rec["size_bytes"]
    assert hashlib.sha256(p.read_bytes()).hexdigest() == rec["sha256"]
