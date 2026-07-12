"""Tests for the Miller T1 downloader + reconstruction provenance (hermetic: localhost only, no external net).

Covers: ODP URL shape, deterministic trio derivation from the PUBLIC runinfo (no labels), resumable byte-
exact download incl. a real partial->resume against a localhost Range server, sha256 + size verification,
and the machine-actionable reconstruction stage map / manifest.
"""

from __future__ import annotations

import hashlib
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer
from pathlib import Path

import pytest

from benchmark import miller_download as md

FIXTURE = Path(__file__).resolve().parent / "fixtures" / "miller_ipv_sra_runinfo.csv"


# ---- pure derivation (no network, no labels) ---------------------------------------------------------
def test_odp_url_shape():
    assert md.odp_url("SRR24836184") == "https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR24836184/SRR24836184"


def test_patient_targets_hu287_trio_from_public_runinfo():
    t = md.patient_targets(FIXTURE, "Hu_287")
    assert [x["role"] for x in t] == ["normal_exome", "tumor_exome", "tumor_rna"]
    runs = {x["role"]: x["run"] for x in t}
    assert runs == {"normal_exome": "SRR24836184", "tumor_exome": "SRR24836169", "tumor_rna": "SRR24836183"}
    assert all(x["url"] == md.odp_url(x["run"]) for x in t)
    assert all(x["runinfo_size_mib"] and x["runinfo_size_mib"] > 0 for x in t)


def test_patient_targets_unknown_patient_raises():
    with pytest.raises(ValueError):
        md.patient_targets(FIXTURE, "Hu_000")


# ---- checksum / offset primitives ------------------------------------------------------------------
def test_offset_verify_sha(tmp_path):
    p = tmp_path / "f.bin"
    p.write_bytes(b"abcdef")
    assert md.partial_offset(p) == 6
    assert md.partial_offset(tmp_path / "missing") == 0
    assert md.verify_size(p, 6) and not md.verify_size(p, 7)
    assert md.sha256_file(p) == hashlib.sha256(b"abcdef").hexdigest()


# ---- resumable download against a localhost Range-capable server ------------------------------------
class _RangeHandler(BaseHTTPRequestHandler):
    payload = b""

    def log_message(self, *a):
        pass

    def do_GET(self):
        data = type(self).payload
        rng = self.headers.get("Range")
        start = 0
        if rng and rng.startswith("bytes="):
            start = int(rng.split("=")[1].split("-")[0])
        body = data[start:]
        self.send_response(206 if start else 200)
        self.send_header("Content-Length", str(len(body)))
        if start:
            self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
        self.end_headers()
        self.wfile.write(body)


def _serve(payload):
    _RangeHandler.payload = payload
    srv = HTTPServer(("127.0.0.1", 0), _RangeHandler)
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    return srv, f"http://127.0.0.1:{srv.server_address[1]}/obj"


def test_resumable_get_full_then_noop(tmp_path):
    payload = bytes(range(256)) * 400          # 102400 bytes
    srv, url = _serve(payload)
    try:
        dest = tmp_path / "SRRx.sra"
        rec = md.resumable_get(url, dest, expected_size=len(payload), chunk=4096)
        assert rec["complete"] and rec["bytes"] == len(payload) and rec["resumed_from"] == 0
        assert rec["sha256"] == hashlib.sha256(payload).hexdigest()
        assert dest.read_bytes() == payload
        # second call: already complete -> no-op, byte-identical
        rec2 = md.resumable_get(url, dest, expected_size=len(payload), chunk=4096)
        assert rec2["skipped_already_complete"] and rec2["sha256"] == rec["sha256"]
    finally:
        srv.shutdown()


def test_resumable_get_resumes_from_partial(tmp_path):
    payload = bytes(range(256)) * 400
    srv, url = _serve(payload)
    try:
        dest = tmp_path / "SRRy.sra"
        dest.write_bytes(payload[:50000])       # simulate an interrupted download
        rec = md.resumable_get(url, dest, expected_size=len(payload), chunk=4096)
        assert rec["resumed_from"] == 50000 and rec["complete"]
        assert dest.read_bytes() == payload     # appended tail reconstructs the exact object
        assert rec["sha256"] == hashlib.sha256(payload).hexdigest()
    finally:
        srv.shutdown()


def test_resumable_get_oversized_partial_restarts_clean(tmp_path):
    payload = b"X" * 1000
    srv, url = _serve(payload)
    try:
        dest = tmp_path / "SRRz.sra"
        dest.write_bytes(b"Y" * 1500)           # corrupt/too-big partial -> must restart, not trust it
        rec = md.resumable_get(url, dest, expected_size=len(payload), chunk=256)
        assert rec["bytes"] == 1000 and rec["complete"] and dest.read_bytes() == payload
    finally:
        srv.shutdown()


# ---- machine-actionable stage map + manifest -------------------------------------------------------
def test_reconstruction_stages_method_specific_tools_and_reference_sentinels():
    present = {"fasterq-dump": "/x/fd", "bwa": "/x/bwa", "samtools": "/x/st", "salmon": "/x/sal", "gatk": "/x/gatk"}
    def ref(p):                                   # only the salmon index exists; GRCh38/BWA sentinels absent
        return "salmon_index" in p
    by = {s["stage"]: s for s in md.reconstruction_stages(resolve=lambda t: present.get(t), ref_exists=ref)}
    assert by["sra_to_fastq"]["status"] == "RUNNABLE" and by["sra_to_fastq"]["runnable_method"] == "fasterq-dump"
    assert by["rna_quant"]["status"] == "RUNNABLE"                          # salmon + index present
    # bwa+samtools present but the full BWA sentinel set is absent -> NOT_EVALUABLE, method shows missing_refs
    wes = by["wes_alignment"]
    assert wes["status"] == "NOT_EVALUABLE"
    m = wes["method_status"][0]
    assert m["missing_tools"] == [] and m["missing_refs"]                   # distinct: refs missing, tools ok
    # Mutect2: gatk present but FASTA/.fai/.dict absent -> missing_refs (method-specific)
    mut = [x for x in by["somatic_calling"]["method_status"] if x["method"] == "Mutect2"][0]
    assert mut["missing_tools"] == [] and any(".dict" in r for r in mut["missing_refs"])
    # HLA: no OptiType/arcasHLA/T1K -> every method missing tools
    assert by["hla_typing_classI"]["status"] == "NOT_EVALUABLE"
    assert all(x["missing_tools"] for x in by["hla_typing_classI"]["method_status"])
    # mutanome: no vep/pvacseq -> missing tools
    assert by["mutanome_enumeration"]["status"] == "NOT_EVALUABLE"
    assert all(x["missing_tools"] for x in by["mutanome_enumeration"]["method_status"])


def test_scoring_stage_recognizes_local_prime_but_stays_upstream_blocked():
    # Fix 4: scoring must recognize the on-disk PRIME/MixMHCpred (not require PATH), yet stay NOT_EVALUABLE
    by = {s["stage"]: s for s in md.reconstruction_stages(resolve=lambda t: None, ref_exists=lambda p: False)}
    sc = by["scoring_prime_epicurus"]
    assert sc["status"] == "NOT_EVALUABLE"
    assert "PRIME" in sc["resolved_tools"]                                  # local path recognized
    if md.PRIME_BIN.exists():
        assert "UPSTREAM-BLOCKED" in sc["reason"] and sc["resolved_tools"]["PRIME"]


def test_build_manifest_records_provenance_and_isolation():
    targets = md.patient_targets(FIXTURE, "Hu_287")
    results = {t["run"]: {"bytes": 100, "expected_size_bytes": 100, "sha256": "deadbeef",
                          "complete": True, "resumed_from": 0} for t in targets}
    man = md.build_manifest("Hu_287", targets, results, which=lambda t: None, ref_exists=lambda p: False)
    assert man["patient_id"] == "Hu_287" and man["bioproject"] == "PRJNA980652" and man["tranche"] == "T1"
    assert man["download_complete"] is True and len(man["runs"]) == 3
    assert all(r["size_verified"] and r["expected_size_bytes"] == r["downloaded_bytes"] for r in man["runs"])
    assert "labels never consulted" in man["isolation"]
    assert all(s["status"] == "NOT_EVALUABLE" for s in man["reconstruction_stages"])
    assert "no released processed" in man["note"].lower()


def test_build_manifest_size_mismatch_fails_verification(tmp_path):
    # Fix 1: expected_size_bytes persisted; a byte-size mismatch => size_verified False => not complete
    targets = md.patient_targets(FIXTURE, "Hu_287")
    results = {t["run"]: {"bytes": 100, "expected_size_bytes": 100, "complete": True} for t in targets}
    results[targets[0]["run"]] = {"bytes": 99, "expected_size_bytes": 100, "complete": True}  # short read
    man = md.build_manifest("Hu_287", targets, results, which=lambda t: None, ref_exists=lambda p: False)
    bad = [r for r in man["runs"] if r["run"] == targets[0]["run"]][0]
    assert bad["size_verified"] is False and man["download_complete"] is False
