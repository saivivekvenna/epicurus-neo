"""Miller IPV (PRJNA980652) — T1 raw-input downloader + reconstruction provenance (LOCKED_TEST isolation).

Downloads the PUBLIC SRA Open-Data-Program (.sra) objects for one patient's trio (normal exome + tumor
exome + tumor RNA) **resumably** (HTTP Range) and **checksummed** (sha256 + byte-size verification against
the server Content-Length), deriving the run list deterministically from the public SRA runinfo — it NEVER
reads the recognition-label table. Labels define the pre-registered metric only; they must not influence
download, HLA, expression, calling, candidate generation, thresholds, or ranking.

It also emits a machine-actionable RECONSTRUCTION provenance manifest: for every downstream stage
(sra->fastq, HLA typing, WES alignment, somatic calling, RNA quant, mutanome enumeration, PRIME/Epicurus
scoring) it records the exact tool(s) + reference(s) required and whether they are RUNNABLE or
NOT_EVALUABLE right now — a concrete, actionable reason, never a vague blocker.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import urllib.request
from pathlib import Path

from benchmark.miller_ingest import parse_sra_runinfo

_ROOT = Path(__file__).resolve().parents[2]
ODP_BASE = "https://sra-pub-run-odp.s3.amazonaws.com/sra"
_ROLE = {"normal_exome": "normal_exome", "tumor_exome": "tumor_exome", "tumor_rna": "tumor_rna"}
# local (gitignored) scoring binaries — recognized by PATH-independent existence, not shutil.which
PRIME_BIN = _ROOT / "data/raw/tools/PRIME/PRIME"
MIXMHC_BIN = _ROOT / "data/raw/tools/MixMHCpred/MixMHCpred"


def _default_ref_exists(rel_path: str) -> bool:
    return (_ROOT / rel_path).exists()


def odp_url(run: str) -> str:
    """Public SRA Open-Data-Program object URL for an accession (anonymous S3, HTTP Range-capable)."""
    return f"{ODP_BASE}/{run}/{run}"


def patient_targets(runinfo_path, patient_id: str) -> list[dict]:
    """Deterministic download targets for one patient's trio, derived from the PUBLIC runinfo only.

    Returns [{run, role, url, runinfo_size_mib}] sorted by role. No label table is consulted."""
    runs = parse_sra_runinfo(runinfo_path)
    g = runs[runs["patient_id"] == patient_id]
    if g.empty:
        raise ValueError(f"no SRA runs for patient {patient_id}")
    out = []
    for r in g.itertuples():
        role = _ROLE.get(r.assay_kind)
        if role is None:
            continue
        out.append({"run": r.run, "role": role, "url": odp_url(r.run),
                    "runinfo_size_mib": None if r.size_mb != r.size_mb else float(r.size_mb)})
    return sorted(out, key=lambda d: d["role"])


# ---- primitive, hermetically testable pieces ---------------------------------------------------------
def partial_offset(dest) -> int:
    """Bytes already on disk for a (possibly partial) destination file; 0 if absent."""
    p = Path(dest)
    return p.stat().st_size if p.exists() else 0


def verify_size(dest, expected_bytes: int) -> bool:
    return Path(dest).exists() and Path(dest).stat().st_size == int(expected_bytes)


def sha256_file(path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def probe_size(url: str, timeout: float = 30.0, opener=urllib.request.urlopen) -> int:
    """Server Content-Length via a ranged 1-byte GET (HEAD is sometimes blocked on S3 anon)."""
    req = urllib.request.Request(url, headers={"Range": "bytes=0-0"})
    with opener(req, timeout=timeout) as resp:
        cr = resp.headers.get("Content-Range")            # 'bytes 0-0/TOTAL'
        if cr and "/" in cr:
            return int(cr.rsplit("/", 1)[1])
        cl = resp.headers.get("Content-Length")
        return int(cl) if cl is not None else -1


def resumable_get(url: str, dest, *, expected_size: int | None = None, chunk: int = 1 << 20,
                  timeout: float = 300.0, opener=urllib.request.urlopen) -> dict:
    """Resumable byte-exact download. If a partial file exists, continue from its size via HTTP Range and
    APPEND; if it already equals expected_size, it is a no-op. Verifies the final size when known and
    returns a provenance record {run-agnostic}: bytes, sha256, resumed_from, complete."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    offset = partial_offset(dest)
    if expected_size is not None and offset == int(expected_size):
        return {"bytes": offset, "sha256": sha256_file(dest), "resumed_from": offset, "complete": True,
                "skipped_already_complete": True}
    if expected_size is not None and offset > int(expected_size):
        # corrupt/oversized partial -> restart clean (fail closed, never trust a too-big partial)
        dest.unlink()
        offset = 0
    headers = {"Range": f"bytes={offset}-"} if offset else {}
    req = urllib.request.Request(url, headers=headers)
    mode = "ab" if offset else "wb"
    with opener(req, timeout=timeout) as resp, open(dest, mode) as fh:
        shutil.copyfileobj(resp, fh, length=chunk)
    final = partial_offset(dest)
    complete = (expected_size is None) or (final == int(expected_size))
    return {"bytes": final, "sha256": sha256_file(dest), "resumed_from": offset, "complete": bool(complete),
            "skipped_already_complete": False}


# ---- machine-actionable reconstruction stage map -----------------------------------------------------
def reconstruction_stages(which=shutil.which, ref_exists=_default_ref_exists) -> list[dict]:
    """Every downstream reconstruction stage with its required tool(s) AND reference/index path(s), and its
    CURRENT status: RUNNABLE only if the tools resolve AND every required reference/index is present on disk;
    otherwise NOT_EVALUABLE with a DISTINCT reason (missing tool vs missing reference). The scoring stage
    recognizes the local (gitignored) PRIME/MixMHCpred binaries by path (not PATH) but stays upstream-blocked
    until the candidate universe + HLA exist. Machine-actionable, not a prose 'blocked'."""
    stages = [
        {"stage": "sra_to_fastq", "tools_any": [["fasterq-dump"], ["fastq-dump"]], "ref_paths": [],
         "references": [], "produces": "paired FASTQ per run"},
        {"stage": "hla_typing_classI", "tools_any": [["OptiType"], ["arcasHLA"], ["hla-la"]],
         "ref_paths": ["data/raw/refs/hla/hla_reference_dna.fasta"],
         "references": ["IMGT/HLA DNA reference (hla_reference_dna.fasta)"],
         "produces": "4-digit class-I HLA (A/B/C) from normal exome — REQUIRED for candidate-HLA pairing"},
        {"stage": "wes_alignment", "tools_any": [["bwa", "samtools"], ["bwa-mem2", "samtools"]],
         "ref_paths": ["data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
                       "data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bwt"],
         "references": ["GRCh38 primary assembly FASTA + .fai + BWA index"],
         "produces": "coordinate-sorted tumor & normal BAM"},
        {"stage": "somatic_calling", "tools_any": [["gatk"], ["strelka"]],
         "ref_paths": ["data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"],
         "references": ["GRCh38 FASTA", "germline resource + PoN (Mutect2)"],
         "produces": "somatic SNV/indel VCF with tumor VAF + depth",
         "caveat": "bcftools alone is NOT a valid somatic caller (no tumor-vs-normal model); requires "
                   "Mutect2 (gatk4) or Strelka2"},
        {"stage": "rna_quant", "tools_any": [["salmon"], ["kallisto"]],
         "ref_paths": ["data/raw/refs/gencode/salmon_index/info.json"],
         "references": ["GENCODE transcriptome FASTA + salmon/kallisto index"],
         "produces": "per-gene TPM + mutant-allele RNA evidence"},
        {"stage": "mutanome_enumeration", "tools_any": [["pvacseq"], ["pvactools"]],
         "ref_paths": ["data/raw/refs/vep/homo_sapiens"],
         "references": ["GRCh38 + GENCODE annotation + Ensembl VEP cache"],
         "produces": "full class-I 8-11mer candidate universe (the SHARED denominator for both arms)"},
        {"stage": "scoring_prime_epicurus", "tools_any": [], "ref_paths": [],
         "references": ["genuine PRIME 2.1 + MixMHCpred 3.0 (gitignored local)",
                        "configs/frozen/epicurus_v0_1.json", "configs/frozen/expression_policy_v1.json"],
         "produces": "genuine PRIME AND frozen Epicurus ranks over the IDENTICAL candidate universe"},
    ]
    for s in stages:
        if s["stage"] == "scoring_prime_epicurus":
            local = {"PRIME": str(PRIME_BIN) if PRIME_BIN.exists() else None,
                     "MixMHCpred": str(MIXMHC_BIN) if MIXMHC_BIN.exists() else None}
            s["resolved_tools"] = local
            s["status"] = "NOT_EVALUABLE"
            s["reason"] = ("local PRIME/MixMHCpred present but UPSTREAM-BLOCKED: needs the re-enumerated "
                           "candidate universe + class-I HLA before any peptide can be scored"
                           if (local["PRIME"] or local["MixMHCpred"])
                           else "PRIME/MixMHCpred binaries not found on disk")
            s.pop("tools_any", None)
            continue
        satisfied = None
        for combo in s["tools_any"]:
            if all(which(t) is not None for t in combo):
                satisfied = combo
                break
        missing_refs = [p for p in s["ref_paths"] if not ref_exists(p)]
        if satisfied is None:
            need = " OR ".join("+".join(c) for c in s["tools_any"])
            s["status"], s["reason"] = "NOT_EVALUABLE", f"missing tool(s): need [{need}] on PATH"
        elif missing_refs:
            s["status"] = "NOT_EVALUABLE"
            s["reason"] = f"tools present ({'+'.join(satisfied)}) but missing reference/index: {missing_refs}"
            s["resolved_tools"] = {t: which(t) for t in satisfied}
        else:
            s["status"] = "RUNNABLE"
            s["resolved_tools"] = {t: which(t) for t in satisfied}
        s.pop("tools_any", None)
    return stages


def build_manifest(patient_id: str, targets: list[dict], results: dict, *, which=shutil.which,
                   ref_exists=_default_ref_exists) -> dict:
    """Provenance manifest: per-run download records (incl. expected vs downloaded byte sizes + a verified
    flag) + the machine-actionable downstream stage map."""
    runs = []
    for t in targets:
        r = results.get(t["run"], {})
        exp = r.get("expected_size_bytes")
        got = r.get("bytes")
        runs.append({**{k: t[k] for k in ("run", "role", "url", "runinfo_size_mib")},
                     "expected_size_bytes": exp, "downloaded_bytes": got, "sha256": r.get("sha256"),
                     "size_verified": bool(exp is not None and exp > 0 and got == exp),
                     "complete": r.get("complete"), "resumed_from": r.get("resumed_from")})
    all_complete = bool(runs) and all(r["complete"] and r["size_verified"] for r in runs)
    return {
        "cohort": "miller_ipv", "bioproject": "PRJNA980652", "patient_id": patient_id,
        "isolation": "LOCKED_TEST: labels never consulted for download/HLA/expression/generation/ranking",
        "tranche": "T1", "runs": runs, "download_complete": all_complete,
        "reconstruction_stages": reconstruction_stages(which=which, ref_exists=ref_exists),
        "note": "web search (2026-07-12) found NO released processed Hu_287 VCF/TPM/HLA; raw reconstruction "
                "is necessary. Candidate universe MUST be re-enumerated from WES (not the IPV-prefiltered "
                "tested set); genuine PRIME and frozen Epicurus share that identical universe.",
    }


def download_patient(patient_id: str, dest_dir, runinfo_path, *, chunk: int = 1 << 20) -> dict:
    """Download a patient's full trio (network), verifying each run's byte-size and recording sha256, then
    write DOWNLOAD_MANIFEST.json. Resumable across invocations."""
    dest_dir = Path(dest_dir)
    targets = patient_targets(runinfo_path, patient_id)
    results = {}
    for t in targets:
        dest = dest_dir / f"{t['run']}.sra"
        expected = probe_size(t["url"])
        rec = resumable_get(t["url"], dest, expected_size=expected if expected > 0 else None, chunk=chunk)
        rec["expected_size_bytes"] = expected
        results[t["run"]] = rec
    manifest = build_manifest(patient_id, targets, results)
    (dest_dir / "DOWNLOAD_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest
