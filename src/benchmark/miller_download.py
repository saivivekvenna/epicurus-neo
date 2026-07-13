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
# tools we may install locally (not on PATH) — resolved by existence in addition to shutil.which
LOCAL_TOOLS = {"gatk": _ROOT / "data/raw/tools/gatk-4.5.0.0/gatk"}
# pinned, reproducible micromamba environments — the EXACT envs that reconstructed Hu_287. Tools inside
# them are not on PATH; they are invoked as `micromamba run -n <env> <tool>`. We resolve them by the
# existence of the env's bin/<tool>, so a Hu_287-identical toolchain is recognized without ad-hoc installs.
_MICROMAMBA = _ROOT / "data/raw/tools/bin/micromamba"
_MAMBA_ENVS = _ROOT / "data/raw/tools/micromamba/envs"
PINNED_ENV_TOOLS = {"OptiTypePipeline.py": "hla", "razers3": "hla", "vep": "vep"}

# reference sentinel paths (repo-relative)
_G = "data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
_G_DICT = "data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.dict"
_BWA_SENTINELS = [_G, _G + ".fai", _G + ".amb", _G + ".ann", _G + ".bwt", _G + ".pac", _G + ".sa"]


def _default_ref_exists(rel_path: str) -> bool:
    return (_ROOT / rel_path).exists()


def pinned_env_tool(name: str) -> str | None:
    """Resolve a tool provided by a pinned micromamba env (env bin/<tool>), else None.

    Fail-closed: returns None unless BOTH the pinned micromamba launcher and the env's bin/<tool> exist,
    so a partially-provisioned toolchain is never reported as runnable."""
    env = PINNED_ENV_TOOLS.get(name)
    if env is None or not _MICROMAMBA.exists():
        return None
    cand = _MAMBA_ENVS / env / "bin" / name
    return str(cand) if cand.exists() else None


def resolve_tool(name: str, which=shutil.which) -> str | None:
    """Resolve a tool by PATH first, then a known local install path, then a pinned micromamba env."""
    hit = which(name)
    if hit:
        return hit
    local = LOCAL_TOOLS.get(name)
    if local and local.exists():
        return str(local)
    return pinned_env_tool(name)


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
def reconstruction_stages(which=shutil.which, ref_exists=_default_ref_exists, resolve=None) -> list[dict]:
    """Every downstream stage as one or more METHODS, each with its specific tool(s) AND reference/index
    SENTINELS. A method is RUNNABLE only if all its tools resolve (PATH or local) AND every sentinel exists;
    a stage is RUNNABLE if ANY method is, else NOT_EVALUABLE listing per-method missing tools/references
    distinctly. Method-specific (Fixes): WES needs the full BWA sentinel set (.amb/.ann/.bwt/.pac/.sa + .fai);
    Mutect2 needs FASTA+.fai+.dict (germline resource is an optional documented strategy) while Strelka2 has
    its own set; mutanome needs GRCh38 + GENCODE GTF + VEP cache, not merely a cache dir; HLA is
    tool-specific. Scoring recognizes local PRIME/MixMHCpred but stays upstream-blocked."""
    stages = [
        {"stage": "sra_to_fastq", "produces": "paired FASTQ per run",
         "methods": [{"method": "fasterq-dump", "tools": ["fasterq-dump"], "refs": []},
                     {"method": "fastq-dump", "tools": ["fastq-dump"], "refs": []}]},
        {"stage": "hla_typing_classI",
         "produces": "4-digit class-I HLA (A/B/C) from normal exome (needed for PRIME/EL ranking)",
         "methods": [{"method": "OptiType", "tools": ["OptiTypePipeline.py", "razers3"], "refs": [],
                      "note": "pinned micromamba `hla` env (OptiType+razers3+GLPK); OptiType bundles its own "
                              "HLA reference so no external FASTA sentinel is required — the exact provider "
                              "that typed Hu_287 via scripts/miller_hu287_hla.sh"},
                     {"method": "arcasHLA", "tools": ["arcasHLA"], "refs": ["data/raw/refs/hla/IMGTHLA"]},
                     {"method": "T1K", "tools": ["run-t1k"], "refs": ["data/raw/refs/hla/t1k_hlaidx"]}]},
        {"stage": "wes_alignment", "produces": "coordinate-sorted, dup-marked tumor & normal BAM",
         "methods": [{"method": "bwa-mem+samtools", "tools": ["bwa", "samtools"], "refs": _BWA_SENTINELS}]},
        {"stage": "somatic_calling", "produces": "PASS somatic SNV/indel VCF with tumor VAF + depth",
         "methods": [{"method": "Mutect2", "tools": ["gatk"], "refs": [_G, _G + ".fai", _G_DICT],
                      "note": "matched tumor-vs-normal; gnomAD af-only germline resource optional "
                              "(documented sensitivity deviation if absent)"},
                     {"method": "Strelka2", "tools": ["configureStrelkaSomaticWorkflow.py"],
                      "refs": [_G, _G + ".fai"], "note": "requires bgzipped+tabixed inputs"}]},
        {"stage": "rna_quant", "produces": "per-gene TPM (+ mutant-allele RNA evidence via genome align)",
         "methods": [{"method": "salmon", "tools": ["salmon"], "refs": ["data/raw/refs/gencode/salmon_index/info.json"]}]},
        {"stage": "mutanome_enumeration",
         "produces": "full class-I 8-11mer lossless peptide universe (shared by the lossless arms; the "
                     "pvac arm generates its own set from the same base variants)",
         "methods": [{"method": "VEP-REST+lossless", "tools": ["bcftools"], "refs": [_G, _G + ".fai"],
                      "note": "the EXACT provider that produced the frozen Hu_287 universe: Ensembl VEP REST "
                              "consequence + reference-protein enumeration (event_b.lossless_peptide_generation, "
                              "responses cached per-patient under ensembl_cache/); bcftools left-aligns/splits "
                              "indels against GRCh38 before enumeration. Needs network, no local VEP cache/GTF."},
                     {"method": "VEP+lossless", "tools": ["vep"],
                      "refs": ["data/raw/refs/vep/homo_sapiens", _G, "data/raw/refs/gencode/gencode.v44.annotation.gtf"]},
                     {"method": "pvacseq", "tools": ["pvacseq", "vep"],
                      "refs": ["data/raw/refs/vep/homo_sapiens", _G]}]},
    ]
    if resolve is None:
        def resolve(name):
            return resolve_tool(name, which)
    for s in stages:
        runnable_method, per_method = None, []
        for m in s["methods"]:
            missing_tools = [t for t in m["tools"] if resolve(t) is None]
            missing_refs = [p for p in m["refs"] if not ref_exists(p)]
            per_method.append({"method": m["method"], "missing_tools": missing_tools,
                               "missing_refs": missing_refs, **({"note": m["note"]} if "note" in m else {})})
            if not missing_tools and not missing_refs and runnable_method is None:
                runnable_method = m["method"]
        s["method_status"] = per_method
        if runnable_method is not None:
            s["status"], s["runnable_method"] = "RUNNABLE", runnable_method
        else:
            s["status"] = "NOT_EVALUABLE"
            s["reason"] = "; ".join(
                f"{p['method']}: " + ", ".join(
                    ([f"missing tools {p['missing_tools']}"] if p["missing_tools"] else [])
                    + ([f"missing refs {p['missing_refs']}"] if p["missing_refs"] else []))
                for p in per_method)
    # scoring is special: recognize local PRIME/MixMHCpred by path, but stay upstream-blocked
    local = {"PRIME": str(PRIME_BIN) if PRIME_BIN.exists() else None,
             "MixMHCpred": str(MIXMHC_BIN) if MIXMHC_BIN.exists() else None}
    stages.append({
        "stage": "scoring_prime_epicurus", "produces": "genuine PRIME AND frozen Epicurus over the IDENTICAL universe",
        "resolved_tools": local, "status": "NOT_EVALUABLE",
        "reason": ("local PRIME/MixMHCpred present but UPSTREAM-BLOCKED: needs the re-enumerated candidate "
                   "universe + class-I HLA before any peptide can be scored" if (local["PRIME"] or local["MixMHCpred"])
                   else "PRIME/MixMHCpred binaries not found on disk")})
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
