"""Miller IPV Hu_287 — raw-input RECONSTRUCTION driver (LOCKED_TEST; labels never read here).

Turns the downloaded SRA trio into north-star inputs as far as the installed toolchain honestly allows,
and records a machine-actionable dependency/reference manifest for everything that cannot yet run. Nothing
is fabricated: a stage with a missing tool/reference emits an explicit NOT_EVALUABLE with the exact install/
acquire command, never a vague blocker.

Subcommands:
    deps      write DEPENDENCY_MANIFEST.{json,md} (tool versions + references + per-stage status + hints)
    convert   fasterq-dump each downloaded .sra -> paired FASTQ (records read counts + provenance)

The candidate universe MUST be re-enumerated from WES (not the IPV-prefiltered tested set); genuine PRIME
and frozen Epicurus share that identical universe. HLA typing, somatic calling, and mutanome enumeration are
the hard blockers on the current toolchain (see the manifest); expression (salmon) is the one north-star
input reachable now once its reference index is built.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from benchmark.miller_download import reconstruction_stages

ROOT = Path(__file__).resolve().parents[1]
HU287 = ROOT / "data/raw/miller_ipv/hu_287"
ART = ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction"
PRIME_BIN = ROOT / "data/raw/tools/PRIME/PRIME"
MIXMHC_BIN = ROOT / "data/raw/tools/MixMHCpred/MixMHCpred"

# Reference data required by downstream stages (not auto-fetched; acquire commands are machine-actionable).
REFERENCES = [
    {"name": "GRCh38_primary_assembly", "purpose": "WES alignment + somatic calling",
     "dest": "data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
     "acquire": "curl -L https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/"
                "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | gunzip > <dest>; bwa index <dest>; "
                "samtools faidx <dest>"},
    {"name": "GENCODE_v44_transcripts", "purpose": "RNA quantification (salmon index)",
     "dest": "data/raw/refs/gencode/gencode.v44.transcripts.fa",
     "acquire": "curl -L https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/"
                "gencode.v44.transcripts.fa.gz | gunzip > <dest>; salmon index -t <dest> -i "
                "data/raw/refs/gencode/salmon_index -k 31"},
    {"name": "VEP_cache_GRCh38", "purpose": "variant annotation for mutanome enumeration (pVACtools)",
     "dest": "data/raw/refs/vep/homo_sapiens",
     "acquire": "curl -L https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/"
                "homo_sapiens_vep_110_GRCh38.tar.gz | tar xz -C data/raw/refs/vep"},
    {"name": "IMGT_HLA_dna", "purpose": "class-I HLA typing (OptiType/arcasHLA reference)",
     "dest": "data/raw/refs/hla/hla_reference_dna.fasta",
     "acquire": "OptiType ships hla_reference_dna.fasta; arcasHLA: 'arcasHLA reference --version latest'"},
]

# Absent tools -> exact reproducible install command (bioconda is the reproducible native path; no Docker).
INSTALL_HINTS = {
    "hla_typing_classI": "conda install -c bioconda optitype razers3   # or: conda install -c bioconda arcas-hla",
    "somatic_calling": "conda install -c bioconda gatk4   # Mutect2 tumor-vs-normal; or bioconda strelka",
    "mutanome_enumeration": "pip install pvactools && conda install -c bioconda ensembl-vep   # + VEP cache",
}

_VERSION_CMD = {"fasterq-dump": ["fasterq-dump", "--version"], "bwa": ["bwa"],
                "samtools": ["samtools", "--version"], "bcftools": ["bcftools", "--version"],
                "salmon": ["salmon", "--version"]}


def _probe_version(tool: str) -> str | None:
    exe = shutil.which(tool)
    if exe is None:
        return None
    try:
        out = subprocess.run(_VERSION_CMD.get(tool, [tool, "--version"]), capture_output=True, text=True,
                             timeout=30)
        text = (out.stdout + out.stderr).strip().splitlines()
        return next((ln.strip() for ln in text if ln.strip()), exe)
    except Exception:
        return exe


def tool_versions(which=shutil.which) -> dict:
    tv = {t: _probe_version(t) for t in _VERSION_CMD}
    tv["PRIME"] = str(PRIME_BIN) if PRIME_BIN.exists() else None
    tv["MixMHCpred"] = str(MIXMHC_BIN) if MIXMHC_BIN.exists() else None
    return tv


def reference_status(exists=lambda p: (ROOT / p).exists()) -> list[dict]:
    return [{**r, "present": bool(exists(r["dest"]))} for r in REFERENCES]


def dependency_manifest(which=shutil.which, exists=None) -> dict:
    exists = (lambda p: (ROOT / p).exists()) if exists is None else exists
    stages = reconstruction_stages(which=which)
    for s in stages:                                   # attach machine-actionable install hint where blocked
        if s["status"] == "NOT_EVALUABLE" and s["stage"] in INSTALL_HINTS:
            s["install_hint"] = INSTALL_HINTS[s["stage"]]
    return {
        "cohort": "miller_ipv", "patient_id": "Hu_287", "tranche": "T1",
        "isolation": "LOCKED_TEST: labels define the metric only; never used for HLA/expr/calling/generation/ranking",
        "tool_versions": tool_versions(which=which),
        "references": reference_status(exists=exists),
        "reconstruction_stages": stages,
        "note": "PRIME/MixMHCpred binaries are present on disk (gitignored) but scoring is upstream-blocked "
                "until the candidate universe is re-enumerated. HLA typing, somatic calling, and mutanome "
                "enumeration are the hard blockers on this toolchain; expression via salmon is reachable once "
                "the GENCODE index is built. No released processed Hu_287 VCF/TPM/HLA exists (web 2026-07-12).",
    }


def write_deps() -> dict:
    ART.mkdir(parents=True, exist_ok=True)
    man = dependency_manifest()
    (ART / "DEPENDENCY_MANIFEST.json").write_text(json.dumps(man, indent=2) + "\n")
    (ART / "DEPENDENCY_MANIFEST.md").write_text(_deps_md(man))
    return man


def _deps_md(m: dict) -> str:
    L = ["# Miller Hu_287 — reconstruction dependency/reference manifest\n",
         f"_LOCKED_TEST. {m['note']}_\n", "\n## Installed tools\n"]
    for t, v in m["tool_versions"].items():
        L.append(f"- **{t}**: {v or 'ABSENT'}")
    L.append("\n## Per-stage status (machine-actionable)\n")
    L.append("| stage | status | produces | missing / reason | install hint |")
    L.append("|---|---|---|---|---|")
    for s in m["reconstruction_stages"]:
        L.append(f"| {s['stage']} | **{s['status']}** | {s['produces']} | "
                 f"{s.get('reason', s.get('resolved_tools', ''))} | {s.get('install_hint', '')} |")
    L.append("\n## References required (not auto-fetched)\n")
    L.append("| name | purpose | present | dest |")
    L.append("|---|---|:--:|---|")
    for r in m["references"]:
        L.append(f"| {r['name']} | {r['purpose']} | {'y' if r['present'] else 'n'} | `{r['dest']}` |")
    L.append("\n### Acquire commands\n")
    for r in m["references"]:
        L.append(f"- **{r['name']}**: `{r['acquire']}`")
    return "\n".join(L) + "\n"


def convert_sra(sra: Path, outdir: Path) -> dict:
    """fasterq-dump one .sra -> paired FASTQ; record read counts. Requires fasterq-dump on PATH."""
    if shutil.which("fasterq-dump") is None:
        return {"sra": str(sra), "status": "NOT_EVALUABLE", "reason": "fasterq-dump absent"}
    outdir.mkdir(parents=True, exist_ok=True)
    run = sra.stem
    cmd = ["fasterq-dump", "--split-files", "--skip-technical", "-O", str(outdir), "-t", str(outdir),
           str(sra)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        return {"sra": str(sra), "status": "ERROR", "returncode": proc.returncode,
                "stderr": proc.stderr[-2000:]}
    fastqs = sorted(str(p.name) for p in outdir.glob(f"{run}*.fastq"))
    reads = {}
    for fq in fastqs:
        n = sum(1 for _ in open(outdir / fq)) // 4
        reads[fq] = n
    return {"sra": str(sra), "run": run, "status": "OK", "fastqs": fastqs, "reads_per_file": reads}


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    cmd = argv[0] if argv else "deps"
    if cmd == "deps":
        man = write_deps()
        runnable = [s["stage"] for s in man["reconstruction_stages"] if s["status"] == "RUNNABLE"]
        blocked = [s["stage"] for s in man["reconstruction_stages"] if s["status"] != "RUNNABLE"]
        print("RUNNABLE:", runnable)
        print("NOT_EVALUABLE:", blocked)
        print("wrote", ART / "DEPENDENCY_MANIFEST.json")
        return 0
    if cmd == "convert":
        ART.mkdir(parents=True, exist_ok=True)
        fqdir = HU287 / "fastq"
        results = []
        for sra in sorted(HU287.glob("*.sra")):
            print("converting", sra.name, flush=True)
            results.append(convert_sra(sra, fqdir))
        (ART / "CONVERT_PROVENANCE.json").write_text(json.dumps(results, indent=2) + "\n")
        for r in results:
            print(r.get("run"), r["status"], r.get("reads_per_file"))
        return 0
    print("usage: miller_reconstruct.py [deps|convert]")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
