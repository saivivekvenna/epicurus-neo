#!/usr/bin/env bash
# Miller patient — class-I HLA typing (OptiType) from the NORMAL exome (frozen prereg §2).
# LOCKED_TEST: no recognition label read. Waits for the normal BAM, extracts MHC-region reads, runs OptiType.
#
# Usage:
#   miller_hu287_hla.sh              # full run: extract MHC reads -> OptiType -> write provenance
#   miller_hu287_hla.sh provenance   # regenerate ONLY HLA_PROVENANCE.json from the EXISTING result+files
#                                    #   (does NOT rerun alignment/OptiType; alleles+objective stay identical)
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
export MAMBA_ROOT_PREFIX="$ROOT/data/raw/tools/micromamba"
MM="$ROOT/data/raw/tools/bin/micromamba"
PATIENT_ID="${PATIENT_ID:-Hu_287}"
PATIENT_SLUG="$(printf '%s' "$PATIENT_ID" | tr '[:upper:]' '[:lower:]')"
OUT="$ROOT/data/raw/miller_ipv/$PATIENT_SLUG/hla"; mkdir -p "$OUT"
if [[ "$PATIENT_ID" = "Hu_287" ]]; then
  PROV="$ROOT/artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction"
else
  PROV="$ROOT/artifacts/milestone_8_generalization/patients/$PATIENT_ID"
fi
NBAM="$ROOT/data/raw/miller_ipv/$PATIENT_SLUG/somatic/${PATIENT_ID}_N.md.bam"
HLA1="$OUT/hla_1.fq"; HLA2="$OUT/hla_2.fq"
REGION="6:29800000-33600000"                     # GRCh38 (Ensembl contig '6') MHC region — HYPHEN (samtools)
LOG="$OUT/run.log"
say(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

# Regenerate HLA_PROVENANCE.json from the EXACT completed files + OptiType result. Records per-file
# relative path/sha256/size for the chain (normal MD BAM+index, extracted FASTQs, result TSV), the
# micromamba/OptiType/GLPK/razers3 environment identity, and the result objective/read count/region/path.
write_provenance(){
  local RES="$1" GLPK="$2"
  local MM_VERSION OPTITYPE_PKG GLPK_PKG RAZERS_PKG
  MM_VERSION="$("$MM" --version 2>/dev/null | head -1)"
  OPTITYPE_PKG="$("$MM" list -n hla 2>/dev/null | awk 'tolower($1)=="optitype"{print $2" "$3; exit}')"
  GLPK_PKG="$("$MM" list -n hla 2>/dev/null | awk 'tolower($1)=="glpk"{print $2" "$3; exit}')"
  RAZERS_PKG="$("$MM" list -n hla 2>/dev/null | awk 'tolower($1)=="razers3"{print $2" "$3; exit}')"
  "$ROOT/.venv/bin/python" - "$RES" "$PROV/HLA_PROVENANCE.json" "$ROOT" "$REGION" \
      "$NBAM" "$HLA1" "$HLA2" "$GLPK" "$MM_VERSION" "$OPTITYPE_PKG" "$GLPK_PKG" "$RAZERS_PKG" "$PATIENT_ID" <<'PY'
import sys, json, csv, os, hashlib
from pathlib import Path
(res, out, root, region, nbam, hla1, hla2, glpk_solver,
 mm_version, optitype_pkg, glpk_pkg, razers_pkg, patient_id) = sys.argv[1:14]

def sha256(p, chunk=1 << 20):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(chunk), b""):
            h.update(b)
    return h.hexdigest()

def rec(path):
    p = Path(path)
    if not p.is_file():
        raise SystemExit(f"FAIL-CLOSED: provenance file missing: {path}")
    return {"path": os.path.relpath(p, root), "sha256": sha256(p), "size_bytes": p.stat().st_size}

row = list(csv.DictReader(open(res), delimiter="\t"))[0]

def norm(a):
    a = str(a).strip()
    return None if not a or a in {"nan", ""} else "HLA-" + a           # OptiType gives e.g. A*02:01
alleles = sorted({x for x in (norm(row.get(k)) for k in ("A1", "A2", "B1", "B2", "C1", "C2")) if x})
if len(alleles) != 6:                                                  # class-I A/B/C, two each
    raise SystemExit(f"FAIL-CLOSED: expected 6 class-I alleles, got {len(alleles)}: {alleles}")

files = {"normal_md_bam": rec(nbam), "normal_md_bam_index": rec(nbam + ".bai"),
         "hla_fastq_1": rec(hla1), "hla_fastq_2": rec(hla2), "optitype_result_tsv": rec(res)}

doc = {
    "patient_id": patient_id,
    "isolation": "LOCKED_TEST: no label read",
    "tool": "OptiType (bioconda osx-64 via micromamba/Rosetta)",
    "solver": glpk_solver,
    "source": "normal exome MHC-region reads -> razers3 -> OptiType",
    "extraction_region_grch38": region,
    "class_i_alleles": alleles,
    "objective": float(row["Objective"]),
    "reads_used": float(row["Reads"]),
    "result_path": os.path.relpath(Path(res), root),
    "environment": {
        "micromamba_version": mm_version,
        "optitype_pkg": optitype_pkg,
        "glpk_pkg": glpk_pkg,
        "razers3_pkg": razers_pkg,
        "glpk_solver": glpk_solver,
        "mamba_env": "hla",
        "mamba_root_prefix": os.path.relpath(Path(os.environ["MAMBA_ROOT_PREFIX"]), root),
    },
    "provenance_files": files,
    "raw": row,
}
json.dump(doc, open(out, "w"), indent=2)
print("HLA class-I:", alleles, "| objective:", doc["objective"], "| reads:", doc["reads_used"])
print("provenance files:", ", ".join(f"{k}={v['sha256'][:12]}…" for k, v in files.items()))
PY
}

# ---- provenance-only mode: rebuild the JSON from existing artifacts; never reruns HLA ----------------
if [ "${1:-}" = "provenance" ]; then
  GLPK_VERSION="$("$MM" run -n hla glpsol --version | head -1)"
  RES="$(find "$OUT/optitype" -name "*_result.tsv" | head -1)"
  [ -n "$RES" ] || { echo "FAIL-CLOSED: no OptiType *_result.tsv under $OUT/optitype"; exit 1; }
  echo "provenance-only: result $RES"
  write_provenance "$RES" "$GLPK_VERSION"
  echo "wrote $PROV/HLA_PROVENANCE.json"
  exit 0
fi

# ---- full run --------------------------------------------------------------------------------------
: > "$LOG"
say "waiting for normal BAM $NBAM"
until [[ -f "$NBAM" && -f "${NBAM}.bai" ]]; do sleep 30; done
say "extracting MHC-region reads ($REGION) -> paired FASTQ"
samtools view -b "$NBAM" "$REGION" 2>>"$LOG" \
  | samtools collate -Oun128 - 2>>"$LOG" \
  | samtools fastq -1 "$HLA1" -2 "$HLA2" -0 /dev/null -s /dev/null -n 2>>"$LOG"
NPAIR=$(( $(wc -l < "$HLA1") / 4 ))
if [ "$NPAIR" -le 0 ]; then say "FAIL-CLOSED: 0 read pairs extracted from $REGION (aborting; not running OptiType)"; exit 1; fi
say "OptiType (DNA) on $NPAIR read pairs"
GLPK_VERSION="$("$MM" run -n hla glpsol --version | head -1)"
"$ROOT/.venv/bin/python" - "$GLPK_VERSION" <<'PY'
import re, sys
m = re.search(r"(?:v|Solver )([0-9]+)\.([0-9]+)", sys.argv[1])
if not m or (int(m.group(1)), int(m.group(2))) < (4, 58):
    raise SystemExit(f"FAIL-CLOSED: OptiType requires GLPK >=4.58; found {sys.argv[1]!r}")
PY
say "GLPK solver: $GLPK_VERSION"
rm -rf "$OUT/optitype"; mkdir -p "$OUT/optitype"
"$MM" run -n hla OptiTypePipeline.py --dna -i "$HLA1" "$HLA2" \
     --outdir "$OUT/optitype" 2>>"$LOG"
RES="$(find "$OUT/optitype" -name "*_result.tsv" | head -1)"
say "OptiType result: $RES"
write_provenance "$RES" "$GLPK_VERSION"
say "wrote $PROV/HLA_PROVENANCE.json"
say "=== HLA DONE ==="
