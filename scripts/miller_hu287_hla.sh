#!/usr/bin/env bash
# Miller Hu_287 — class-I HLA typing (OptiType) from the NORMAL exome (frozen prereg §2).
# LOCKED_TEST: no recognition label read. Waits for the normal BAM, extracts MHC-region reads, runs OptiType.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
export MAMBA_ROOT_PREFIX="$ROOT/data/raw/tools/micromamba"
MM="$ROOT/data/raw/tools/bin/micromamba"
OUT="$ROOT/data/raw/miller_ipv/hu_287/hla"; mkdir -p "$OUT"
PROV="$ROOT/artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction"
NBAM="$ROOT/data/raw/miller_ipv/hu_287/somatic/Hu_287_N.md.bam"
REGION="6:29800000-33600000"                     # GRCh38 (Ensembl contig '6') MHC region — HYPHEN (samtools)
LOG="$OUT/run.log"; : > "$LOG"
say(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

say "waiting for normal BAM $NBAM"
until [[ -f "$NBAM" && -f "${NBAM}.bai" ]]; do sleep 30; done
say "extracting MHC-region reads ($REGION) -> paired FASTQ"
samtools view -b "$NBAM" "$REGION" 2>>"$LOG" \
  | samtools collate -Oun128 - 2>>"$LOG" \
  | samtools fastq -1 "$OUT/hla_1.fq" -2 "$OUT/hla_2.fq" -0 /dev/null -s /dev/null -n 2>>"$LOG"
NPAIR=$(( $(wc -l < "$OUT/hla_1.fq") / 4 ))
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
"$MM" run -n hla OptiTypePipeline.py --dna -i "$OUT/hla_1.fq" "$OUT/hla_2.fq" \
     --outdir "$OUT/optitype" 2>>"$LOG"
RES="$(find "$OUT/optitype" -name "*_result.tsv" | head -1)"
say "OptiType result: $RES"
"$ROOT/.venv/bin/python" - "$RES" "$PROV/HLA_PROVENANCE.json" "$GLPK_VERSION" <<'PY'
import sys, json, csv
res, out, glpk_version = sys.argv[1], sys.argv[2], sys.argv[3]
row = list(csv.DictReader(open(res), delimiter="\t"))[0]
def norm(a):
    a = str(a).strip()
    return None if not a or a in {"nan", ""} else "HLA-" + a  # OptiType gives e.g. A*02:01
alleles = sorted({x for x in (norm(row.get(k)) for k in ("A1","A2","B1","B2","C1","C2")) if x})
json.dump({"patient_id": "Hu_287", "isolation": "LOCKED_TEST: no label read",
           "tool": "OptiType (bioconda osx-64 via micromamba/Rosetta)",
           "solver": glpk_version,
           "source": "normal exome MHC-region reads -> razers3 -> OptiType",
           "class_i_alleles": alleles, "raw": row}, open(out, "w"), indent=2)
print("HLA class-I:", alleles)
PY
say "wrote $PROV/HLA_PROVENANCE.json"
say "=== HLA DONE ==="
