#!/usr/bin/env bash
# Miller patient — label-blind tumor-RNA spliced alignment (HISAT2) for RNA-alt evidence (evidence-only).
# Never a hard filter. Produces a coordinate-sorted RNA BAM; the universe script pileups it at somatic sites.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
export MAMBA_ROOT_PREFIX="$ROOT/data/raw/tools/micromamba"
MM="$ROOT/data/raw/tools/bin/micromamba"
REF="$ROOT/data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
IDX="$ROOT/data/raw/refs/GRCh38/hisat2_index/grch38"
PATIENT_ID="${PATIENT_ID:-Hu_287}"
PATIENT_SLUG="$(printf '%s' "$PATIENT_ID" | tr '[:upper:]' '[:lower:]')"
RNA_RUN="${TUMOR_RNA_RUN:-SRR24836183}"
FQ="$ROOT/data/raw/miller_ipv/$PATIENT_SLUG/fastq"
OUT="$ROOT/data/raw/miller_ipv/$PATIENT_SLUG/rna"; mkdir -p "$OUT" "$(dirname "$IDX")"
LOG="$OUT/run.log"; : > "$LOG"
say(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
THREADS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"

if [[ ! -f "${IDX}.1.ht2" ]]; then
  say "hisat2-build GRCh38 index (~1-2h)"
  $MM run -p "$MAMBA_ROOT_PREFIX/envs/rna" hisat2-build -p "$THREADS" "$REF" "$IDX" >>"$LOG" 2>&1
fi
BAM="$OUT/${PATIENT_ID}_tumor_rna.sorted.bam"
if [[ ! -f "$BAM" ]]; then
  say "hisat2 spliced align tumor RNA ($RNA_RUN)"
  $MM run -p "$MAMBA_ROOT_PREFIX/envs/rna" hisat2 -p "$THREADS" -x "$IDX" \
      -1 "$FQ/${RNA_RUN}_1.fastq" -2 "$FQ/${RNA_RUN}_2.fastq" --no-unal 2>>"$LOG" \
    | $MM run -p "$MAMBA_ROOT_PREFIX/envs/rna" samtools sort -@ "$THREADS" -o "$BAM" - 2>>"$LOG"
  $MM run -p "$MAMBA_ROOT_PREFIX/envs/rna" samtools index "$BAM"
fi
say "RNA BAM ready: $BAM"
say "=== RNA DONE ==="
