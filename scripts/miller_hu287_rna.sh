#!/usr/bin/env bash
# Miller Hu_287 — label-blind tumor-RNA spliced alignment (HISAT2) for RNA-alt evidence (evidence-only).
# Never a hard filter. Produces a coordinate-sorted RNA BAM; the universe script pileups it at somatic sites.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
export MAMBA_ROOT_PREFIX="$ROOT/data/raw/tools/micromamba"
MM="$ROOT/data/raw/tools/bin/micromamba"
REF="$ROOT/data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
IDX="$ROOT/data/raw/refs/GRCh38/hisat2_index/grch38"
FQ="$ROOT/data/raw/miller_ipv/hu_287/fastq"
OUT="$ROOT/data/raw/miller_ipv/hu_287/rna"; mkdir -p "$OUT" "$(dirname "$IDX")"
LOG="$OUT/run.log"; : > "$LOG"
say(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }
THREADS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"

if [[ ! -f "${IDX}.1.ht2" ]]; then
  say "hisat2-build GRCh38 index (~1-2h)"
  $MM run -p "$MAMBA_ROOT_PREFIX/envs/rna" hisat2-build -p "$THREADS" "$REF" "$IDX" >>"$LOG" 2>&1
fi
BAM="$OUT/Hu_287_tumor_rna.sorted.bam"
if [[ ! -f "$BAM" ]]; then
  say "hisat2 spliced align tumor RNA (SRR24836183)"
  $MM run -p "$MAMBA_ROOT_PREFIX/envs/rna" hisat2 -p "$THREADS" -x "$IDX" \
      -1 "$FQ/SRR24836183_1.fastq" -2 "$FQ/SRR24836183_2.fastq" --no-unal 2>>"$LOG" \
    | $MM run -p "$MAMBA_ROOT_PREFIX/envs/rna" samtools sort -@ "$THREADS" -o "$BAM" - 2>>"$LOG"
  $MM run -p "$MAMBA_ROOT_PREFIX/envs/rna" samtools index "$BAM"
fi
say "RNA BAM ready: $BAM"
say "=== RNA DONE ==="
