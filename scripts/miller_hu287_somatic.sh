#!/usr/bin/env bash
# Miller Hu_287 — WES align -> Mutect2 matched somatic (frozen per RECONSTRUCTION_METHOD_PREREG.md).
# LOCKED_TEST: reads NO recognition label. Idempotent-ish (skips steps whose outputs exist).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT"
export PATH="$ROOT/data/raw/tools/bin:$PATH"          # provides the `python` shim gatk needs
GATK="$ROOT/data/raw/tools/gatk-4.5.0.0/gatk"
REF="$ROOT/data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
DICT="${REF%.fa}.dict"
BED="$ROOT/data/raw/refs/intervals/cds.merged.bed"
FQ="$ROOT/data/raw/miller_ipv/hu_287/fastq"
OUT="$ROOT/data/raw/miller_ipv/hu_287/somatic"
PROV="$ROOT/artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction"
mkdir -p "$OUT" "$PROV"
THREADS="$(sysctl -n hw.ncpu 2>/dev/null || echo 4)"
LOG="$OUT/run.log"; : > "$LOG"
say(){ echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

# frozen sample names / read groups (prereg §2) — plain vars (macOS bash 3.2 has no assoc arrays)
RUN_N=SRR24836184; SM_N=Hu_287_N       # normal exome
RUN_T=SRR24836169; SM_T=Hu_287_T       # tumor exome

align(){  # $1=run  $2=sample_name
  local run="$1" sm="$2"
  local md="$OUT/${sm}.md.bam"
  if [[ -f "$md" ]]; then say "align: $md exists, skip"; return; fi
  say "align bwa-mem $run SM=$sm (RG)"
  bwa mem -t "$THREADS" -R "@RG\tID:${run}\tSM:${sm}\tPL:ILLUMINA\tLB:${sm}_xE" \
      "$REF" "$FQ/${run}_1.fastq" "$FQ/${run}_2.fastq" 2>>"$LOG" \
    | samtools sort -@ "$THREADS" -o "$OUT/${sm}.sorted.bam" - 2>>"$LOG"
  samtools index "$OUT/${sm}.sorted.bam"
  say "align MarkDuplicates $sm"
  "$GATK" MarkDuplicates -I "$OUT/${sm}.sorted.bam" -O "$md" -M "$OUT/${sm}.dupmetrics.txt" 2>>"$LOG"
  samtools index "$md"
  rm -f "$OUT/${sm}.sorted.bam" "$OUT/${sm}.sorted.bam.bai"
}

say "=== Hu_287 somatic reconstruction (threads=$THREADS) ==="
say "waiting for BWA index sentinels (.bwt/.sa)"
until [[ -f "${REF}.bwt" && -f "${REF}.sa" ]]; do sleep 30; done
say "BWA index ready"
[[ -f "$DICT" ]] || { say "CreateSequenceDictionary"; "$GATK" CreateSequenceDictionary -R "$REF" 2>>"$LOG"; }
align "$RUN_N" "$SM_N"
align "$RUN_T" "$SM_T"

VCF="$OUT/Hu_287.mutect2.vcf.gz"
if [[ ! -f "$VCF" ]]; then
  say "Mutect2 matched (--normal $SM_N) over CDS intervals (label-blind)"
  "$GATK" Mutect2 -R "$REF" \
     -I "$OUT/${SM_T}.md.bam" -I "$OUT/${SM_N}.md.bam" --normal "$SM_N" \
     -L "$BED" --interval-padding 0 \
     -O "$VCF" 2>>"$LOG"
fi
FILT="$OUT/Hu_287.mutect2.filtered.vcf.gz"
if [[ ! -f "$FILT" ]]; then
  say "FilterMutectCalls"
  "$GATK" FilterMutectCalls -R "$REF" -V "$VCF" -O "$FILT" 2>>"$LOG"
fi
PASS="$OUT/Hu_287.somatic.pass.vcf.gz"
say "select PASS"
bcftools view -f PASS "$FILT" -Oz -o "$PASS" 2>>"$LOG"; bcftools index -t "$PASS"

NPASS="$(bcftools view -H "$PASS" | wc -l | tr -d ' ')"
NMISS="$(bcftools view -H -v snps "$PASS" | wc -l | tr -d ' ')"
say "PASS somatic variants: $NPASS (snv rows ~$NMISS)"

cat > "$PROV/SOMATIC_PROVENANCE.json" <<EOF
{
  "patient_id": "Hu_287", "isolation": "LOCKED_TEST: no label read",
  "reference": "Ensembl GRCh38 primary assembly r110",
  "aligner": "bwa mem + samtools sort + GATK MarkDuplicates",
  "read_groups": {"tumor": "SM=${SM_T}", "normal": "SM=${SM_N}"},
  "caller": "GATK4 Mutect2 matched-normal (--normal ${SM_N}) + FilterMutectCalls",
  "germline_resource": "NONE (matched-normal-only; documented specificity deviation, prereg §2/§4)",
  "calling_region": "label-blind Ensembl CDS merged BED ($BED)",
  "n_pass_variants": ${NPASS:-0},
  "outputs": {"filtered_vcf": "$FILT", "pass_vcf": "$PASS"},
  "tool_versions": {"bwa": "$(bwa 2>&1 | grep -i version | head -1)", "samtools": "$(samtools --version | head -1)", "gatk": "4.5.0.0"}
}
EOF
say "wrote $PROV/SOMATIC_PROVENANCE.json"
say "=== DONE ==="
