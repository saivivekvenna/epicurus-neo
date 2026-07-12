# Sid identical-input end-to-end benchmark — PROTOCOL (frozen before competitor results)

_Falsification-first. An honest LOSS or NOT_EVALUABLE is an acceptable outcome. Post-hoc n=1 patient / 3
recognized positives — always labelled as such; never a powered superiority claim._

## 0. Leakage correction that motivated this protocol (do NOT preserve the old L3 claim)

The prior "lossless 1/3 → 3/3" result is **INVALID as an end-to-end benchmark**. `scripts/osteosarc_peptide_recovery.py`
hard-codes `TARGETS = {ASPM, MAP2, DYNC1H1}` — the exact recognized positives — and generates candidates
**only for those**. Target selection leaks the answer even though Hudson labels are joined afterwards.
Verified: the existing recovered-candidate set covers the 3 positives but only **10.2%** of the eligible
variant universe → `assert_generation_label_blind` (src/event_b/sid_benchmark.py) **correctly fails** it.
Also `artifacts/osteosarc_audit/somatic.vcf.gz` (1,213 records) contains **none** of the 3 target
coordinates; the targets live in the website VAF table. The presentation-only top-20 "recovery" is only a
**target-conditioned peptide-reconstruction sensitivity test**, NOT end-to-end proof.

## 1. Guardrail (hard, test-enforced)

Candidate-generation input MUST be the **complete, label-blind variant universe**, never a TARGETS list.
`tests/test_sid_benchmark.py` fails if a generation's variant set is a positive-selected subset or covers
< 95% of the eligible universe.

## 2. Primary patient + inputs

Public osteosarc / Sid. Complete label-blind universe = **all 200 unique `variant_id`** in
`data/raw/osteosarc/site_cache/variant_vafs_long.tsv` (longitudinal, multi-pipeline). **Consequence
eligibility, declared BEFORE any label join:** a variant is class-I-eligible iff its consequence ∈
{missense, frameshift, stop_gained, inframe_deletion/insertion/indel} AND it has somatic tumor read
support (max tumor `alt_reads` > 0) ⇒ **147 eligible** variants (127 missense + 10 frameshift + 7
stop-gain + 3 in-frame deletions;
excludes synonymous/intron/splice/UTR/blank). All 3 exact recognized mutations are eligible: ASPM
`chr1:197102716`, DYNC1H1 `chr14:101980529`, and MAP2 `chr2:209694772`.

Those exact Hudson recognition labels (IFNγ/TCR) are **EVALUATION-ONLY**, joined only after every pipeline
output is frozen. Other mutations in the same genes are not positives. Labels never gate the generation
universe.

## 3. Metrics (mutation identity)

Normalize every pipeline's output to stable `variant_id` (GENE-chrom-pos); collapse peptide/HLA rows to the
best-scoring row per variant BEFORE top-20 so duplicates cannot inflate mutation hits.
- **Primary:** mutation-level recognized **hits@20** from common inputs.
- Also: candidate-generation recall (positives reaching the candidate set), rankable recall (positives with
  a scoreable peptide+HLA), conditional ranking (rank of recognized among generated), **stage of first
  loss** (variant-calling / generation / HLA / scoring / ranking), runtime, reproducibility, input
  completeness.

## 4. Arms (Epicurus configs FROZEN before competitor results)

Epicurus (all consume the complete 147-variant universe, label-blind):
1. **lossless generation + presentation-only MixMHCpred** (NO PRIME immunogenicity).
2. **lossless generation + genuine PRIME** (`data/raw/tools/PRIME`).
3. **full frozen Epicurus** (`configs/frozen/epicurus_v0_1.json`).

Controls: **pVAC-style + same scorer** (binding-first selection then the same MixMHCpred/PRIME) — a
documented approximation, never labelled as the pVACtools binary.

## 5. Competitors (audited; run from primary official code/docs where feasible)

Peptide scorers are NOT end-to-end pipelines and are excluded as competitors.

| pipeline | boundary needed | status here |
|---|---|---|
| pVACtools/pVACseq (griffithlab) | VEP-annotated VCF + HLA (+expr) | **PARTIAL / boundary-mismatch**: an existing `data/raw/osteosarc/pvactools_all_epitopes.tsv` (**2025.01**) exists, but it is a DIFFERENT input boundary from the 200-variant longitudinal VAF universe (single-pipeline VCF vs multi-pipeline table) ⇒ **not identical inputs**; reported separately, not as a matched arm. A fresh Docker pVACseq run at the 147-variant boundary is the matched comparator (feasibility TBD). |
| OpenVax Vaxrank | annotated VCF + tumor RNA BAM + isovar/varcode | needs an RNA BAM we do not have at this boundary → **NOT_EVALUABLE** (input boundary) |
| nextNEOpi (nextflow) | raw tumor/normal FASTQ + many tools | raw FASTQ unavailable; nextflow not installed → **NOT_EVALUABLE** |
| NeoDisc | registration/licensed distribution | **NOT_EVALUABLE** — license/registration blocker (explicit) |
| CNNeoPP | research code + peptide inputs | audit; likely peptide-scorer-like, not end-to-end |

Environment audited: `docker` present; `mhcflurry 2.2.1`; genuine `PRIME`+`MixMHCpred` under
`data/raw/tools` (gitignored); **no** pvactools/netMHCpan/vaxrank/nextflow/Rscript in PATH; network up.

## 6. Rules

- Freeze Epicurus configs before any competitor result. No per-tool option tuning on Sid labels — use
  documented defaults / preregistered common choices. Record every tool version/commit/license.
- Same raw/processed starting boundary explicit per arm. If a competitor cannot accept the 147-variant
  processed boundary (or needs unavailable raw FASTQ), mark **NOT_EVALUABLE** — never invent a proxy.
- Generic `expected=None` generation (no per-target verification), recording per-variant success/failure;
  never silently exclude a variant.
- If full all-variant Epicurus generation is infeasible in this environment, mark the Epicurus end-to-end
  arm **NOT_EVALUABLE** and report exactly why + the runnable path.

## 7. Verdict vocabulary

WIN / TIE / LOSS / NOT_EVALUABLE per arm×competitor, at the declared boundary, mutation-level hits@20 with
the stage-of-first-loss for every missed positive. n=1/3 post-hoc — descriptive only.
