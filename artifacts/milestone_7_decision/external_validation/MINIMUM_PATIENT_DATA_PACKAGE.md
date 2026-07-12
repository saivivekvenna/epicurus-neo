# Minimum patient data package — exact de-identified schema to request

The canonical, minimal, **de-identified** package that makes a patient usable for the north-star benchmark:
*from identical WES/RNA/HLA inputs, does Epicurus place more experimentally recognized neoantigens in the
final top-20 than standard pVAC-style generation + genuine PRIME (patient-level paired Δ recognized
hits@20)?* Request DATA, never PHI. Controlled access is acceptable and expected. Nothing below requires a
patient identity, date of birth, or free-text clinical note.

This schema is what we ask authors/DACs to provide (or to confirm is already deposited). It maps onto the
repo's unified candidate schema consumed by `event_b.prime_transfer.external_validate` and the four-arm
harness (`src/benchmark/four_arm.py`).

---

## Tier 0 — the reachability inputs (enable raw → candidate GENERATION; the Level-1/3 requirement)

Without these, a cohort can only support **conditional ranking (Level 2)** on the peptides the authors
already chose — it cannot test the generation stage, which is our only demonstrated lever.

| field | requirement | notes |
|---|---|---|
| `patient_pseudo_id` | pseudonymous, stable across all files | e.g. `P01`; the ONLY linkage key; no PHI |
| `genome_build` | exact (GRCh38 / hg19) + reference FASTA name | needed for lossless generation from the raw allele |
| tumor somatic variants | **tumor/normal VCF** (preferred) OR tumor+normal **BAM/FASTQ** | multi-caller callset preferred; keep per-variant caller support |
| per-variant evidence | chrom, pos, ref, alt, **DNA VAF**, read depth | VAF/depth drive confidence, never a hard filter |
| RNA quantification | gene/transcript **counts + TPM** (e.g. RSEM/Salmon) | matched tumor sample; timepoint labelled |
| mutant-allele RNA evidence | RNA reads supporting the mutant allele (or RNA VAF) | distinguishes expressed vs silent mutations |
| HLA genotype | **4-digit or higher** class I (and class II if assayed) | from WES (OptiType/arcasHLA/HLA-HD) or direct typing; per patient |

## Tier 1 — the candidate universe + labels (the denominator and the ground truth)

| field | requirement | notes |
|---|---|---|
| **full preselection candidate universe** | EVERY neoantigen candidate the selection pipeline enumerated, BEFORE down-selection | this is the denominator; without it top-20 is undefined. If only the tested subset exists, mark `preselection_universe = tested_only` |
| `selection_algorithm` + version | the tool/model + version that ranked/chose tested peptides | e.g. pVACtools x.y, NetMHCpan-4.1, proprietary vX; lets us reproduce or bypass their selection |
| tested peptide sequence | exact assayed peptide (and its length/register) | 8–14mer class I; for long peptides record the minimal epitope if deconvolved |
| **peptide↔HLA pairing** | the restricting allele per tested peptide-HLA pair | required for PRIME/EL scoring and per-pair labels |
| source mutation link | tested peptide → originating variant (gene, chrom:pos, ref/alt, transcript) | connects labels back to Tier-0 generation |
| **outcome label** | one of `POSITIVE` / `TESTED_NEGATIVE` / `UNTESTED` | NEVER coerce untested to negative; UNTESTED is not a negative |
| assay | ELISpot / ICS / tetramer-multimer / TIL reactivity / genetic screen | functional recognition, not predicted/MS-presented-only |
| `assay_timepoint` | pre/post-treatment, cycle, or sample date-offset | keep as a column; do not collapse |
| `assay_threshold` | the positivity rule used (e.g. SFC cutoff, tetramer gate) | so we can re-derive labels if thresholds differ |
| effector detail (if any) | CD8/CD4, cytokine, TCR clonotype | optional; helps interpret class I vs II |

## Tier 2 — pooled / barcoded screens (deconvolution provenance)

Required whenever peptides were assayed in POOLS or via DNA-barcoded multimers (Miller IPV, EVX-01,
tetramer baskets). Without deconvolution, a pool "hit" is not a per-peptide label.

| field | requirement | notes |
|---|---|---|
| `pool_id` / `barcode_id` | membership of each tested peptide in each pool/barcode | the tested denominator lives here, not just the deconvolved hits |
| deconvolution mapping | pool/barcode → individual peptide-HLA outcome | how a pool positive was resolved to a peptide |
| pool composition manifest | full list of peptides per pool (incl. those never resolved) | the **true tested denominator**; public TCR/assay files alone are NOT this |

---

## Longitudinal integrity rule (do not collapse)

Preserve **one row per (patient, peptide, HLA, assay, timepoint)**. A peptide POSITIVE at one timepoint and
NEGATIVE at another is TWO rows, not a contradiction to be averaged away. The benchmark decides how to
aggregate downstream; the source package must retain the granularity. Collapsing contradictory longitudinal
outcomes silently manufactures or destroys positives.

## What we return / commit to (state this in every request)

- Only **aggregate, patient-level** benchmark statistics are published (paired Δ hits@20, CIs) — never
  individual-level records; controlled data stays in its controlled environment.
- Split assignment (TRAIN / DEV / LOCKED_TEST) is fixed **before** any label touches the model; locked-test
  cohorts are never used for development.
- Exact and near-peptide leakage controls are applied (patient- and study-level holdout; k-mer near-dup
  screen); the cohort's leakage status is recorded in `COHORT_ACQUISITION_TRACKER.csv`.

## Mapping to the unified benchmark schema (one row per candidate)

`patient_id, mutation_id, mutant_peptide, hla_allele, label∈{POSITIVE,TESTED_NEGATIVE,UNTESTED},
candidate_source∈{pvac,lossless_recovery}, assay, assay_timepoint, prime(%rank), el(EL %rank),
expr(TPM/decile), n_callers, tumor_vaf` (+ optional `cmp_*` comparators). Tier-0 fields regenerate
`prime/el/expr` and the full candidate universe; Tier-1/2 fields supply `label` and provenance.
