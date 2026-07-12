# Miller IPV (PRJNA980652) — selective download plan

Derived deterministically from the **public** SRA runinfo (`SRA_RUNINFO.csv`, 39 runs) via
`benchmark.miller_ingest`. **Do not bulk-download the ~0.22 TB blindly** — the label table (S1/S2) that
makes the cohort usable is still blocked (see status), so the first tranche exists only to prove the
raw→generation→PRIME/Epicurus loop once labels arrive.

## Cohort inputs (verified public)
- **13 patients**, sample ids `Hu_<NNN>`; every patient has a **complete trio**: normal exome (WXS) +
  tumor exome (WXS) + tumor RNA-seq. 26 WXS + 13 RNA-Seq = 39 runs, 26 BioSamples.
- **Total ~549.7 Gbases / ~0.217 TB** SRA archive. Per-patient crosswalk + sizes: `INPUT_CROSSWALK.csv`.
- Access: **OPEN** (no dbGaP/DAR). Cloud mirrors on S3/GS.

## Tranche policy (smallest scientifically valid first)
| tranche | scope | runs | size |
|---|---|---|---|
| **T1 (pilot)** | **1 patient trio** (smallest = `Hu_287`) | 3 (normal exome + tumor exome + tumor RNA) | **~7.2 GB** |
| T2 | +2 more small patients | 6 | ~23 GB cumulative |
| T3…​ | remaining patients in ascending size | … | up to ~0.217 TB total |

**Rule:** download **T1 only** first, run the full loop end-to-end on one patient (generation → genuine
PRIME + frozen Epicurus → per-patient recognized hits@20) to prove reproducibility, THEN scale. Never pull
the full 0.217 TB before T1 validates and the label table is in hand.

## Per-patient pipeline once a trio lands (destination `data/raw/miller_ipv/`, gitignored)
1. `prefetch`/`fasterq-dump` the 3 runs (or stream from the S3/GS mirror).
2. HLA-type from the normal+tumor exome (OptiType/arcasHLA) → 4-digit class I (no published HLA table).
3. Somatic calling tumor-vs-normal exome → per-mutation VAF/depth.
4. Quantify tumor RNA-seq → TPM + mutant-allele RNA evidence.
5. Re-enumerate the **full** class-I 8–11mer mutanome (lossless/pVACtools) = the denominator (the paper's
   349 tested is IPV-prefiltered; both arms must share the re-enumerated universe).
6. Score genuine PRIME + frozen Epicurus features; join S1/S2 labels **after** ranking.

## Rough runtime/storage (order-of-magnitude, one patient trio)
- Storage: ~7 GB compressed SRA + ~20–40 GB working (BAMs, intermediates).
- Compute: HLA typing minutes; alignment+calling ~1–3 h/patient on a workstation; RNA quant ~20–40 min;
  mutanome enumeration + PRIME/MHCflurry scoring minutes. Full 13-patient cohort ≈ a day of wall-clock on a
  single machine, dominated by alignment.

## Status
- Inputs: **DONE (verified, open)** — this plan is executable now.
- The download itself is deferred: **RUNNABLE BUT BLOCKED ON FILE** (S1/S2 labels), so a T1 pull would
  only exercise the generation half. Pull T1 when you decide to validate the pipeline mechanics, or wait
  for the label table to run the full north-star loop.
