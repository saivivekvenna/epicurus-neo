# Miller IPV — pre-registered end-to-end benchmark protocol (frozen BEFORE headline)

Registered 2026-07-12, before any Miller label is seen (S1/S2 are still paywalled). This fixes the
analysis so the headline cannot be chosen after looking at results. Miller IPV is a **LOCKED_TEST** cohort:
it is never used for development or tuning, and no constant is fit to it.

## Hypothesis (the north star)
On the Miller patients, from **identical** WES/RNA/HLA inputs, does the frozen Epicurus pipeline place
**more experimentally recognized neoantigens in the final top-20** than standard pVAC-style generation +
genuine PRIME? **Primary metric: patient-level paired Δ in recognized hits@20** (Wilcoxon signed-rank +
paired bootstrap CI across the ≤13 patients). Reachability (L1) and conditional ranking (L2) are reported
as **diagnostics only**, never as the headline.

## Arms (four-arm harness `src/benchmark/four_arm.py`, frozen)
1. `pvac_prime` — standard pVAC-style candidates + genuine PRIME (incumbent baseline).
2. `lossless_prime` — lossless-generation union + genuine PRIME (protected incumbent; isolates generation).
3. `lossless_epicurus` — lossless union + frozen Epicurus v0.1 (`configs/frozen/epicurus_v0_1.json`).
4. `full_epicurus` — lossless union + Epicurus + route-aware selection.
Expression enters **confidence-only** per the frozen policy (`configs/frozen/expression_policy_v1.json`);
no expression rank penalty. No model retuning.

## Pre-registered exclusions (decided now)
- Peptides with non-standard residues or class-I-invalid length are `NOT_EVALUABLE`, not negatives.
- Patients whose inputs do not reconstruct (failed HLA typing, no somatic calls, or missing RNA) are
  `NOT_EVALUABLE` and excluded from the paired test, reported explicitly.
- Labels are three-state: `POSITIVE` / `TESTED_NEGATIVE` / `UNTESTED`. **UNTESTED is never a negative.**
  Only `POSITIVE` counts as a recognized hit; the denominator for recall is the patient's `POSITIVE` set.
- Contradictory longitudinal rows (same peptide, different timepoint/outcome) are **preserved**; the
  primary analysis uses the pre-registered timepoint (baseline/pre-stimulation ELISpot). A peptide is
  `POSITIVE` for the primary metric iff it is `POSITIVE` at that timepoint. Alternate-timepoint sensitivity
  is reported separately, not as the headline.

## Pre-registered leakage controls (decided now)
- **Split by patient** (and this whole cohort is a locked study-level holdout); never a random peptide split.
- Exact + near-peptide (k-mer, threshold 0.8) de-duplication of Miller peptides against the Epicurus
  training peptides (cd8_multimer) and PRIME's training set; overlaps are **reported**, and Epicurus is
  frozen (trained only on multimer), so no refit is possible regardless.
- Report, but do not filter, whether Miller's 754 assayed peptides appear in IEDB (would contaminate
  NetMHCpan-EL/BigMHC comparators, not genuine PRIME 2.1).

## Candidate universe (the denominator)
The 349 tested variants are **IPV-prefiltered**; the primary run **re-enumerates the full class-I mutanome**
from each patient's WES so both arms share one candidate universe. A secondary "tested-set-only"
conditional-ranking (L2) analysis on the 349/754 is reported as a diagnostic.

## GO / NO-GO (headline gate)
No north-star headline is published unless, for each evaluable patient, the full loop is reproducible:
raw input → candidate generation → the **same** candidate universe → genuine PRIME **and** frozen Epicurus
→ paired hits@20, with `NOT_EVALUABLE` explicit wherever a requirement is missing. If only the tested
peptides can be reconstructed (no full universe), the cohort yields an **L2 diagnostic only**.

## Current status
`RUNNABLE BUT BLOCKED ON FILE`: the ingestion contract, run manifest, download tranches, and this protocol
are done; the run is blocked on the paywalled S1/S2 label table (`EXTERNAL ACTION REQUIRED`). Inputs are
open and the pilot tranche (7.2 GB) is executable now to validate the generation half.
