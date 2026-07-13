# Canonical Epicurus product-path end-to-end audit

**Locked after prior Hu_287 and Sid outcomes were known.** This is a reproducible
product integration audit, not prospective validation and not a model-selection
experiment. No threshold, weight, or portfolio constraint is changed here.

## Deliverable being tested

The only headline arm is the actual Epicurus-owned production path:

```text
raw patient WES/RNA evidence
  -> somatic variants + patient HLA + expression/RNA support
  -> label-blind lossless peptide×HLA generation
  -> production input normalization
  -> production deterministic validity gate
  -> production translated/presented/recognized/coverage evidence score
  -> production eligibility policy
  -> production capped 20-slot portfolio
```

The genomics tools are replaceable upstream adapters, consistent with the public
Epicurus product boundary. This audit nevertheless requires evidence that the
candidate table was reconstructed from that patient's WES/RNA/HLA; a standalone
assay peptide table is not end-to-end eligible.

## Frozen product configuration

- `InferenceConfig()` defaults;
- `k=20`;
- `max_per_mutation=2`;
- `max_per_gene=4`;
- no HLA cap;
- deterministic validity gate enabled;
- core threshold 0.55;
- supporting threshold 0.35;
- measured zero expression and zero mutant RNA reads at depth >=10 retain the
  currently shipped exclusion semantics for this audit. They may be criticized,
  but cannot be silently removed to improve the result.

## Patients

- **Hu_287:** complete local raw WES/RNA/HLA reconstruction and lossless class-I
  universe; discovery patient, not independent.
- **Sid:** raw longitudinal genomic/RNA reconstruction with patient HLA and
  label-blind generation; 137/147 eligible mutations are generated and the
  remaining 10 have documented non-enumerable consequences; repeatedly
  inspected, not blind.

No other local patient currently has both a reconstructed mutation denominator
and a runnable patient-level candidate universe. Peptide-table-only cohorts are
excluded from the headline.

## Primary endpoint and funnel

- unique experimentally recognized mutations represented in the actual product
  top 20;
- generated -> deterministic-valid -> product-eligible -> selected counts;
- recognized-mutation reachability at every stage after selections are frozen;
- exact exclusion reason for every lost recognized mutation.

## Comparators

Genuine PRIME on the same generated universe is a diagnostic comparator, not the
deliverable. Report ordinary PRIME top 20 and PRIME with the same mutation cap so
generation quality and set selection are not conflated with product behavior.

## Label barrier

Every product and comparator selection, stage membership set, input hash, policy
hash, and code hash is written to `FROZEN_PIPELINE.json` before either patient's
recognition labels are loaded. Results are then joined mechanically.

## Interpretation

The product can be called runnable end to end if raw-derived inputs traverse all
stages and produce a valid portfolio. It can be called successful on a patient
only by its actual final hits@20. Two known patients cannot establish general
superiority. Any disagreement between this product result and earlier research
arms invalidates using the research-arm number as the shipped-product claim.
