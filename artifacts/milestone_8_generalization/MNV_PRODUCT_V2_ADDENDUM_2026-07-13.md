# Label-blind MNV product-v2 addendum — 2026-07-13

## Trigger

The pre-incident-pinned Hu_315 freeze failed closed before universe generation
because its PASS VCF contains two equal-length multi-nucleotide substitutions
(MNVs). The label-blind variant-class audit found 624 supported SNVs, 20
supported simple deletions, two unsupported MNVs, and no other unsupported
classes among 646 PASS records. No recognition outcome was read or joined.

The legacy result is permanently retained as `NOT_EVALUABLE` in
`patients/Hu_315/LEGACY_FREEZE_FAILURE.json`. It is not repaired, reclassified,
or presented as a successful execution of the pre-incident lane.

## Product-v2 rule

A separately versioned, patient-agnostic generation rule adds atomic MNV
support:

1. An equal-length genomic block substitution is represented as one HGVS
   `g.<start>_<end>delins<ALT>` event. It is never atomized into independent
   SNVs.
2. Ensembl VEP must accept the event and route it to an already enumerable
   protein-altering consequence. Unsupported consequences remain
   `NOT_ENUMERABLE`.
3. If VEP reports a multi-residue `amino_acids` block, the complete reference
   segment must match the fetched Ensembl protein.
4. The block is replaced atomically and peptide windows are anchored only on
   residue positions whose amino acid actually changes.
5. Empty, unequal-length, reference-mismatching, non-standard, or no-op protein
   substitutions fail closed.

This rule contains no patient ID, gene, coordinate, peptide, outcome, or known
positive. It applies identically to Hu_287, Sid, every Miller calibration
patient, every final-held-out patient, and the identical-raw-input nextNEOpi
comparison.

## Isolation and interpretation

- The change is motivated solely by an unlabeled raw variant class and an
  explicit fail-closed execution error.
- The final-held-out label seal remains intact.
- The pre-incident calibration claim remains failed; product-v2 is a new
  development lane registered before final portfolios and final label access.
- Product-v2 must rerun development/calibration inputs uniformly and must be
  frozen, hashed, and committed before any final-held-out unseal.
- PRIME remains component-level only. nextNEOpi remains the primary end-to-end
  comparator under `END_TO_END_COMPARATOR_ADDENDUM_2026-07-13.md`.

## Registered implementation evidence

- `src/event_b/lossless_peptide_generation.py` SHA-256:
  `d49008b414f7aa69654ade7f1168a18eb57920018d2d3e860559d06cdc1b749c`
- Targeted generator/universe/product suite: 89 tests passed.
- The observed genomic delins form was accepted by Ensembl VEP and routed to a
  missense consequence before implementation; this check used no recognition
  data.
