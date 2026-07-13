# Pre-calibration selection arms

**Registered before opening any additional Miller outcomes.** These are controls,
not six opportunities to cherry-pick a final claim. All arms operate on the same
patient-specific, label-blind, lossless candidate universe.

## Why the Hu_287 result needs a harder test

Hu_287 has 14 eligible mutations. With a one-route-per-mutation cap, any ranker
can expose all 14 in a 20-slot portfolio; consequently both genuine PRIME and
Epicurus recover all three measured-recognized mutations. This remains a valid
end-to-end reachability result, but it does not identify the better recognition
score. Sid has 137 eligible mutations and therefore separates the problems:
cap-one Epicurus recovers 1/3 while cap-one PRIME recovers 2/3.

The additional Miller patients must distinguish three possible sources of value:

1. lossless generation and safe validity handling;
2. mutation-level portfolio coverage rather than duplicate peptide/HLA routes;
3. recognition scoring beyond genuine PRIME.

## Mandatory controls

- `prime_plain`: genuine PRIME route-level top 20.
- `prime_mutation_cap1`: genuine PRIME, at most one route per mutation.
- `epicurus_plain`: the shipped product's frozen
  `epicurus_lower_evidence_score`, route-level top 20.
- `epicurus_mutation_cap1`: the same shipped product score, at most one route
  per mutation.

This score definition was clarified before any additional Miller outcomes were
opened. The legacy research-only `epicurus`/`frozen_epicurus_score` column remains
available for provenance and diagnostics, but is not an Epicurus product arm and
does not participate in fusion or evidence-lane selection.

## Pre-registered candidate policies

- `rank_fusion_cap1`: mutation-level, within-patient percentile-rank fusion of
  genuine PRIME and the shipped Epicurus product score; missing evidence
  contributes no advantage.
- `evidence_lane_portfolio`: mutation-level round-robin union of independent
  PRIME, the shipped Epicurus product score, and presentation-evidence lanes,
  deduplicated by mutation.
  A lane that is unavailable for a patient is skipped and its slots are shared
  by the remaining lanes. This is dynamic with respect to evidence availability,
  never patient identity.

Both candidates use one selected route per mutation until every eligible
mutation has one route; only then may a second route be admitted. The best route
for a mutation is selected by the arm's own score, so route choice and mutation
choice are not conflated.

## Calibration decision

After all six calibration portfolios are frozen, outcomes may be joined once.
The universal final policy is selected by the locked lexicographic objective in
`PROTOCOL.md`. Candidate-arm comparison and any one permitted policy revision are
development only. The six final patients receive exactly one frozen policy and
one final outcome join.
