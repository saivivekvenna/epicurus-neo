# Frozen-policy portfolio generalization stress test

**Locked 2026-07-12 after the Hu_287 result was known.** This is a post-result
generalization audit, not a preregistered prospective validation. No parameter is
fit in this experiment. The primary policy is the already-frozen
`configs/frozen/evidence_router_v1.json` policy that produced the Hu_287 result:
`k=20`, `max_per_mutation=2`, `max_per_gene=4`, no HLA cap, at most three
one-slot non-CORE route reserves, and deterministic
`md5(mutant_peptide|hla_allele)` tie-breaking.

## Question

Does the frozen route-aware portfolio improve **unique recognized mutations in
20 slots** beyond ordinary scalar top-20 selection, and is any improvement due to
the selector itself or specifically to the Epicurus score?

## Primary arms

All arms use the same router-valid/rankable candidate rows within a patient.

1. `prime_plain`: genuine PRIME score, ordinary top 20.
2. `prime_route_aware`: genuine PRIME score, frozen Epicurus route-aware selector.
3. `epicurus_plain`: frozen Epicurus v0.1 score, ordinary top 20.
4. `epicurus_route_aware`: frozen Epicurus v0.1 score, frozen route-aware selector.

This crossed 2×2 comparison prevents attributing a generic diversification gain
to the Epicurus recognition score.

## Eligible local patients

- **Hu_287:** mutation-resolved lossless class-I universe, genuine PRIME,
  Epicurus v0.1, patient HLA, and three IFN-g-recognized mutations. This is the
  discovery/replay patient; it is not independent evidence.
- **Sid/osteosarcoma:** mutation-resolved lossless class-I universe, genuine
  PRIME, Epicurus v0.1, patient HLA, and three recognized mutations. This patient
  has been inspected repeatedly and is a post-hoc stress test, not a blind test.

Gartner, IMPROVE, multimer, Zhao, CEDAR, and the Event-B label backbone are
excluded from the primary mechanism test because their locally scored tables do
not preserve a reliable underlying mutation identity across multiple generated
peptide×HLA routes. Gene or peptide identity will not be substituted for mutation
identity.

## Endpoints

Primary endpoint per arm and patient:

- number of unique recognized mutations represented in the selected 20 slots.

Required diagnostics:

- selected slot count and saturation;
- number of unique selected mutations;
- duplicate-slot burden (`selected slots - unique selected mutations`);
- recognized mutation identities;
- paired selector delta within each score (`route_aware - plain`);
- paired score delta within each selection mode (`Epicurus - PRIME`).

No pooled significance claim is permitted for two post-hoc patients. Aggregate
totals are descriptive only.

## Sensitivity analyses

- Repeat the crossed comparison at `k=10` and `k=30` using the same cap values.
- Descriptively test mutation caps 1, 2, 3, 5, and no cap at `k=20`. These are
  explicitly post-hoc sensitivity analyses and cannot replace the primary cap-2
  result.
- Compare the full route-aware selector with a cap-only selector to distinguish
  diversity from route-reserve effects.

## Label barrier and reproducibility

Each arm's selected mutation IDs must be serialized in a freeze payload before
recognition labels are joined for scoring. Input file SHA-256 hashes, frozen
policy hash, code hash, and exact source columns are recorded. Missing inputs or
unreliable mutation identity produce `NOT_EVALUABLE`, never an inferred result.

## Interpretation guardrail

- A Hu_287 win alone is replication of the discovery result.
- Improvement on Sid is supportive but not independent validation.
- Failure or regression on Sid falsifies broad claims that diversification is
  already proven to generalize.
- If `prime_route_aware` matches `epicurus_route_aware`, the supported invention
  is a scorer-agnostic vaccine-portfolio layer, not a superior recognition model.
- A publishable superiority claim still requires additional mutation-resolved,
  multi-patient, frozen end-to-end reconstructions.
