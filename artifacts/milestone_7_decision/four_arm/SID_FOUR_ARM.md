# osteosarc.com / Sid — the one END-TO-END (Level-3) evaluable patient

> **Status: post_hoc_diagnostic_n3_single_patient_not_blinded_not_powered.** Post-hoc, n=3 measured positives, single patient. Strict label isolation: the 3 recognized mutations are joined only AFTER each arm's ranking. NOT a blinded, powered, or prospective superiority test — a reachability/attribution DIAGNOSTIC.

- k = 20; evaluation-only positives: ASPM, DYNC1H1, MAP2
- available inputs: epicurus_features, genuine_prime, lossless_generation, measured_labels, pvac_candidates, router_features

This single patient instantiates all three benchmark levels; each is read separately:

- **L1 reachability** — how many recognized mutations survive raw→generation.
  Frozen: generation recall pVAC 1/3 -> lossless union 3/3 (recovered 2 recognized mutation(s) at the generation stage).
- **L2 conditional ranking** — ordering among the generated/rankable candidates (within this patient's denominator only); see the per-arm hits@20 below.
- **L3 end-to-end patient utility (PRIMARY)** — recognized mutations in the final top-20 from common raw inputs vs standard pVAC + genuine PRIME; see `total` in the stage attribution. `lossless_prime` (lossless generation + genuine PRIME) is the protected incumbent.

## Primary — frozen Epicurus v0.1

> EL feature: NetMHCpan-EL MT %rank (frozen Epicurus feature); MISSING on lossless-recovered candidates -> NaN -> 0.5 percentile (frozen policy)

| arm | gen recall | rankable recall | hits@20 | recall@20 | covered |
|---|---|---|---:|---:|---|
| `pvac_prime` | 1/3 | 1/3 | 1 | 0.3333 | DYNC1H1 |
| `lossless_prime` | 3/3 | 3/3 | 3 | 1.0 | ASPM, DYNC1H1, MAP2 |
| `lossless_epicurus` | 3/3 | 3/3 | 1 | 0.3333 | DYNC1H1 |
| `full_epicurus` | 3/3 | 3/3 | 1 | 0.3333 | DYNC1H1 |

Stage attribution: generation **+2** · scorer **-2** · selection **+0** · total **+0** (top-20 hits)

## Sensitivity — MixMHCpred EL populated on recovered candidates (NON-FROZEN)

> MixMHCpred %rank for ALL candidates (present on recovered rows too). NON-FROZEN sensitivity: isolates the missing-EL-feature confound on recovered rows.

| arm | gen recall | rankable recall | hits@20 | recall@20 | covered |
|---|---|---|---:|---:|---|
| `pvac_prime` | 1/3 | 1/3 | 1 | 0.3333 | DYNC1H1 |
| `lossless_prime` | 3/3 | 3/3 | 3 | 1.0 | ASPM, DYNC1H1, MAP2 |
| `lossless_epicurus` | 3/3 | 3/3 | 2 | 0.6667 | ASPM, DYNC1H1 |
| `full_epicurus` | 3/3 | 3/3 | 2 | 0.6667 | ASPM, DYNC1H1 |

Stage attribution: generation **+2** · scorer **-1** · selection **+0** · total **+1** (top-20 hits)

## Leakage controls

- Frozen Epicurus training cohort: `cd8_multimer`; Sid out-of-sample: True
- Candidate peptides: 1471; near/exact PRIME-training overlap: 1 (0.0007)
- Label isolation: measured positives joined only AFTER each arm's ranking (harness-enforced)

## Interpretation

Generation recovers all 3 recognized mutations (recall 1/3 -> 3/3; +2 top-20 hits under genuine PRIME). Under the FROZEN Epicurus scorer the recovered ASPM+MAP2 fall back out of the top-20 (scorer stage -2, net 0) — but this is partly a missing-EL-feature artifact: the sensitivity arm (EL populated for recovered rows) recovers ASPM (scorer -1, net +1). Low-expression MAP2 still drops under Epicurus (a real expression-reweighting effect, consistent with prior cohorts where a learned recognition score on top of presentation hurts). Net: the RELIABLE, reproducible win here is candidate GENERATION feeding genuine PRIME; the Epicurus scorer neither clearly helps nor is fairly testable until presentation features are computed on recovered candidates. n=3, post-hoc — not a gate.

