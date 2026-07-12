# osteosarc.com / Sid — the one END-TO-END (Level-3) evaluable patient

> **Status: post_hoc_diagnostic_n3_single_patient_not_blinded_not_powered.** Post-hoc, n=3 measured positives, single patient. Strict label isolation: the 3 recognized mutations are joined only AFTER each arm's ranking. NOT a blinded, powered, or prospective superiority test — a reachability/attribution DIAGNOSTIC.

- k = 20; evaluation-only positives: ASPM, DYNC1H1, MAP2
- available inputs: epicurus_features, genuine_prime, lossless_generation, measured_labels, pvac_candidates, router_features

This single patient instantiates all three benchmark levels; each is read separately:

- **L1 reachability** — how many recognized mutations survive raw→generation.
  generation recall pVAC 1/3 -> lossless union 3/3 (recovered 2 recognized mutation(s) at the generation stage).
- **L2 conditional ranking** — ordering among the generated/rankable candidates (within this patient's denominator only); see the per-arm hits@20 below.
- **L3 end-to-end patient utility (PRIMARY)** — recognized mutations in the final top-20 from common raw inputs vs standard pVAC + genuine PRIME; see `total` in the stage attribution. `lossless_prime` (lossless generation + genuine PRIME) is the protected incumbent.

> Epicurus feature provenance: prime = genuine PRIME %rank; el = presentation %rank (see each variant); expr = RSEM gene TPM. Frozen formula prime+el+expr only. NetMHCpan runnable locally: False; MHCflurry vs NetMHCpan-EL Spearman on 5710 pVAC rows = 0.518.

## Primary (FAIR) — frozen Epicurus, genuine MHCflurry EL on all candidates

> EL feature: GENUINE MHCflurry presentation %rank for ALL candidates (independent learned predictor; recovered candidates get real presentation evidence, no 0.5 impute). Fair four-arm attribution.

| arm | gen recall | rankable recall | hits@20 | recall@20 | covered |
|---|---|---|---:|---:|---|
| `pvac_prime` | 1/3 | 1/3 | 1 | 0.3333 | DYNC1H1 |
| `lossless_prime` | 3/3 | 3/3 | 3 | 1.0 | ASPM, DYNC1H1, MAP2 |
| `lossless_epicurus` | 3/3 | 3/3 | 2 | 0.6667 | ASPM, DYNC1H1 |
| `full_epicurus` | 3/3 | 3/3 | 2 | 0.6667 | ASPM, DYNC1H1 |

Stage attribution: generation **+2** · scorer **-1** · selection **+0** · total **+1** (top-20 hits)

## Reference — literal frozen NetMHCpan-EL (recovered candidates imputed to 0.5)

> EL feature: Literal frozen NetMHCpan-EL for pVAC candidates; MISSING on recovered -> 0.5 impute (shows the imputation artifact the fair run removes).

| arm | gen recall | rankable recall | hits@20 | recall@20 | covered |
|---|---|---|---:|---:|---|
| `pvac_prime` | 1/3 | 1/3 | 1 | 0.3333 | DYNC1H1 |
| `lossless_prime` | 3/3 | 3/3 | 3 | 1.0 | ASPM, DYNC1H1, MAP2 |
| `lossless_epicurus` | 3/3 | 3/3 | 1 | 0.3333 | DYNC1H1 |
| `full_epicurus` | 3/3 | 3/3 | 1 | 0.3333 | DYNC1H1 |

Stage attribution: generation **+2** · scorer **-2** · selection **+0** · total **+0** (top-20 hits)

## Sensitivity — MixMHCpred EL (PRIME backbone) on all candidates

> MixMHCpred %rank (PRIME's backbone) for all candidates — secondary sanity.

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

L1 reachability: generation recovers all 3 recognized mutations (recall 1/3 -> 3/3; +2 top-20 hits under genuine PRIME = the protected lossless_prime incumbent). L3 end-to-end: with GENUINE presentation features computed on recovered candidates (MHCflurry, no impute), the frozen Epicurus scorer stage is -1 and the full stack nets 1 vs pVAC+PRIME. The earlier -2 frozen scorer loss was substantially a 0.5-impute artifact on recovered rows (reference vs primary). Any residual drop is the Epicurus expression/EL reweighting demoting a low-expression true positive — consistent with prior cohorts where a learned recognition score on top of presentation does not help. NetMHCpan is not locally runnable, so el uses MHCflurry (independent predictor); its moderate agreement with NetMHCpan-EL is disclosed in feature_provenance. n=3, post-hoc, descriptive — NOT a gate, no constant tuned to Sid.

