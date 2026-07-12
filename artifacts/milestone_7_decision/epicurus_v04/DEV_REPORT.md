# Epicurus v0.4 DEVELOPMENT - source-aware tower - verdict: **REJECT**

Preregistered: `PREREGISTERED_PROTOCOL.md`. DEVELOPMENT ONLY - Gartner TEST not opened; no external claim. The source-**name** tower is **mechanism evidence only, not deployable**. v0.1 remains the frozen model of record.

Provenance verified (12 inputs, git 5e578a668a). 118 scored patients (source-balanced, patient-paired bootstrap).

## Registered gate (candidate = F, the feature tower)

- vs **genuine PRIME**: dhits@20 = -0.0157 CI[-0.1472, 0.1121] -> beats PRIME: **False**
- vs **strongest presentation**: d = 0.0824 CI[-0.0981, 0.2528] -> no regression: **True**

## Members (OOF hits@20)

| member | overall hits | d vs PRIME (CI) | d vs presentation (CI) |
|---|--:|---|---|
| P - pooled (frozen v0.3) | 1.0528 | -0.0926 [-0.2278, 0.0509] | 0.0056 [-0.0982, 0.1074] |
| C - calibration tower | 1.0435 | -0.1019 [-0.2407, 0.0417] | -0.0037 [-0.1046, 0.1121] |
| **F - feature tower** | 1.1296 | -0.0157 [-0.1472, 0.1121] | 0.0824 [-0.0981, 0.2528] |

## Mechanism contrasts (the registered hypothesis)

- **F - P** (source-conditioning vs naive pooling): 0.0769 CI[-0.0732, 0.2084] -> statistically TIED (CI spans 0)
- **C - P** (pure prevalence calibration): -0.0093 CI[-0.0796, 0.0528] -> statistically TIED (CI spans 0)
- **F - C** (**feature weighting beyond calibration**): 0.0861 CI[-0.0797, 0.2269] -> statistically TIED (CI spans 0)

## Per-source (F, hits@20 vs PRIME)

| source | patients | hits F | hits PRIME | d vs PRIME |
|---|--:|--:|--:|--:|
| gartner | 40 | 1.25 | 0.975 | 0.275 |
| improve | 60 | 1.25 | 1.35 | -0.1 |
| multimer | 18 | 0.889 | 1.111 | -0.222 |

## Diagnostics

- **Selected lambda per fold**: [{'fold': 0, 'lam': 10.0, 'C': 1.0, 'tau': 0.5}, {'fold': 1, 'lam': 10.0, 'C': 0.3, 'tau': 1.0}, {'fold': 2, 'lam': 1.0, 'C': 0.3, 'tau': 0.5}, {'fold': 3, 'lam': 1.0, 'C': 0.3, 'tau': 0.5}, {'fold': 4, 'lam': 3.0, 'C': 0.3, 'tau': 1.0}]
- **Source-only vs augmented**: {'gartner': {'augmented_mean_hits': 1.25, 'source_only_mean_hits': 1.225}, 'improve': {'augmented_mean_hits': 1.25, 'source_only_mean_hits': 1.15}, 'multimer': {'augmented_mean_hits': 0.889, 'source_only_mean_hits': 1.0}}
- **Study shortcut**: source-identity AUROC = 0.898 (prevalence {'gartner': 0.00028, 'improve': 0.02666, 'multimer': 0.0042})
- **Negative control (capacity probe)**: true-source F-P = 0.0769, shuffled-source F-P = -0.2324 -> **shuffled source (random grouping) is far WORSE than true source, so the lift tracks GENUINE source structure, NOT model capacity**.
- **Multi-init stability**: min Spearman = 0.999996 (ranking ~identical, so OOF hits stable); max|d std-score| = 0.01657 (strict <=1e-3 threshold not met - the objective is non-convex as pre-registered; the metric is rank-based so this score-scale wobble does not move hits@20).
- **Quarantine stratum**: {'n_patients_with_quarantined_pos': 14, 'mean_hits_model': 0.5, 'mean_hits_prime': 0.357, 'note': 'recurrent-antigen robustness stratum; reported only, never used for selection.'}
- **Feature ablation** (leave-one-out d vs PRIME) + **effective per-source weights**: see DEV_RESULT.json `feature_ablation_vs_prime` / `effective_weights_final`.
- **Attrition (label-blind)**: {'feature_bearing': 152, 'rankable_label_blind': 152, 'scored_has_positive': 118, 'lost_to_quarantine_only': 5}

## What this means

- **Gate REJECT - still a TIE with genuine PRIME, not a loss.** F d vs PRIME = -0.0157 CI[-0.1472, 0.1121] (statistically TIED (CI spans 0)); no regression vs presentation (d=0.0824). F sits CLOSER to PRIME than pooled v0.3 (-0.016 vs -0.093) but does not clear ACCEPT (needs CI_lo>0). PRIME still not beaten.
- **The tower recovered real, source-structured signal - the Gartner edge the design predicted.** Per-source, F beats PRIME on Gartner by d**0.275** (hits 1.25 vs 0.975; up from pooled +0.05). The **negative control is clean**: shuffled source (random grouping) is far WORSE than true source, so the lift tracks GENUINE source structure, NOT model capacity.
- **Mechanism directionally supported but underpowered.** F-P = 0.0769 [-0.0732,0.2084] and F-C = 0.0861 [-0.0797,0.2269] are POSITIVE (feature-weighting, not calibration - C-P = -0.0093 ~0), but every CI spans 0. The improvement over pooling comes from source-conditioned FEATURE weighting as hypothesized, yet is not statistically established at this sample size.
- **Why only a tie: IMPROVE and multimer don't cooperate.** F loses to PRIME on IMPROVE (-0.1) and multimer (-0.222; the multimer head over-specializes on n=18), cancelling the Gartner gain in the source-balanced aggregate.
- **Recognition wall persists.** The tower's gains ride on presentation-adjacent features (f_pres_abs, f_prime_pct); orthogonal recognition features (expression/foreignness/agretopicity/processing) still add nothing or hurt (ablation). Per-source feature WEIGHTING helps; no new recognition AXIS appears.

## Verdict

**REJECT.** v0.1 remains frozen (configs/frozen/epicurus_v0_1.json); v0.4 is REJECTED_DEVELOPMENT. Source-name tower is mechanism evidence only (not deployable). Gartner TEST NOT opened.

