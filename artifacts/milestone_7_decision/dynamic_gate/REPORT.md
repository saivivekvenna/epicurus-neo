# Dynamic upstream gate — evaluation

`python -m scripts.dynamic_gate` · frozen config `t=0.25` (target 0.95, LOCO-calibrated on gartner, improve, multimer).

Layered safe-rejection gate (Layer0 deterministic impossibility + Layer1 AND-of-core-vetoes on el+prime + Layer2 expression rescue + Layer3 patient-adaptive rails). Survivors are handed UNCHANGED to genuine PRIME / frozen Epicurus v0.1. Label-blind; cohort identity never an input.


## Gate metrics — frozen config

| cohort | pos retention (CP-LB) | worst-patient | neg removal | pts losing a pos | EL-gate matched retention |
|---|--:|--:|--:|--:|--:|
| gartner | 1.000 (0.937) | 1.000 | 0.045 | 0 | 1.000 |
| improve | 0.987 (0.975) | 0.800 | 0.023 | 6 | 0.976 |
| multimer ⚠️in-sample | 0.971 (0.868) | 0.667 | 0.031 | 1 | 1.000 |
| **checkmate153 (LOCKED)** | 0.957 (0.920) | 0.571 | 0.104 | 3 | 0.963 |

_CheckMate expression coverage 0.19 → gate acts on core el+prime there._


## Downstream consequence (gate → UNCHANGED ranker; paired per-patient, gated vs ungated)

| cohort | ranker | ungated hits@20 | gated hits@20 | Δ [CI] | no-regression |
|---|---|--:|--:|--:|:--:|
| gartner | genuine_prime | 0.692 | 0.692 | 0.000 [0.0, 0.0] | yes |
| gartner | frozen_epicurus | 0.808 | 0.808 | 0.000 [0.0, 0.0] | yes |
| improve | genuine_prime | 1.361 | 1.361 | 0.000 [0.0, 0.0] | yes |
| improve | frozen_epicurus | 1.230 | 1.230 | 0.000 [0.0, 0.0] | yes |
| multimer | genuine_prime | 1.105 | 1.105 | 0.000 [0.0, 0.0] | yes |
| multimer | frozen_epicurus | 1.263 | 1.263 | 0.000 [0.0, 0.0] | yes |
| checkmate153 | genuine_prime | 5.000 | 5.000 | 0.000 [0.0, 0.0] | yes |
| checkmate153 | frozen_epicurus | 5.071 | 5.071 | 0.000 [0.0, 0.0] | yes |

## Pareto frontier — LOCO calibration (negative removal achievable at each retention target)

| target | cohort | LOCO t | pos retention (CP-LB) | worst-pt | neg removal | meets CP-LB≥0.95 |
|--:|---|--:|--:|--:|--:|:--:|
| 0.9 | gartner | 0.4 | 1.000 (0.937) | 1.000 | 0.111 | no |
| 0.9 | improve | 0.65 | 0.730 (0.694) | 0.000 | 0.332 | no |
| 0.9 | multimer* | 0.4 | 0.971 (0.868) | 0.667 | 0.093 | no |
| 0.95 | gartner | 0.25 | 1.000 (0.937) | 1.000 | 0.045 | no |
| 0.95 | improve | 0.15 | 0.991 (0.981) | 0.800 | 0.005 | yes |
| 0.95 | multimer* | 0.25 | 0.971 (0.868) | 0.667 | 0.031 | no |
| 0.975 | gartner | 0.2 | 1.000 (0.937) | 1.000 | 0.020 | no |
| 0.975 | improve | 0.0 | 1.000 (0.994) | 1.000 | 0.000 | yes |
| 0.975 | multimer* | 0.25 | 0.971 (0.868) | 0.667 | 0.031 | no |
| 0.99 | gartner | 0.05 | 1.000 (0.937) | 1.000 | 0.000 | no |
| 0.99 | improve | 0.0 | 1.000 (0.994) | 1.000 | 0.000 | yes |
| 0.99 | multimer* | 0.05 | 1.000 (0.916) | 1.000 | 0.001 | no |

_`*` = multimer, frozen-Epicurus in-sample (optimistic)._


## Decisive diagnostic — does the gate EVER improve top-20 at a SAFE retention?

Threshold sweep (rails off). Watch retention and the paired downstream Δhits@20 together: where retention stays ≈1.0 the downstream Δ is ≈0 (the gate removes only sub-top-20 negatives); any positive Δ appears only after retention has already fallen (denominator effect bought by dropping positives — not a safe win).


**gartner** (ungated hits@20: epicurus 0.808, prime 0.692)

| t | neg removal | pos retention (CP-LB) | Δhits@20 epicurus | Δhits@20 prime |
|--:|--:|--:|--:|--:|
| 0.3 | 0.064 | 1.000 (0.937) | 0.000 [0.0, 0.0] | 0.000 [0.0, 0.0] |
| 0.5 | 0.182 | 1.000 (0.937) | 0.000 [0.0, 0.0] | 0.000 [0.0, 0.0] |
| 0.6 | 0.277 | 1.000 (0.937) | -0.038 [-0.1154, 0.0] | 0.000 [0.0, 0.0] |
| 0.7 | 0.410 | 0.957 (0.869) | 0.000 [-0.1154, 0.1154] | 0.038 [0.0, 0.1154] |
| 0.75 | 0.478 | 0.913 (0.812) | 0.038 [-0.0769, 0.1923] | 0.038 [0.0, 0.1154] |
| 0.85 | 0.655 | 0.848 (0.733) | 0.038 [-0.1154, 0.2308] | 0.077 [0.0, 0.1923] |
| 0.95 | 0.856 | 0.435 (0.310) | -0.077 [-0.3846, 0.2308] | 0.038 [-0.2308, 0.3077] |

_Where the oracle lift lives (t=0.6): the gate removes **0.000** of HIGH-presentation negatives (n=927) vs **0.540** of LOW-presentation ones (n=925). The decoys that outrank positives are the top-EL negatives a label-blind presentation gate cannot touch._


**improve** (ungated hits@20: epicurus 1.230, prime 1.361)

| t | neg removal | pos retention (CP-LB) | Δhits@20 epicurus | Δhits@20 prime |
|--:|--:|--:|--:|--:|
| 0.3 | 0.041 | 0.962 (0.943) | 0.000 [0.0, 0.0] | 0.000 [0.0, 0.0] |
| 0.5 | 0.161 | 0.880 (0.853) | 0.016 [-0.0328, 0.082] | 0.000 [0.0, 0.0] |
| 0.6 | 0.270 | 0.797 (0.763) | 0.033 [-0.0492, 0.1148] | 0.000 [0.0, 0.0] |
| 0.7 | 0.402 | 0.692 (0.655) | 0.033 [-0.0656, 0.1311] | 0.000 [0.0, 0.0] |
| 0.75 | 0.475 | 0.625 (0.587) | 0.098 [0.0, 0.1967] | 0.016 [0.0, 0.0492] |
| 0.85 | 0.652 | 0.452 (0.413) | 0.147 [0.0164, 0.2787] | 0.016 [0.0, 0.0492] |
| 0.95 | 0.869 | 0.214 (0.183) | -0.049 [-0.2459, 0.1475] | -0.147 [-0.3443, 0.0492] |

_Where the oracle lift lives (t=0.6): the gate removes **0.000** of HIGH-presentation negatives (n=3837) vs **0.485** of LOW-presentation ones (n=3819). The decoys that outrank positives are the top-EL negatives a label-blind presentation gate cannot touch._


**multimer** (ungated hits@20: epicurus 1.263, prime 1.105)

| t | neg removal | pos retention (CP-LB) | Δhits@20 epicurus | Δhits@20 prime |
|--:|--:|--:|--:|--:|
| 0.3 | 0.047 | 0.971 (0.868) | 0.000 [0.0, 0.0] | 0.000 [0.0, 0.0] |
| 0.5 | 0.163 | 0.971 (0.868) | 0.000 [0.0, 0.0] | 0.000 [0.0, 0.0] |
| 0.6 | 0.262 | 0.971 (0.868) | 0.000 [0.0, 0.0] | 0.000 [0.0, 0.0] |
| 0.7 | 0.394 | 0.882 (0.751) | 0.000 [0.0, 0.0] | 0.000 [0.0, 0.0] |
| 0.75 | 0.471 | 0.853 (0.715) | 0.000 [0.0, 0.0] | 0.000 [0.0, 0.0] |
| 0.85 | 0.659 | 0.824 (0.681) | -0.053 [-0.1579, 0.0] | 0.000 [0.0, 0.0] |
| 0.95 | 0.882 | 0.706 (0.552) | -0.158 [-0.3158, 0.0] | 0.000 [-0.2105, 0.2632] |

_Where the oracle lift lives (t=0.6): the gate removes **0.000** of HIGH-presentation negatives (n=1766) vs **0.478** of LOW-presentation ones (n=1770). The decoys that outrank positives are the top-EL negatives a label-blind presentation gate cannot touch._


**checkmate153** (ungated hits@20: epicurus 5.071, prime 5.000)

| t | neg removal | pos retention (CP-LB) | Δhits@20 epicurus | Δhits@20 prime |
|--:|--:|--:|--:|--:|
| 0.3 | 0.147 | 0.957 (0.920) | 0.000 [0.0, 0.0] | 0.000 [0.0, 0.0] |
| 0.5 | 0.336 | 0.870 (0.819) | 0.000 [-0.2143, 0.2143] | -0.071 [-0.2143, 0.0] |
| 0.6 | 0.444 | 0.815 (0.757) | -0.071 [-0.2857, 0.1429] | -0.071 [-0.2143, 0.0] |
| 0.7 | 0.569 | 0.741 (0.678) | -0.143 [-0.4286, 0.1429] | -0.071 [-0.2143, 0.0] |
| 0.75 | 0.642 | 0.698 (0.633) | -0.071 [-0.3571, 0.2143] | 0.071 [-0.2143, 0.4286] |
| 0.85 | 0.778 | 0.488 (0.421) | -0.286 [-0.7143, 0.1429] | -0.214 [-0.7143, 0.2143] |
| 0.95 | 0.918 | 0.235 (0.181) | -2.357 [-2.9286, -1.7857] | -2.286 [-3.1429, -1.4286] |

_Where the oracle lift lives (t=0.6): the gate removes **0.000** of HIGH-presentation negatives (n=227) vs **0.803** of LOW-presentation ones (n=279). The decoys that outrank positives are the top-EL negatives a label-blind presentation gate cannot touch._


## Verdict
**B for the gate as a safe recall-preserving pruner. C is SCOPED: falsified only for same-feature presentation/PRIME-derived monotone gating with current peptide features. The general orthogonal-feature dynamic-gate hypothesis is UNTESTED / data-blocked — NOT falsified.**

The layered gate is real, safe, and Pareto-dominates the incumbent EL-percentile gate on positive retention at matched removal (FEASIBILITY.md). Its honest value is **shrinking the candidate universe handed to expensive downstream steps (genuine-PRIME scoring, wet-lab) without losing recognized positives** — NOT lifting top-20.

**Circularity caveat (CIRCULARITY_AUDIT.md).** The v1 gate's veto axes {EL, PRIME} ARE the dominant inputs to the downstream rankers (genuine PRIME uses PRIME; frozen Epicurus is PRIME-dominated over prime/el/expr). So the gate removes exactly the candidates those same rankers already bury: 0 gate-removed candidates ever sat in a ranker top-20, and it removes **0%** of the high-EL/high-PRIME decoy stratum by construction. The zero downstream Δhits@20 at safe retention is therefore a **structural tautology of feature overlap**, and proves only that a same-feature monotone presentation gate cannot move a same-feature top-20 — it does NOT prove a label-blind gate in general cannot help. The aggressive-gating premise is falsified for this feature class only.

CP retention lower bounds are sample-size-capped (Gartner 46 positives ⇒ max CP-LB 0.937), so small cohorts cannot certify ≥0.95 regardless of the gate — an underpowering limit, not an unsafety.

**What could close the gap (v2, data-blocked; V2_CONTRACT.md / V2_PREREGISTRATION.md):** a gate whose veto uses features ORTHOGONAL to the downstream rank — mutant-allele RNA VAF + read support, tumor DNA VAF/depth/CCF/clonality, proteasomal processing/stability, agretopicity/WT-differential, transcript confidence — as cross-fitted residual signals after conditioning on PRIME/EL (which define only the hard-decoy stratum + rescue/floor, never veto axes). Missing⇒KEEP. **Constraints learned:** (i) a SEQUENCE-ONLY residual gate is FALSIFIED — an independent train-Gartner→test-IMPROVE run retained 1.5% of positives (peptide/HLA motifs encode severe study/assay shift); (ii) leave-one-STUDY-out calibration + explicit OOD abstention are mandatory. An in-sample multimer probe shows only a faint orthogonal signal (AUROC ~0.59), unvalidatable cross-study here. Concrete path = open WES+RNA of Miller IPV (PRJNA980652) + Gartner reconstruction.


> DEV features are peptide-/presentation-only (prime/el/expr). No raw WES depth, mutant-allele RNA VAF, purity/CCF, or processing available cross-cohort -> this gate does NOT learn genomics.
> multimer is frozen Epicurus' training cohort -> in-sample; excluded from the safety headline.
> CheckMate expr coverage is sparse; the gate operates on the core el+prime axes there.
> Oracle retention (100%) is a ceiling, never a validation.
> Downstream hits@20/recall@20 rise partly because gating shrinks the top-20 denominator; the no-regression check (gated >= ungated) is the honest safety statement, not the raw lift.
