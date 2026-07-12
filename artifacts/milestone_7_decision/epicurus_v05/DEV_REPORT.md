# Epicurus v0.5 DEVELOPMENT — context-conditioned pairwise challenger — verdict: **REJECT**

Preregistered: `PREREGISTERED_PROTOCOL.md`. DEVELOPMENT ONLY — Gartner TEST not opened; no external claim. The ONLY gate is R vs GENUINE PRIME (raw unmasked prime_rank).

Provenance verified (18 inputs, git 03488b896e). 118 scored patients (source-balanced, patient-paired bootstrap).

## Frozen comparator reproduction (§2.1)

- P (convex): {'max_abs_coef_diff': 4.9e-05, 'overall_hits_refit': 1.0528, 'overall_hits_frozen': 1.0528, 'abs_hits_diff': 0.0, 'tol_coef': 0.002, 'tol_hits': 1e-06, 'reproduced': True}
- A (convex): {'max_abs_coef_diff': 4.9e-05, 'overall_hits_refit': 0.9343, 'overall_hits_frozen': 0.9343, 'abs_hits_diff': 0.0, 'tol_coef': 0.002, 'tol_hits': 1e-06, 'reproduced': True}
- F (nonconvex, honest tolerance): {'overall_hits_refit': 1.1296, 'overall_hits_frozen': 1.1296, 'abs_hits_diff': 0.0, 'tol_hits': 0.005, 'nonconvex': True, 'reproduced': True, 'note': 'F is nonconvex (v0.4 multi-init wobble); hits verified to v0.4 tolerance, residual reported.'}

## Registered gate (candidate = R)

- vs **genuine PRIME**: Δhits@20 = -0.012 CI[-0.125, 0.1084] → beats PRIME: **False** (statistically TIED (CI spans 0))
- vs **strongest presentation** (mix): Δ = 0.0324 CI[-0.0722, 0.125] → no regression: **True**

## Members (OOF hits@20)

| member | overall hits | Δ vs PRIME (CI) |
|---|--:|---|
| P — pooled (frozen v0.3) | 1.0528 | -0.0926 [-0.2278, 0.0509] |
| Q — shared pairwise | 1.1222 | -0.0231 [-0.1472, 0.1111] |
| **R — context pairwise** | 1.1333 | -0.012 [-0.125, 0.1084] |
| F — source-name tower (frozen v0.4) | 1.1296 | -0.0157 [-0.1472, 0.1121] |

## Contrasts (paired; descriptive — the only gate is R vs PRIME)

- **R − Q** (isolates context): 0.0111 CI[-0.0574, 0.0843] → statistically TIED (CI spans 0)
- **Q − P** (objective + exact-witness supervision, not a pure isolation): 0.0694 CI[-0.0454, 0.1796] → statistically TIED (CI spans 0)
- **A − Q** (DESCRIPTIVE — objective form AND Gartner bag discipline, not a pure objective isolation): -0.188 CI[-0.2981, -0.0861]
- **A − P** (DESCRIPTIVE — supervision granularity): -0.1185 CI[-0.2565, 0.013]
- **R − F** (portable context vs source-name ceiling): 0.0037 CI[-0.1389, 0.1565] → statistically TIED (CI spans 0)

## Verdict

**REJECT.** v0.1 remains the frozen model of record; v0.5 is REJECTED_DEVELOPMENT. Gartner TEST NOT opened; no external claim.

