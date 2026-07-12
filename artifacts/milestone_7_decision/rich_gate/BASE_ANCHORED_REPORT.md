# Rich-gate v2 — base-anchored residual gate (IMPROVE nested)

`python -m scripts.rich_gate_base_anchored`

Utility: `U = base_percentile + alpha*(feature_percentile-0.5); Epicurus UNCHANGED (base anchor)`.
_DEVELOPMENT/LOCAL only. rich-gate v1 (orthogonal-only utility) remains NULL/falsified. The IMPROVE hydrophobic gain does NOT transfer externally (see external_transfer) => regime-aware abstention required; DO NOT freeze hydrophobicity into product._

Base-only no-op sanity (alpha=0 == Epicurus): **PASS** (max|Δ|=0.0).


| feature | class | **Δ hits@20** [CI] | p_better | imp/tie/harm | random-matched Δ | beats random | pos/neg swapped-in | chosen alphas |
|---|---|--:|--:|--:|--:|:--:|--:|--:|
| PropHydroAro | primary | **+0.371** [0.1143, 0.6286] | 0.998 | 23/40/7 | +0.055 | yes | 51/712 | [2.0, 1.5, 1.5, 1.5, 2.0] |
| HydroCore | primary | **+0.257** [0.0429, 0.4714] | 0.992 | 23/36/11 | +0.032 | yes | 41/671 | [1.0, 0.75, 1.0, 1.0, 0.75] |
| HydroAll | attribution | **+0.229** [0.0, 0.4571] | 0.968 | 21/38/11 | -0.038 | yes | 45/734 | [2.0, 0.5, 1.5, 2.0, 2.0] |
| VarAlFreq | attribution | **+0.086** [-0.0429, 0.2143] | 0.888 | 13/51/6 | +0.062 | yes | 19/369 | [0.2, 0.2, 0.2, 0.3, 0.2] |
| SelfSim | attribution | **-0.057** [-0.1857, 0.0714] | 0.152 | 4/59/7 | +0.001 | NO | 8/193 | [0.2, 0.2, 0.0, 0.0, 0.2] |
| DAI | attribution | **-0.200** [-0.3714, -0.0429] | 0.004 | 4/55/11 | +0.013 | NO | 8/301 | [0.2, 0.1, 0.0, 2.0, 0.0] |

## Verdict (local IMPROVE)
DEVELOPMENT-POSITIVE (single study): base-anchored gate on PropHydroAro lifts NET hits@20 by +0.371 (CI [0.1143, 0.6286], beats matched-random +0.055) on nested patient-disjoint IMPROVE folds. Feature choice followed inspection => development discovery, NOT external proof. Requires an UNTOUCHED external rich cohort to confirm before any freeze.


## Transport 1 — leave-source-out WITHIN IMPROVE (PropHydroAro)

Pooled Δ **+0.371**; all sources positive: **True**. Per source: Basket +0.500 (α=1.5), bladder +0.083 (α=1.5), melanoma +0.538 (α=2.0).


## Transport 2 — EXTERNAL transfer (frozen on IMPROVE → Gartner/multimer)

Frozen α=1.5 on IMPROVE, applied label-blind with a peptide-derived hydro/aromatic proxy:

| external cohort | Δ hits@20 | imp/tie/harm |
|---|--:|--:|
| gartner | **+0.000** | 3/20/3 |
| multimer ⚠️in-sample | **-0.421** | 1/10/8 |

_External transfer FALSIFIES a universal hard gate: the IMPROVE hydrophobic gain does not transport (harms Gartner/multimer). Consistent with the user's independent external run (Gartner −0.154, multimer −0.577).


## Why IMPROVE differs (hydrophobicity of positives vs negatives by source)

| source | n | pos | mean pep len | hydro(pos) | hydro(neg) |
|---|--:|--:|--:|--:|--:|
| Basket | 5362 | 116 | 9.62 | 0.545 | 0.44 |
| bladder | 6237 | 147 | 9.36 | 0.494 | 0.462 |
| melanoma | 5921 | 204 | 9.35 | 0.516 | 0.444 |

_If positives are systematically MORE hydrophobic than negatives only in IMPROVE (a TIL/screened set), the gain reflects IMPROVE's candidate-selection/assay regime, not a universal recognition rule._


## Regime-aware transport verdict

- local significant (IMPROVE nested): **True**
- source-invariant within IMPROVE: **True**
- external non-harmful: **False**

**Decision:** REGIME-LOCAL: significant within IMPROVE and across its tissue sources, but FALSIFIED externally (Gartner/multimer harmed) => a regime-aware gate must ABSTAIN (no-op) off the IMPROVE-like regime. NOT frozen into product.

_Rule: activate only if local_significant AND source_invariant AND external_nonharmful; else abstain._


> Epicurus is UNCHANGED (the base anchor); the feature only nudges. alpha=0 reproduces Epicurus (no-op sanity). rich-gate v1 (orthogonal-only utility) stays NULL/falsified. The base-anchored architecture is validated LOCALLY on IMPROVE but hydrophobicity is regime-specific and MUST NOT be frozen into product; a regime-aware gate abstains off-regime. Next lever = an untouched external RICH cohort + leave-source-out-supported features.
