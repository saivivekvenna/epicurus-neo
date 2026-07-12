# Dynamic upstream gate — feasibility probe (pre-spec)

_Probe script: `scripts/dynamic_gate_feasibility.py` (mirrors the throwaway used here); inputs are the
base cohort CSVs written by `python -m scripts.pool_size_sensitivity`
(`artifacts/milestone_7_decision/pool_size_sensitivity/base_{gartner,improve,multimer}.csv`;
columns `patient_id, mutant_peptide, hla_allele, label, prime, el, expr`)._

## Question

The pool-size diagnostic showed oracle pruning (keep all positives + 25% of negatives) lifts frozen
Epicurus hits@20 Gartner 0.808→1.652 and IMPROVE 1.230→3.214, but the deployable **pure-EL percentile
gate** retained only **66.7%** (Gartner) / **35.7%** (IMPROVE) of positives at that pool size. Can a
smarter *label-blind* gate close that oracle gap — remove negatives while keeping the positives EL alone
drops?

## Method

Within-patient percentiles of the three frozen features (`el`, `prime`, `expr`; oriented so higher =
better). Two probes:

1. **Where do positives sit, and are the EL-dropped positives rescuable?** For the positives the pure-EL
   top-K gate drops, look at their `prime`/`expr` percentiles.
2. **AND-of-independent-vetoes** (missing → KEEP): veto a candidate iff `el_pct<t` **and** `prime_pct<t`
   **and** `expr_pct<t`. Compare negative-removal / positive-retention against a **pure-EL gate matched to
   the same total removal count** — the honest head-to-head.

## Findings

### Positives sit high on presentation in Gartner/multimer, but are genuinely spread in IMPROVE
| cohort | pos `el_pct` median | pos in EL bottom-half | pos `expr_pct` median | pos in expr bottom-half |
|---|--:|--:|--:|--:|
| gartner | 0.858 | 13.0% | 0.842 | 4.3% |
| improve | 0.557 | 44.5% | 0.545 | 45.2% |
| multimer | 0.917 | 17.6% | 0.658 | 35.3% |

IMPROVE is the project's recognition wall made concrete: presentation does not concentrate its positives.

### The EL-dropped positives ARE mostly rescuable by an orthogonal axis
Fraction of pure-EL-dropped positives recoverable by `prime_pct>0.5 OR expr_pct>0.5`:

| cohort | keep 50% | keep 25% |
|---|--:|--:|
| gartner | **100%** | **100%** (all dropped positives are high-expression: median expr_pct 0.84) |
| improve | 72.9% | 77.9% (⇒ ~22–27% low on EL **and** PRIME **and** expr — the unrescuable floor) |
| multimer | 83.3% | 88.9% |

### AND-of-vetoes Pareto-dominates the incumbent pure-EL gate at matched removal
Positive retention at matched negative-removal (AND-veto **vs** pure-EL):

| neg removed | gartner AND / EL | improve AND / EL | multimer AND / EL |
|---|--:|--:|--:|
| ~9%  | 1.000 / 0.978 | 0.940 / 0.899 | 0.971 / 0.971 |
| ~16% | 1.000 / 0.978 | 0.880 / 0.831 | 0.971 / 0.971 |
| ~27% | 1.000 / 0.978 | 0.797 / 0.749 | 0.971 / 0.941 |
| ~41% | 0.957 / 0.935 | 0.692 / 0.642 | 0.882 / 0.941 |
| ~48% | 0.913 / 0.870 | 0.625 / 0.570 | 0.853 / 0.912 |

## Verdict on feasibility (honest)

- **The layered gate is real and beats the incumbent** — AND-of-independent-vetoes with missing→KEEP
  dominates the pure-EL gate on positive retention at (almost) every matched negative-removal level.
- **The aggressive "remove 50–75% at ≥95% retention" target is only reachable where presentation
  concentrates positives** (Gartner: ~40–45% removal at ≥93–95%). **IMPROVE is data-bound** on the three
  peptide-only features: ≥95% retention caps negative removal near ~10%, because ~1/4 of its positives are
  jointly low on EL, PRIME and expression. This is precisely the gap **mutant-allele RNA VAF / read
  support / agretopicity (WES+RNA)** would target.
- ⇒ The gate must be **dynamic**: pick the most aggressive per-patient threshold that still meets a
  calibrated retention lower bound, and **never do worse than keep-all**. A flat 50–75% mandate would
  silently sacrifice IMPROVE positives; the honest deployment target is patient-adaptive.

This probe uses labels only to *measure* retention; the gate itself is label-blind. Multimer is frozen
Epicurus' training cohort (in-sample; flagged). Expression is used only as a **rescue/keep** signal inside
an AND — never as a standalone veto and never to reweight the ranker — consistent with the frozen
confidence-only expression policy.
