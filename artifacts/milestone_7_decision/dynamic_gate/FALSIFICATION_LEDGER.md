# Dynamic upstream gate — falsification ledger

Chronological record of what was tried and what the data said. Nothing here is tuned to the LOCKED
CheckMate test; entries below are DEV-only unless stated.

| # | variant / hypothesis | result | kept? |
|---|---|---|---|
| 1 | **Pure-EL percentile gate** (incumbent from the pool-size diagnostic) | Retains only 66.7% (Gartner) / 35.7% (IMPROVE) of positives at the 25%-negative pool. Baseline to beat. | baseline |
| 2 | **AND-of-independent-vetoes** (el ∧ prime ∧ expr, missing→KEEP) | **Pareto-dominates** the pure-EL gate on positive retention at every matched negative-removal level on Gartner/IMPROVE; ties/edges on multimer (FEASIBILITY.md). Adopted as Layer 1. | **YES** |
| 3 | **expr as a required veto axis** (3-way AND) | Identical to #2 on full-coverage dev cohorts, but disables the gate entirely on CheckMate (81% expr-missing → every candidate auto-kept). Rejected in favor of **core={el,prime} + expr rescue-only**, which is byte-identical on dev and still operates on the external cohort. | replaced by #4 |
| 4 | **core veto {el,prime} + expression rescue-only** | Adopted. Preserves the frozen confidence-only expression policy (expr never vetoes alone, never reweights the ranker) and keeps the fail-open guarantee on the core axes. | **YES** |
| 5 | **Single global threshold `t` calibrated on pooled dev positives (target 0.95)** | Dragged to `t=0.25` by IMPROVE's spread positives → removes only 2–10% of negatives. Safe but modest. This is the frozen deployment config. | frozen |
| 6 | **Aggressive per-cohort operating points to hit 50–75% removal** | To remove ~50% of negatives, retention collapses (IMPROVE 0.45–0.62, Gartner 0.85–0.91). Violates the ≥95% safety bar. The 50–75%@≥95% target is **unreachable** on peptide-only features. | **FALSIFIED** |
| 7 | **Premise: a safe gate improves downstream top-20 (closes the oracle gap)** | **FALSIFIED.** Threshold sweep: at any safe retention (≈1.0) the paired downstream Δhits@20 ≈ 0 — the gate removes only sub-top-20 negatives. It removes **0%** of the high-presentation decoys that outrank positives (the ones the oracle's random pruning deletes). Apparent lift appears only after retention has fallen (denominator effect). | **FALSIFIED** |
| 8 | **CP retention lower bound ≥ 0.95 as a universal bar** | Sample-size-capped: `CP-LB(46/46) = 0.937`, so Gartner (46 positives) can never certify ≥0.95 regardless of the gate. The bar is only meaningful on the larger cohorts (IMPROVE 467). Reported as an underpowering limit, not an unsafety. | noted |

## Not yet tried (would need new data / next iteration)
- Predictor-disagreement rescue (Layer 2) on Gartner's 5-predictor set — implemented + unit-tested but OFF
  in the frozen config; its effect is a rescue (raises retention, lowers removal), so it cannot change the
  falsification of #6/#7, only make the gate safer. Evaluate once a cohort needs it.
- WES/RNA rescue axes (mutant-allele RNA VAF, DNA VAF/depth/purity/CCF, processing, agretopicity) — the
  ONLY lever that can separate positives from high-presentation decoys and thus close the oracle gap.
  Blocked on data (absent across all current eval cohorts). Contract in SPEC §7.
