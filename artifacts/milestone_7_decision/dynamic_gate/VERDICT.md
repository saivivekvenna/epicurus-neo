# Dynamic upstream gate — verdict

**B for the gate as a safe recall-preserving pruner. C for the aggressive premise (close the oracle gap /
improve downstream top-20 with a label-blind presentation gate).**

## What was built
A layered, label-blind safe-rejection gate (`src/event_b/dynamic_gate.py`, frozen
`configs/frozen/dynamic_gate_v1.json`, CLI `epicurus-neo dynamic-gate`):
- **Layer 0** deterministic impossibility gate (reused `apply_deterministic_gate`; inert on research
  cohorts, active on WES/RNA product inputs);
- **Layer 1** AND-of-independent-vetoes over core axes `{el, prime}` (missing → KEEP; strong on any axis → KEEP);
- **Layer 2** expression rescue-only + optional predictor-disagreement / near-boundary rescue;
- **Layer 3** patient-adaptive rails (top-M-by-EL floor, removal cap, pool/coverage floors → keep-all fallback);
- calibration by Clopper–Pearson lower bound on calibration positives (Neyman–Pearson / conformal),
  leave-one-cohort-out, frozen before the LOCKED CheckMate 153 test.

16 safety-invariant unit tests lock the guarantees (`tests/test_dynamic_gate.py`).

## What the data said
1. **The gate is real and beats the incumbent.** The AND-of-vetoes Pareto-dominates the pure-EL percentile
   gate on positive retention at matched negative-removal on Gartner/IMPROVE (`FEASIBILITY.md`).
2. **The 50–75%-removal-at-≥95%-retention target is unreachable on peptide-only features.** To remove ~50%
   of negatives, retention collapses (IMPROVE 0.45–0.62). At a universally-safe operating point removal is
   only ~2–10%.
3. **A safe gate does NOT improve downstream top-20.** At any safe retention (≈1.0) the paired downstream
   Δhits@20 is ≈0 — the gate removes only negatives already below top-20. It removes **0%** of the
   high-presentation decoys that outrank positives, which is exactly what the oracle's random pruning
   deletes to earn its lift (Gartner 0.808→1.652). Apparent lift appears only after positives are being
   dropped (denominator effect). **The oracle gap is unclosable by any label-blind presentation gate** —
   this re-derives the project's recognition wall from the gate/oracle decomposition.
4. **CP lower bounds are sample-size-capped** (`CP-LB(46/46)=0.937`): Gartner cannot certify ≥0.95
   retention regardless of the gate — an underpowering limit, not an unsafety.

## Honest deployment value
Not top-20 lift. The gate's real use is **shrinking the candidate universe handed to expensive downstream
steps** (genuine-PRIME scoring, wet-lab validation) at a calibrated recall floor, and doing so strictly
better than the incumbent EL gate. On the LOCKED CheckMate 153 test the frozen config removed 10.4% of
negatives at 95.7% retention with no downstream regression.

## What unlocks the next level (WES/RNA extension contract, SPEC §7)
An orthogonal signal that separates true positives from high-presentation decoys — **mutant-allele RNA VAF
+ read support, tumor DNA VAF/depth/purity/CCF, proteasomal processing, agretopicity** — consumed as
KEEP-only rescue axes (never imputed to a vetoing value). Absent across all current eval cohorts. Concrete
path: open WES+RNA of Miller IPV (`PRJNA980652`) and Gartner reconstruction.
