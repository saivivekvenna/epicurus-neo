# Dynamic upstream gate — verdict

**B for the gate as a safe recall-preserving pruner. C is SCOPED — falsified only for same-feature
presentation/PRIME-derived monotone gating with current peptide features. The general orthogonal-feature
dynamic-gate hypothesis is UNTESTED / data-blocked, NOT falsified.**

> **Correction (post-review).** An earlier draft wrote "C: any label-blind gate cannot improve top-20."
> That over-generalized. The v1 gate's veto axes {EL, PRIME} ARE the dominant downstream-ranker inputs, so
> its zero top-20 effect is a **structural tautology of feature overlap** (CIRCULARITY_AUDIT.md), proving
> only that a same-feature monotone presentation gate cannot move a same-feature top-20. A gate built on
> features orthogonal to the rank (WES/RNA/clonality) is a different, untested hypothesis — see
> V2_CONTRACT.md / V2_PREREGISTRATION.md. Verdict wording narrowed accordingly.

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
3. **A same-feature presentation gate does NOT improve downstream top-20 (circular by construction).** At
   any safe retention (≈1.0) the paired downstream Δhits@20 is ≈0 — the gate removes only negatives
   already below top-20, and **0%** of the high-EL/high-PRIME decoys that outrank positives. But this is a
   **tautology**: the veto axes {EL, PRIME} are the dominant downstream-ranker inputs, so the gate cannot
   move a top-20 those same features produce (CIRCULARITY_AUDIT.md — 0 removed candidates in any ranker
   top-20; Spearman(keep-margin, ranker)≈+1). The oracle's lift (Gartner 0.808→1.652) comes from deleting
   the high-presentation decoys — the very stratum a same-feature gate cannot touch. **Falsified for
   presentation/PRIME-derived gating only.** Whether an ORTHOGONAL-feature gate can delete those decoys is
   untested (v2, data-blocked).
4. **CP lower bounds are sample-size-capped** (`CP-LB(46/46)=0.937`): Gartner cannot certify ≥0.95
   retention regardless of the gate — an underpowering limit, not an unsafety.

## Honest deployment value
Not top-20 lift. The gate's real use is **shrinking the candidate universe handed to expensive downstream
steps** (genuine-PRIME scoring, wet-lab validation) at a calibrated recall floor, and doing so strictly
better than the incumbent EL gate. On the LOCKED CheckMate 153 test the frozen config removed 10.4% of
negatives at 95.7% retention with no downstream regression.

## What unlocks the next level (v2, data-blocked — V2_CONTRACT.md / V2_PREREGISTRATION.md)
A gate whose veto uses features **orthogonal to the downstream rank** — mutant-allele RNA VAF + read
support, tumor DNA VAF/depth/CCF/clonality, proteasomal processing/stability, agretopicity/WT-differential,
transcript confidence — as **cross-fitted residuals after conditioning on PRIME/EL** (EL/PRIME define only
the hard-decoy stratum + rescue/floor, never veto axes). Missing⇒KEEP.

**Constraints already learned (do not repeat):**
- A **sequence-only** residual gate is FALSIFIED — an independent train-Gartner→test-IMPROVE run retained
  **1.5%** of positives; peptide/HLA motifs encode severe study/assay shift. v2 forbids sequence-motif
  vetoes and requires **leave-one-study-out** calibration + **explicit OOD abstention**.
- An in-sample multimer probe (orthogonal features, hard-decoy stratum) shows only a **faint** signal
  (AUROC ~0.59, ~9% decoy removal at 95% retention) — in-sample, underpowered, not cross-study validatable
  here (orthogonal_probe.json).

Concrete path: open WES+RNA of Miller IPV (`PRJNA980652`) + Gartner reconstruction. **v2 is correctly
blocked on that data.**
