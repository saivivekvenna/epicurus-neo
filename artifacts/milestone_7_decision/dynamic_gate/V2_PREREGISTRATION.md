# Dynamic gate v2 — pre-registration of the primary diagnostic (Miller / Gartner reconstruction)

_Pre-registered BEFORE the WES/RNA data exists, so the analysis has no post-hoc degrees of freedom. Fixed
2026-07-12. Runs only once the reconstructed cohorts (Miller IPV `PRJNA980652`, Gartner WES/RNA) provide
the orthogonal features in V2_CONTRACT.md. Verdict rule = the project's standard (ACCEPT iff paired CI
lower bound > 0 vs both baselines; else CONSISTENT_WITH_NO_EFFECT / REJECT)._

## Hypothesis

An orthogonal-residual gate (V2_CONTRACT.md) can remove a meaningful fraction of **high-EL/high-PRIME
TESTED_NEGATIVE decoys** while retaining recognized positives, and thereby lift downstream recognized
hits@20 — something a same-feature presentation gate structurally cannot do (CIRCULARITY_AUDIT.md).

## Population, units, roles (fixed before fitting)

- **Stratum S** = candidates with within-patient EL percentile > 0.75 AND PRIME percentile > 0.75 (the
  hard-decoy region). All primary metrics are computed **within S**.
- **TRAIN/calibration:** ≥2 studies with reconstructed orthogonal features, mixed (never single-study).
- **DEV:** leave-one-STUDY-out across the reconstructed cohorts.
- **LOCKED TEST:** the last-arriving reconstructed study, frozen config, scored once. CheckMate 153 is NOT
  used (consumed locked evidence for v1).
- Only explicit POSITIVE / TESTED_NEGATIVE; PU cohorts excluded from removal claims.

## Primary endpoint (the one number that decides v2)

Within stratum S, under leave-one-study-out:

1. **Negative removal at ≥95% positive retention** — the largest fraction of stratum tested-negatives the
   gate removes while the Clopper–Pearson lower bound on stratum positive retention stays ≥ 0.95. Reported
   per held-out study and worst-study.
2. **Paired downstream hits@20** — gate→UNCHANGED genuine PRIME and frozen Epicurus vs (a) the v1 AND-gate
   and (b) ungated, paired per patient, 20k bootstrap, `pre_registered_verdict`.

**Pre-declared success (v2 ACCEPT):** worst-study stratum negative removal ≥ 20% at CP-LB retention ≥ 0.95,
**and** downstream hits@20 delta CI lower bound > 0 vs BOTH the v1 gate and ungated, **and** OOD abstention
never had to fire to achieve it. Anything less is CONSISTENT_WITH_NO_EFFECT or REJECT (no superiority
claim).

## Secondary / diagnostic

- Orthogonal-residual AUROC within S (leave-one-study-out) — the mechanism check.
- OOD-abstention rate per held-out study and the retention it protected.
- Ablation: which orthogonal feature carries the signal (mutant-RNA VAF vs DNA VAF/CCF vs processing vs
  agretopicity), leave-one-feature-out.
- Negative control: shuffle the orthogonal features across candidates within study → removal at 95%
  retention must collapse to ~random (else the "signal" is leakage).

## Falsification conditions (pre-committed)

- Cross-study transfer collapses (retention ≪ 0.95 out-of-study without abstention) → REJECT; the OOD
  detector must have caught it. (This already happened for the **sequence-only** variant — retained 1.5%
  of positives train-Gartner→test-IMPROVE — which is why sequence-motif vetoes are excluded a priori.)
- Orthogonal AUROC within S ≈ 0.5 leave-one-study-out → mechanism absent → REJECT; the recognition wall
  holds even with genomics.
- Any downstream lift that appears only after retention < 0.95 → denominator artifact, not a win.

## Guardrails

- EL/PRIME define S and the rescue/floor ONLY; never v2 veto inputs (no circularity).
- No sequence-motif / n-gram / anchor-token veto features (falsified; encode study shift).
- Mixed-study calibration + leave-one-study-out + OOD abstention are mandatory.
- Missing orthogonal evidence ⇒ KEEP. No field imputed to a value that could enable a veto.
- Freeze the config before the LOCKED study; no peeking, no re-tuning on consumed cohorts.
