# End-to-end benchmark synthesis — three levels, genuine features, frozen expression policy

This ties together the milestone-7 benchmark work into one deployable recommendation and one honest
negative result. It is NOT a superiority claim; the only end-to-end-evaluable patient is post-hoc, n=3.

## The benchmark is three levels, never pooled

| level | task | eligible cohorts | status |
|---|---|---|---|
| **L1 reachability** | raw variants → peptide generation, stage-loss attribution | osteosarc/Sid only (raw callset + labels) | lossless generation recovers **1/3 → 3/3** recognized on Sid |
| **L2 conditional ranking** | rank among generated candidates, within each denominator | multimer, Gartner, IMPROVE, CheckMate (+Sid) | each stands alone; presentation is the ceiling |
| **L3 end-to-end patient utility** (north star) | recognized in final top-20 from common raw inputs vs pVAC+PRIME | osteosarc/Sid only (post-hoc, n=3) | protected incumbent lossless+PRIME = **3/3** |

Cohort roles are fixed and never merged (CEDAR/Zhao = recognition priors; multimer = presentation/
T-cell compatibility + Epicurus training; Gartner = broad-denominator; IMPROVE = prefiltered subset;
Sid = end-to-end; RTTP = deployment). See `COHORT_ELIGIBILITY.md`.

## What we can and cannot beat

- **Reachability (generation) is the real, reproducible lever.** The input-only lossless generator
  regenerates the missed recognized windows from the raw allele + Ensembl alone; feeding genuine PRIME
  it lifts Sid recognized-in-top-20 from 1/3 to 3/3. This is a Level-1/Level-3 gain and the protected
  incumbent.
- **Ranking beyond lossless+PRIME is data-bound.** With GENUINE presentation features computed on the
  recovered candidates (MHCflurry, no imputation), the Epicurus scorer nets at best +1 on Sid but still
  demotes low-expression MAP2; and on the development cohorts EVERY lever that moves the rank on
  expression or reserves slots regresses at least one conditional-ranking cohort within a fixed top-20
  budget (expression rank penalty → −multimer/−IMPROVE; expression-stratum reserve → −Gartner;
  predictor-disagreement reserve → −Gartner). See `SID_FOUR_ARM.md` and `../expression_policy/`.

## Frozen deployable policy

**Deployable end-to-end = lossless generation + genuine PRIME + confidence-only expression annotation.**
RNA expression is frozen as confidence-only (`configs/frozen/expression_policy_v1.json`); it never moves
the rank. lossless+PRIME is the protected incumbent (already keeps high-presentation candidates
irrespective of expression, preserving reachability of low-expression recognized ones). The
soft-saturating guard and the portfolio reserve are retained but off-by-default (no-regression-equivalent
or worse on current cohorts).

## The one lever that remains: an untouched end-to-end patient

Every ranking result above is either post-hoc (Sid) or on a semi-consumed / range-restricted cohort. A
genuine superiority claim requires a Level-3-eligible, PRIME-untouched patient: **raw multi-caller
callset + WES/RNA + class-I HLA + measured recognition, labels joined only after ranking.** Best routes
(documented, no external action taken here):

1. **Miller IPV `PRJNA980652`** — the only fully-OPEN end-to-end (L3) build prospect: raw WES+RNA public,
   per-peptide ELISpot negatives recorded (n=13). Verify the S1/S2 label supplement, then build.
2. **CheckMate 153 raw WES/RNA — EGA `EGAD00001011302` (controlled, BMS discretion), NOT dbGaP** —
   upgrades the already-run L2 external cohort to L1/L3.
3. **More Hudson-lab / RTTP patients** with the IFNγ peptide-expansion assay **plus the stimulation-pool
   composition** (the true recognition denominator) — converts n=3 into a real denominated benchmark.

Canonical acquisition plan + tracker + data-request schema + outreach drafts:
`../external_validation/ACQUISITION_EXECUTION_PLAN.md`. Until an L3 run clears the go/no-go gate, this
benchmark is honest infrastructure + a post-hoc diagnostic, not a gate. See `../../../NORTH_STAR_HISTORY.md`
and `memory/benchmark-three-level-hierarchy.md`, `memory/expression-ranking-policy.md`.
