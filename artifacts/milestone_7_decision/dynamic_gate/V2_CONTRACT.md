# Dynamic gate v2 — orthogonal residual gate contract (data-blocked)

_Design-only. v2 is BLOCKED on WES/RNA reconstruction (Miller IPV `PRJNA980652`, Gartner). This contract
fixes the design so that when the data lands, no degrees of freedom are chosen post-hoc. Supersedes
nothing in v1; v1 stays frozen (`configs/frozen/dynamic_gate_v1.json`)._

## Why v2 exists (the v1 circularity finding)

v1's core veto axes {EL, PRIME} are the dominant inputs to the downstream rankers (genuine PRIME; frozen
Epicurus is PRIME-dominated over prime/el/expr). So v1 provably cannot move a top-20 those same features
produce (CIRCULARITY_AUDIT.md: 0 removed candidates in any ranker top-20; 0% hard-decoy removal). v1's
falsification is therefore **scoped to same-feature presentation gating** — it says nothing about a gate
built on features ORTHOGONAL to the downstream rank. v2 tests exactly that.

## The only candidates that matter: the hard-decoy stratum

The oracle's top-20 lift comes from deleting **high-EL / high-PRIME TESTED_NEGATIVE decoys** that outrank
positives. v2 must target THAT stratum specifically:

- **Stratum S** = candidates with within-patient EL percentile > 0.75 **and** PRIME percentile > 0.75.
- EL/PRIME are used ONLY to (a) define S and (b) provide the v1 rescue/floor safety layer. They are
  **NEVER v2 veto inputs** (that would reintroduce the circularity).

## v2 veto features (orthogonal to the downstream rank; all KEEP-biased)

Each is a per-candidate biological quantity NOT consumed by PRIME/EL/frozen-Epicurus. Missing ⇒ KEEP.

| feature | orientation (veto-supporting = "candidate is a dead decoy") | source |
|---|---|---|
| mutant-allele RNA VAF | ~0 mutant RNA reads ⇒ not transcribed | RNA-seq (reconstruction) |
| RNA mutant read support / depth | no/low mutant read coverage | RNA-seq |
| allele-specific / mutant-transcript expression | mutant isoform not expressed | RNA-seq |
| tumor DNA VAF + depth | very low VAF / low depth ⇒ subclonal or artifact | WES |
| CCF / clonality | subclonal (low CCF) | WES + purity |
| proteasomal processing / cleavage / stability | not processed / unstable | predictors (NetChop/NetMHCstab) |
| agretopicity / WT-differential | mutant ≈ WT (no neo-signal; self-tolerance) | mut-vs-WT binding |
| transcript / annotation confidence | low-confidence transcript / artifact | VEP/gene model |

A candidate is veto-eligible only if it is in S **and** an orthogonal-residual model (below) confidently
rates it a dead decoy **and** no rescue fires. No single missing field can cause a veto.

## Learning: cross-fitted residual after conditioning on PRIME/EL

1. Restrict to stratum S (this conditions on presentation by stratification — EL/PRIME are ~flat within S).
2. Fit the orthogonal veto model on the **residual** label (POSITIVE vs TESTED_NEGATIVE within S), i.e. the
   part of recognition NOT explained by presentation. Optionally regress out residual EL/PRIME variation
   inside S before fitting so the model cannot smuggle presentation back in.
3. **Cross-fitting is leave-one-STUDY-out**, never within-study CV and never leave-one-patient-out only.
   Rationale — falsified variant #S1 below.

## Mandatory OOD abstention (falsified variant #S1 forced this)

An independent **sequence-only** hard-decoy gate (peptide n-grams + allele/anchor + length/chemistry;
PRIME/EL only for the stratum) collapsed cross-study: train-Gartner→test-IMPROVE retained **1.5%** of
positives (removed 198/201); train-IMPROVE→test-Gartner retained **12.8%** (removed 34/39). Peptide/HLA
recognition motifs encode severe study/assay/domain shift. Therefore v2:

- **forbids sequence-motif / n-gram / anchor-token veto features** (they memorize study, not recognition);
- requires an **OOD detector**: if a held-out study's stratum distribution (on the orthogonal features) is
  far from the calibration studies (e.g. energy/MMD distance or per-feature coverage below a floor), the
  gate **abstains → KEEP-all** for that study/patient;
- requires **mixed-study calibration** (≥2 studies) and reports leave-one-study-out transfer, not a single
  train/test direction.

## Safety bar (unchanged from v1 §6, plus)

- v2 veto acts ONLY inside stratum S; outside S the v1 rescue/floor governs (nothing new removed).
- worst-study positive-retention CP lower bound ≥ 0.95 under leave-one-study-out;
- OOD abstention must trigger before retention can fall (verified on a held-out study);
- deployment candidate frozen before any untouched test; no tuning on consumed cohorts (CheckMate is
  consumed locked evidence for v1 and is off-limits).

## What v2 will report (primary diagnostic in V2_PREREGISTRATION.md)

Within stratum S: negative removal at ≥95% positive retention (leave-one-study-out, CP-bounded), then the
**paired downstream hits@20** of gate→unchanged PRIME/Epicurus vs the v1 AND-gate and vs ungated. A
positive, OOD-robust result here is the first thing in the project that could move top-20 without leaking
presentation.
