# Recognition-transfer — POSTMORTEM (no policy changes)

_Written after the final one-shot Stage-2 result. This document explains the outcome and names the next
falsifiable work. It changes NO model, threshold, or artifact. Sid is now LOCKED from all further
selection._

## What was tested

A recognition gate frozen **entirely on non-Sid data** (Stage 1, `configs/frozen/sid_recognition_gate_v1.json`,
SHA `e44d3d1a…`) — the combined family = anchored logistic `{PRIME, EL, Expression, WES-VAF}` (C=0.5),
`score = prime_pct + 0.1·pred_pct`, plus a q=1 promote-side mutant-RNA (`rna_af`) reserve — was applied
**once**, apply-only (verified SHAs, no refit), to Sid's full accounted denominator (137 generated + 10
documented-unrepresentable = 147).

## Result (final; tie-aware mutation-level hits@20)

| arm | guaranteed hits@20 |
|---|--:|
| genuine PRIME (baseline) | **2/3** (DYNC1H1, MAP2) |
| frozen Epicurus v0.1 | 1/3 (DYNC1H1) |
| **frozen non-Sid gate** | **2/3 — identical to PRIME** (DYNC1H1, MAP2) |

ASPM stays at PRIME score-rank **41** ([41,41], no boundary tie). The label-blind q=1 reserve slot went to
**NBPF1** (`rna_af` > ASPM's 0.5), not to a recognized positive. So the gate preserves PRIME's two hits and
does **not** reach 3/3.

## Why it did not transfer

- Stage 1's non-Sid evidence was **weak**: nested per-family Δ = **+3 hits** (83 vs 80), transport-positive
  and beating matched-random (77), **but the patient-level paired bootstrap Δ = +0.043, CI [−0.043, +0.129],
  p>0 = 0.80 — spanning zero**. It was directional, not significant. A signal that underpowered on IMPROVE
  should not be expected to convert to a discrete top-20 gain on a single osteosarcoma patient.
- The one candidate the gate *could* rescue (ASPM, strong mutant RNA VAF 0.5) is not the top-`rna_af`
  non-protected mutation on Sid, and the anchored logistic coefficients are tiny
  (`[0.0157, 0.0055, 0.0007, 0.0022]`), so the presentation anchor dominates and the residual barely moves
  ranks. This is consistent with the whole milestone's finding: **presentation is the ceiling; a
  recognition signal that generalizes across cohorts/regimes has not been found.**
- Sid is n=1 and was previously inspected, so even the tie (2/3 = 2/3) is **exploratory confirmation, not
  pristine external validation**.

## Honest verdict

**The non-Sid-trained recognition gate ties genuine PRIME at 2/3 and does not achieve a defensible 3/3.**
Stage 1's non-Sid signal did not transfer to a Sid top-20 gain. No overfit 3/3 was manufactured — the
firewall held (Sid informed no selection), and the honest answer is a tie.

## Process disclosure (not retroactively claimed)

The hard one-shot guard/sentinel (refuse if a result exists; atomic receipt before reading Sid) and the
exhaustive brute-force tie-resolution tests requested in the final audit were **NOT** present in the
Stage-2 code commit `19f505b` and were not added before the run. The gate's guaranteed-tie accounting used
an analytic (not brute-forced) rule; it happens not to be decision-relevant here (all three positives have
degenerate score intervals `[b,b]` and ASPM is far outside top-20), but the exhaustive verification was not
performed. These are acknowledged gaps in the Stage-2 tooling, not in the frozen result.

## Next falsifiable, non-overfit work (Sid LOCKED from all further selection)

1. **Move evaluation off n=1.** Build/evaluate on **multi-patient held-out FULL denominators** (not the
   pre-screened tested subsets) — the honest metric is per-patient recognized hits@20 / recall@20 on the
   complete somatic universe. Miller IPV (`PRJNA980652`, open WES+RNA, ELISpot negatives) is the
   pre-registered target; the same leakage-clean, patient/study/peptide-grouped, nested protocol applies.
2. **Reframe the lever as a high-RECALL upstream reducer, not a top-20 reranker.** The consistent wall is
   that recognition doesn't reorder the top of the list; the tractable, biologically-grounded win is
   removing candidates that are impossible neoantigens (unexpressed mutant allele, lost HLA, subclonal/
   low-VAF, unprocessed) to shrink the denominator at guaranteed recall — derived from WES/RNA/presentation
   biology, validated cross-cohort with OOD abstention, never tuned on any single patient's labels.
3. **Keep Sid as locked descriptive evidence only.** It may illustrate a deployed pipeline end-to-end, but
   it must never again inform feature/threshold/model/arm selection.
