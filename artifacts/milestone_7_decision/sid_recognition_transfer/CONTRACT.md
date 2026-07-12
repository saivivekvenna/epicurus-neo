# Non-Sid recognition-gate → single-shot Sid evaluation — PREREGISTERED CONTRACT

_Frozen before any Sid label/rank is consulted for selection. Goal: a DEFENSIBLE 3/3 on Sid, or an honest
report that a non-Sid-trained recognition gate does not reach it. Falsification-first; n=1 patient / 3
recognized positives on the Sid side is descriptive, never a powered claim._

## 0. Hard firewall (Sid is locked)

Sid's recognition labels (ASPM/DYNC1H1/MAP2) and Sid ranks are FORBIDDEN inputs to: feature selection,
threshold selection, model fitting, calibration, hyperparameter/alpha selection, AND arm/pipeline
selection. Every such choice is made **solely** on non-Sid cohorts. Sid is touched exactly **once**, at
the end, with a single frozen pipeline whose config + SHA-256 are serialized before that touch. Any
post-Sid change voids the result.

## 1. Training/selection cohorts (non-Sid only)

- IMPROVE (rich 88-col; real response labels), CD8 multimer (POSITIVE/TESTED_NEGATIVE; in-sample to frozen
  Epicurus — flagged), Gartner NCI (POSITIVE/TESTED_NEGATIVE). Zhao/CEDAR: used only if a valid comparable
  feature + tested-negative subset exists; otherwise documented as excluded.
- All are patient/study-disjoint from Sid (osteosarcoma, single patient).

## 2. Comparable feature set (intersection with Sid's 137 generated mutations)

Verified computable for BOTH the training cohorts AND Sid's 137 without consulting Sid labels:
- **presentation**: genuine PRIME %rank, MixMHCpred %rank (both, via the installed tool);
- **peptide physicochemical**: hydrophobic/aromatic fraction + core hydrophobicity + length (from the
  mutant peptide sequence);
- **expression**: within-patient percentile of expression (IMPROVE Expression / multimer TPM / Gartner
  decile / Sid gene TPM).
Partial bridges (declared, NOT in the core set): agretopicity/WT-differential (needs WT-peptide scoring for
Sid's lossless peptides) and mutant-RNA VAF / DNA VAF (IMPROVE↔Sid only). If the core set fails, these are
the documented missing bridges, not a post-hoc rescue.

## 3. Architecture (base-anchored; frozen ranker unchanged)

Base-anchored residual (validated form): `U = base_pct + α·residual_pct`, where base = the frozen ranker's
own within-patient percentile (genuine PRIME as primary; MixMHCpred as the presentation-only secondary),
and `residual` = a recognition model over the comparable NON-presentation features (physchem, expression).
`α=0` reproduces the base exactly (no-op). The base ranker is never replaced; the gate only reselects.

## 4. Model families attempted (declared up front — full multiplicity reported)

1. **null / no-op** (α=0) — the honest control.
2. **base-anchored physchem residual** (logistic on hydrophobicity+length).
3. **base-anchored expression residual**.
4. **base-anchored physchem+expression residual**.
5. **absence-safe expression/unexpressed gate** (biological, non-fit).
Base ∈ {genuine PRIME, MixMHCpred}. α ∈ {0, 0.25, 0.5, 1, 1.5, 2}. All combinations are the multiplicity;
the count is reported and a single winner is chosen by §6.

## 5. Splits (grouping to prevent leakage)

Leave-one-STUDY-out (cohort) OUTER across {IMPROVE, multimer, Gartner}; within the training studies,
patient-group inner folds for α/feature selection. Exact + near-peptide (≥0.8) groups kept together;
recurrent peptides quarantined. Study identity is never a model feature.

## 6. Objective + freeze criterion (non-Sid only)

Clinically-aligned objective: mean per-patient recognized **hits@20** (primary) and **recall@20**
(secondary) after gate→unchanged base ranker, with a catastrophic-loss constraint (worst-study mean Δ must
not be < −0.10 vs the base). Select the ONE (family, base, α) maximizing the outer leave-one-study-out mean
Δhits@20 subject to the constraint. If no candidate beats the null by a positive outer mean that also
satisfies the constraint, **freeze the null (α=0)** — i.e. do not act. Serialize the winner's config +
SHA-256 to `configs/frozen/sid_recognition_gate_v1.json` BEFORE touching Sid.

## 7. The single Sid evaluation

Apply the frozen pipeline to Sid's full accounted denominator (137 generated; 10 documented-
unrepresentable). Report tie-aware (exact-score interval) nominal / guaranteed hits@20 for: genuine PRIME,
frozen Epicurus v0.1, and the frozen gate on each base — vs the same baselines. One run, no iteration.

## 8. Reporting

All attempted (family × base × α) arms with their non-Sid outer scores; the multiplicity count; the frozen
choice + hash; the single Sid result; and an explicit statement of whether a DEFENSIBLE 3/3 was achieved
(guaranteed, not tiebreak-nominal) WITHOUT any Sid-informed choice. If the frozen gate ties PRIME / fails to
reach guaranteed 3/3, that is the honest verdict, consistent with the recognition-transfer wall; document
the missing feature bridge (agretopicity / mutant-RNA VAF) as the next lever.
