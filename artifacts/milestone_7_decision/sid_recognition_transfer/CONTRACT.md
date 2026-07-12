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

## 3b. Two-stage protocol with a hard pre-Sid checkpoint

**Stage 1 (non-Sid only):** select and FREEZE one arm/null using only IMPROVE (+ multimer/Gartner for the
core arm). Serialize config+SHA-256. Commit. Audit that no Sid file/label was accessed. **STOP.**
**Stage 2 (later, separate):** load Sid and run exactly one tie-aware evaluation with the frozen config.
This file is committed at the Stage-1 checkpoint; Stage-2 code does not exist yet.

## 4. Model families attempted (declared up front — full multiplicity reported)

1. **null / no-op** (α=0) — the honest control.
2. **base-anchored physchem residual** (logistic on hydrophobicity+length).
3. **base-anchored expression residual**.
4. **base-anchored physchem+expression residual**.
5. **absence-safe expression/unexpressed gate** (biological, non-fit).
6. **promote-side mutant-RNA portfolio-reserve** (predeclared; the key new family). Protect the top
   `20−q` presentation/PRIME slots; allocate `q` of the 20 to the candidates with the strongest
   mutant-allele RNA evidence (IMPROVE `rna_af`/`rna_var`/`ValMutRNACoef` ↔ Sid `variant_vafs_long`
   tumor-RNA mutant reads/VAF — a real comparable bridge, 70 pts / 467 pos on IMPROVE), deduplicating
   overlap with the protected lane and backfilling any unused reserve by PRIME. `q ∈ {0 (incumbent), 1, 2,
   3, 4}` chosen ONLY by nested patient-grouped CV inside non-Sid IMPROVE, with per-cancer-cohort
   (bladder/melanoma/Basket) transport + matched-random-reserve controls and a no-catastrophic-regression
   requirement. This replaces the previously-falsified DEMOTION gate with the PROMOTE-side prior its own
   audit recommended.

7. **IMPROVE-rich anchored logistic (partial-bridge)** — patient-balanced official-partition OOF logistic
   on within-patient percentiles of {PRIME, EL, Expression, rna_var, rna_af, VarAlFreq, CelPrev};
   `score = prime_pct + α·OOF_pred_pct`. IMPROVE-only rich features (VarAlFreq/CelPrev/RNA); Sid mapping is
   a partial bridge (missing CelPrev → neutral). Per-cancer transport (bladder/melanoma/Basket).
8. **CORE-DEPLOYABLE anchored logistic** — same anchored form on {PRIME, EL/Mix, Expression, DNA VAF}
   (NO CelPrev, NO mutant-RNA). Maps cleanly to Sid (prime/mix/gene-TPM/WES-VAF), Gartner
   (expr_decile/vaf_decile) and multimer (VAF neutral-abstain). Transport tested on Gartner **and**
   multimer, not just within IMPROVE. This is the preferred deployable arm if it wins.

α ∈ {0, 0.10, 0.20, 0.25, 0.30, 0.50}; logistic C ∈ {0.5, 1, 2}; q ∈ {0,1,2,3,4}; base ∈ {PRIME, MixMHCpred}.
The full (arm × α × C × base) grid is the declared multiplicity; the count is reported and ONE winner
(or the null) is chosen by §6 on non-Sid data only, with bootstrap CIs over patients.

**Non-Sid preliminary (author-reported, no Sid consulted; to be reproduced under §6):** the pure q-reserve
family FAILS (nested folds select q=0 — rejected). The core-deployable arm is strongest: OOF hits (70 pts)
α 0=83 / .25=87 / .30=87, with α=.25 giving Basket +1, bladder +1, melanoma +2 and all cohort means ≥0.
The IMPROVE-rich arm is +3 at α=.25 (86). These are reproduced independently in Stage 1 before any freeze.

**Disclosure.** Sid has been inspected in prior turns of this session (its per-mutation ranks are known to
the author). This single-shot evaluation is therefore **exploratory confirmation, NOT pristine external
validation** — the firewall (§0) prevents Sid from informing any selection, but Sid is not a naïve
holdout. A genuinely pristine test requires an untouched external rich cohort.

## 5. Splits (grouping to prevent leakage)

Leave-one-STUDY-out (cohort) OUTER across {IMPROVE, multimer, Gartner}; within the training studies,
patient-group inner folds for α/feature selection. Exact + near-peptide (≥0.8) groups kept together;
recurrent peptides quarantined. Study identity is never a model feature.

## 6. Objective + freeze criterion (non-Sid only)

Clinically-aligned objective: total recognized **hits@20** after gate→unchanged base ranker. A config is
(arm, cols, base, C, α, q). Every config is scored by **true nested evaluation**: outer = the 5 official
patient-disjoint Partitions; the anchored logistic is fit OOF (train outer-4, predict held-out), and the
q-reserve/combined policy is applied on the held-out score — so the reported total is out-of-fold. Inner
patient-grouped CV on the outer-train selects C/α/q for the nested-CV validation number (reported to show
the SELECTION PROCEDURE generalizes and is not overfit to the multiplicity).

Freeze rule (registered; NO significance/CI gate): among ALL declared configs, select the one with the
largest nested-outer hits total that satisfies **(i)** positive Δ vs the null (α=0, q=0), **(ii)** no
catastrophic regression (worst cancer-cohort Δ ≥ −0.10), **(iii)** transport (every cancer-cohort Δ ≥ 0),
and **(iv)** beats the matched-random reserve control (mean over ≥20 seeds). If no config qualifies, freeze
the null. The bootstrap patient CI is REPORTED as a disclosed limitation, never used as a gate. Serialize
the winner + SHA-256 to `configs/frozen/sid_recognition_gate_v1.json` BEFORE touching Sid.

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
