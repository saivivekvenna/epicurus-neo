# Risk-controlled negative reducer — PREREGISTRATION CONTRACT (frozen before any experiment)

_Written and committed **before** loading a single feature/label. Sid and Miller are LOCKED and are not
read anywhere in the design, dev, or freeze of this task. This document is the binding protocol; the runner
and tests implement it exactly. Any deviation must be recorded as an explicit "PROTOCOL CORRECTION" section
here, committed separately, before it takes effect._

## 0. Motivation and honesty framing (Requirement 9)

The prior recognition-transfer gate optimized a **positive score** plus one crude RNA reserve slot and did
not transfer (tie at 2/3 on Sid). Oracle experiments in this milestone showed patient top-20 rises sharply
when **negatives are removed while positives survive**. This task therefore builds a **RISK-CONTROLLED
NEGATIVE REDUCER**: a label-blind model that confidently *removes* tested-negatives at guaranteed high
positive retention, hands survivors unchanged to genuine PRIME, and never reranks toward positives.

**Disclosure:** this hypothesis was motivated *after* seeing the Sid failure. Therefore even a later Sid
replay is exploratory. Only **Miller IPV (PRJNA980652, LOCKED_TEST)** can serve as pristine external
validation, and it is not touched in this task. Sid remains **LOCKED_DESCRIPTIVE**.

## 1. Data allocation (roles are fixed and named honestly)

| cohort | role | patients | positives | negatives | features present |
|---|---|--:|--:|--:|---|
| IMPROVE | **consumed DEV** | 70 | 467 | 17,053 | portable {PRIME, EL, expr} + rich {WES-VAF, RNA-VAF/coef, CelPrev} |
| Gartner NCI | **consumed DEV** | 26 | 46 | 3,722 | portable only |
| CD8 multimer | **consumed DEV** | 26 | 34 | 8,069 | portable only |
| **DEV total** | | **122** | **547** | **28,844** | |
| Miller IPV | **LOCKED_TEST** | ~13 | — | — | not read |
| Sid / osteosarc | **LOCKED_DESCRIPTIVE** | 1 | — | — | not read |

These three DEV cohorts are **consumed** — they are used for model/threshold/family selection and are never
again described as "untouched" or "held-out external". The only genuinely external validation targets are
Miller (locked here) and, exploratorily, Sid (locked here).

## 2. Primary target and the retention guarantee (Requirement 2)

**Primary target:** maximize the fraction of TESTED_NEGATIVES the gate confidently removes, **subject to a
high-positive-retention guarantee**, NOT AUROC. AUROC/AUPRC may be reported as diagnostics only; selection
never optimizes them.

**Retention guarantee.** Retention = (positives surviving the gate) / (positives eligible). The gate is only
allowed to *claim* removal on a fold when the **one-sided Clopper–Pearson 95% lower bound** on retention is
**≥ 0.95**. Power reality (exact, `0.05^(1/n)` at 100% retention):

| stratum | positives n | CP-95% LB at 100% retention | CP-claim powered? |
|---|--:|--:|:--:|
| DEV aggregate | 547 | 0.9945 | **yes** |
| IMPROVE | 467 | 0.9936 | **yes** |
| Gartner | 46 | 0.9370 | no |
| multimer | 34 | 0.9157 | no |
| (min powered n) | 59 | 0.9505 | threshold |

A CP-95%-LB≥0.95 retention claim requires **n ≥ 59 positives**, so it is **only powered at the DEV aggregate
and on IMPROVE**. On Gartner and multimer (n=46, 34) the CP claim is **impossible even at 100% retention** →
those folds fall to the abstain clause: **report raw retention, require ZERO positive loss (catastrophic
threshold below), and make no CP-backed removal claim**. Raw retention is always reported for every stratum.

## 3. Dynamic, label-blind inference (Requirement 3)

At inference the gate is a **recognition-risk model + calibrated removal threshold**:

1. Model outputs a per-candidate **keep-score** (higher = more likely a genuine recognizable positive).
2. Candidates with keep-score **below the calibrated threshold τ are REMOVED** (predicted confident
   negatives). All others are KEPT.
3. **Abstention (no removal) is triggered by** (a) patient/study **OOD**: the patient's feature distribution
   is out-of-support vs the training studies (per-feature min/max envelope + a coverage check); or (b)
   **missing features**: any model input absent/NaN for that candidate or patient → that candidate is KEPT
   (never removed on missing evidence). Missing evidence always defaults to KEEP.
4. **Top-20 PRIME safety lane:** the 20 best genuine-PRIME candidates per patient are **never removable**,
   regardless of keep-score. The reducer may only prune outside that lane.
5. After removal, **genuine PRIME ranks the survivors**. If a patient has < 20 survivors, **backfill** the
   top-20 by re-admitting the highest-PRIME removed candidates until 20 (or the pool) are present — removal
   never shrinks a patient below 20 rankable candidates.

The gate is a *pre-ranker denominator reducer*. It does not reorder toward positives; PRIME ranks whatever
survives.

## 4. Bounded, predeclared families and models (Requirement 4)

**Feature families (fixed; no post-hoc additions):**
- **portable** = within-patient oriented percentiles of {PRIME (lower better), EL (lower better),
  expression (higher better)}. Available in all three studies.
- **rich** = portable + {WES tumor VAF (`VarAlFreq`), mutant-RNA VAF (`rna_af`), mutant-RNA coefficient
  (`ValMutRNACoef`), clonality proxy `CelPrev`} — present **only in IMPROVE**. Under leave-one-study-out the
  rich family is **structurally non-transferable** (train studies lack the columns when IMPROVE is held out;
  test studies lack them when IMPROVE is train), so it can be reported **only as an IMPROVE-internal
  patient-CV diagnostic** and is **NOT eligible** for the transferable freeze. This is stated up front so
  rich is not silently promoted.

**Models (fixed):** (a) L2-regularized **logistic regression** (interpretable, monotone in oriented
percentiles by construction of the balanced fit); (b) **monotonic shallow gradient-boosted trees** (max
depth ≤ 3, monotone constraints in the recognition-favoring direction) *iff the installed sklearn exposes
`HistGradientBoostingClassifier` with `monotonic_cst`*; otherwise this model is recorded as UNAVAILABLE and
skipped, not substituted. (c) **NULL** = no gate (remove nothing). No other model families may be added
after seeing results.

**Threshold grid** (removal aggressiveness, calibrated on inner CV, not free-tuned on outer): keep-score
quantile cut `τ ∈ {none, and the largest cut whose inner CP-95%-LB retention ≥ 0.95}` — i.e. τ is chosen as
the **most aggressive inner cut that still guarantees inner retention**, per model, per outer-train.

## 5. Truly nested evaluation (Requirement 5)

- **Outer:** leave-one-**study**-out across {IMPROVE, Gartner, multimer} → 3 outer folds. The held-out study
  is the outer test; the model, hyperparameters, family, AND removal threshold are selected **without it**.
- **Inner:** patient/group CV **within the outer-train studies** for model + threshold selection. Inner CP
  retention is computed on inner-validation positives.
- **NULL / no-gate is always a candidate.** Conservative tie-break at equal inner objective prefers
  **NULL > portable > simpler model (logistic before trees) > less aggressive threshold**.
- **Reported per outer study:** raw positive retention, CP-95%-LB retention (or "underpowered/abstain"),
  negative-removal fraction, **patient-macro hits@20 delta** (primary utility), pooled hits@20 (secondary),
  per-positive rank changes, **patient-level paired-bootstrap CI** (fixed seed) of the hits@20 delta, a
  **matched-random removal control** (remove the same count of non-safe-lane candidates uniformly at random,
  averaged over seeds), and the **worst-study** result.
- **Eligibility for a non-null freeze (ALL must hold):** (i) **every** outer study has **zero catastrophic
  hit loss** (per-study patient-macro Δhits@20 ≥ `CATASTROPHIC = -0.02`); (ii) **aggregate** patient-macro
  Δhits@20 > 0 **and beats the matched-random control**; (iii) on the powered strata (aggregate/IMPROVE) the
  retention guarantee CP-95%-LB ≥ 0.95 holds; (iv) negative-removal fraction > 0 (the gate actually does
  something). If any fail → **freeze NULL** and say so explicitly.

Primary utility is **patient-macro hits@20**; pooled counts are secondary.

## 6. Leakage & determinism (Requirement 6)

- **Peptide quarantine recomputed inside each split.** A held-out candidate whose mutant peptide
  exact- or near-(≥0.8 normalized similarity)-matches any training-side peptide is **quarantined from
  hit-counting only** (it still competes for the 20 slots; the pool is never shrunk). The quarantine mask is
  **recomputed within each inner split's own train** so outer-test/inner-val peptides cannot leak into
  selection. Cross-study peptide matches are included in the quarantine.
- **No study-identity features.** Study/cohort membership is never a model input.
- **Stable seeds/hashes.** Per-patient seeds derive from `sha256(patient_id)` (process-independent), never
  Python `hash()`. All randomness (matched-random control, any tree subsampling) is seeded.

## 7. Freeze semantics & reproducibility (Requirement 7)

- If §5 eligibility passes, freeze a **fully fitted apply-only payload**: feature order, coefficients/
  intercept (full-precision JSON floats) or the serialized tree model, the family, the calibrated τ
  (refit/recalibrated on **all DEV** using the selected recipe), the OOD envelope, the percentile/direction
  policy, sklearn version, a model-payload SHA-256, a config SHA-256, and `stage2_must_not_refit`.
- If eligibility fails, freeze **NULL** (`{"frozen": "NULL", ...}`) and state the honest negative result.
- Record code commit, config/data-file SHA-256s. Assert **reproducibility across `PYTHONHASHSEED`** (run the
  selection twice in separate processes; frozen SHA must be identical).

## 8. Later locked application — SAFEGUARDS (deferred; NOT run in this task) (Requirement 8)

Before the frozen gate is ever applied to Miller or Sid, a **separate future task** must first commit:
standalone **apply-only** code (AST/import guard proving no refit / no sklearn fit call), a **hard atomic
one-shot sentinel** that refuses if a result already exists and writes an execution receipt before reading
any locked input, **exhaustive small-N brute-force tie tests** (score-boundary and any reserve/quarantine
ties), and an **explicit human approval checkpoint**. **This task does NOT run Sid or Miller.**

## 9. Explicit disclosure (Requirement 9)

Restated: hypothesis motivated post-Sid-failure → later Sid replay is exploratory; only Miller can validate
externally; both are locked here. The three DEV cohorts are consumed, not pristine. A NULL result is a
fully acceptable, reportable outcome and will be frozen honestly if the evidence does not pass §5.

---

### Declared multiplicity (fixed enumeration)
outer folds = 3 studies × {logistic, monotonic-trees(if available), NULL} × inner-calibrated τ, portable
family only for the transferable freeze; rich = IMPROVE-internal diagnostic only. Aggregate matched-random
control seeded over 20 draws. Nothing outside this enumeration is selected.
