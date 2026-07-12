# Risk-controlled negative reducer — PREREGISTRATION CONTRACT (frozen before any experiment)

_Written and committed **before** loading a single feature/label. Sid and Miller are LOCKED and are not
read anywhere in the design, dev, or freeze of this task. This document is the binding protocol; the runner
and tests implement it exactly. Any deviation must be recorded as an explicit "PROTOCOL CORRECTION" section
here, committed separately, before it takes effect._

## PROTOCOL CORRECTION 2 (recorded and committed BEFORE the runner reads any cohort data — audit fixes)

Audit of the (unrun) runner surfaced leakage/serialization issues; fixed here before any data read:

1. **No in-fold weight leakage.** Inner-OOF fold models are fit with weights recomputed on that fold's own
   train subset (`balanced_weights(train.iloc[tr])`), never with weights derived from the full outer-train
   (which would use validation-fold class totals).
2. **τ from OUT-OF-FOLD scores, never in-sample.** Both (a) the outer-test application and (b) the final
   full-DEV freeze calibrate τ on **patient-grouped OOF scores** for the selected (model, C, m); the full
   model is then refit only to *score* the untouched target (outer-test / — for freeze — to serialize
   coefficients). In-sample τ calibration (optimistic retention) is prohibited.
3. **HGB is EVALUATION-ONLY and INELIGIBLE for the frozen payload.** A monotonic HGB cannot be serialized to
   applyable JSON coefficients under this task's apply-only constraint, so it is reported as a *diagnostic*
   (does a monotone tree ensemble beat the nonnegative logistic in nested LOSO?) but is **excluded from
   `inner_select`'s selectable/freezable set** — the deployable recipe space is **{NULL, nonnegative-
   logistic(C-grid)} × m**. This prevents freezing an unusable payload. (If a future task needs trees, it
   must first commit deterministic tree serialization + an apply-equivalence test.)
4. **Aggregate CP eligibility is enforced.** §5.iii is applied to BOTH the pooled DEV aggregate
   (Σ n_pos, Σ pos_removed over all outer-test folds) AND the IMPROVE outer-test stratum: both CP-95%-LB
   retention values must be ≥ 0.95 and both are recorded.
5. **Deployment recipe is chosen by the SAME preregistered `inner_select` run once on ALL DEV** (not a
   post-hoc majority vote across outer folds); the full-DEV candidate table is recorded, then τ is
   OOF-calibrated and the model refit on all DEV for the payload.
6. **Provenance.** The runner records the git commit, SHA-256 of every data file actually read, config +
   model-payload SHA-256, sklearn/scipy versions, and a cross-`PYTHONHASHSEED` reproducibility check
   (identical frozen SHA in a second process). Per-fold **full inner candidate tables** are stored in JSON.
7. **`fit_nnlogistic` fails closed.** If the optimizer does not converge or returns non-finite parameters, it
   falls back to a constant keep-score (zero coefficients) ⇒ the gate removes nothing (KEEP-all), never
   garbage removals.

## PROTOCOL CORRECTION 1 (recorded and committed BEFORE any code/experiment — supersedes §3.4, §3.3, §4, §5)

**Fatal flaw fixed.** Protecting the *full* original PRIME top-20 as an unremovable safety lane forces the
final PRIME top-20 to equal the original PRIME top-20, so **patient-macro Δhits@20 is mathematically pinned
to 0** and a 3/3 improvement is impossible. The safety lane is therefore **removed** and replaced by:

1. **Bounded protected-core grid `m ∈ {0, 5, 10}`.** The top-`m` genuine-PRIME candidates per patient are
   unremovable; everything ranked below the core competes for removal. `m` is a hyperparameter **selected
   inside nested CV**, with a **conservative tie-break preferring the LARGER `m`** at equal inner objective.
   **Exact-PRIME-score ties at the `m` boundary are all protected** (the core may thus be larger than `m`).
   With `m < 20`, removing a negative ranked between the core and rank 20 can promote a below-20 positive
   into the top-20, so **Δhits@20 is now live** and a genuine gain is achievable.
2. **CP-95% retention is the real safety mechanism** (not the core). The core only bounds how much of the
   PRIME head is exempt from pruning; positive protection comes from the retention-guaranteed threshold.
3. **Matched-random control uses the SAME `m`/protected-core and the SAME removal count** as the gate, and
   the **full pool always competes** for the 20 slots.

**Monotonicity honesty (supersedes §4).** The earlier claim that ordinary L2 logistic is "monotone by
construction" because inputs are oriented is **false** — orientation does not constrain coefficient signs.
Corrected: the primary model is a **nonnegative-coefficient-constrained logistic** (coefficients bounded
`≥ 0` via L-BFGS-B on the balanced, L2-penalized weighted log-loss; intercept free), which *does* make the
keep-score monotone nondecreasing in each recognition-favoring percentile — implemented with a fixed method
and unit tests asserting `coef ≥ 0` and equivalence of the linear→sigmoid apply path. The monotonic shallow
HGB (if available) keeps its `monotonic_cst = +1` constraints. If the nonnegative logistic proves not
robustly implementable it falls back to **unconstrained** logistic, **disclosed as non-monotone** (never
described as monotone). No "monotone by construction" claim is made for any unconstrained fit.

**Small-study abstention (supersedes §3.3 wording).** Two distinct abstentions: (a) **CP-claim abstention** —
on underpowered strata (Gartner n=46, multimer n=34, where CP-95%-LB≥0.95 is impossible even at 100%
retention) the gate **still applies and still removes**; we only **abstain from the CP-backed retention
*claim*** there and report raw retention + require zero catastrophic hits@20 loss. (b) **OOD inference
abstention** — a genuinely out-of-support patient/feature causes **no removal** (KEEP). These are not the
same: small n ⇒ abstain from the *claim*, not from *removal*; OOD ⇒ abstain from *removal*.

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
3. **OOD / missing-feature abstention (no removal → KEEP).** (a) patient/study **OOD**: the patient's
   feature distribution is out-of-support vs the training studies (per-feature min/max envelope + coverage
   check) → KEEP all of that patient. (b) **missing features**: any model input absent/NaN for a candidate →
   that candidate is KEPT. Missing evidence always defaults to KEEP. (This is distinct from CP-claim
   abstention on small studies — see PROTOCOL CORRECTION 1: there the gate still removes.)
4. **Protected PRIME core `m ∈ {0,5,10}` (per PROTOCOL CORRECTION 1, replacing the top-20 safety lane).**
   The top-`m` genuine-PRIME candidates per patient are unremovable; `m` is nested-CV-selected (tie-break
   prefers larger `m`); exact-PRIME-score ties at the `m` boundary are all protected (core may exceed `m`).
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

**Models (fixed; see PROTOCOL CORRECTION 1 for the monotonicity honesty fix):** (a) **nonnegative-
coefficient-constrained logistic** — coefficients bounded `≥ 0` (L-BFGS-B on balanced, L2-penalized weighted
log-loss; intercept free), which makes the keep-score monotone nondecreasing in each recognition-favoring
percentile; falls back to **unconstrained logistic disclosed as non-monotone** only if the constrained fit
is not robustly implementable. (b) **monotonic shallow gradient-boosted trees** (max depth ≤ 3,
`monotonic_cst = +1`) *iff sklearn exposes `HistGradientBoostingClassifier` with `monotonic_cst`* (verified
present: sklearn 1.9.0); else recorded UNAVAILABLE and skipped, not substituted. (c) **NULL** = no gate. No
other model families may be added after seeing results.

**Protected-core + threshold grid** (calibrated on inner CV, not free-tuned on outer): protected PRIME core
`m ∈ {0, 5, 10}`; keep-score cut `τ` chosen as the **most aggressive inner cut whose inner CP-95%-LB
retention ≥ 0.95** (τ = none is always a candidate), per model, per `m`, per outer-train. Conservative
tie-break at equal inner objective prefers **NULL > larger m > simpler model (logistic before trees) > less
aggressive τ**.

## 5. Truly nested evaluation (Requirement 5)

- **Outer:** leave-one-**study**-out across {IMPROVE, Gartner, multimer} → 3 outer folds. The held-out study
  is the outer test; the model, hyperparameters, family, AND removal threshold are selected **without it**.
- **Inner:** patient/group CV **within the outer-train studies** selects **(model, protected-core `m`,
  τ)**. τ is the most aggressive inner cut whose pooled inner CP-95%-LB retention ≥ 0.95; among (model, `m`)
  the **inner objective is inner-OOF patient-macro Δhits@20** (now non-degenerate under the core grid).
- **NULL / no-gate is always a candidate.** Conservative tie-break at equal inner objective prefers
  **NULL > larger `m` > simpler model (logistic before trees) > less aggressive τ**.
- **Reported per outer study:** raw positive retention, CP-95%-LB retention (or "CP-underpowered → claim
  abstained; gate still applied"), negative-removal fraction, **patient-macro hits@20 delta** (primary
  utility), pooled hits@20 (secondary), per-positive rank changes, **patient-level paired-bootstrap CI**
  (fixed seed) of the hits@20 delta, a **matched-random removal control** (same selected `m`/protected-core,
  same removal count, uniform random over non-core candidates, averaged over seeds), and the **worst-study**
  result.
- **Eligibility for a non-null freeze (ALL must hold):** (i) **every** outer study has **zero catastrophic
  hit loss** (per-study patient-macro Δhits@20 ≥ `CATASTROPHIC = -0.02`); (ii) **aggregate** patient-macro
  Δhits@20 > 0 **and beats the matched-random control**; (iii) on the powered strata (aggregate/IMPROVE) the
  retention guarantee CP-95%-LB ≥ 0.95 holds (small studies abstain from the CP *claim* but the gate still
  applied); (iv) negative-removal fraction > 0. If any fail → **freeze NULL** and say so explicitly.

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
outer folds = 3 studies; inner selection over {nonnegative-logistic (C grid), monotonic-trees (if
available), NULL} × protected-core `m ∈ {0,5,10}` × inner-CP-calibrated τ; portable family only for the
transferable freeze; rich = IMPROVE-internal diagnostic only. Matched-random control seeded over 20 draws at
the selected `m`/count. Nothing outside this enumeration is selected.
