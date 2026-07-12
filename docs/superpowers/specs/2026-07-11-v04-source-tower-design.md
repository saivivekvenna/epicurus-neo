# Epicurus v0.4 — source-aware tower: design spec

_Date: 2026-07-11. Status: design approved (partial pooling), pending spec sign-off before implementation._

Development-only experiment on the frozen `configs/frozen/mil_dev_split_v1.json`. The Gartner TEST holdout is
never loaded/scored. The frozen split is used verbatim. No external-superiority claim is produced. Nothing is
committed by the experiment; v0.1 stays the frozen model of record unless v0.4 clears the registered gate.

---

## 1. Motivation and the registered hypothesis

The v0.3 MIL ranker ties genuine PRIME on the source-balanced OOF metric (Δhits@20 = −0.093, CI[−0.228,+0.051]).
Its diagnostics show the tie is not a power problem — it is a **pooling** problem:

- Source-only OOF models are much stronger where naive pooling dilutes them: Gartner source-only **1.30** vs
  pooled 1.025 (PRIME 0.975 → **+0.325 latent**), IMPROVE source-only 1.167 vs pooled 0.967.
- Multimer goes the **other** way — it *prefers* pooling: pooled 1.167 > source-only 1.056.
- Study-shortcut source-identity AUROC = **0.898**; positive prevalence differs ~100× (gartner 0.0003 /
  improve 0.027 / multimer 0.004). One shared decision boundary cannot fit all three regimes.

These source-only numbers are honest OOF (training is restricted to one source; the held-out fold is still
scored), so the latent Gartner edge is real, not in-sample.

**Registered hypothesis (narrowed, per review):**

> Does **source-conditioned feature weighting** recover ranking signal erased by naive pooling, **beyond what
> source-prevalence calibration alone can explain**?

This is a *mechanism test*. v0.4 is **not** proposed as the deployable universal Epicurus ranker (its heads are
indexed by dataset name; see §8). It asks only whether known source heterogeneity is suppressing development
performance and, if so, whether the recoverable part is feature-weight heterogeneity rather than base-rate
calibration.

## 2. Model family (three nested members; one is the gated candidate)

All three share v0.3's linear instance scorer, MIL log-sum-exp bag aggregation (temperature τ), fail-closed core
features, orthogonal-feature masks, PRIME-leakage guard, and the source×patient×bag-balanced training objective.
They differ only in how source identity enters the scorer.

Let `x` be the standardized feature vector (§4), `s(i)` the source of instance `i`.

| member | instance scorer | free source params | role |
|---|---|---|---|
| **P — pooled** (= v0.3) | `f(x) = w₀·x + b₀` | none | **frozen corrected-v0.3 OOF**, loaded not retuned (§5.6) |
| **C — calibration tower** | `f(x) = w₀·x + b₀ + c_{s(i)}` | per-source intercept `c_s` | isolates prevalence-calibration; **independently nested-selected** |
| **F — feature tower** *(gated candidate)* | `f(x) = (w₀ + v_{s(i)})·x + b₀ + c_{s(i)}` | per-source head `v_s` + intercept `c_s` | source-conditioned feature weighting; independently nested-selected |

**Why C is the right middle comparator.** A patient belongs to exactly one source, so a per-source intercept
`c_s` shifts that patient's entire candidate pool equally and therefore **cannot change within-patient ranking
(hits@20) at all**. Its only effect is during fitting: it absorbs each source's base rate so the shared weights
`w₀` are not distorted by trying to fit three different prevalences at once. Thus any ranking lift of **C over P**
comes purely from *freeing `w₀`*, and the lift of **F over C** is exactly the registered "feature weighting
beyond calibration" signal.

## 3. Objective, shrinkage, and optimization

### 3.1 Penalized loss

Minimize, over `(w₀, b₀, {v_s}, {c_s})`:

```
L(θ) = L_MIL(θ)  +  α‖w₀‖²  +  α·λ·Σ_s ‖v_s‖²  +  α_b·λ_b·Σ_s c_s²
```

- `L_MIL` is the source-balanced bag negative-log-likelihood (§3.3), identical in form to v0.3.
- `α = 1/(C·n_bags)` — the same complexity control as v0.3 (`C` = inverse regularization strength).
- **`λ` controls *relative* pooling of the feature heads**: large λ → heads shrink to 0 (pooled); small λ →
  heads free (toward independent). λ is defined on **standardized** features so its scale is comparable across
  features and across folds.
- `c_s` are regularized by a fixed `α_b·λ_b := α` (the same ridge strength as the shared weights `w₀`, **not**
  tuned) — they are nuisance/calibration parameters, kept from ballooning and reported separately from
  feature-head effects.

### 3.2 Identifiability (no fake λ=0 claim)

The map `(w₀,{v_s}) → {w_s = w₀+v_s}` is invariant under `w₀→w₀+δ, v_s→v_s−δ`, but the ridge penalties
`α‖w₀‖² + αλΣ‖v_s‖²` are **not** — they select a unique split for every finite `λ > 0`. So:

- Predictions (functions of the **effective** weights `w_s = w₀ + v_s`) are well-defined for all finite λ.
- We **report `w_s` and the deviation norms `‖v_s‖`**, never an identifiability-ambiguous `w₀`/`v_s` split as a
  scientific quantity.
- **`λ = ∞`** is implemented as an **explicit pooled branch** (`v_s = c_s = 0`), not a large float, and is
  verified numerically equal to frozen v0.3 (§7 test 1).
- **Independent per-source heads** are a **separate branch** (each source fit alone), used only for the
  source-only diagnostic and the "opposite feature weight" synthetic test — **not** a gated config and **not**
  claimed to be the `λ→0` limit.

### 3.3 Source-balanced objective (patients as the unit)

`L_MIL = Σ_bags w_b · bce(y_b, σ(lse_b))`, with per-bag weight `w_b = 1/(S · n_patients_in_source(s) ·
n_bags_of_patient(p))`. Each patient's bags sum to `1/(S·N_s)`; each **source** contributes total weight `1/S`;
each **patient** within a source contributes equally regardless of how many candidate rows or bags it has. This
is v0.3's `_bag_targets` weighting, preserved unchanged, so the shared backbone `w₀` cannot be dominated by the
largest source. (Gartner negatives are counted at the **bag** level, never the ~285k child instances.)

### 3.4 Optimization — deterministic, not asserted-convex

`lse_b = τ·log(mean_i exp(f(x_i)/τ))` is convex in the scores, and the **negative-bag** term `softplus(+lse)` is
convex; the **positive-bag** term `softplus(−lse)` is convex-non-increasing ∘ convex and is **not** guaranteed
convex. We therefore **do not** claim a convex objective. We preregister:

> deterministic smooth optimization (L-BFGS-B) from a **fixed zero initialization**, `maxiter=500`,
> `ftol=1e-9`; the analytic gradient is supplied.

**Multi-init robustness is a diagnostic, not a tuning dimension:** each outer-fold candidate is additionally
refit from ≥2 further fixed-seed initializations (small Gaussian perturbations); we report that the OOF
predictions are materially identical — pre-registered pass threshold: Spearman rank-correlation ≥ 0.999 **and**
max |Δ standardized score| ≤ 1e-3 across every pair of inits. If not met, the result is flagged unstable.
Initialization is never selected to improve the metric.

## 4. Features, scaling, and leakage (unchanged from v0.3, made explicit)

- Same CORE (`f_prime_pct, f_mix_pct, f_el_pct, f_pres_abs`) + ORTHO (`f_expr, f_agreto, f_foreign, f_bindstab,
  f_proc`), same source-correct EL orientation (Gartner `Score_EL` higher-better vs IMPROVE/multimer %rank
  lower-better), same fail-closed drop of any unit missing a core score, same masked-neutral orthogonals, same
  PRIME-leakage neutralization.
- **All fitted preprocessing is fit inside outer-train only.** Within-patient percentiles are patient-local
  (depend only on that patient's own pool → no cross-fold leakage). `el_strength` is row-local. The only
  train-fitted transform is **feature standardization** (z-score), whose mean/std come from the fold's training
  rows exclusively; held-out rows are transformed with the stored stats. λ is defined on these standardized
  features. A test asserts standardization stats are a function of train rows only (§7 test 7).

## 5. Evaluation, nested selection, and the gate

### 5.1 Nested OOF (unchanged frame)

Outer loop = the 5 frozen folds; train on 4, score the held-out fold; concatenate → one OOF hits@20 per patient.
Bootstrap = source-balanced, patient-paired, 2000 iterations, 95% percentile CI on the weighted-mean Δ. No
patient is ever scored by a model that saw it (patient-blocked folds).

### 5.2 Pre-registered grid

Standardized features make λ comparable, so:

- `λ ∈ {0.03, 0.1, 0.3, 1, 3, 10, 30, ∞}`  (∞ = explicit pooled branch)
- `C ∈ {0.1, 0.3, 1.0}` , `τ ∈ {0.5, 1.0}`  (unchanged from v0.3)
- `α_b·λ_b` fixed (intercept regularization; not tuned).

### 5.3 Nested selection rule (frozen, deterministic)

Inside outer-train, by inner grouped OOF on the **same source-balanced patient-level hits@20** used as the north
star, select the configuration by this lexicographic rule:

1. **Maximize** inner source-balanced mean hits@20.
2. Among configs within a preregistered tolerance `ε = 0.01` weighted hits@20 of the best, **prefer the more
   pooled** model (larger λ; `∞` beats any finite).
3. Then **prefer stronger regularization** (smaller `C`).
4. Then a deterministic lexical tie-break on `(λ, C, τ)`.

This biases toward pooling and simplicity, so noisy inner folds cannot conjure over-specialized heads.

### 5.4 Registered gate (unchanged from v0.3)

The **gated candidate is F (feature tower)**. ACCEPT iff, on the OOF result:

1. `Δ vs genuine PRIME`: CI_lower > 0, **and**
2. `Δ vs strongest presentation`: no regression (point ≥ 0 and 95% CI not entirely below 0).

Otherwise **REJECT** = "failed to establish superiority" (not "worse"). C and P are context/diagnostic, never a
second gated candidate (avoids multiplicity). ACCEPT licenses only *one* pre-registered look at the
SEMI-CONSUMED Gartner TEST; an external claim still requires an untouched cohort.

### 5.5 Mechanism contrast (the registered hypothesis, separate from the gate)

Report, with source-balanced patient-paired CIs:

- **F − P** (does source-conditioning help at all vs naive pooling?)
- **C − P** (how much is pure prevalence calibration?)
- **F − C** (**the registered mechanism signal**: feature-weight heterogeneity beyond calibration).

Pre-registered readings: **(a)** F ACCEPTs the PRIME gate; **(b)** F ties PRIME but F−P CI_lo>0 and F−C CI_lo>0
= mechanism confirmed, scientific progress, still a gate REJECT; **(c)** F−C CI spans 0 = heterogeneity is mostly
calibration/none; **(d)** no lift. No outcome is re-spun after the fact.

### 5.6 Comparator fairness (frozen P; honest best members)

- **P is the exact frozen corrected-v0.3 rung-3 OOF**: per-fold `(C, τ)` + coefficients are **loaded from
  `epicurus_v03/DEV_RESULT.json` (`ladder.rung3_MIL.folds`) and used verbatim — never re-tuned**. The λ=∞ tower
  branch is verified numerically equal to it (§7 test 1) but P's *reported* numbers come from the frozen configs.
- **C and F are each independently nested-selected** inside outer-train (C picks its own `(C,τ)`; F its own
  `(C,τ,λ)`). F−C therefore compares two honest best members, not F's hyperparameters borrowed by C.
- **F−P, F−C, C−P are paired on the identical set of eligible patients** (§5.7), with the identical per-patient
  candidate pools used for every model and baseline.

### 5.7 Eligibility and attrition (label-blind; explains 152→118)

- **Eligibility is label-blind:** a patient is *rankable* iff, after fail-closed core-feature/mask filtering, it
  has ≥1 candidate with all core scores present. This uses **no outcome labels** (a test asserts the eligibility
  function never references `eval_positive`/`bag_label`).
- **Scored ⊂ eligible:** hits@20 additionally requires ≥1 in-pool positive (intrinsic to a capture metric — you
  cannot measure capture with no target). The 152→118 gap is reported as an explicit by-source attrition ladder:
  feature-bearing patients → rankable (label-blind) → scored (has in-pool positive), plus core-feature/mask
  availability rates.
- Every baseline and model is scored on **exactly the same per-patient candidate pool** (asserted).

## 6. Fallback for unseen sources (defined now)

At score time: **known source → `w₀ + v_s`; unseen source → shared backbone `w₀`** (heads default off). This is
tested (§7 test 4). It is a stopgap: the intended future direction is to replace dataset-name heads with
**assay/regime heads** keyed on observable properties (EL-score semantics, candidate-generation protocol, assay
type, peptide-length regime) so the mechanism transfers to new cohorts. That generalization is explicitly **out
of scope** for v0.4 and noted as future work.

## 7. Tests (superset of v0.3; the review's required additions included)

1. **λ=∞ ≡ pooled v0.3** — the pooled branch's coefficients and OOF predictions equal frozen `MILRanker` at the
   same `(C, τ)` to optimizer tolerance.
2. **Source-label-ordering invariance** — permuting the integer source encoding leaves OOF predictions unchanged.
3. **Equal source contribution** — objective weights sum to `1/S` per source despite unequal patient/candidate
   counts (assert on a constructed unequal frame).
4. **Unseen-source fallback** — a row with an unknown source scores with `w₀` only (equals the pooled score).
5. **Multi-init stability** — ≥3 fixed-seed initializations give materially identical OOF predictions (corr≈1).
6. **Recovers opposite per-source weight** — synthetic frame where source A needs `+feature` and source B needs
   `−feature`; the feature tower (small λ) fits opposite-sign effective weights and beats pooled P; the pooled
   model cannot separate them.
7. **Preprocessing fit inside train only** — standardization mean/std are a pure function of train rows; a change
   to held-out rows does not change the fitted scaler.
8. **Intercept-only effects are rank-inert** — a calibration tower with `v_s=0` produces identical within-patient
   hits@20 to a pooled model with matched `w₀` (per-source intercepts cannot masquerade as ranking lift).
9. Plus v0.3-style mechanics: MIL bag aggregation scores a witness high, bag-balanced (not instance-inflated)
   weights, OOF never scores a training patient, gate logic, determinism of the whole run.

## 8. Diagnostics (v0.3 §8 suite + additions)

Per-source hits/recall/AUROC vs each baseline; source-only vs augmented; study shortcut; leave-one-feature-out
ablation; score-orientation Spearman; family-leakage assertions (no patient/genomic-family/non-quarantined
peptide crosses a fold; 0 Gartner TEST patients); recurrent-peptide quarantine stratum (reported, never
selected). **Additions:** Δ vs frozen pooled v0.3; the C/F/P mechanism contrasts (§5.5); the **selected λ per
fold** and per-source deviation norms `‖v_s‖` and effective weights `w_s` (how much, and in what direction, each
source departs from shared); multi-init stability numbers.

## 9. Files

- `src/event_b/epicurus_v04.py` — `TowerMILRanker` (members P/C/F via flags + `λ=∞` branch), unseen-source
  fallback; imports `assemble_frame, run_oof, evaluate_model, per_patient_metrics, baseline_score,
  paired_bootstrap, _fit_weights, MILRanker` from `epicurus_v03`. v0.3 is **not modified** (its frozen negative
  stays immutable).
- `scripts/epicurus_v04_dev.py` — runner: nested C×τ×λ selection, gate, mechanism contrasts, full diagnostics →
  `artifacts/milestone_7_decision/epicurus_v04/{DEV_RESULT.json, DEV_REPORT.md}` and
  `configs/frozen/epicurus_v0_4_dev.json` (`status` per gate, `supersedes_frozen: false` unless ACCEPT).
- `artifacts/milestone_7_decision/epicurus_v04/PREREGISTERED_PROTOCOL.md` — the scientific pre-commitment
  (this spec's §1–§6 in protocol form), written and committed **before** any fit.
- `tests/test_epicurus_v04.py` — §7 tests.
- `NORTH_STAR_HISTORY.md` + memory updated as the final step.

## 10. Guardrails (unchanged)

Gartner TEST not read/scored/tuned/evaluated. Frozen split verbatim. No external claim; ACCEPT means only
"passed the registered *development* gate." PRIME/MixMHCpred binaries + training data are non-commercial and
gitignored. Nothing committed by the experiment.

## 11. Final review safeguards (pre-fit, minimal additions)

1. **Comparator fairness** — §5.6: P = frozen corrected-v0.3 (loaded, not retuned); C and F each independently
   nested-selected; F−P/F−C/C−P paired on identical eligible patients + pools.
2. **Provenance frozen before the first fit** — `PROVENANCE.json` records SHA256 of: the frozen split, every
   input data file (the three PRIME/Mix caches, IMPROVE zip, multimer corpus, crosswalk inputs), the feature
   code files (`epicurus_v03.py`, `epicurus_v04.py`, `nci_crosswalk.py`) + git HEAD, and the prereg protocol.
   The runner recomputes and **fails fast** on any mismatch. **Fail-fast TEST-absence invariant:** every Gartner
   dev patient must be contained in the TRAIN crosswalk patient set (structurally excludes TEST-only patients),
   asserted label-free. **I/O guard:** the run executes inside a context manager that patches `open` to raise if
   any known Gartner TEST path is opened — TEST is *never opened*, not merely filtered post-load.
3. **Attrition** — §5.7: label-blind eligibility/attrition by source + core-feature/mask availability; test that
   eligibility ignores labels; identical pools for all comparators.
4. **Source-label negative control (diagnostic only)** — refit F with a **deterministically shuffled** per-patient
   pseudo-source (same 3-way marginal), applied consistently at train and score time; report its F−P lift. A
   lift comparable to true-source F would indicate the gain is model *capacity*, not genuine source structure.
   **Never used for selection or gating.**
5. **hits@20 is the sole gate**; **best-positive rank, nDCG@20, recall@20** are reported as **non-selected**
   diagnostics so a coarse capture metric cannot hide whether the actual ordering improved.
6. **PRIME masking/availability** — report, per source, the PRIME-leakage **mask rate** and `prime_rank`
   availability. Preserve the conservative split: the **genuine PRIME comparator** uses raw unmasked `prime_rank`
   (harder bar); the **Epicurus PRIME *feature*** (`f_prime_pct`) is masked-to-neutral on near-training peptides.
7. **Guardrails preserved** — source-balanced paired bootstrap; per-source results; multi-init stability; no
   Gartner TEST; the source-**name** tower is **mechanism evidence only**, explicitly **not** deployable/
   generalizable; any ACCEPT is **development-only**.
8. **Commit discipline** — commit *only* this design + the preregistration + `PROVENANCE.json` before fitting,
   staging those paths explicitly (no sweep of unrelated dirty files). The v0.3 EL bug is recorded transparently
   in the protocol; v0.3's frozen corrected artifacts are **not modified**. Implementation, tests, and results
   are run but not committed (consistent with the standing no-commit guardrail).
