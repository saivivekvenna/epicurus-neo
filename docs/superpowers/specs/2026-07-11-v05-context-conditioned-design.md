# Epicurus v0.5 — deployable context-conditioned pairwise challenger: design spec

**Design source of truth for the v0.5 DEVELOPMENT experiment.** This is the pre-fit design that the
preregistration (`artifacts/milestone_7_decision/epicurus_v05/PREREGISTERED_PROTOCOL.md`) formalizes and
freezes. It fixes the hypothesis, the two new model members (Q, R), the frozen comparators (P, F, and the
descriptive additive rung A), the approved/rejected scoring-time contexts (with the pre-fit **source-alias**
audit that justifies each), the pairwise objective and its exact weights, the grid, the gate, and every
required diagnostic. Written before any Q/R model is fit or inspected. **Gartner TEST is not loaded, not scored,
not opened.** The frozen `mil_dev_split_v1` is used verbatim.

> **This file is the pre-fit CORRECTION of the original v0.5 preregistration (commit `06cb1d1`).** It fixes six
> scientifically necessary flaws found in review: (1) the fitted intercept is removed from Q/R because it is
> unidentifiable under a pairwise objective; (2) the frozen comparators P/F are **refit from their exact frozen
> code**, not "loaded" (their JSON does not persist the coefficients/scaler needed to score); (3) the Q−P
> contrast is stated to mix objective **and** exact-witness supervision, with the frozen additive rung A added
> as a descriptive comparator; (4) each Gartner negative bag is represented by a convex **log-mean-exp over all
> its children**, not one fixed representative; (5) the context feasibility gate is replaced by an **actual
> patient-grouped, source-balanced source-classification alias audit** (range overlap was not an alias test),
> which rejects every pool-enumeration context and keeps only candidate-level contexts; (6) all five prereg
> assets are updated consistently.

---

## 1. Where v0.4 left us, and the v0.5 hypothesis

v0.4's feature tower **F** (`(w₀+v_s)·x`) TIED genuine PRIME (Δhits@20 = −0.0157, CI[−0.147, 0.112]) with no
regression vs presentation, and **recovered the predicted Gartner edge** (+0.275 vs PRIME). Its negative control
was clean (shuffled source far worse), so the lift tracked genuine source structure, not capacity. But F is
**nondeployable**: its heads are indexed by dataset **name**, so it cannot score a brand-new cohort, and the
per-source heads harmed IMPROVE (−0.10) and multimer (−0.222), cancelling the Gartner gain. F was also
**nonconvex** — multi-init score wobble (max|Δstd-score| 0.0166) rode on a non-convex MIL log-sum-exp objective.

> **v0.5 registered hypothesis.** Does **portable, scoring-time context** (variables knowable while ranking a
> brand-new patient's candidate pool) plus a **stable within-patient pairwise ranking objective** capture the
> useful heterogeneity that v0.4's source-**name** tower captured — *without* any study identity — and does the
> resulting deployable ranker **R** beat genuine PRIME on source-balanced hits@20?

Two things change from v0.4, and we isolate each as cleanly as the data allow:
- **Objective** (Q vs P): replace v0.4's MIL log-sum-exp bag NLL with a **strictly convex** weighted
  within-patient **pairwise** ranking loss, so scores/coefficients are identical across initializations
  (eliminating the v0.4 wobble). Q carries **no context**. *(This contrast is not a pure objective isolation —
  see §2; the frozen additive rung A makes the objective-vs-supervision split visible.)*
- **Context** (R vs Q): add a small, **portable** context-interaction block on top of Q. R is the gated
  candidate. R vs F asks whether *deployable* context recovers what the *nondeployable* source-name tower got.

This is still a DEVELOPMENT experiment. No external-superiority claim is made regardless of outcome; ACCEPT
licenses at most one pre-registered look at the SEMI-CONSUMED Gartner TEST, and an external claim would still
need an untouched cohort.

## 2. Candidate ladder (frozen comparators + two new)

`x` = standardized feature vector (§6, unchanged from corrected v0.3/v0.4). `z` = centered context vector (§4).

| member | scorer | context | provenance | role |
|---|---|---|---|---|
| **P** — pooled | `w₀·x + b₀` | none | **frozen** corrected-v0.3 rung-3 MIL, **refit from frozen code** (§2.1) | objective/supervision baseline (bag-MIL) |
| **A** — additive | `w₀·x + b₀` | none | **frozen** v0.3 rung-2 `additive_logistic`, **refit from frozen code** | *descriptive* pointwise exact-witness comparator |
| **F** — feature tower | `(w₀+v_s)·x + b₀ + c_s` | source **name** | **frozen** v0.4 `F_feature_tower`, **refit from frozen code** (§2.1) | nondeployable mechanism ceiling |
| **Q** — shared pairwise | `w₀·x` *(no intercept)* | none | **new**, convex pairwise objective, exact-witness | isolates pairwise objective |
| **R** — context ranker *(gated)* | `w₀·x + Σ_{f∈C4}Σ_c β_{f,c}·x_f·z_c` *(no intercept)* | **portable** context | **new**, R ⊃ Q | the deployable candidate |

**No fitted intercept in Q/R.** A global intercept cancels in every within-patient pairwise difference
`s_a − s_b`, so it is **unidentifiable** under a pairwise objective (its data gradient is exactly zero). We
therefore fit **no intercept** in Q or R (equivalently `b ≡ 0`). This is what makes the objective strictly
convex: the data term is convex and **every remaining fitted coefficient** (`w₀`, and R's `β`) carries a
**positive L2 penalty** (§5.1). The frozen comparators P/A/F keep the intercepts they were fit with, but those
are **rank-inert within a patient** (a global `b₀` and per-source `c_s` shift a patient's whole pool equally), so
they do not affect hits@20 — they are shown only to describe the frozen scorers faithfully.

**Contrasts (honest, not over-claimed).**
- **Q−P is NOT a pure objective isolation.** For Gartner, Q is supervised by **exact positive peptide
  witnesses** (`eval_positive`) while P's MIL is supervised by **positive bag labels**. So Q−P confounds two
  changes: pairwise-vs-pointwise *objective* **and** exact-witness-vs-bag *supervision*.
- **A (frozen v0.3 rung-2 additive_logistic)** is a pointwise logistic on the **same exact-witness positives** as
  Q (instance-level negatives; **not** bag-disciplined — a caveat, so A is descriptive only). It lets the split
  be seen: **A vs Q** ≈ pointwise-vs-pairwise on exact witnesses (objective form); **A vs P** ≈
  exact-witness-pointwise vs bag-MIL (supervision granularity). No causal isolation is claimed from any single
  contrast; the three together describe the mechanism.
- **R−Q** isolates *context* (both new, same objective, same supervision, R ⊃ Q). **R−F** asks whether portable
  context matches the source-name ceiling.

### 2.1 How the frozen comparators are reconstructed (they are REFIT, not "loaded")

The frozen result JSONs **do not persist the state needed to score**: v0.4 `members.F_feature_tower.folds[*]`
records only the selected hyperparameters `(C, tau, lam)` and **no fitted coefficients**; v0.3
`ladder.rung3_MIL.folds[*]` / `ladder.rung2_additive.folds[*]` record coefficients but **no feature
standardizer** (`mean_`, `std_`) and **no intercept**. So P/A/F cannot be "loaded." Each is reconstructed by
**re-running its exact frozen code** (`MILRanker`, `AdditiveLogistic`, `TowerMILRanker`) on **each original
outer-train fold** using the **already-selected per-fold hyperparameters**, with **ZERO retuning**. Because those
fits are deterministic (fixed zero init, deterministic L-BFGS-B), a refit reproduces the frozen model. The runner
**verifies reproduction to numerical tolerance** before any comparison: P and A (convex) must reproduce the
frozen per-fold coefficients and per-patient/aggregate hits to tight tolerance; **F is nonconvex**, so it is
verified only to v0.4's own documented multi-init tolerance (its wobble), and any residual is reported. If a
comparator fails to reproduce, the run **fails fast** — we do not silently substitute a differently-fit model.

## 3. Pre-fit context **source-alias** audit (the gate on *what may enter R*)

Ran `scripts/epicurus_v05_context_audit.py` (read-only, OUTCOME-label-blind, **no Q/R model, no ranking
metric**; source labels used only descriptively to quantify aliasing; outcome labels never read) on the
non-quarantined eval pool (**315,050 rows / 152 patients** — gartner 56, improve 70, multimer 26 — over
{gartner, improve, multimer}). Full numbers:
`artifacts/milestone_7_decision/epicurus_v05/CONTEXT_FEASIBILITY_AUDIT.json`.

**Range overlap is not an alias test** — it is retained in the JSON as a descriptive statistic only. The **gate**
is an **actual source-classification audit**: predict source from each candidate context with **group-held-out
patients** (`StratifiedGroupKFold(5, seed=0)`, stratified on source, grouped on patient) and **source-balanced ×
patient-balanced** row weights (each source totals weight 1; patients equal within source; rows within a patient
share — so patients are the effective independent unit and Gartner's ~290k rows cannot inflate the score). We
report **macro balanced accuracy** (source-balanced chance = 1/3) and **macro one-vs-rest AUROC** (chance = 0.5).

**FROZEN decision rule (fixed before any Q/R fit).** A context is **APPROVED** only if
(a) it is **candidate-level** — a property of the single candidate (its peptide sequence, its HLA allele, its
leakage-safe already-computed predictor values) — and **not derived from source-specific pool enumeration or a
source-specific denominator**; and
(b) it is **not a source alias**: standalone **macro OVR AUROC ≤ 0.70**.
Any context needing pool enumeration / a shared-but-unavailable denominator, or that is a deterministic transform
of already-modeled features, is **rejected on structural grounds regardless of its alias number**.

### APPROVED (candidate-level, weak alias; the interaction-block contexts)

| context | formula (candidate `i` in patient `p`) | macro bal-acc | **macro AUROC** | within-patient | why deployable / not a proxy |
|---|---|--:|--:|---|---|
| `ctx_pep_len` | `len(peptide_i)` | 0.425 | **0.609** | varies (std 0.96, ~3.9 lengths/pt) | from the sequence alone; only weakly predicts source |
| `ctx_pred_disagree` | `SD([f_prime_pct, f_mix_pct, f_el_pct]_i)` | 0.366 | **0.564** | varies (std 0.074) | leakage-safe predictor **disagreement** from the already-masked percentiles; lets the ensemble learn *when predictors disagree*; near-chance alias |
| `ctx_hla_locus` | class-I locus **A/B/C** of `HLA_i` (A=ref → `loc_B, loc_C`) | 0.403 | **0.584** | varies (std 0.45) | biological locus of the candidate's own allele; portable across the three allele-string formats; weak alias |

All three are per-candidate (vary within patient → carry genuine within-patient ranking signal) and inference-
computable while scoring a brand-new pool. Contexts are **centered/scaled on outer-train only**, so R with all
β = 0 collapses **exactly** to Q (nested). Because each enters only through the interaction block `x_f · z_c`
(reweighting core feature `f`), **R − Q is exactly "context reweighting of the four presentation features."**

### REJECTED (with reasons; alias number corroborates the structural call)

| candidate | macro AUROC | reason for rejection |
|---|--:|---|
| `ctx_log_pool` (log candidate-pool size) | **0.830** | pool size has **no source-invariant denominator** — Gartner is combinatorially peptide/HLA-expanded, IMPROVE/multimer denominators are built differently; **strong source alias**. |
| `ctx_hla_multiplicity` (#distinct HLA in pool) | **0.777** | #distinct HLA observed in a pool is **enumeration-dependent**, not documented genotype breadth (no completeness flag across sources); **strong source alias**. |
| `ctx_hla_comp_frac` (HLA competition / pool) | **0.706** | denominator is the **source-specific pool size** → still pool-enumeration-dependent; alias also above threshold. |
| `ctx_hla_competition_raw` (unnormalized count) | **0.804** | pool-size-driven; near-disjoint across sources. |
| `ctx_routes_per_bag` (bag cardinality) | **0.831** | routes/mutation = a **categorical bag-vs-instance ascertainment proxy** (Gartner 1–74 vs instance sources ≈1); a perfect source/ascertainment proxy. |
| `ctx_pred_consensus` (mean of the 3 percentiles) | 0.500 | *not* an alias, but **redundant** — a deterministic transform of features already in `x`; carries no portable information beyond `x`. |
| `ctx_prime_masked` (PRIME leakage-mask flag) | 0.543 | *not* an alias, but **already modeled** (f_prime_pct→0.5 on masked rows) and its **rate** is source-skewed (multimer 0.197 vs ~0.01); reusing as context risks re-introducing leaked-PRIME structure. Kept strictly as the existing leakage-safe feature semantics, never as context. |

**Residual joint aliasing of the approved block** (`[pep_len, pred_disagree, locus]`): macro AUROC **0.682**
(< 0.70), balanced accuracy 0.485 — reported as a diagnostic, not a per-context gate. One honest caveat: the
multimer cohort has **no C-restricted candidates** in this pool, so `loc_C` is mildly source-informative; this is
a real cohort-composition fact (a per-candidate allele property, not a pool-enumeration artifact) and is
absorbed into the block number, which still clears the threshold. Diagnostic §5.4-C re-runs this alias audit at
fit time and forbids any transferability claim for any context that turns out near-perfectly source-aliasing.

## 4. The context vector `z`

For candidate `i` in patient `p`: `z_i = [ctx_pep_len, ctx_pred_disagree, ctx_locus_B, ctx_locus_C]`, each
**centered and scaled** using the outer-train mean/std (fit on train rows only, exactly like features — no
test/holdout leakage). `ctx_hla_locus` uses the **preregistered identifiable encoding**: locus ∈ {A, B, C}
parsed from the normalized allele string (`A0201`, `HLA-A01:01`, `HLA-A*02:01` all → the same locus), with **A as
the reference level** and `loc_B, loc_C` indicators; anything unparseable → reference. The interaction block
multiplies each of these four centered context columns by each of the **four core presentation features**
`C4 = [f_el_pct, f_pres_abs, f_prime_pct, f_mix_pct]` → up to **16** interaction terms `β_{f,c}`. **Orthogonal
recognition features (`f_expr, f_agreto, f_foreign, f_bindstab, f_proc`) remain shared** (no interactions), so no
recognition axis is source-specialized. (`ctx_pred_disagree` is a deterministic function of three of the core
features; interacting it with them adds label-blind curvature — "trust these features differently when they
disagree" — not a new recognition axis; diagnostic §5.4-J inspects which interaction terms carry weight.)

## 5. Objective, optimization, evaluation, gate

### 5.1 Pairwise ranking loss with bag-aware LME negatives (Q and R), strictly convex

Score `s(x,z) = θ·φ(x,z)`, where `φ = x` for Q and `φ = [x ; vec(x_{C4} ⊗ z)]` for R — **linear in θ, no
intercept**. For each within-patient (positive peptide `a`, **negative bag** `b`) pair, the negative bag is
represented by a convex **log-mean-exp over ALL its children's linear scores** at a **fixed** temperature `τ = 1`
(in score/logit units; not tuned):

```
s_j    = θ·φ_j                                             (linear in θ)
LME_b  = τ · log( (1/|b|) · Σ_{j∈b} exp(s_j / τ) )          (convex in θ; τ = 1 fixed)
L(θ)   = Σ_p Σ_{a∈pos_p} Σ_{b∈negbag_p} w_{p,a,b} · softplus( LME_b − s_a )
         + ½ λ_w ‖w₀‖² + ½ λ_ctx Σ_{f,c} β_{f,c}²
```

**Strict convexity (verified).** `s_j` is linear in θ ⇒ `exp(s_j/τ)` is convex ⇒ the mean is convex ⇒
`LME_b = LSE(s/τ)·τ − τ·log|b|` is convex (log-sum-exp is convex; subtracting the constant `τ·log|b|` preserves
it). `−s_a` is affine, so `(LME_b − s_a)` is convex. `softplus` is convex **and nondecreasing**, so
`softplus(convex)` is convex; a nonnegative-weighted sum of convex terms is convex. Adding `½λ_w‖w₀‖² +
½λ_ctx Σβ²` — a **positive-definite** quadratic over **every** fitted coefficient (`λ_w ≥ 0.1 > 0`, and either
`λ_ctx ≥ 3 > 0` or `λ_ctx = ∞` which constrains `β = 0`) — makes **L strictly convex with a unique global
minimum**. This is why there is no fitted intercept: with one, `b`'s data gradient is exactly zero (it cancels in
every difference) and it would need its own penalty to be identified; dropping it is cleaner and keeps the "every
fitted coefficient is L2-penalized" guarantee exact. The **interaction block gets stronger group shrinkage**
`λ_ctx ≥ λ_w`. Solved by deterministic L-BFGS-B from fixed zero init with analytic gradient. **Singleton bags**
(all instance-source negatives, and any Gartner bag with one child) satisfy `LME_b = s_b`, so the loss reduces
to **ordinary pairwise logistic** there.

**Negative bags (bag discipline; no model-selected mining).** Positives = `eval_positive` peptides.
Tested-negatives: instance sources = each response-0 instance (a singleton bag); **Gartner = each NEGATIVE bag**,
now represented by the **LME over all of its children** — *not* one fixed representative. This is the correction:
the LME rises whenever the model scores **any** child of a negative bag highly (via that child's `exp(s_j)`
weight), so children Q/R rank highly are **not** ignored, while the aggregation stays **convex** and **fixed**
(it is a smooth function of the current scores, not an adaptively selected hard negative → not model-selected
mining, and `s_b` never leaves the convex regime). Resource check (feasible): the eval pool has **6,685 Gartner
negative bags** with **≤ 74 children each** (mean 42), and `Σ_p P_p·B_p ≈ 1.7×10⁵` positive×negative-bag pairs;
per gradient evaluation costs one pass over the ~3×10⁵ instance scores plus a per-bag segmented softmax
(`numpy.add.reduceat` over bag-contiguous rows) plus ~1.7×10⁵ O(d) pair terms — trivial on CPU, fully
vectorized, no subsampling, no stochasticity.

**Weights (source-balanced, bag-disciplined).** For a pair (pos `a`, negative bag `b`) in patient `p` of
source `s`:

```
w_{p,a,b} = 1/S · 1/N_s · 1/P_p · 1/B_p
```

`S` = #sources (3), `N_s` = #patients in source `s`, `P_p` = #positive units in `p`, `B_p` = **#negative bags**
in `p`. Then each **source** contributes total pair-weight `1/S`; each **patient** equal within source; each
**positive unit** equal within patient; each **negative bag** equal within positive unit. Because a bag's
children are collapsed by the LME **before** weighting and each bag weighs `1/B_p`, **a 74-child bag cannot
dominate a 1-child bag and Gartner's ~290k child rows cannot dominate**. (Derivation: sum over `b` =
`1/(S N_s P_p)`; over `a` = `1/(S N_s)`; over `p∈s` = `1/S`; over `s` = 1.) All pairs are **fixed** (enumerated
a-priori), not model-selected.

### 5.2 Grid, nested selection, tie-breaks

- Grid: `λ_w ∈ {0.1, 0.3, 1, 3, 10}`, `λ_ctx ∈ {3, 10, 30, 100, ∞}` with the constraint `λ_ctx ≥ λ_w`.
  `λ_ctx = ∞` sets the whole interaction block to 0 → **R collapses to Q** (nested check). `τ = 1` is fixed, not
  on the grid.
- **Q** is selected over `λ_w` only (no interaction block). **R** is selected over `(λ_w, λ_ctx)`. Both selected
  **entirely inside outer-train** by inner source-balanced mean-patient hits@20 (same selection metric and
  source-balanced weighting as v0.4). Outer frozen folds unchanged; no patient scored by a model that saw it.
- Selection rule (deterministic): (1) maximize inner source-balanced hits@20; (2) within ε = 0.01 prefer **more
  shrinkage** (larger `λ_ctx`, then larger `λ_w`); (3) deterministic lexical tie-break on `(λ_ctx, λ_w)`.

### 5.3 Evaluation

Nested OOF on the 5 frozen folds; **source-balanced patient-paired bootstrap** (2000×, 95% percentile CI), the
identical harness as v0.3/v0.4 (`per_patient_metrics`, `paired_bootstrap`, `_source_balanced_weight`), reused
verbatim. Primary metric = **source-balanced mean patient hits@20**. Baselines: **genuine GfellerLab PRIME 2.1**
`%rank` (raw, unmasked — the harder bar), MixMHCpred 3.0 `%rank`, and the strongest per-source presentation
(`el_oriented`), orientation fixed a-priori and verified on outer-train only.

### 5.4 Gate and required diagnostics (none used for selection)

- **Primary gate (candidate = R vs genuine PRIME).** ACCEPT iff paired-bootstrap CI **lower bound > 0** on
  source-balanced mean hits@20. Otherwise an honest **REJECT/TIE** ("failed to establish superiority," not
  necessarily "worse"). No victory language unless this clears.
- **A. Non-regression** vs the strongest presentation (point ≥ 0, 95% CI not entirely < 0); **R vs Q** mechanism
  contrast; **Q vs P** (objective **+ exact-witness supervision** — not a pure isolation); **A vs Q** and **A vs
  P** (the objective-vs-supervision split, descriptive); **R vs F** (portable vs source-name ceiling).
- **B. Per-source** deltas + best-positive rank + nDCG@20 + recall@20, R vs each baseline; **PRIME
  availability/masking** report per source.
- **C. Context alias audit (the transferability guard).** Re-run the §3 source-classification alias audit at fit
  time (macro balanced-accuracy / macro OVR AUROC, patient-grouped, source-balanced), report each context's
  within-source variation, and **leave-one-context-out** ablation of R. Any context whose alias AUROC is
  near-perfect ⇒ **no transferability claim** permitted for it.
- **D. Shuffled-context negative control.** Refit R with contexts deterministically **shuffled within
  source/patient** (per-candidate contexts shuffled within patient). Comparable lift ⇒ capacity artifact, not
  context signal. Diagnostic.
- **E. Leave-one-source-out transfer stress test.** R vs Q trained without each source, scored on the held-out
  source. **Descriptive only (n = 3 sources).**
- **F. Convexity / determinism.** ≥2 extra perturbed fixed-seed inits; **required pass**: Spearman = 1.0 and
  max|Δcoefficient| ≤ 1e-6 (strict — the objective is provably convex, so this is a *hard gate*, unlike v0.4).
- **G. Label-blind attrition** reproduces **152 rankable** and **118 evaluated** positive-bearing patients
  (feature-bearing 152 → rankable-label-blind 152 → scored-has-positive 118), else explained.
- **H. Leakage assertions**: no patient/genomic-family/non-quarantined peptide crosses a fold; **0 Gartner TEST
  patients**; PRIME masking + recurrent-peptide quarantine **exactly as v0.4**.
- **I. Provenance** hashes re-verified (fail-fast on mismatch); Gartner-TEST I/O guard active; frozen comparators
  P/A/F **reproduction-verified** before any comparison (§2.1).
- **J. Interpretation.** State explicitly whether any gain is **recognition signal or presentation engineering**
  (inspect which interaction terms carry weight — context × presentation features is presentation reweighting,
  not a new recognition axis; flag if `pred_disagree` interactions dominate).

## 6. Features, scaling, leakage (unchanged from corrected v0.3 / v0.4)

CORE `f_prime_pct, f_mix_pct, f_el_pct, f_pres_abs` + ORTHO `f_expr, f_agreto, f_foreign, f_bindstab, f_proc`;
source-correct EL orientation; **fail-closed** drop of any unit missing a core score; masked-neutral orthogonals;
PRIME-leakage neutralization (f_prime_pct → 0.5 on near-PRIME-training peptides). **All fitted preprocessing
(feature std, context center/scale) fit inside outer-train only.** The frozen split and quarantine are used
verbatim.

## 7. Files (this commit is pre-fit; only design/prereg/provenance/audit)

- `docs/superpowers/specs/2026-07-11-v05-context-conditioned-design.md` — this design (source of truth).
- `artifacts/milestone_7_decision/epicurus_v05/PREREGISTERED_PROTOCOL.md` — the frozen preregistration.
- `artifacts/milestone_7_decision/epicurus_v05/CONTEXT_FEASIBILITY_AUDIT.json` — the pre-fit alias-audit numbers.
- `scripts/epicurus_v05_context_audit.py` — the read-only, no-Q/R-model alias audit that produced the JSON.
- `artifacts/milestone_7_decision/epicurus_v05/PROVENANCE.json` — pinned SHA256 of split/inputs/feature code/
  frozen comparators/this design/the prereg; re-verified by the future runner (fail-fast on mismatch).

Model implementation (`src/event_b/epicurus_v05.py`), tests, runner, and results are **not** part of this commit
and are **not** produced in this invocation.

## 8. Guardrails (unchanged from v0.4)

Gartner TEST not read/scored/tuned/evaluated (I/O guard patches `open`/`read_csv`/`read_excel`/`ZipFile` to raise
on any known TEST path). Frozen `mil_dev_split_v1` verbatim. v0.3/v0.4 frozen artifacts not modified. No
external-superiority claim; ACCEPT is development-only. PRIME/MixMHCpred binaries + training data are
non-commercial and gitignored. Only design + prereg + provenance + the audit (script/JSON) are committed before
fitting.
