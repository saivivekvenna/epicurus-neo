# Epicurus v0.5 — PREREGISTERED development protocol (deployable context-conditioned pairwise challenger)

**Registered before any Q/R model is fit or inspected.** This protocol fixes the hypothesis, units, model
members, baselines, evaluation, success gate, selection rule, diagnostics, and guardrails for the v0.5
DEVELOPMENT experiment on the frozen `mil_dev_split_v1`. It is falsifiable and cannot be re-specified after seeing
outcomes. The Gartner TEST holdout is **not loaded, not scored, not opened**. The frozen split is **not altered**.
No external-superiority claim will be made regardless of outcome.

> **PRE-FIT CORRECTION of the original v0.5 preregistration (commit `06cb1d1`).** Six review-mandated fixes,
> all pre-fit: (1) **no fitted intercept** in Q/R (a global intercept cancels under a pairwise objective and is
> unidentifiable); (2) frozen comparators P/A/F are **refit from their exact frozen code**, not "loaded" (the
> JSONs do not persist scoreable state); (3) Q−P is stated to mix objective **and** exact-witness supervision,
> and the frozen additive rung **A** is added as a descriptive comparator; (4) each Gartner negative bag is a
> convex **log-mean-exp over all children**, not one fixed representative; (5) the context gate is an **actual
> patient-grouped, source-balanced source-classification alias audit**, which rejects all pool-enumeration
> contexts and keeps only candidate-level ones; (6) design, protocol, audit script+JSON, and provenance updated
> together.

Design source of truth: `docs/superpowers/specs/2026-07-11-v05-context-conditioned-design.md`. Split:
`configs/frozen/mil_dev_split_v1.json` (k=5). Pre-fit context **source-alias** audit (justifies the
approved/rejected contexts): `CONTEXT_FEASIBILITY_AUDIT.json` from `scripts/epicurus_v05_context_audit.py`
(read-only, OUTCOME-label-blind, no Q/R model). Provenance pinned in `PROVENANCE.json` (SHA256 of split, inputs,
feature code, frozen comparators, the design, this protocol, the audit; git HEAD) and re-verified by the runner
(fail-fast on mismatch).

---

## 1. Registered hypothesis

> Does **portable, scoring-time context** (variables knowable while ranking a brand-new patient's candidate pool,
> with **no study identity**) plus a **strictly convex within-patient pairwise ranking objective** capture the
> heterogeneity that v0.4's nondeployable source-**name** tower captured, and does the deployable ranker **R**
> beat genuine PRIME on source-balanced mean patient hits@20?

Two changes vs v0.4: (i) **objective** — Q replaces the MIL log-sum-exp bag NLL with a strictly convex pairwise
loss (kills the v0.4 multi-init score wobble); (ii) **context** — R adds a small portable context-interaction
block on top of Q. **R−Q isolates context.** **Q−P is NOT a pure objective isolation** — for Gartner, Q uses
exact-witness positives while P uses bag labels, so Q−P mixes pairwise objective **and** exact-witness
supervision; the frozen additive rung **A** (pointwise, exact-witness) makes the split visible (A−Q objective
form; A−P supervision granularity). **R−F** compares portable context to the source-name ceiling.

## 2. Candidate ladder (§2 of design)

| member | scorer | context | provenance | role |
|---|---|---|---|---|
| **P** | `w₀·x + b₀` | none | **frozen** corrected-v0.3 rung-3 MIL (**refit from code**, §2.1) | bag-MIL objective/supervision baseline |
| **A** | `w₀·x + b₀` | none | **frozen** v0.3 rung-2 `additive_logistic` (**refit from code**) | *descriptive* pointwise exact-witness comparator |
| **F** | `(w₀+v_s)·x + b₀ + c_s` | source **name** | **frozen** v0.4 `F_feature_tower` (**refit from code**, §2.1) | nondeployable mechanism ceiling |
| **Q** | `w₀·x` *(no intercept)* | none | new, convex pairwise objective, exact-witness | isolates the pairwise objective |
| **R** *(gated)* | `w₀·x + Σ_{f∈C4}Σ_c β_{f,c}·x_f·z_c` *(no intercept)* | **portable** context | new, R ⊃ Q (β=0 ⇒ Q) | the deployable candidate |

`x` = standardized features (§6); `z` = centered context (§4). **Q/R fit no intercept** — a global intercept
cancels in every within-patient pairwise difference (its data gradient is exactly zero), so it is
unidentifiable; dropping it keeps the "every fitted coefficient is L2-penalized" guarantee that makes the
objective strictly convex. P/A/F keep the intercepts they were fit with, but a global `b₀` / per-source `c_s` is
**rank-inert within a patient** and does not affect hits@20.

### 2.1 Frozen comparators are REFIT (not "loaded"), then reproduction-verified

The result JSONs do not persist scoreable state: v0.4 `members.F_feature_tower.folds[*]` stores only `(C, tau,
lam)` and **no coefficients**; v0.3 `ladder.rung3_MIL.folds[*]` / `ladder.rung2_additive.folds[*]` store
coefficients but **no feature standardizer** (`mean_`/`std_`) and **no intercept**. So P/A/F are reconstructed by
re-running their **exact frozen code** (`MILRanker`, `AdditiveLogistic`, `TowerMILRanker`) on **each original
outer-train fold** at the **already-selected per-fold hyperparameters**, with **ZERO retuning**. Fits are
deterministic (fixed zero init, deterministic L-BFGS-B), so a refit reproduces the frozen model. The runner
**verifies reproduction to numerical tolerance before any comparison**: P and A (convex) reproduce the frozen
coefficients + per-patient/aggregate hits to tight tolerance; **F is nonconvex**, so it is verified only to
v0.4's documented multi-init tolerance and any residual is reported. Reproduction failure ⇒ **fail fast** (no
silent substitution).

## 3. Decision units and labels (unchanged from v0.3/v0.4; never blended)

| source | decision unit | supervision | evaluation label |
|---|---|---|---|
| Gartner TRAIN × Müller | **mutation BAG** (25mer key) | bag CD8+ / screened-neg (UNTESTED/AMBIGUOUS excluded) | peptide-level exact positives (`VALIDATED=1`) in the patient's Müller pool |
| IMPROVE | instance (peptide-HLA) | response 0/1 | response 0/1 |
| CD8 multimer | instance (peptide-HLA) | response 0/1 | response 0/1 |

One shared instance scorer; metrics per source then aggregated; absolute hits never pooled across unit types.
Gartner negatives count at the **bag** level (a bag is one weighted unit; its children enter only through the
convex LME of §5), never as ~285k independent child instances.

## 4. Approved contexts (§3 of design; pre-fit **source-alias** audit is the gate on what may enter R)

`z = [ctx_pep_len, ctx_pred_disagree, ctx_locus_B, ctx_locus_C]`, each **centered/scaled on outer-train rows
only**. All are OUTCOME-label-blind and **candidate-level** (a property of the single candidate — its peptide
sequence, its HLA allele, its leakage-safe already-computed predictor percentiles), computable while scoring a
brand-new patient's pool. **The gate is an actual source-classification audit** (patient-grouped
`StratifiedGroupKFold(5, seed=0)`, source-balanced × patient-balanced weights, macro balanced-accuracy / macro
OVR AUROC); **range overlap is descriptive only, not the gate**.

**FROZEN rule (before any Q/R fit):** APPROVE iff candidate-level & not pool-enumeration/denominator-dependent
**AND** standalone macro OVR AUROC ≤ 0.70.

| context | formula | macro bal-acc | **macro AUROC** | verdict |
|---|---|--:|--:|---|
| `ctx_pep_len` | `len(peptide_i)` | 0.425 | **0.609** | APPROVE |
| `ctx_pred_disagree` | `SD([f_prime_pct,f_mix_pct,f_el_pct]_i)` (leakage-safe masked pcts) | 0.366 | **0.564** | APPROVE |
| `ctx_hla_locus` | class-I locus A/B/C of `HLA_i` (A=ref → `loc_B,loc_C`) | 0.403 | **0.584** | APPROVE |

Approved-block joint alias: macro AUROC **0.682** (< 0.70) — diagnostic. Interaction block = these four context
columns × the **four core presentation features** `C4 = [f_el_pct, f_pres_abs, f_prime_pct, f_mix_pct]` (≤16
`β_{f,c}`), with stronger group shrinkage than the shared weights. **Orthogonal recognition features remain
shared** (no interactions). Contexts centered ⇒ β = 0 makes R ≡ Q.

**REJECTED contexts (audited, will not enter any model):**
- `ctx_log_pool` (macro AUROC 0.830): pool size has **no source-invariant denominator** (Gartner is
  combinatorially peptide/HLA-expanded); strong source alias.
- `ctx_hla_multiplicity` (0.777): #distinct HLA in a pool is **enumeration-dependent**, not documented genotype
  breadth (no completeness flag across sources); strong source alias.
- `ctx_hla_comp_frac` (0.706): denominator is the **source-specific pool size** → pool-enumeration-dependent;
  alias also above threshold.
- `ctx_hla_competition_raw` (0.804) and `ctx_routes_per_bag` (0.831, a categorical bag-vs-instance ascertainment
  proxy): pool-size / ascertainment source proxies.
- `ctx_pred_consensus` (0.500): not an alias but **redundant** — a deterministic mean of already-modeled
  features.
- `ctx_prime_masked` (0.543): not an alias but **already modeled** (f_prime_pct→0.5) with a source-skewed rate;
  kept only as the existing leakage-safe feature semantics.

No source/study/corpus name, assay modality, label prevalence, fold ID, outcome-derived statistic, or proxy for
these enters any model.

## 5. Objective, negatives, weights, optimization (§5.1 of design)

Score **linear in θ, no intercept** (`φ = x` for Q; `φ = [x ; vec(x_{C4}⊗z)]` for R). Weighted within-patient
pairwise logistic loss with **bag-aware log-mean-exp negatives** (fixed `τ = 1`):

```
s_j    = θ·φ_j
LME_b  = τ · log( (1/|b|) Σ_{j∈b} exp(s_j/τ) )              # convex; over ALL children of negative bag b
L(θ)   = Σ_p Σ_{a∈pos_p} Σ_{b∈negbag_p} w_{p,a,b}·softplus(LME_b − s_a) + ½λ_w‖w₀‖² + ½λ_ctx Σ_{f,c} β_{f,c}²
w_{p,a,b} = 1/S · 1/N_s · 1/P_p · 1/B_p          # B_p = #negative BAGS in patient p
```

**Strictly convex** (verified): `s_j` linear ⇒ `LME_b` convex (log-sum-exp minus a constant); `−s_a` affine ⇒
`(LME_b−s_a)` convex; `softplus` convex **and nondecreasing** ⇒ `softplus(·)` convex; nonnegative-weighted sum
convex; plus a **positive-definite L2 over every fitted coefficient** (`λ_w ≥ 0.1`; `λ_ctx ≥ 3` or `∞`⇒β=0) ⇒
**unique global minimum** ⇒ init-invariant coefficients. Solved deterministically (L-BFGS-B, fixed zero init,
analytic gradient, `maxiter=500, ftol=1e-9`). **Singleton bags ⇒ `LME_b = s_b` ⇒ ordinary pairwise logistic.**

**Negatives (bag discipline; no model-selected mining).** Positives = `eval_positive`. Tested-negatives: instance
sources = each response-0 instance (a singleton bag); **Gartner = each NEGATIVE bag**, represented by the convex
**LME over all its children** (not one fixed representative). The LME rises whenever the model scores **any**
child of a negative bag highly, so children R ranks highly are not ignored, yet the aggregation is **fixed and
convex** (a smooth function of current scores — not adaptive hard-negative mining). Feasibility: 6,685 Gartner
negative bags, ≤74 children each (mean 42), `Σ_p P_p·B_p ≈ 1.7×10⁵` pairs — a segmented softmax over
bag-contiguous rows plus O(d) per pair, fully vectorized, no subsampling. The weights make each **source**
contribute total pair-weight `1/S`, each **patient** equal within source, each **positive unit** equal within
patient, each **negative bag** equal within positive unit — so a 74-child bag cannot dominate a 1-child bag and
Gartner's child rows cannot dominate. All pairs **fixed** (enumerated a-priori).

## 6. Features, scaling, leakage (unchanged from corrected v0.3/v0.4)

CORE + ORTHO exactly as v0.3; source-correct EL orientation; fail-closed core-score drop; masked-neutral
orthogonals; PRIME-leakage neutralization (f_prime_pct → 0.5 on near-PRIME-training peptides — **PRIME masking
exactly as v0.4; masked Gartner rows are never anchored to leaked PRIME**). Feature std **and** context
center/scale fit **inside outer-train only**; within-patient/patient-local quantities computed from the
candidate/pool without labels.

## 7. Grid, selection, gate (§5.2–5.4 of design)

- Grid: `λ_w ∈ {0.1,0.3,1,3,10}`, `λ_ctx ∈ {3,10,30,100,∞}` with `λ_ctx ≥ λ_w` (`∞` ⇒ R≡Q, nested check).
  `τ = 1` fixed (not on the grid). Q selected over `λ_w` only; R over `(λ_w, λ_ctx)`. **Both selected entirely
  inside outer-train** by inner source-balanced mean-patient hits@20. **Outer frozen folds unchanged**; no
  patient scored by a model that saw it.
- Selection rule (frozen, deterministic): (1) maximize inner source-balanced hits@20; (2) within ε=0.01 prefer
  more shrinkage (larger λ_ctx, then larger λ_w); (3) lexical tie-break on `(λ_ctx, λ_w)`.
- **Primary gate — candidate = R vs genuine PRIME.** ACCEPT iff paired-bootstrap **CI_lower > 0** on
  source-balanced mean patient hits@20. Else **REJECT/TIE** (honest; "failed to establish superiority").
  No victory language unless the gate clears. ACCEPT licenses **one** pre-registered look at the SEMI-CONSUMED
  Gartner TEST; an external claim still needs an untouched cohort.
- Evaluation: nested OOF on the 5 frozen folds; source-balanced patient-paired bootstrap (2000×, 95% percentile
  CI) — the identical harness as v0.3/v0.4, reused verbatim. Genuine PRIME uses **raw unmasked `prime_rank`**
  (harder bar). Orientations fixed a-priori, verified on outer-train only.

## 8. Eligibility / attrition (label-blind; must reproduce 152 → 118)

Rankable (label-blind): ≥1 candidate with all core scores present after fail-closed filtering — **no labels**.
Scored ⊂ rankable: additionally ≥1 in-pool positive. **Required:** feature-bearing 152 → rankable-label-blind 152
→ scored-has-positive 118, else explained. Every baseline and model scored on **exactly the same** per-patient
pool.

## 9. Required diagnostics (all reported; none used for selection)

1. **Contrasts** with CIs: R vs PRIME (gate), R vs strongest presentation (non-regression), **R−Q** (context),
   **Q−P** (objective **+ exact-witness supervision**, not a pure isolation), **A−Q** and **A−P** (the
   objective-vs-supervision split, descriptive), **R−F** (portable vs source-name ceiling).
2. **Per-source** hits@20 / recall@20 / best-positive-rank / nDCG@20, R vs each baseline; **PRIME
   availability/masking rate per source**.
3. **Context alias audit** (transferability guard): re-run the §4 patient-grouped source-balanced
   source-classification audit at fit time (macro bal-acc / macro OVR AUROC); within-source variation;
   **leave-one-context-out** ablation of R. Near-perfect aliasing ⇒ **no transferability claim** for that
   context.
4. **Shuffled-context negative control**: refit R with contexts shuffled within source/patient (per-candidate
   within patient). Comparable lift ⇒ capacity artifact.
5. **Leave-one-source-out transfer stress test**: R vs Q trained without each source, scored on it —
   **descriptive (n=3)**.
6. **Convexity/determinism (hard gate)**: ≥2 perturbed fixed-seed inits ⇒ Spearman = 1.0 and max|Δcoef| ≤ 1e-6.
7. **Leakage assertions**: no patient/genomic-family/non-quarantined peptide crosses a fold; **0 Gartner TEST
   patients**; recurrent-peptide quarantine + PRIME masking exactly as v0.4; frozen comparators
   **reproduction-verified** (§2.1).
8. **Selected (λ_w, λ_ctx) per fold**; effective coefficients incl. which interaction terms carry weight.
9. **Interpretation**: state whether any gain is **recognition signal or presentation engineering** (context ×
   presentation-feature interactions are presentation reweighting, not a new recognition axis; flag if
   `pred_disagree` interactions dominate).

## 10. Guardrails

- Gartner TEST outcomes not read/scored/tuned/evaluated. **I/O guard** patches `open`/`read_csv`/`read_excel`/
  `ZipFile` to raise on any known TEST path (never opened, not merely filtered). Provenance re-verified
  (fail-fast on mismatch).
- Frozen `mil_dev_split_v1` used verbatim; not regenerated/altered. v0.3/v0.4 frozen artifacts not modified.
- No external-superiority claim; any ACCEPT is **development-only**. PRIME/MixMHCpred binaries + training data are
  non-commercial and gitignored. **Only the design, this protocol, `PROVENANCE.json`, and the read-only context
  audit (script + JSON) are committed before fitting**; the v0.5 model implementation, tests, runner, and results
  are **not** produced or committed in this invocation.
