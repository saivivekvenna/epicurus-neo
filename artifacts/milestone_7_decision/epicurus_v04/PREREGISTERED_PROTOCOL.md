# Epicurus v0.4 — PREREGISTERED development protocol (source-aware tower)

**Registered before any model is fit.** This protocol fixes the hypothesis, units, model members, baselines,
evaluation, success gate, selection rule, diagnostics, and guardrails for the v0.4 DEVELOPMENT experiment on the
frozen `mil_dev_split_v1`. It is written so the result is falsifiable and cannot be re-specified after seeing
outcomes. The Gartner TEST holdout is **not loaded, not scored, not opened**. The frozen split is **not
altered**. No external-superiority claim will be made from this experiment regardless of outcome.

Design source of truth: `docs/superpowers/specs/2026-07-11-v04-source-tower-design.md`. Split:
`configs/frozen/mil_dev_split_v1.json` (k=5). Provenance pinned in `PROVENANCE.json` (SHA256 of split, inputs,
feature code, this protocol; git HEAD) and re-verified by the runner (fail-fast on mismatch).

**Prior-run note (transparency).** The v0.3 first run was invalidated by a source-EL-semantics bug (Gartner
`Score_EL` is a 0-1 likelihood, higher-better, not a lower-better %rank). It was fixed and v0.3 rerun; **v0.4's
pooled comparator P is the frozen, corrected v0.3 result.** v0.3's corrected artifacts are **not modified** here.

---

## 1. Registered hypothesis (narrowed)

> Does **source-conditioned feature weighting** recover ranking signal erased by naive pooling, **beyond what
> source-prevalence calibration alone can explain**?

v0.4 is a **mechanism test**, not a deployable universal ranker: its heads are indexed by dataset **name**, so it
cannot natively score a brand-new cohort (see §9 fallback + future work). Motivation (v0.3 diagnostics): source-
only OOF beats pooled on Gartner (1.30 vs 1.025; PRIME 0.975 → +0.325 latent) and IMPROVE (1.167 vs 0.967) but
multimer *prefers* pooling (1.167 vs 1.056); study-shortcut source-identity AUROC 0.898, prevalence differs ~100×.

## 2. Decision units and labels (unchanged from v0.3; never silently blended)

| source | decision unit | supervision label | evaluation label |
|---|---|---|---|
| Gartner TRAIN × Müller | **mutation BAG** (transcript 25mer key) | bag CD8+ / screened-neg (UNTESTED/AMBIGUOUS excluded) | **peptide-level** exact positives (`VALIDATED=1`) within the patient's Müller pool |
| IMPROVE | **instance** (peptide-HLA) | response 0/1 | response 0/1 |
| CD8 multimer | **instance** (peptide-HLA) | response 0/1 | response 0/1 |

One shared instance scorer; metrics computed **per source, then aggregated**; absolute hit counts never pooled
across unit types. Gartner negatives count at the **bag** level, never the ~285k child instances.

## 3. Model members (three nested; one gated)

`x` = standardized feature vector (§6). `s(i)` = source of instance `i`.

| member | instance scorer | source params | role |
|---|---|---|---|
| **P — pooled** (= corrected v0.3) | `w₀·x + b₀` | none | frozen comparator, **loaded not retuned** |
| **C — calibration tower** | `w₀·x + b₀ + c_{s}` | intercept `c_s` | isolates prevalence calibration |
| **F — feature tower** *(gated)* | `(w₀+v_{s})·x + b₀ + c_{s}` | head `v_s` + intercept `c_s` | source-conditioned feature weighting |

**Per-source intercepts are rank-inert within patient** (a patient is entirely in one source; `c_s` shifts its
whole pool equally). So any ranking lift of **C over P** comes purely from freeing `w₀` of base-rate compromise,
and the lift of **F over C** is the registered "feature weighting beyond calibration" signal.

## 4. Objective and optimization

Minimize over `(w₀, b₀, {v_s}, {c_s})`:

```
L(θ) = L_MIL(θ) + α‖w₀‖² + α·λ·Σ_s‖v_s‖² + α·Σ_s c_s²      with  α = 1/(C·n_bags)
```

- `L_MIL` = source-balanced bag NLL: `Σ_bags w_b·bce(y_b, σ(lse_b))`, `lse_b = τ·log(mean_i exp(f(x_i)/τ))`,
  `w_b = 1/(S·N_s·n_bags_p)` → **each source contributes total weight 1/S; each patient equal within source**
  (patients, not candidate rows, are the weighting unit). Identical to v0.3's weighting.
- `C` (inverse reg) and `τ` control complexity/aggregation as in v0.3; `λ` controls **relative pooling of the
  feature heads**, defined on standardized features. Per-source intercepts share the `w₀` ridge (`α`), fixed.
- **Not asserted convex:** the positive-bag term `softplus(−lse)` is convex-non-increasing ∘ convex, not
  provably convex. Optimization is **deterministic** L-BFGS-B from fixed zero init (`maxiter=500, ftol=1e-9`)
  with analytic gradient. **Multi-init robustness** is a *diagnostic* (≥2 extra fixed-seed perturbed inits; pass
  = Spearman ≥0.999 and max|Δ std-score| ≤1e-3), **never** a tuning dimension.
- **Identifiability:** the ridge penalties select a unique split for every finite `λ>0`. We report the
  **effective weights `w_s = w₀+v_s`** and deviation norms `‖v_s‖`, never an ambiguous split. `λ=∞` is an
  **explicit pooled branch** (`v=c=0`), verified equal to P; independent per-source heads are a **separate
  diagnostic branch**, not claimed as the `λ→0` limit.

## 5. Baselines (orientations fixed a priori, verified not fit)

Genuine **GfellerLab PRIME 2.1** %rank (lower better; the incumbent), **MixMHCpred 3.0** %rank, and the
**strongest presentation** per source (Gartner `Score_EL` higher-better; IMPROVE `RankEL` / multimer
`EL (%Rank)` lower-better). The **genuine PRIME comparator uses raw unmasked `prime_rank`** (harder bar); the
Epicurus PRIME *feature* is masked-to-neutral on near-training peptides (§6). Orientation is checked on
outer-train only, never chosen from the outcome.

## 6. Features, scaling, leakage (unchanged from corrected v0.3)

CORE `f_prime_pct, f_mix_pct, f_el_pct, f_pres_abs` + ORTHO `f_expr, f_agreto, f_foreign, f_bindstab, f_proc`;
source-correct EL orientation; **fail-closed** drop of any unit missing a core score (never imputed);
masked-neutral orthogonals; PRIME-leakage neutralization (label-correlated flags never features). **All fitted
preprocessing is fit inside outer-train only** — within-patient percentiles are patient-local, `el_strength` is
row-local, and feature standardization mean/std come from the fold's train rows exclusively (λ defined on these).

## 7. Evaluation, selection, gate

- **Nested OOF** on the 5 frozen folds; source-balanced **patient-paired bootstrap** (2000×, 95% percentile CI).
  No patient scored by a model that saw it.
- **Grid:** `λ ∈ {0.03,0.1,0.3,1,3,10,30,∞}`, `C ∈ {0.1,0.3,1.0}`, `τ ∈ {0.5,1.0}`. **C and F are each
  independently nested-selected** inside outer-train.
- **Selection rule (frozen, deterministic), by inner source-balanced hits@20:** (1) maximize; (2) within
  tolerance ε=0.01 prefer **more pooled** (larger λ; ∞ beats finite); (3) then stronger reg (smaller C);
  (4) deterministic lexical tie-break on `(λ,C,τ)`.
- **Registered gate — candidate = F.** ACCEPT iff (i) Δ vs genuine PRIME CI_lower > 0 **and** (ii) no regression
  vs strongest presentation (point ≥0, 95% CI not entirely <0). Else **REJECT** = "failed to establish
  superiority" (not "worse"). C and P are context, never a second gated candidate. ACCEPT licenses **one**
  pre-registered look at the SEMI-CONSUMED Gartner TEST; an external claim still needs an untouched cohort.
- **Comparator fairness:** P = frozen corrected-v0.3 rung-3 (per-fold `(C,τ)`+coefficients loaded verbatim, not
  retuned). **F−P, F−C, C−P are paired on the identical eligible patients and identical per-patient pools.**
- **Mechanism contrast (registered, separate from the gate):** report **F−P, C−P, F−C** with CIs. Readings:
  (a) F ACCEPTs PRIME gate; (b) F ties PRIME but F−P and F−C CI_lo>0 = mechanism confirmed, still gate REJECT;
  (c) F−C spans 0 = mostly calibration/none; (d) no lift. No outcome re-spun.

## 8. Eligibility / attrition (label-blind; explains 152→118)

- **Rankable (label-blind):** ≥1 candidate with all core scores present after fail-closed filtering — uses **no
  labels** (asserted by test). **Scored ⊂ rankable:** additionally ≥1 in-pool positive (intrinsic to capture).
- Report the by-source attrition ladder feature-bearing → rankable → scored, plus core-feature/mask availability.
- Every baseline and model scored on **exactly the same** per-patient candidate pool (asserted).

## 9. Unseen-source fallback + future work

Score-time: known source → `w₀+v_s`; **unseen source → shared backbone `w₀`**. This is a stopgap; the intended
generalization is **assay/regime heads** keyed on observable properties (EL semantics, candidate-generation
protocol, assay type, peptide-length regime), out of scope for v0.4.

## 10. Required diagnostics (all reported; none used for selection)

1. **Per-source** hits@20 / recall@20 / best-positive-rank / nDCG@20, model vs each baseline.
2. **Source-only vs augmented** per source (honest OOF, training restricted to one source).
3. **Study shortcut** — source-identity AUROC + per-source prevalence.
4. **Coefficient/feature ablation** — leave-one-feature-out OOF Δ vs PRIME; per-source **effective weights
   `w_s`** + deviation norms `‖v_s‖` (how much / which direction each source departs from shared).
5. **Selected λ per fold** (how much specialization the data actually wanted).
6. **Score-orientation** Spearman vs label (confirms a-priori orientation, outer-train only).
7. **Family-leakage assertions** — no patient / genomic-family / non-quarantined peptide crosses a fold; **0
   Gartner TEST patients** (structural containment in TRAIN crosswalk).
8. **Recurrent-peptide quarantine stratum** — reported separately; **no selection on it**.
9. **Source-label negative control** — refit F with a deterministically shuffled per-patient pseudo-source;
   report its F−P lift. Comparable lift ⇒ capacity artifact, not source structure. **Diagnostic only.**
10. **Multi-init stability** numbers (§4). **PRIME mask/availability rate per source** (§5/§6).
11. **Mechanism contrasts** F−P / C−P / F−C (§7).

## 11. Guardrails

- Gartner TEST outcomes not read/scored/tuned/evaluated. **I/O guard** patches `open` to raise if any known TEST
  path is opened (never opened, not merely filtered). **Fail-fast** TEST-absence invariant (label-free
  containment). Provenance re-verified (fail-fast on mismatch).
- Frozen `mil_dev_split_v1` used verbatim; not regenerated/altered. v0.3 corrected artifacts not modified.
- No external-superiority claim; the source-**name** tower is **mechanism evidence only, not deployable/
  generalizable**; any ACCEPT is **development-only**. PRIME/MixMHCpred binaries + training data are
  non-commercial and gitignored. Only the design + this protocol + `PROVENANCE.json` are committed before
  fitting; implementation/tests/results are run but not committed.
