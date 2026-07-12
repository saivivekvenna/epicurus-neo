# Epicurus v0.5 — deployable context-conditioned pairwise challenger: design spec

**Design source of truth for the v0.5 DEVELOPMENT experiment.** This is the pre-fit design that the
preregistration (`artifacts/milestone_7_decision/epicurus_v05/PREREGISTERED_PROTOCOL.md`) formalizes and
freezes. It fixes the hypothesis, the two new model members (Q, R), the two frozen comparators (P, F), the
approved/rejected scoring-time contexts (with the pre-fit feasibility audit that justifies each), the pairwise
objective and its exact weights, the grid, the gate, and every required diagnostic. Written before any Q/R model
is fit or inspected. **Gartner TEST is not loaded, not scored, not opened.** The frozen `mil_dev_split_v1` is
used verbatim.

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

Two things change from v0.4, and we isolate each:
- **Objective** (Q vs P): replace v0.4's MIL log-sum-exp bag NLL with a **strictly convex** weighted
  within-patient **pairwise** ranking loss, so scores/coefficients are identical across initializations
  (eliminating the v0.4 wobble). Q carries **no context**.
- **Context** (R vs Q): add a small, **portable** context-interaction block on top of Q. R is the gated
  candidate. R vs F asks whether *deployable* context recovers what the *nondeployable* source-name tower got.

This is still a DEVELOPMENT experiment. No external-superiority claim is made regardless of outcome; ACCEPT
licenses at most one pre-registered look at the SEMI-CONSUMED Gartner TEST, and an external claim would still
need an untouched cohort.

## 2. Candidate ladder (two frozen, two new)

`x` = standardized feature vector (§6, unchanged from corrected v0.3/v0.4). `z` = centered context vector (§4).

| member | scorer | context | provenance | role |
|---|---|---|---|---|
| **P** — pooled | `w₀·x + b₀` | none | **frozen** corrected-v0.3 rung-3, loaded not retuned | objective baseline |
| **F** — feature tower | `(w₀+v_s)·x + b₀ + c_s` | source **name** | **frozen** v0.4 F, loaded not retuned | nondeployable mechanism ceiling |
| **Q** — shared pairwise | `w₀·x + b₀` | none | **new**, convex pairwise objective | isolates the objective |
| **R** — context ranker *(gated)* | `w₀·x + b₀ + Σ_{f∈C4}Σ_c β_{f,c}·x_f·z_c` | **portable** context | **new**, R ⊃ Q | the deployable candidate |

**Contrasts.** Q−P isolates the *objective* (pairwise vs MIL-bag-NLL); R−Q isolates *context*; R−F asks whether
portable context matches the source-name ceiling. P and F are loaded verbatim from their frozen DEV_RESULT.json
(no retuning) and paired on the identical eligible patients / identical per-patient pools.

## 3. Pre-fit context feasibility audit (the gate on *what may enter R*)

Ran `scripts/epicurus_v05_context_audit.py` (read-only, label-blind, **no model, no ranking metric**; source
labels used only descriptively to quantify aliasing; outcome labels never read) on the non-quarantined eval pool
(315,050 rows / 152 patients / {gartner, improve, multimer}). Full numbers:
`artifacts/milestone_7_decision/epicurus_v05/CONTEXT_FEASIBILITY_AUDIT.json`. A context is **APPROVED** only if it
is (a) computable while scoring a brand-new patient's pool from peptide sequence + HLA typing + the pool itself,
with **no labels**, and (b) **not a near-perfect source proxy** (pairwise range overlap materially > 0 and
within-source variation present). "Pairwise range overlap" = shared span / union span of two sources' value ranges.

### APPROVED (≤4, exactly the interaction-block contexts)

| context | formula (per candidate `i` in patient `p`) | min pairwise overlap | within-patient | why deployable / not a proxy |
|---|---|--:|---|---|
| `ctx_pep_len` | `len(peptide_i)` | **0.50** | varies (std 0.96, ~3.9 lengths/pt) | from the sequence alone; lengths overlap heavily across sources (8–12 everywhere) |
| `ctx_log_pool` | `log1p(|pool_p|)` | **0.37** | patient-constant → interaction-only | you know your candidate-pool size when scoring; ranges overlap (Gartner 4.5–10.3, IMPROVE 4.1–7.0) |
| `ctx_hla_multiplicity` | `#distinct HLA alleles in pool_p` | **0.25** | patient-constant → interaction-only | patient HLA-typing breadth = real product biology; IMPROVE spans the full 1–6, Gartner clusters at full 6 |
| `ctx_hla_comp_frac` | `(#cand in pool_p sharing HLA_i) / |pool_p|` | **0.48** | varies (std 0.055) | scale-free share of the pool competing for this allele; computable from pool + HLA typing |

Contexts are centered (mean subtracted, fit on outer-train only) so R with all β = 0 collapses **exactly** to Q
(nested). Patient-constant contexts (`ctx_log_pool`, `ctx_hla_multiplicity`) have **no within-patient main-effect
ranking value** — a constant shifts a patient's whole pool equally — so they enter **only** through the
interaction block `x_f · z_c`, which reweights feature `f` per patient. Per-candidate contexts (`ctx_pep_len`,
`ctx_hla_comp_frac`) also enter only through the interaction block, so **R − Q is exactly "context reweighting of
the four presentation features,"** a clean mechanism contrast.

### REJECTED (with reasons)

| candidate | min overlap | reason for rejection |
|---|--:|---|
| `ctx_routes_per_bag` (bag cardinality = routes/mutation) | **0.00** | Gartner 1–74 vs instance sources ≡ 1 (disjoint). This is exactly the **categorical bag-vs-instance benchmark ascertainment** the prohibitions forbid calling deployable biology — it is a perfect source/ascertainment proxy in this corpus. |
| `ctx_hla_competition_raw` (unnormalized count) | **0.05** | pool-size-driven; Gartner (1–7006) near-disjoint from IMPROVE/multimer. A source proxy. Superseded by its scale-free `ctx_hla_comp_frac` (overlap 0.48). |
| ORTHO feature-availability mask | — (source-deterministic) | which orthogonal features exist is fixed *per source* (Gartner=agreto+bindstab; multimer=expr+agreto+foreign+proc; IMPROVE=expr) → a **perfect source name**. Forbidden. |
| `ctx_prime_masked` (PRIME leakage-mask rate) | binary (range-overlap uninformative) | mask **rate** is source-skewed (multimer 0.197 vs ~0.01 elsewhere) AND the mask is **already modeled** (f_prime_pct→0.5 neutral on masked rows). Re-using it as a context risks re-introducing leaked-PRIME structure. Kept strictly as the existing leakage-safe feature semantics, never as context. |

**Alias caveat carried into the gate.** None of the four approved contexts is a *near-perfect* proxy, but
`ctx_log_pool` and `ctx_hla_multiplicity` are the most source-correlated (patient-constant, moderate overlap).
Gate diagnostic §5.4-D quantifies how well the approved context block predicts source and ablates each context;
**near-perfect aliasing (there is none here, but if a re-run shows it) forbids any transferability claim.**

## 4. The context vector `z`

For candidate `i` in patient `p`: `z_i = [ctx_pep_len, ctx_log_pool, ctx_hla_multiplicity, ctx_hla_comp_frac]`,
each **centered** using the outer-train mean and scaled by the outer-train std (fit on train rows only, exactly
like features — no test/holdout leakage; patient-local quantities computed within each patient's own pool). The
interaction block multiplies each centered context by each of the **four core presentation features**
`C4 = [f_el_pct, f_pres_abs, f_prime_pct, f_mix_pct]` → up to 16 interaction terms `β_{f,c}`. **Orthogonal
recognition features (`f_expr, f_agreto, f_foreign, f_bindstab, f_proc`) remain shared** (no interactions), so no
recognition axis is source-specialized.

## 5. Objective, optimization, evaluation, gate

### 5.1 Pairwise ranking loss (Q and R), strictly convex

Score `s(x,z) = θ·φ(x,z)`, where `φ = x` for Q and `φ = [x ; vec(x_{C4} ⊗ z)]` for R — **linear in θ**. For each
within-patient ordered pair (positive unit `a`, tested-negative unit `b`):

```
L(θ) = Σ_{p} Σ_{a∈pos_p} Σ_{b∈neg_p} w_{p,a,b} · softplus( −(s_a − s_b) )
       + ½ λ_w ‖w₀‖² + ½ λ_ctx Σ_{f,c} β_{f,c}²
```

`softplus(−Δ) = log(1+exp(−Δ))` is convex in Δ, Δ is **linear** in θ, and the L2 terms are strictly convex →
**L is strictly convex with a unique global minimum**. Orthogonal-feature and intercept weights share the `λ_w`
ridge; the **interaction block gets its own, stronger shrinkage `λ_ctx ≥ λ_w`** ("stronger group shrinkage than
shared coefficients"). Solved by deterministic L-BFGS-B from fixed zero init with analytic gradient.

**Pairs and negative units (bag discipline; no model-selected mining).** Positives = `eval_positive` peptides.
Tested-negatives: instance sources = each response-0 instance; **Gartner = each NEGATIVE bag** represented by
**one fixed label-blind, model-blind instance** — its highest a-priori presentation strength `f_pres_abs`
(deterministic peptide-string tie-break), so a bag counts **once** (never its ~285k children) and the ranker is
forced to beat the *hardest-presenting decoy per mutation*. The representative is chosen **before** optimization
from a fixed prior, so it is **not** model-selected hard-negative mining, and `s_b` stays linear in θ (convex).

**Weights (source-balanced, bag-disciplined).** For a pair (pos `a`, neg unit `b`) in patient `p` of source `s`:

```
w_{p,a,b} = 1/S · 1/N_s · 1/P_p · 1/B_p
```

`S` = #sources (3), `N_s` = #patients in source `s`, `P_p` = #positive units in `p`, `B_p` = #tested-negative
units in `p`. Then: each **source** contributes total pair-weight `1/S`; each **patient** equal within source;
each **positive unit** equal within patient; each **negative unit** equal within positive unit. Gartner's child
rows collapse to bag units before weighting, so **292k children cannot dominate**. (Derivation: sum over `b` =
`1/(S N_s P_p)`; over `a` = `1/(S N_s)`; over `p∈s` = `1/S`; over `s` = 1.) All pairs are **fixed** (enumerated
a-priori), not model-selected. Resource bound: pair count ≈ `Σ_p P_p·B_p` (~10⁵–10⁶ with Gartner negatives at
bag granularity) — vectorized; if it exceeds the budget, negatives are represented by the same one-per-bag units
already specified (no subsampling of positives, no stochasticity).

### 5.2 Grid, nested selection, tie-breaks

- Grid: `λ_w ∈ {0.1, 0.3, 1, 3, 10}`, `λ_ctx ∈ {3, 10, 30, 100, ∞}` with the constraint `λ_ctx ≥ λ_w`.
  `λ_ctx = ∞` sets the whole interaction block to 0 → **R collapses to Q** (nested check).
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
  contrast; **Q vs P** objective contrast; **R vs F** (portable vs source-name ceiling).
- **B. Per-source** deltas + best-positive rank + nDCG@20 + recall@20, R vs each baseline; **PRIME
  availability/masking** report per source.
- **C. Context alias audit (the transferability guard).** Report how well the approved context block predicts
  source (descriptive AUROC/accuracy), each context's within-source variation, and **leave-one-context-out**
  ablation of R. Near-perfect aliasing ⇒ **no transferability claim** permitted for that context.
- **D. Shuffled-context negative control.** Refit R with contexts deterministically **shuffled within
  source/patient** where meaningful (per-candidate contexts shuffled within patient; patient-constant contexts
  permuted across patients within source). Comparable lift ⇒ capacity artifact, not context signal. Diagnostic.
- **E. Leave-one-source-out transfer stress test.** R vs Q trained without each source, scored on the held-out
  source. **Descriptive only (n = 3 sources).**
- **F. Convexity / determinism.** ≥2 extra perturbed fixed-seed inits; **required pass**: Spearman = 1.0 and
  max|Δcoefficient| ≤ 1e-6 (strict — the objective is provably convex, so this is a *hard gate*, unlike v0.4).
- **G. Label-blind attrition** reproduces **152 rankable** and **118 evaluated** positive-bearing patients
  (feature-bearing 152 → rankable-label-blind 152 → scored-has-positive 118), else explained.
- **H. Leakage assertions**: no patient/genomic-family/non-quarantined peptide crosses a fold; **0 Gartner TEST
  patients**; PRIME masking + recurrent-peptide quarantine **exactly as v0.4**.
- **I. Provenance** hashes re-verified (fail-fast on mismatch); Gartner-TEST I/O guard active.
- **J. Interpretation.** State explicitly whether any gain is **recognition signal or presentation
  engineering** (inspect which interaction terms carry weight — context × presentation features is presentation
  reweighting, not a new recognition axis).

## 6. Features, scaling, leakage (unchanged from corrected v0.3 / v0.4)

CORE `f_prime_pct, f_mix_pct, f_el_pct, f_pres_abs` + ORTHO `f_expr, f_agreto, f_foreign, f_bindstab, f_proc`;
source-correct EL orientation; **fail-closed** drop of any unit missing a core score; masked-neutral orthogonals;
PRIME-leakage neutralization (f_prime_pct → 0.5 on near-PRIME-training peptides). **All fitted preprocessing
(feature std, context center/scale) fit inside outer-train only.** The frozen split and quarantine are used
verbatim.

## 7. Files (this commit is pre-fit; only design/prereg/provenance)

- `docs/superpowers/specs/2026-07-11-v05-context-conditioned-design.md` — this design (source of truth).
- `artifacts/milestone_7_decision/epicurus_v05/PREREGISTERED_PROTOCOL.md` — the frozen preregistration.
- `artifacts/milestone_7_decision/epicurus_v05/CONTEXT_FEASIBILITY_AUDIT.json` — the pre-fit audit numbers.
- `scripts/epicurus_v05_context_audit.py` — the read-only, no-model audit that produced the JSON.
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
