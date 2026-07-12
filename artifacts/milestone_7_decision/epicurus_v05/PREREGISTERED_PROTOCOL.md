# Epicurus v0.5 — PREREGISTERED development protocol (deployable context-conditioned pairwise challenger)

**Registered before any Q/R model is fit or inspected.** This protocol fixes the hypothesis, units, model
members, baselines, evaluation, success gate, selection rule, diagnostics, and guardrails for the v0.5
DEVELOPMENT experiment on the frozen `mil_dev_split_v1`. It is falsifiable and cannot be re-specified after seeing
outcomes. The Gartner TEST holdout is **not loaded, not scored, not opened**. The frozen split is **not altered**.
No external-superiority claim will be made regardless of outcome.

Design source of truth: `docs/superpowers/specs/2026-07-11-v05-context-conditioned-design.md`. Split:
`configs/frozen/mil_dev_split_v1.json` (k=5). Pre-fit context feasibility audit (justifies the approved/rejected
contexts): `CONTEXT_FEASIBILITY_AUDIT.json` from `scripts/epicurus_v05_context_audit.py` (read-only, label-blind,
no model). Provenance pinned in `PROVENANCE.json` (SHA256 of split, inputs, feature code, frozen comparators, the
design, this protocol, the audit; git HEAD) and re-verified by the runner (fail-fast on mismatch).

---

## 1. Registered hypothesis

> Does **portable, scoring-time context** (variables knowable while ranking a brand-new patient's candidate pool,
> with **no study identity**) plus a **strictly convex within-patient pairwise ranking objective** capture the
> heterogeneity that v0.4's nondeployable source-**name** tower captured, and does the deployable ranker **R**
> beat genuine PRIME on source-balanced mean patient hits@20?

Two isolated changes vs v0.4: (i) **objective** — Q replaces the MIL log-sum-exp bag NLL with a strictly convex
pairwise loss (kills the v0.4 multi-init score wobble); (ii) **context** — R adds a small portable
context-interaction block on top of Q. **Q−P isolates the objective; R−Q isolates context; R−F compares portable
context to the source-name ceiling.**

## 2. Candidate ladder (§2 of design)

| member | scorer | context | provenance | role |
|---|---|---|---|---|
| **P** | `w₀·x + b₀` | none | **frozen** corrected-v0.3 rung-3 (loaded, not retuned) | objective baseline |
| **F** | `(w₀+v_s)·x + b₀ + c_s` | source **name** | **frozen** v0.4 F (loaded, not retuned) | nondeployable mechanism ceiling |
| **Q** | `w₀·x + b₀` | none | new, convex pairwise objective | isolates objective |
| **R** *(gated)* | `w₀·x + b₀ + Σ_{f∈C4}Σ_c β_{f,c}·x_f·z_c` | **portable** context | new, R ⊃ Q (β=0 ⇒ Q) | the deployable candidate |

`x` = standardized features (§6); `z` = centered context (§4). P, F loaded verbatim from their frozen
DEV_RESULT.json, paired on the identical eligible patients / identical per-patient pools as Q, R.

## 3. Decision units and labels (unchanged from v0.3/v0.4; never blended)

| source | decision unit | supervision | evaluation label |
|---|---|---|---|
| Gartner TRAIN × Müller | **mutation BAG** (25mer key) | bag CD8+ / screened-neg (UNTESTED/AMBIGUOUS excluded) | peptide-level exact positives (`VALIDATED=1`) in the patient's Müller pool |
| IMPROVE | instance (peptide-HLA) | response 0/1 | response 0/1 |
| CD8 multimer | instance (peptide-HLA) | response 0/1 | response 0/1 |

One shared instance scorer; metrics per source then aggregated; absolute hits never pooled across unit types.
Gartner negatives count at the **bag** level (one fixed label-blind representative per negative bag — §5), never
the ~285k child instances.

## 4. Approved contexts (§3 of design; pre-fit audit is the gate on what may enter R)

`z = [ctx_pep_len, ctx_log_pool, ctx_hla_multiplicity, ctx_hla_comp_frac]`, each **centered/scaled on
outer-train rows only**. All four are label-blind and computable while scoring a brand-new patient's pool from
peptide sequence + HLA typing + the pool itself:

| context | formula | min pairwise source overlap | within-patient |
|---|---|--:|---|
| `ctx_pep_len` | `len(peptide_i)` | 0.50 | varies |
| `ctx_log_pool` | `log1p(|pool_p|)` | 0.37 | patient-constant (interaction-only) |
| `ctx_hla_multiplicity` | `#distinct HLA in pool_p` | 0.25 | patient-constant (interaction-only) |
| `ctx_hla_comp_frac` | `(#cand sharing HLA_i)/|pool_p|` | 0.48 | varies |

Interaction block = these four × the **four core presentation features** `C4 = [f_el_pct, f_pres_abs,
f_prime_pct, f_mix_pct]` (≤16 `β_{f,c}`), with stronger group shrinkage than the shared weights.
**Orthogonal recognition features remain shared** (no interactions). Contexts centered ⇒ β = 0 makes R ≡ Q.

**REJECTED contexts (audited, will not enter any model):** `ctx_routes_per_bag` (bag cardinality; overlap 0.00,
categorical bag-vs-instance ascertainment proxy — prohibited); `ctx_hla_competition_raw` (overlap 0.05,
pool-size-driven source proxy — superseded by its fraction); ORTHO feature-availability mask (source-deterministic
= a source name); `ctx_prime_masked` as a context (source-skewed rate + already modeled — kept only as the
existing leakage-safe feature semantics). No source/study/corpus name, assay modality, label prevalence, fold ID,
outcome-derived statistic, or proxy for these enters any model.

## 5. Objective, negatives, weights, optimization (§5.1 of design)

Score linear in θ (`φ = x` for Q; `φ = [x ; vec(x_{C4}⊗z)]` for R). Weighted within-patient pairwise logistic loss:

```
L(θ) = Σ_p Σ_{a∈pos_p} Σ_{b∈neg_p} w_{p,a,b}·softplus(−(s_a − s_b)) + ½λ_w‖w₀‖² + ½λ_ctx Σ_{f,c} β_{f,c}²
w_{p,a,b} = 1/S · 1/N_s · 1/P_p · 1/B_p
```

**Strictly convex** (convex softplus ∘ linear + strictly-convex L2) ⇒ unique global minimum ⇒ init-invariant
coefficients. `λ_ctx ≥ λ_w` (stronger interaction shrinkage). Solved deterministically (L-BFGS-B, fixed zero
init, analytic gradient, `maxiter=500, ftol=1e-9`).

**Negatives (bag discipline; no model-selected mining).** Positives = `eval_positive`. Tested-negatives: instance
sources = each response-0 instance; **Gartner = each NEGATIVE bag** represented by **one fixed, label-blind,
model-blind** instance (highest a-priori `f_pres_abs`, deterministic peptide-string tie-break). One unit per bag
(never children); representative chosen before optimization from a fixed presentation prior ⇒ not model-selected
and `s_b` linear in θ. All pairs **fixed** (enumerated a-priori). The weights make each **source** contribute
total pair-weight `1/S`, each **patient** equal within source, each **positive unit** equal within patient, each
**negative unit** equal within positive unit — so Gartner's child rows cannot dominate.

## 6. Features, scaling, leakage (unchanged from corrected v0.3/v0.4)

CORE + ORTHO exactly as v0.3; source-correct EL orientation; fail-closed core-score drop; masked-neutral
orthogonals; PRIME-leakage neutralization (f_prime_pct → 0.5 on near-PRIME-training peptides — **PRIME masking
exactly as v0.4; masked Gartner rows are never anchored to leaked PRIME**). Feature std **and** context
center/scale fit **inside outer-train only**; within-patient quantities patient-local.

## 7. Grid, selection, gate (§5.2–5.4 of design)

- Grid: `λ_w ∈ {0.1,0.3,1,3,10}`, `λ_ctx ∈ {3,10,30,100,∞}` with `λ_ctx ≥ λ_w` (`∞` ⇒ R≡Q, nested check).
  Q selected over `λ_w` only; R over `(λ_w, λ_ctx)`. **Both selected entirely inside outer-train** by inner
  source-balanced mean-patient hits@20. **Outer frozen folds unchanged**; no patient scored by a model that saw
  it.
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
   **Q−P** (objective), **R−F** (portable vs source-name ceiling).
2. **Per-source** hits@20 / recall@20 / best-positive-rank / nDCG@20, R vs each baseline; **PRIME
   availability/masking rate per source**.
3. **Context alias audit** (transferability guard): descriptive prediction of source from the context block;
   within-source variation; **leave-one-context-out** ablation of R. Near-perfect aliasing ⇒ **no transferability
   claim**.
4. **Shuffled-context negative control**: refit R with contexts shuffled within source/patient (per-candidate
   within patient; patient-constant permuted across patients within source). Comparable lift ⇒ capacity artifact.
5. **Leave-one-source-out transfer stress test**: R vs Q trained without each source, scored on it —
   **descriptive (n=3)**.
6. **Convexity/determinism (hard gate)**: ≥2 perturbed fixed-seed inits ⇒ Spearman = 1.0 and max|Δcoef| ≤ 1e-6.
7. **Leakage assertions**: no patient/genomic-family/non-quarantined peptide crosses a fold; **0 Gartner TEST
   patients**; recurrent-peptide quarantine + PRIME masking exactly as v0.4.
8. **Selected (λ_w, λ_ctx) per fold**; effective coefficients incl. which interaction terms carry weight.
9. **Interpretation**: state whether any gain is **recognition signal or presentation engineering**
   (context × presentation-feature interactions are presentation reweighting, not a new recognition axis).

## 10. Guardrails

- Gartner TEST outcomes not read/scored/tuned/evaluated. **I/O guard** patches `open`/`read_csv`/`read_excel`/
  `ZipFile` to raise on any known TEST path (never opened, not merely filtered). Provenance re-verified
  (fail-fast on mismatch).
- Frozen `mil_dev_split_v1` used verbatim; not regenerated/altered. v0.3/v0.4 frozen artifacts not modified.
- No external-superiority claim; any ACCEPT is **development-only**. PRIME/MixMHCpred binaries + training data are
  non-commercial and gitignored. **Only the design, this protocol, `PROVENANCE.json`, and the read-only context
  audit (script + JSON) are committed before fitting**; the v0.5 model implementation, tests, runner, and results
  are **not** produced or committed in this invocation.
