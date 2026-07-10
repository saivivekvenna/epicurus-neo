# M6: First Recognition Swing — Design

**Status:** approved design, pre-registration. No model trained yet.
**Date:** 2026-07-10.
**Standing verdict:** `INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA` (45 candidate-resolved
patients < 100). M6 is a **diagnostic swing under that verdict**, not a headline or clinical claim.
The gate exists to prevent strong claims, not to prevent learning.

## The question

> Does a learned recognition model improve patient-level top-*k* vaccine-target selection on
> *entirely unseen Event-B studies*, beyond prevalence and — where a valid presentation baseline
> exists — beyond presentation alone?

The first swing must be **boring, hard to game, and capable of hurting our feelings.** A wide
interval, `CONSISTENT_WITH_NO_EFFECT`, or `REJECT` are all valid, informative outcomes. Success is a
leakage-safe answer plus confound diagnostics plus a precise statement of the missing variable — not
necessarily a green `ACCEPT`.

## Sequencing (locked)

Three deliverables, built and frozen **in order**. Later stages must not alter earlier stages'
feature registration, model selection, or interpretation.

| Stage | Scope | Enters primary gate? |
|---|---|---|
| **M6A** | Event-B-only models + LOSO. The clean scientific test. | Yes |
| **M6B** | Event-A → Event-B transfer (IMPROVE). Exploratory extension. | No (own secondary gate) |
| **M6C** | Osteosarc full-patient replay. **Contract document only, no run.** | No |

> Note on renaming: an earlier draft used "M6B" for both the transfer arm and the Osteosarc
> contract. Resolved here: **M6A** = Event-B-only, **M6B** = Event-A transfer, **M6C** = Osteosarc
> contract.

M6A is committed and frozen before M6B is touched. If transfer helps, good; if it hurts, that is
biologically informative (IMPROVE measures Event-A, not Event-B). It must not blur the first clean
answer.

---

## The corpus, as it actually is

Every design choice below is forced by these grounded facts (pinned corpus
`outputs/event_b_backbone/combined/`, Event-B only — **not** the legacy IMPROVE-mixed combined).

- **45 candidate-resolved patients** across **4 studies**; Nous-209's 37 patient-level-only patients
  are **excluded from the peptide classifier entirely** (no peptide labels).
- Per-study candidate-resolved patients: **pdac 16, mKRAS 12, braun 9, hu 8.**
- **974 primary candidate labels** — 272 POSITIVE, 693 TESTED_NEGATIVE, **9 UNTESTED dropped.**
- Class mix (primary): 207 class-I, 334 class-II, **433 unknown/both-class.**
- **Per-patient candidate count is small and study-shaped:** median 14; only **11/45 patients have
  ≥20 candidates**; 15/45 have <10. Per study: hu 37–128 (rich), pdac 7–20, braun 8–19, **mKRAS
  exactly 6 for every patient** (shared KRAS panel).

**Consequence that reshapes the framing:** for any patient where `n_candidates ≤ k`, "select the top
*k*" is **degenerate** — you select all of them, ranking is irrelevant, and the metric collapses to
prevalence × n. A non-degenerate *selection* question exists for only ~11 patients, heavily hu. So
this corpus **barely supports a "top-20 from the full universe" claim**; it primarily supports
(a) reranking among *reported/tested* candidates and (b) per-peptide recognition *classification*.
M6 will name each result for what it actually is.

Additional standing blockers (unchanged): no presentation score exists in the corpus yet; braun's
HLA is not public; study_id is near-perfectly confounded with cancer type, platform, and label
prevalence.

---

## Data substrate & leakage discipline

- **Unit:** one label per candidate (dedup assays to the primary candidate label). POSITIVE=1 vs
  TESTED_NEGATIVE=0. UNTESTED dropped. Nous excluded.
- **Features are pre-vaccine only.** Allowed: peptide sequence, biochemistry, mutant-vs-WT
  contrasts, HLA/class, presentation score (WP1), and — conditionally — out-of-fold external
  cross-reactivity and PLM embeddings. **Never features:** `response_label`, `timepoint`,
  `relative_to_vaccine`, `event_type`, `qualitative_result`, `quantitative_result`, any assay-result
  or recognition-evidence table, and — critically — `study_id`, `cancer_type`, `vaccine_platform`
  (see WP2/WP5). Enforce `information_timing = PRE_SELECTION`.
- **Grouping key for all splits and bootstraps is `patient_id`.** No peptide-row-level resampling.

---

## Work packages

### WP0 — Candidate-universe completeness gate *(new; runs before any metric)*

Arguably more important than another model. Before top-*k* metrics, characterize what each patient's
candidate table actually is. Per patient emit: `n_candidates`, `n_tested`, `n_positive`,
`denominator_type ∈ {COMPLETE_TESTED_SET, PARTIAL_CANDIDATE_SET, POSITIVE_ENRICHED,
UNKNOWN_DENOMINATOR}`, `k_patient = min(20, n_eligible)`, and
`ranking_informative = (n_eligible > k_patient)`.

`n_eligible` is defined once and used everywhere: **the count of that patient's candidates that are
label-resolved (POSITIVE or TESTED_NEGATIVE, not UNTESTED) and scoreable under the model being
evaluated (all required features present).** Ranking is over exactly those candidates.

- Patients whose table is `POSITIVE_ENRICHED` or `UNKNOWN_DENOMINATOR` **cannot enter primary top-*k*
  metrics.** (Grounded: all 45 candidate-resolved patients carry ≥1 tested negative, so none are
  trivially positive-only; the gate's real work here is *naming the denominator type per study* and
  flagging ranking-degenerate patients, not culling.)

  > **Erratum (post-implementation, 2026-07-10):** the parenthetical grounding above is wrong.
  > 38 of 45 candidate-resolved patients carry ≥1 tested negative; 7 mKRAS 6/6-responders carry
  > none, so the shipped gate *does* cull those 7 from primary top-*k*. The implemented
  > `denominator_type` is binary and encodes only presence/absence of a tested negative
  > (rankability) — not the four categories listed above, and not denominator completeness or
  > selection bias. Shipped labels: `HAS_TESTED_NEGATIVE` / `NO_TESTED_NEGATIVE`. See
  > `artifacts/milestone_6/m6a_audit.md` and the `completeness_report` docstring.
- The report states plainly whether M6 is evaluating **"reranking among reported vaccine candidates"**
  (the honest current name) versus **"selecting the best 20 from the full generated candidate
  universe"** (not yet supported). Do not use the second name for the first experiment.

### WP1 — HLA resolution → presentation baseline *("resolve HLA first", scoped honestly)*

Per-study feasibility (from grounding):

| Study | HLA status | Presentation baseline |
|---|---|---|
| **hu** | candidate-level HLA + class present | ✅ ready now |
| **pdac** | per-epitope best-prediction MHC-I/II in `antigens.parquet`, not joined to candidates | ✅ bounded: antigen→candidate join; class-I minimal epitopes → MHCflurry |
| **mKRAS** | patient genotypes in source paper, untranscribed; shared *long* SLPs, class UNKNOWN | ⚠️ **stretch only** |
| **braun** | not public | ❌ **permanently prevalence-only** |

- Compute an MHCflurry `presentation_score` per candidate where HLA is resolved. Emit a per-study
  presentation-availability table.
- **hu + pdac form the primary presentation-comparable subset.** braun is HLA-agnostic (prevalence
  only). Any study whose HLA is not resolved within the WP1 time-box **falls back to prevalence-only,
  stated explicitly** — no uniform-baseline pretense.
- **mKRAS presentation is secondary and conditional.** Its long SLPs must be tiled into minimal
  epitopes, which introduces a multiple-testing bias (one long peptide → many short candidates →
  best predicted binder selected → longer peptides look artificially better). mKRAS presentation runs
  **only if its tiling + per-peptide aggregation rule is registered before results**, and is reported
  as a separately labelled sensitivity analysis, never folded into the primary presentation track.

### WP2 — Pre-vaccine feature matrix (shared core; missingness is not a study label)

Cross-study missingness can leak study identity. The **headline model uses a core feature set
available under equivalent definitions across all included studies.** Extended tiers are separate,
labelled analyses.

- **Core:** sequence + universally available biochemistry (length, composition, hydrophobicity,
  charge). Available for all 4 studies.
- **Core + contrastive:** mutant-vs-WT deltas (`hydrophobicity_delta`, `charge_delta`, anchor/TCR-face
  counts). See WP2 class-gating below.
- **Core + presentation:** MHCflurry score, on the HLA-resolved compatible subset only.
- **Extended (secondary only):** study-dependent expression/clonality, source-native scores, PLM
  embeddings.

Rules: use `infer_numeric_feature_columns` with the non-feature exclusion list; drop zero-variance
columns (e.g. constant assay_type); **`study_id`/`cancer_type`/`platform` are never predictive
features.** Missingness indicators are **diagnostic first**, predictive only in a clearly labelled
secondary analysis — a field being absent must not act as a hidden study label.

**Class-gated anchor features (change 5):** class-I anchor vs TCR-facing positions are derived from
the binding register only where the register is defined. **In the required M6A ladder, class-II
anchor/TCR-face features are marked `unavailable` (not zero-filled)** — the class-II binding core can
shift and no frozen class-II core predictor is wired in. Wiring one (e.g. a NetMHCIIpan-style core
assignment) is a conditional extension, not a required deliverable. Never apply the class-I position
scheme to class-II peptides blindly. `unavailable` is represented so that missingness cannot act as a
hidden class/study label (see missingness rule above).

### WP3 — Model ladder (small required set; complexity is opt-in after freeze)

**Required ladder (all of M6A):**

- **B0 — prevalence:** training positive rate.
- **B1 — presentation-only:** MHCflurry score, on the presentation-compatible subset (WP1).
- **M1 — regularized logistic regression** on core (+ contrastive/presentation tiers).
- **M2 — shallow gradient boosting** (low-depth HistGradientBoosting / XGBClassifier).

**Conditional, only after M1/M2 are frozen and only if addable without model-selection leakage:**
ESM paired-delta features, a calibrated two-head ensemble (M3), and richer combinations. With 45
patients, architecture shopping is a bigger risk than underfitting — extensions must earn their way
in against a frozen baseline, not be assumed deliverables.

Class I vs II: **one shared model with MHC class as a feature** (n is too small to split cleanly);
class-stratified performance is **mandatory** in reporting. Calibration policy is WP6.

### WP4 — Evaluation (two explicit headline tracks; no mixed baseline)

Do **not** pool "strongest available baseline per fold" — a fold-varying opponent (hu vs
presentation, braun vs prevalence) is uninterpretable. Two separate, explicitly named tracks:

1. **Universal track** — all 4 candidate-resolved studies, **learned model vs prevalence.** Answers
   "does the model learn anything useful across the whole corpus?" This is *not* "beating
   presentation."
2. **Presentation track** — **hu + pdac only**, **learned model vs frozen presentation-only.** The
   real recognition-vs-presentation test. mKRAS added only as a labelled sensitivity analysis if WP1's
   long-peptide baseline passes validation.

Source-native rankings (TESLA/PRIME/IMPROVE-RF/NetMHCpan where compatible) are reported as **secondary
reference points**, with compatibility stated explicitly — never part of a primary mixed baseline.

**Cross-validation:** leave-one-study-out (LOSO), 4 folds via a rotation wrapper over
`STUDY_HOLDOUT`. LOSO is the headline precisely because it punishes study-memorization. Patient-grouped
CV within training studies is secondary, never headline.

**Top-*k* metric definition (change 3):** `k_patient = min(20, n_eligible)`. Report `hits@k_patient`,
`precision` with denominator `k_patient`, capture fraction, P(≥1 hit) — never divide by 20 for a
patient with fewer candidates. Retain a **literal fixed-@20** metric on the ~11 patients with ≥20
candidates as a labelled sensitivity analysis (noting its hu-heavy composition). Because selection is
degenerate where `n ≤ k`, also report ranking metrics split by `ranking_informative`, and run a
per-peptide **recognition-classification** analysis (AUROC / average precision / Brier / calibration,
LOSO) on all candidate-resolved patients as the honest home for signal that the thin selection
question cannot show. **AUROC is explicitly not the headline.**

**Reporting geometry (change 9):**
- **macro-LOSO** (equal weight per held-out study) — **primary interpretation.**
- **micro-LOSO** (all held-out patients pooled) — reported beside it.
- **per-study** — every fold shown separately; name which study drives any result.
- **Bootstrap at the patient level** (20,000 resamples), resampling patients within each fold; for
  pooled analyses preserve study structure / report study-stratified patient bootstraps. With only 4
  studies, do not pretend study-level uncertainty is precisely estimated.

**Metric set to add (honoring the spec; harness `scorecard` lacks these):** P(≥2), P(≥4), AUROC,
Brier, calibration/reliability — as diagnostics alongside the existing hits@k / precision@k / P(≥1) /
capture / MRR.

### WP5 — Study-confound diagnostics (the highest-risk item)

The biggest risk is that the model learns study design, not recognition. Required:
- **Study-only classifier:** can the feature matrix predict `study_id`? Quantifies confound strength.
- **Feature distribution + label prevalence by study.**
- **Study-ablation** and **study-balanced training weights.**
- LOSO itself is the structural defense; these diagnostics quantify how much to trust it.

### WP6 — Calibration policy (restraint)

Separate class-I/class-II calibration is valid **only** if each class has enough training patients and
positives *inside each LOSO fold*, fit entirely within training data via nested grouped splits.
Otherwise use **one global training-only calibrator**, or make **no calibrated claim** for that class.
Never fit isotonic on tiny subsets. Class-stratified *performance* remains mandatory regardless.

### WP7 — Registered success gate (locked before scores are generated)

No vague "meaningful regression." Explicit, pre-registered criteria:

Let Δ = macro-study mean `hits@k_patient` (learned − opponent), with a 90% patient-level bootstrap
interval `[Δ_lo, Δ_hi]`. Proposed default rule (confirm or override at kickoff, before any score):

- **`ACCEPT`** — universal track: `Δ > 0` **and `Δ_lo > 0`** (interval excludes zero), **and** no harm
  margin crossed. A positive point estimate with an interval spanning zero is explicitly **not** a win.
- **`CONSISTENT_WITH_NO_EFFECT`** — interval spans zero and no harm margin crossed.
- **`REJECT`** — `Δ_hi < 0`, **or** a harm margin is crossed.
- **Harm margins (proposed default):** at the 90% patient-level lower bound, learned **P(≥1)** is not
  more than **0.05 absolute** below the opponent, and **capture fraction** not more than **0.05
  absolute** below. Crossing either → `REJECT` regardless of Δ.
- **Presentation claim (separate verdict):** the same rule on the **hu+pdac** subset, learned vs
  presentation-only. Reported independently; it never substitutes for the universal verdict.

Report every fold + macro + micro + CIs + which study drives, plus an explicit "what evidence would
change this answer" section.

### M6B — Event-A → Event-B transfer (auxiliary, own gate)

Built and reported **after M6A is frozen.** Valid mechanisms only: representation pretraining,
multitask learning with **separate Event-A / Event-B heads**, or a **frozen Event-A teacher as one
pre-vaccine feature** (via `transfer_ranker`). **Do not merge Event-A and Event-B labels.** The arm
answers one question with its own delta table:

> Does Event-A information improve Event-B transfer **beyond the Event-B-only model?**

It enters no primary success gate; if it earns a secondary gate, that gate is declared in advance.

### M6C — Osteosarc replay (contract document only)

A tiny contract doc, no run (no data/harness exists — the adapter is a declaration-only stub). Define
only: timeline freeze points, source mappings, vaccine targets, ELISPOT outcomes, TCR evidence, and
replay metrics. It must **not** become another week of infrastructure before M6A runs.

---

## Registered decisions (frozen before training)

| Decision | Value |
|---|---|
| Corpus | `outputs/event_b_backbone/combined/` (Event-B only) |
| Population | 45 candidate-resolved patients, 4 studies; Nous excluded; 9 UNTESTED dropped |
| Label | POSITIVE=1 vs TESTED_NEGATIVE=0 |
| Primary CV | leave-one-study-out (4 folds) |
| Primary metric | macro-study mean `hits@k_patient`, `k_patient=min(20,n_eligible)` |
| Primary opponent (universal) | prevalence (B0) |
| Presentation opponent | presentation-only (B1) on hu+pdac only |
| Required models | B0, B1, M1 (logistic), M2 (shallow boosting) |
| Bootstrap | 20,000 resamples, patient-level grouping |
| Verdict rule | `ACCEPT` iff `Δ>0 ∧ Δ_lo>0 ∧ no harm`; `REJECT` iff `Δ_hi<0 ∨ harm`; else `CONSISTENT_WITH_NO_EFFECT` |
| Harm margins (P≥1, capture) | proposed default −0.05 absolute at 90% lower bound; **confirm/override at kickoff** |
| Seed / sort / tie-break | 17 / mergesort / md5 |

Every threshold has a concrete proposed default above; the only value flagged for explicit
confirmation at implementation kickoff (before any model is scored) is the harm-margin magnitude.

---

## Artifacts, conventions, testing

- **Reusable code** → `src/epicurus_neo/`: LOSO rotation wrapper, the Event-B feature-matrix builder,
  the completeness-gate + k_patient scorer, the added metrics (P≥2/≥4, AUROC/Brier/calibration), a
  regularized-logistic model. **Thin runner** → `experiments/m6_recognition_swing.py`.
- **Outputs** → `outputs/` (scored CSV + JSON sidecar per track/model/fold). **Headline audit** →
  `artifacts/milestone_6/` (JSON + Markdown), carrying the `INSUFFICIENT_CANDIDATE_RESOLVED` caveat.
- **Determinism:** seed 17, mergesort, md5 tie-break; reproduce byte-identical outputs on re-run.
- **Tests:** leakage guards (banned columns never reach features; `information_timing` enforced),
  LOSO-rotation correctness (each study held out exactly once), k_patient / completeness-gate
  correctness on a toy, added-metric correctness on a toy, deterministic output hashes.

## Out of scope (deferred to M7–M10)

TCR-aware recognition, portfolio/selection optimizer, abstention, prospective validation, and any new
external ingestion beyond IMPROVE (already present, used only in M6B).

## Expectation, stated plainly

At n=45, with study fully confounded with cancer type / platform / prevalence, per-patient candidate
counts this small, and a genuine top-*k* selection question supported for only ~11 patients, honest
LOSO most likely returns wide CIs / `CONSISTENT_WITH_NO_EFFECT`. The deliverable is a leakage-safe,
hard-to-game answer, the confound diagnostics that qualify it, and a precise statement of the missing
variable — under the standing insufficiency verdict.
