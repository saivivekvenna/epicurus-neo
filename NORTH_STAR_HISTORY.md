# NORTH_STAR_HISTORY — canonical project ledger

Single source of truth for *where the project stands*: the north star, what has been tried, what was
falsified, what is frozen/preserved, which datasets are consumed, and what the next gate is. This file **links
to** the artifacts that hold the actual results; it does not restate or re-derive them. Keep it updated as the
last step of every milestone.

_Last updated: 2026-07-11 (v0.5 deployable context-conditioned pairwise — REJECT/TIE with PRIME (Δ −0.012); closes the assay/regime-head lever: deployable context matches the v0.4 tower's aggregate with NO study identity, but recovers only ~27% of the Gartner edge and adds ~nothing over the pairwise objective → five straight parity results prove the wall is DATA, not model form)._

---

## North star (unchanged)

Ship an open-source tool that takes a patient's WES + RNA-seq and returns a ranked list of their
biologically-eligible neoantigens that **ranks better than genuine GfellerLab PRIME**, proven on **untouched
external patients**. Inputs expandable only to WES/RNA-accessibility parity. The wet-lab loop is a label
factory *behind* the tool. Presentation is largely solved; **recognition is the hard wall.**

**Verdict rule (pre-registered, applies to every comparison):** no superiority claim until a *frozen* model
wins **ACCEPT** (bootstrap Δ CI lower bound > 0) vs genuine PRIME **and** the strongest presentation baseline,
on patients no part of the pipeline has touched. `CONSISTENT_WITH_NO_EFFECT` (positive point estimate, CI spans
0) and `REJECT` both mean *superiority not established* — never "proven equal", never "significantly worse".

---

## Chronological ledger

| when | milestone | outcome | verdict | artifact |
|---|---|---|---|---|
| M6 | Recognition swings (M6A/M6B) | Both honest negatives; binding constraint = data/independent studies, not model richness | REJECT | memory `m6-recognition-swings-outcome` |
| — | LLM biology-only triage (20 Sonnet agents) | Rank neoantigens at ≈chance, worse than PRIME; quantifies the recognition wall | REJECT | `experiments/llm_neoantigen_triage/` |
| — | Data-quality measurement | Per-patient reconstructability scorecard; corpus was a *label* corpus not a *decision-problem* corpus | diagnostic | `docs/data_quality_measurement.md`, `audit-data-quality` CLI |
| — | NeoVax reconstruction scaffold | WES/RNA behind dbGaP phs001451 (controlled); strict gate = all NOT_DECISION_READY | BLOCKED | `docs/neovax_reconstruction.md`, `audit-reconstruction` CLI |
| — | Zhao DC 2026 ingest | 352 pts / 2317 SNV; frozen PRIME-residual benchmark = honest NULL vs MixMHCpred | REJECT | memory `zhao-dc-2026-asset` |
| — | CEDAR recognition asset | 32,465 cancer T-cell rows normalized; 716-peptide backbone overlap = leakage flag | asset (not trained) | memory `cedar-recognition-asset` |
| — | Naive recognition-row augmentation | Falsified (Δhits@2 −0.026 REJECT; CEDAR cross-study transfer AUROC 0.478 below chance) | REJECT | memory `prime-augmentation-outcome` |
| M7 | Real within-patient decision benchmark (Müller NCI ELISpot + 2025 multimer) | Presentation near-solved on top-20; NO orthogonal recognition feature beats it, several hurt | REJECT (all features) | `artifacts/milestone_7_decision/decision_benchmark_report.md`, `decision-benchmark` CLI |
| M7 | Ascertainment correction | "presentation 0.985 near-solved" was a PU artifact; Müller VALIDATED=0 = UNTESTED not tested-neg | correction | memory `m7-ascertainment-correction` |
| M7 | Pre-registered frozen transfer residual (multimer-trained) on Gartner + IMPROVE | Does NOT beat genuine PRIME on either untouched cohort; PRIME edge is cohort-dependent, not universal | REJECT / CONSISTENT | `artifacts/milestone_7_decision/prime_transfer/`, memory `transfer-and-frozen-outcome` |
| M7 | Genuine PRIME 2.1 + MixMHCpred 3.0 head-to-head | PRIME wins on multimer; on frozen Gartner NCI genuine PRIME 0.729 does NOT beat presentation (EL 0.776) | honest negative | `artifacts/milestone_7_decision/prime_headtohead/`, memory `prime-genuine-headtohead` |
| — | **Epicurus v0.2** development ranker | Family-weighted group-balanced pairwise ranker REJECTED in development; v0.1 stays frozen; provenance-leakage guard required | REJECT (dev) | `artifacts/milestone_7_decision/epicurus_v02/MODEL_FAILURES_AND_V02.md`, `configs/frozen/epicurus_v0_2_dev.json` |
| — | **NCI crosswalk data-first audit** | "24 positives / data-bound" was an ARTIFACT of omitting Gartner TRAIN; true dev picture ~166 bag-pools / ~640 positives; constraint = fragmentation + MIL granularity + selection bias + no pristine external proof | data unlock | `artifacts/milestone_7_decision/nci_crosswalk/DATA_DIAGNOSIS.md`, memory `nci-crosswalk-data-unlock` |
| — | **TEST lineage fix + MIL dev split** | Reconciled Mmps propagation to raw exact-key ascertainment (47→46 CD8 bags, 1 conflict); froze a patient/family/study-blocked development split | groundwork | `artifacts/milestone_7_decision/nci_crosswalk/MIL_DEV_SPLIT.md`, `configs/frozen/mil_dev_split_v1.json` |
| — | **Epicurus v0.3** MIL development experiment | Preregistered nested-OOF MIL ranker on the frozen split. **REJECT** (no ACCEPT) but the registered candidate is **statistically TIED** with genuine PRIME (Δhits@20 −0.093, CI[−0.228,+0.051] spans 0) — NOT worse. Presentation alone also ties PRIME; modeling adds **nothing over presentation** (stop rule: no promotion past rung1). Per-source: ties/edges PRIME on Gartner (+0.05) & multimer (+0.056), loses on IMPROVE (−0.38). Orthogonal recognition features add nothing; augmentation hurts. A source-EL-semantics bug invalidated the FIRST run (its "−0.25 significantly worse" was an artifact). | REJECT (tie) | `artifacts/milestone_7_decision/epicurus_v03/{PREREGISTERED_PROTOCOL.md,DEV_REPORT.md}`, `configs/frozen/epicurus_v0_3_dev.json` |
| — | **Epicurus v0.4** source-aware tower | Preregistered partial-pooling MIL tower (shared backbone `w₀` + shrunk per-source feature heads `v_s` + rank-inert intercepts `c_s`; λ nested-selected) on the SAME frozen split/gate. **REJECT** — F still **TIED** with genuine PRIME (Δ −0.016, CI[−0.147,+0.112]), though CLOSER than pooled v0.3 (−0.093). **New:** the tower **recovers the predicted Gartner edge** (F beats PRIME on Gartner Δ**+0.275**, up from pooled +0.05); the **negative control is clean** (shuffled-source F−P −0.232 vs true +0.077 → the lift is genuine source structure, NOT capacity). Mechanism directionally right but **underpowered**: F−P +0.077 [−0.073,+0.208] and F−C +0.086 [−0.080,+0.227] positive (feature-weighting, not calibration: C−P −0.009), all CIs span 0. IMPROVE (−0.10) and multimer (−0.222, head over-specializes on n=18) cancel Gartner in the aggregate. Recognition wall persists (gains ride on presentation-adjacent features). | REJECT (tie; Gartner edge real) | `artifacts/milestone_7_decision/epicurus_v04/{PREREGISTERED_PROTOCOL.md,DEV_REPORT.md,PROVENANCE.json}`, `configs/frozen/epicurus_v0_4_dev.json` |
| — | **Epicurus v0.5** deployable context-conditioned pairwise | Preregistered convex within-patient pairwise ranker; **no study/assay identity** — observable-only context (peptide-length, HLA-locus, predictor-disagreement) × presentation features, nested-λ, same frozen split/gate. **REJECT** — R still **TIED** with genuine PRIME (Δ −0.012, CI[−0.125,+0.108]); closest yet. **Three findings:** (1) the *pairwise objective* is the workhorse (Q−P **+0.069**, biggest single move) — *context adds ~nothing* (R−Q **+0.011**, CI spans 0); (2) *deployable ≈ non-deployable* v0.4 tower in aggregate (R−F **+0.004**) → **no deployability penalty**, BUT the deployable Gartner edge is only **+0.075** (1.05 v 0.975) ≈27% of v0.4's +0.275 → **most of the v0.4 Gartner win was parasitic on study identity**; (3) context×presentation interactions all **|β|≤0.02 = presentation reweighting, NOT a new recognition axis**. Per-source R vs PRIME: Gartner +0.075, IMPROVE −0.167, multimer +0.056. Multi-init stability now **PASSES** (coef Δ 2.1e-7, Spearman 1.0 — v0.4's failure fixed). **Closes the assay/regime-head lever.** | REJECT (tie; deployable ceiling = tower) | `artifacts/milestone_7_decision/epicurus_v05/{PREREGISTERED_PROTOCOL.md,DEV_REPORT.md,PROVENANCE.json}`, `configs/frozen/epicurus_v0_5_dev.json` |

---

## Frozen / preserved models

- **Epicurus v0.1** — the model of record. Serialized coefficients + data checksums,
  `configs/frozen/epicurus_v0_1.json`, applied via `score_with_frozen`. Untouched.
- **Epicurus v0.2** — `configs/frozen/epicurus_v0_2_dev.json` (`status: REJECTED_DEVELOPMENT`,
  `supersedes_frozen: false`). Kept as a documented negative, not shipped.
- **MIL dev split v1** — `configs/frozen/mil_dev_split_v1.json` (`FROZEN_DEVELOPMENT_SPLIT`, no model fitted).
- **Epicurus v0.3** — `configs/frozen/epicurus_v0_3_dev.json` (`status: REJECTED_DEVELOPMENT`,
  `supersedes_frozen: false`). Ties PRIME in development but does not beat it; v0.1 stays the model of record.
- **Epicurus v0.4** — `configs/frozen/epicurus_v0_4_dev.json` (`status: REJECTED_DEVELOPMENT`,
  `supersedes_frozen: false`). Source-aware tower; a **mechanism test with dataset-name heads → NOT deployable/
  generalizable**. Recovers a real Gartner edge but still only ties PRIME in aggregate; v0.1 stays frozen.
- **Epicurus v0.5** — `configs/frozen/epicurus_v0_5_dev.json` (`status: REJECTED_DEVELOPMENT`,
  `supersedes_frozen: false`). Deployable context-conditioned pairwise — **the deployable version of v0.4's
  heads**. Matches the tower's aggregate with observable-only context (no study identity) and fixes the
  stability failure, but recovers only ~27% of the Gartner edge (rest was study-identity-bound) and still ties
  PRIME. v0.1 stays the model of record.

---

## Current data roles (deduplicated)

Source of truth for counts: `artifacts/milestone_7_decision/nci_crosswalk/{CROSSWALK_AUDIT.json,
DATA_ROLE_MATRIX.md}`. The **assay unit** for Gartner is the transcript 25mer `key`; the **leakage-blocking
unit** is `genomic_mutation_family_id` (patient + raw Variant key). Anti-inflation is non-negotiable: negatives
count at the BAG level, never the ~285k/1.09M child rows.

| source | granularity | positives | role |
|---|---|---|---|
| Gartner TRAIN × Müller | mutation bag (transcript key) | 82 exact / 111 CD8 bags w/ features (139 total CD8 transcript / 137 genomic) | **DEVELOPMENT** (in MIL split) |
| IMPROVE | pMHC-multimer instance | 467 | **DEVELOPMENT** (leakage-clean, in MIL split) |
| CD8 multimer | instance (orthogonal features) | 34 | **DEVELOPMENT** (only source with expr/agretopicity/foreignness) |
| Gartner TEST (Mmps, by key) | instance over transcript bag | 27 exact / **46** CD8 transcript bags (raw-reconciled from Mmps' propagated 47) | **HOLDOUT — SEMI-CONSUMED** (aggregate metrics already seen; never fit/tuned) |
| Zhao DC | tiny patient pools (≤20) | 313 | AUXILIARY (AUROC/nDCG only; hits@20 non-informative) |
| CEDAR | PMID recognition buckets | ~7,491 | AUXILIARY recognition only (not patient-level reranking; PRIME-training provenance confound) |
| PRIME/BigMHC TableS4 | training set (in-sample) | 596 | EXCLUDED from fitting a PRIME residual (leakage reference) |
| Hu/Ott NeoVax (dbGaP phs001451) | full WES+RNA denominator | measured POS/TESTED_NEG | **EXTERNAL PROOF — PENDING acquisition** |

---

## Falsified / failed-to-establish (do not retry without a new lever)

- Learned PRIME-residual (multimer-trained, frozen) — no win on Gartner or IMPROVE.
- No **universal** recognition edge over presentation; PRIME-vs-EL direction flips by cohort.
- Naive recognition-row augmentation; CEDAR cross-study transfer (below chance); LLM biology-only triage.
- Epicurus v0.2 ranker in development (data-bound, not model-bound; orthogonal features at/below chance on the
  only denominator cohort).
- Epicurus v0.3 MIL ranker — **ties** genuine PRIME in development but does not beat it (no ACCEPT); modeling
  adds nothing over a source-correct presentation baseline; orthogonal recognition features add nothing;
  cross-source augmentation hurts. (Not "worse than PRIME" — a parity/tie.)
- Epicurus v0.4 source-aware tower — partial pooling **recovers the Gartner edge that naive pooling diluted**
  (Δ+0.275 vs PRIME on Gartner) and the lift is **genuine source structure, not capacity** (clean negative
  control), BUT the aggregate still only **ties** PRIME (IMPROVE/multimer cancel it) → no ACCEPT. The mechanism
  (source-conditioned feature weighting **beyond** prevalence calibration) is directionally confirmed but
  **underpowered** (all mechanism-contrast CIs span 0). Recognition wall unchanged (gains are presentation-side).
  NOTE: heads keyed on dataset **name** — a mechanism result, not a deployable model; do NOT re-run the same
  dataset-name tower expecting a different verdict — the next lever is power/assay-regime heads, not re-tuning.
- Epicurus v0.5 deployable context-conditioned pairwise — **closes the assay/regime-head lever.** Portable,
  inference-time context (peptide-length, HLA-locus, predictor-disagreement) reproduces the v0.4 tower's
  *aggregate* (R−F +0.004) with NO study identity, and the convex **pairwise objective — not the context —**
  carries it (Q−P +0.069 vs R−Q +0.011). But it recovers only +0.075 of the +0.275 Gartner edge (the rest was
  dataset-name-bound), the context×presentation interactions are presentation **reweighting not a new
  recognition axis** (|β|≤0.02), and it still **ties** PRIME (−0.012). → Re-keying heads on observable regime is
  DONE and does not clear the gate. Model form is now exhausted on these 3 cohorts; the remaining levers are
  **DATA** (power + a new denominator cohort + pristine external proof).
- **Source-EL-semantics trap (fixed):** Gartner/Müller `Score_EL` is a 0-1 likelihood (higher-better), NOT a
  %rank like IMPROVE/multimer EL. Conflating them inverted Gartner presentation and invalidated the v0.3 first
  run (see `epicurus_v03/DEV_RESULT.json` `invalidated_first_run`). Always check per-source score orientation.
- **Provenance-leakage trap (fixed):** PRIME-training membership is label-correlated in CEDAR → the reliability
  flag is excluded from features and models fit only on prime-reliable rows.
- **Lineage-propagation trap (fixed):** `MmpsTestingSet_extract` propagates a genomic-family CD8+ onto
  unscreened sibling transcripts (RPS13/4359). Raw `NmersTestingSet` exact-key ascertainment wins.

---

## Where we are now → next gate

- **State — model levers are exhausted; the constraint is data.** v0.5 tested and **closed lever #1**
  (deployable assay/regime conditioning). A convex within-patient pairwise ranker with observable-only context
  **ties** genuine PRIME (Δ −0.012, CI[−0.125,+0.108]) — the closest yet, with **no deployability penalty**: it
  matches the v0.4 dataset-name tower's aggregate (R−F +0.004) using only what a new patient would have. Three
  sharp findings: (a) the **objective is the lever, not the conditioning** — fair within-patient pairwise
  ranking (Q−P +0.069) recovers ~90% of the gain, context adds ~nothing (R−Q +0.011); (b) the deployable Gartner
  edge is only +0.075 vs v0.4's +0.275 → **most of v0.4's Gartner win required knowing the study name**
  (non-portable); (c) context×presentation interactions are pure **presentation reweighting** (|β|≤0.02), **not a
  new recognition axis**. **Net: five straight development candidates (v0.2→v0.5), zero ACCEPTs — every model
  converges to PRIME parity.** The *invariance* of that result across MIL vs pairwise objectives and pooled vs
  source-tower vs context conditioning is itself the finding: the binding constraint is **not model form** — it
  is **data**. What's missing is (i) independent Gartner-like **decision-problem** cohorts (dense per-patient
  candidate denominator + explicit **tested-negatives** + WES/RNA/HLA + a T-cell recognition assay), (ii)
  **statistical power** (n=118 across 3 heterogeneous sources can't establish the directionally-real effects),
  and (iii) a **still-pristine external proof cohort**.
- **Indicated levers (model form now exhausted on these 3 cohorts):**
  1. ~~Assay/regime heads~~ — **DONE (v0.5); closes at a tie.** Do NOT iterate v0.6 model form against the
     current 3 cohorts expecting ACCEPT — v0.3/v0.4/v0.5 prove parity is data-bound, not a tuning miss.
  2. **DATA — acquire a second Gartner-like denominator cohort** (real per-patient candidate ranking + explicit
     immunogenicity negatives + WES/RNA). Now the single highest-value action; see the data-acquisition search
     brief (`artifacts/milestone_7_decision/external_validation/`).
  3. **External proof** — acquire Hu/Ott via dbGaP phs001451
     (`artifacts/milestone_7_decision/external_validation/ACQUISITION_PACKET.md`).
- **Not yet spent:** one pre-registered look at the SEMI-CONSUMED Gartner TEST holdout is still available, but
  only *after* a development candidate clears the gate — none of v0.3/v0.4/v0.5 did, so the holdout stays closed.
  **A superiority claim still requires a win on a genuinely untouched cohort.**

## Maintenance

Update this ledger (row in the chronological table + any role/verdict change) as the final step of each
milestone. Link artifacts; never paste result tables here — they drift. Memory index: `memory/MEMORY.md`.
