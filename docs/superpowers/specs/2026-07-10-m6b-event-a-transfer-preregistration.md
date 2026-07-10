# M6B — Event-A → Event-B transfer: pre-registration (declared gate)

**Status:** frozen before implementation, 2026-07-10. Executes the registered M6B arm of
`2026-07-10-m6-first-recognition-swing-design.md` (§ "M6B — Event-A → Event-B transfer"). M6B is
**auxiliary**: it enters no primary success gate. This document declares its secondary gate *in
advance*, as the parent spec requires ("if it earns a secondary gate, that gate is declared in
advance").

Standing corpus verdict `INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA` still applies; M6B is a
diagnostic swing, not a headline claim.

## The one question

> Does Event-A information improve Event-B transfer **beyond the Event-B-only model?**

## Mechanism (as implemented)

A **frozen Event-A teacher**, trained **only on Event-A labels**, emits one recognition score per
Event-B candidate. That score is added as **one feature** to the Event-B-only model. Labels are never
merged: the teacher sees only Event-A labels, the student sees only Event-B labels.

- **Teacher:** the same sklearn stack and the same **length-agnostic M6A core feature tier** used in
  M6A, fit on the IMPROVE Event-A frame (`outputs/event_b_corpus_combined/`, `study_id == improve`,
  `event_type == EVENT_A_PREEXISTING_REACTIVITY`): **17,082 rows, 458 POSITIVE / 16,624
  TESTED_NEGATIVE, all CLASS_I, peptides 8–11mer**. Trained once, frozen (it never sees any Event-B
  row, so it is identical across all LOSO folds).
- **Student / candidate model:** `logistic(core + event_a_teacher_score)` on Event-B.
- **Baseline:** `logistic(core)` — the frozen M6A universal model. **The teacher feature is the only
  difference between candidate and baseline.**
- **Population:** the M6A universal track — all 4 Event-B studies, leave-one-study-out (4 folds).
  `NO_TESTED_NEGATIVE` patients are excluded from primary top-k (same completeness gate as M6A);
  all candidates are retained for pooled/per-fold classification.

## Registered deviation from the parent spec

The parent spec names `transfer_ranker` as the teacher vehicle. **We deliberately do not use it for
this first swing** (though `xgboost` 3.3.0 and `transfer_ranker` are both available):

1. **Clean isolation.** `transfer_ranker`'s features are mhcflurry- and PLM-based
   (`COMMON_TRANSFER_FEATURES`). Using them, a positive M6B could be "richer features help," not
   "Event-A transfer helps." Using the *same* M6A core features for teacher and student makes the
   added `event_a_teacher_score` a **pure** transfer signal — the sole changed variable.
2. **Universal applicability.** mhcflurry features are class-I / ≤15mer only, which would restrict
   Event-B to a thin subset and reproduce the presentation track's `NOT_VIABLE` collapse. The M6A
   core features are **length-agnostic**, so the teacher scores every Event-B candidate, including
   long SLPs — matching the M6A universal population exactly.
3. **Consistency.** It slots into the existing M6A sklearn LOSO harness with no new model family.

A richer `transfer_ranker`-based teacher is a legitimate *future* M6B variant; it is out of scope for
this first, deliberately-boring swing.

## Declared secondary gate (frozen)

Primary transfer signal is **discrimination (AUROC)**, because at n=45 with only ~8 ranking-informative
patients, macro hits@k is underpowered; AUROC is defined on every candidate in every fold.

- **`ACCEPT_TRANSFER`** iff **all** hold:
  - macro-mean per-fold AUROC delta (candidate − baseline) **> 0**, and
  - AUROC improves on **≥ 3 of 4** held-out studies (a majority-plus of independent studies, guarding
    against a single-study fluke — the failure mode that made one study drive all of M6A), and
  - **no** significant harm on the ranking-informative subset: macro Δ hits@k and Δ P(≥1) each above a
    −0.05 absolute margin at the 90% lower bound (same margin as M6A).
- **`REJECT_TRANSFER`** iff macro AUROC delta < 0 with ≥3/4 folds worse, or significant harm.
- else **`CONSISTENT_WITH_NO_EFFECT_TRANSFER`**.

Anything other than `ACCEPT_TRANSFER` → **freeze M6B and move on.** No retries, no threshold shopping.
macro Δ hits@k is reported for completeness but is **not** required to be positive (underpowered).

## Leakage discipline

- The teacher never sees any Event-B row (frozen, trained on Event-A only).
- **Guard:** assert no shared `(mutant_peptide, hla_allele)` between the Event-A training frame and any
  Event-B evaluation row; raise if found.
- Banned-column discipline is inherited from `m6/features.py` unchanged.

## Data provenance / determinism

- Event-A: `outputs/event_b_corpus_combined/{candidates,assays}.parquet` (processed, <2MB).
- Event-B: `outputs/event_b_backbone/combined/` (the pinned M6A corpus).
- **No raw data is touched** (the terabyte raw-sequencing pipeline is conceptual; nothing GB–TB exists
  in-repo). Seed 17, mergesort, md5 tie-break, 20,000-resample patient-level bootstrap.

## Expected outcome, stated plainly

The teacher is trained on short class-I neoepitopes and applied largely to long class-II SLPs. Honest
transfer is unlikely to clear a ≥3/4-fold AUROC gate; `CONSISTENT_WITH_NO_EFFECT_TRANSFER` is the
most probable result. The deliverable is a leakage-safe, pre-registered answer to whether Event-A
recognition transfers to Event-B at all — not a win.
