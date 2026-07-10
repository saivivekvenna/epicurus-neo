# Milestone 6B audit: Event-A -> Event-B transfer (auxiliary)

**Corpus verdict (standing):** `INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA`

M6B is auxiliary: an Event-A -> Event-B transfer probe that enters no primary success gate. Diagnostic swing under the standing insufficiency verdict; not a headline claim.

**Question:** Does Event-A information improve Event-B transfer beyond the Event-B-only model?

## Declared-gate verdict: **REJECT_TRANSFER**
- Macro per-fold AUROC delta (candidate - baseline): -0.0204
- Folds improved: 1 / 4 (ACCEPT_TRANSFER needs >=3 and a positive macro delta and no harm)
- Macro Δ hits@k (reported, underpowered): 0.0000 CI [-0.0938, 0.0938]
- Ranking-informative patients: 8

## Per-held-out-study AUROC
Baseline = frozen M6A `logistic(core)`; candidate = `logistic(core + event_a_teacher_score)`.

| Study | baseline | candidate | delta | n_eval |
|---|---|---|---|---|
| braun_rcc_2025 | 0.4785 | 0.4682 | -0.0104 | 129 |
| hu_neovax_2021 | 0.4358 | 0.4376 | 0.0019 | 541 |
| mkras_vax_2026 | 0.6750 | 0.6750 | 0.0000 | 72 |
| pdac_neovax_2023 | 0.5576 | 0.4846 | -0.0730 | 223 |

## Teacher
- Frozen Event-A teacher: `logistic` on the `core` tier, trained on 17082 IMPROVE Event-A rows (458 positive). Labels never merged; the teacher never sees an Event-B row.
- Event-A is short class-I (8-11mer); most Event-B is long SLP. The teacher score is added as one feature to the Event-B-only model and is the sole candidate-vs-baseline difference.
