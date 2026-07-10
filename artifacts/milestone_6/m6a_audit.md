# Milestone 6A audit: Event-B-only recognition swing

**Corpus verdict (standing):** `INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA`

Diagnostic swing under the standing insufficiency verdict; not a headline claim.

## Universal track (learned vs prevalence (all 4 studies))
- Verdict: **REJECT**
- Macro-study Δ hits@k_patient: -0.2188 CI [-0.4688, 0.0625]
- Macro AUROC (per-study mean): 0.5367 | Brier: 0.4227
- Per-held-out-study AUROC: braun_rcc_2025=0.4785, hu_neovax_2021=0.4358, mkras_vax_2026=0.6750, pdac_neovax_2023=0.5576
- Pooled OOF AUROC (caveated): 0.3165 — conflates cross-study calibration shift with discrimination
- Ranking-informative patients (n_eligible > k): 8

## Presentation track (learned vs presentation-only (hu + pdac))
- Verdict: **NOT_VIABLE_PRESENTATION_INCOMPATIBLE**
- Macro-study Δ hits@k_patient: nan CI [None, None]
- Macro AUROC (per-study mean): nan | Brier: nan
- Ranking-informative patients (n_eligible > k): 0
- Reason: MHCflurry (class-I, <=15mer) scores only short peptides; long SLPs are incompatible, so fewer than two studies retain a viable presentation-comparable subset (the comparison would collapse to one study).
- Presentation-compatible candidates: hu_neovax_2021=207 (18 pos), pdac_neovax_2023=4 (0 pos)

## Study confound
- Study-only classifier accuracy: 0.9627 (majority rate 0.5606)
- Positive rate by study: braun_rcc_2025=0.4729, hu_neovax_2021=0.2366, mkras_vax_2026=0.8333, pdac_neovax_2023=0.1031
