# Milestone 1 frozen-score re-audit

This table is generated from the scored out-of-fold artifacts that remain on disk from iterations
001–029. Every measured difference uses the same patient rows, the mandatory identity-hash tie-break,
and 20,000 paired bootstrap resamples. No score column was selected and no model was fit.

The old Gartner artifact does not retain the iteration-002 pre-change learned score, so its measured
comparison is against the named Gartner numeric baseline; the differing claimed comparator is shown
explicitly. The development RF/Epicurus blend is reconstructed from its frozen 50/50 rule and audited
against RF, not promoted as a result.

Two build-spec numbers need explicit interpretation. The analytic expectation over uniformly random
rankings is `0.5818`, while the one fixed ordering produced by the mandated ascending MD5 tie-break is
`0.3143`; the leakage canary asserts that the latter is order-independent and does not reproduce the
source-order leak (`2.4714`). The registered current-n MDE (`0.237`) is the two-sided 95% CI exclusion
threshold (`1.96 × paired SE`). Prospective `n_required` calculations separately use the requested
80% power and reproduce the plan's sample-size table.

| iteration | artifact | claimed Δ | measured Δ | 95% CI | verdict |
|---:|---|---:|---:|---:|---|
| 002 | `gartner_patient_oof_scored.csv` | +0.0769 vs prior learned ranker | -0.1154 | [-0.4231, +0.1538] | REJECT |
| 026 | `improve_cv.oof_scored.csv` | +0.2714 vs PRIME | +0.2714 | [+0.0143, +0.5429] | ACCEPT |
| 028 | `results.zip:pred_df_TME_excluded.txt` | +0.2429 vs PRIME | +0.2429 | [+0.0143, +0.4857] | ACCEPT |
| 028 | `improve_xgb_slot_oof.scored.csv` | +0.0143 vs Epicurus | +0.0143 | [-0.1857, +0.2143] | CONSISTENT_WITH_NO_EFFECT |
| 028 | `reconstructed frozen 50/50 RF-Epicurus blend` | +0.0714 vs RF | +0.0714 | [-0.1000, +0.2429] | CONSISTENT_WITH_NO_EFFECT |
| 028 | `improve_none_esm_oof.csv` | +0.2143 vs PRIME | +0.2143 | [-0.0143, +0.4571] | CONSISTENT_WITH_NO_EFFECT |
| 028 | `improve_delta_esm_oof.csv` | -0.0143 vs PRIME | -0.0143 | [-0.2571, +0.2286] | REJECT |
| 028 | `improve_paired_esm_oof.csv` | +0.2143 vs PRIME | +0.2143 | [-0.0286, +0.4571] | CONSISTENT_WITH_NO_EFFECT |
| 028 | `improve_neoguider_official_cv.scored.csv:tme_excluded` | +0.1429 vs PRIME | +0.1429 | [-0.1143, +0.4143] | CONSISTENT_WITH_NO_EFFECT |
| 029 | `improve_neoprecis_approx.scored.csv` | -0.6286 vs PRIME | -0.6286 | [-0.9571, -0.3143] | REJECT |

The surviving `ACCEPT` results are the already-known direct-recognition ranker versus PRIME and the
official IMPROVE RF versus PRIME. Every smaller claimed improvement collapses to
`CONSISTENT_WITH_NO_EFFECT` or `REJECT`.
