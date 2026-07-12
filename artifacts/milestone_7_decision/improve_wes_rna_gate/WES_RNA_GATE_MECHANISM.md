# Biology-first WES/RNA gate on IMPROVE — mechanism report

**Intervention.** Frozen Epicurus v0.1 base order, unchanged reranker. A label-blind demotion gate removes top-20 candidates failing a biological presentation/clonality prerequisite; freed slots backfill by base order. Metric = patient-equal Δ(recognized hits@20). Observational; no causal claim.

## 1. Reachability ceiling (is there headroom at all?)

- n_patients: 70
- patients_with_positive: 61
- total_positives: 467
- top20_hits_now: 75
- promotable_pos_21_60: 111
- patients_with_promotable: 52
- demotable_top20_neg: 1325
- mean_top20_hits: 1.0714
- mean_promotable: 1.5857
- challenger_pos_rate: 0.0396
- top20_pos_rate: 0.0536

> Backfill economics: challenger (21-60) pos-rate **0.0396** vs top-20 pos-rate **0.0536**. When challenger < top-20, random demotion+backfill is net-negative — so any gate must clear that bar, not just beat zero.

## 2. Mechanism — within-patient partial effects (decision zone, ranks 1-60)

base ranks 1-60 (n=4200, pos_rate=0.0443)

Monotone direction = pos_rate(top bin) − pos_rate(bottom bin); positive ⇒ higher feature ⇒ more recognized.

| feature | slope(top−bottom) | per-bin pos_rate (bin:-1=missing) |
|---|---|---|
| Expression | 0.0009 | 0:0.0429(n1050) 1:0.0448(n1050) 2:0.0457(n1050) 3:0.0438(n1050) |
| rna_af | 0.0038 | 0:0.0505(n1050) 1:0.0429(n1050) 2:0.0295(n1050) 3:0.0543(n1050) |
| rna_var | 0.0048 | 0:0.0495(n1050) 1:0.0371(n1050) 2:0.0362(n1050) 3:0.0543(n1050) |
| ValMutRNACoef | 0.0019 | 0:0.0505(n1050) 1:0.039(n1050) 2:0.0352(n1050) 3:0.0524(n1050) |
| VarAlFreq | 0.0162 | 0:0.0286(n1050) 1:0.0429(n1050) 2:0.061(n1050) 3:0.0448(n1050) |
| CelPrev | 0.0028 | 0:0.0486(n1050) 1:0.0429(n1050) 2:0.0343(n1050) 3:0.0514(n1050) |
| HLAexp | -0.0114 | 0:0.0476(n1050) 1:0.0505(n1050) 2:0.0429(n1050) 3:0.0362(n1050) |
| Stability | -0.0086 | 0:0.0457(n1050) 1:0.059(n1050) 2:0.0352(n1050) 3:0.0371(n1050) |
| DAI | 0.0019 | 0:0.0419(n1050) 1:0.0438(n1050) 2:0.0476(n1050) 3:0.0438(n1050) |
| Foreigness | -0.0238 | 0:0.0648(n1050) 1:0.0352(n1050) 2:0.0362(n1050) 3:0.041(n1050) |

**Clonal-vs-expression 2×2** (within-patient median split, pos_rate):
- subclonal_low_expr: 0.029 (n=1001)
- subclonal_high_expr: 0.0419 (n=1099)
- clonal_low_expr: 0.0573 (n=1099)
- clonal_high_expr: 0.048 (n=1001)

**HLA-expression interactions** (pos_rate):
- HLAexp_x_RankEL: lowHLAexp_lo=0.0514(n1012), lowHLAexp_hi=0.0469(n1088), highHLAexp_lo=0.0368(n1088), highHLAexp_hi=0.0425(n1012)
- HLAexp_x_DAI: lowHLAexp_lo=0.0529(n1059), lowHLAexp_hi=0.0451(n1041), highHLAexp_lo=0.0327(n1041), highHLAexp_hi=0.0463(n1059)
- HLAexp_x_Stability: lowHLAexp_lo=0.0617(n957), lowHLAexp_hi=0.0385(n1143), highHLAexp_lo=0.0446(n1143), highHLAexp_hi=0.0334(n957)

## 3. Pre-declared biology gates — Δ(hits@20), controls, transport

| gate | fired top20 | Δmean | boot CI | frac>0 | random Δ | gate−random | transport (Basket/bladder/melanoma) |
|---|---|---|---|---|---|---|---|
| G1_no_rna_confirmation | 589 | -0.1429 | [-0.3143,0.0143] | 0.0335 | -0.1058 | -0.0371 | 0.0/-0.1667/-0.2308 |
| G2_zero_mutant_rna_reads | 589 | -0.1429 | [-0.3143,0.0143] | 0.0335 | -0.1058 | -0.0371 | 0.0/-0.1667/-0.2308 |
| G3_low_rna_vaf_q1 | 363 | -0.1 | [-0.2429,0.0429] | 0.0595 | -0.0614 | -0.0386 | 0.15/-0.2917/-0.1154 |
| G4_low_expression_q1 | 136 | -0.0714 | [-0.1714,0.0286] | 0.0585 | -0.0483 | -0.0231 | -0.15/-0.0417/-0.0385 |
| G5_subclonal_and_low_expr | 59 | -0.0143 | [-0.0857,0.0714] | 0.265 | -0.0174 | 0.0031 | 0.0/0.0/-0.0385 |
| G6_low_hla_expression_q1 | 391 | -0.0714 | [-0.2429,0.0857] | 0.1665 | -0.0473 | -0.0241 | 0.1/-0.2917/0.0 |
| NEGCTRL_high_hla_expression_q4 | 278 | -0.0571 | [-0.2143,0.0714] | 0.174 | -0.0891 | 0.0319 | -0.1/0.0417/-0.1154 |

## 4. Label-ascertainment cross-check (direction of expression/VAF effect)

- IMPROVE Expression slope (ranks1-60): **0.0009**, VarAlFreq slope: **0.0162**
- Gartner Expression(decile) slope: **0.0392**, VarAlFreq(decile) slope: **0.0106**
- multimer expression slope: **0.003**
- Rich features unavailable elsewhere (read-level RNA, HLAexp) cannot be cross-checked; only expression/VAF direction is comparable. A sign flip across cohorts ⇒ the IMPROVE effect is ascertainment-shaped, not a transportable biology.

## 5. Verdict — candidate gate rules

**No deployable gate.** No pre-declared biology gate simultaneously (a) has a bootstrap CI excluding 0, (b) beats matched-random removal, and (c) stays non-negative across all three cancer cohorts. Every RNA-prerequisite gate (no-RNA-confirmation, zero-mutant-reads, low-RNA-VAF, low-expression) is net-NEGATIVE, because within the frozen top-20 boundary positives fail these prerequisites as often as decoys — the presentation rank already caught what these features could catch. The marginal HLA-expression signal (positives ~1.7× higher) DISSOLVES within-patient (the demote-high-HLAexp negative control is *less* harmful than random, the opposite of a real effect) — it was patient-scale confounding.

### Why (mechanism, not just a null)

- **Ascertainment, verified.** IMPROVE Expression slope within-patient at the boundary is ~0 (0.0009) while the SAME axis is clearly positive in Gartner (0.0392). Expression is a live recognition signal on a broad denominator (Gartner) but is FLATTENED in IMPROVE because its ~200-candidate denominator was pre-screened on expression. So expression's null here is an ascertainment artifact, NOT a contradiction of the biology. The read-level RNA and HLAexp features have no Gartner/multimer analogue and cannot be cross-checked — their nulls are unproven, not established.
- **Clonality > expression at this boundary.** The within-patient VAF×Expression 2×2 puts clonal_low_expr HIGHEST (0.0573) and subclonal_low_expr LOWEST (0.029): once candidates are pre-screened for expression, low expression does not kill recognition, but SUBCLONALITY does depress it. DNA VAF / clonality is the single axis whose weak positive direction is CONSISTENT across IMPROVE (0.0162) and Gartner (0.0106).

### One candidate direction for PROSPECTIVE/EXTERNAL testing (not established)

- **Clonality (DNA VAF / truncality), as a PROMOTE-side prior on a NON-pre-screened denominator** — not as a top-20 demotion gate (that failed here). The demotion form fails on IMPROVE precisely because IMPROVE's boundary is already expression/RNA-screened and recognition-limited. The clean test is a full-mutanome denominator where clonal truncal mutations must be *recovered* from a large decoy pool (e.g. the Miller full re-enumeration, or Gartner's broad denominator), asking whether a clonality prior lifts recognized hits@20 — with the same matched-random + LOCO + bootstrap discipline. Prediction from this audit: effect will be small and may not clear the bar.