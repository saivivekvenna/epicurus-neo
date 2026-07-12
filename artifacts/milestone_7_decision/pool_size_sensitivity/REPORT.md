# Pool-size sensitivity diagnostic — frozen Epicurus vs genuine PRIME

`python -m scripts.pool_size_sensitivity` · 25 seeds · frozen scorers only (no model fit, no PRIME binary re-run).

**Question.** Does Epicurus vs genuine PRIME improve when the pool shrinks but retains all positives?

Three arms, all frozen: **genuine PRIME** (`-%rank`, pool-invariant), **frozen Epicurus** (`configs/frozen/epicurus_v0_1.json`; within-patient percentiles → *pool-dependent*), and a tuning-free **pVAC-style+PRIME** proxy (equal-weight EL⊕PRIME percentiles — *not* the pVACtools binary).

Pools per patient: **LARGE** = all tested candidates; **MEDIUM** = all positives + 50% of negatives; **SMALL** = all positives + 25% of negatives. SMALL ⊂ MEDIUM ⊂ LARGE, identical positives/features.


## gartner

Gartner NCI Testing (frozen external transfer test); clean for Epicurus.  
Eligible patients **26** · positives **46** · tested-negatives **3722** · excluded 0 (none).


**Variant A (oracle — reranker stress test; positives 100% retained by construction).** Mean over patients×seeds.

| pool | arm | hits@5 | hits@10 | hits@20 | recall@20 | prec@20 | MRR | med.pos.rank | sat≤20 |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| LARGE | genuine_prime | 0.308 | 0.423 | 0.692 | 0.404 | 0.035 | 0.163 | 35.7 | 0.00 |
| LARGE | frozen_epicurus | 0.538 | 0.615 | 0.808 | 0.474 | 0.040 | 0.269 | 26.9 | 0.00 |
| LARGE | pvac_style_prime | 0.462 | 0.654 | 0.731 | 0.423 | 0.037 | 0.187 | 27.8 | 0.00 |
| MEDIUM | genuine_prime | 0.403 | 0.703 | 1.057 | 0.608 | 0.053 | 0.236 | 18.6 | 0.00 |
| MEDIUM | frozen_epicurus | 0.623 | 0.835 | 1.232 | 0.715 | 0.062 | 0.379 | 14.3 | 0.00 |
| MEDIUM | pvac_style_prime | 0.580 | 0.771 | 1.202 | 0.724 | 0.060 | 0.292 | 14.6 | 0.00 |
| SMALL | genuine_prime | 0.640 | 1.042 | 1.529 | 0.890 | 0.078 | 0.339 | 10.0 | 0.08 |
| SMALL | frozen_epicurus | 0.778 | 1.192 | 1.652 | 0.953 | 0.084 | 0.475 | 8.0 | 0.08 |
| SMALL | pvac_style_prime | 0.769 | 1.198 | 1.652 | 0.952 | 0.084 | 0.422 | 8.0 | 0.08 |

**Epicurus − PRIME hits@20 (paired, mean [95% band over seeds]):**  
LARGE 0.115 [0.1154, 0.1154]  ·  
MEDIUM 0.175 [0.0769, 0.2846]  ·  
SMALL 0.123 [0.0385, 0.2077]  ·  


**Variant B (label-blind EL-percentile gate; size-matched; positives may drop).**

| pool | pos. retention (mean / min) | all retained? | # patients losing a positive | epi−PRIME hits@20 |
|---|--:|:--:|--:|--:|
| LARGE | 1.000 / 1.000 | yes | 0 | 0.115 |
| MEDIUM | 0.917 / 0.333 | NO | 4 | 0.038 |
| SMALL | 0.667 / 0.000 | NO | 12 | 0.115 |

## improve

IMPROVE SRHgroup (frozen external validation); clean for Epicurus; best-powered.  
Eligible patients **61** · positives **467** · tested-negatives **15437** · excluded 9 (no_positive).


**Variant A (oracle — reranker stress test; positives 100% retained by construction).** Mean over patients×seeds.

| pool | arm | hits@5 | hits@10 | hits@20 | recall@20 | prec@20 | MRR | med.pos.rank | sat≤20 |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| LARGE | genuine_prime | 0.475 | 0.836 | 1.361 | 0.237 | 0.068 | 0.252 | 95.3 | 0.00 |
| LARGE | frozen_epicurus | 0.426 | 0.770 | 1.230 | 0.200 | 0.061 | 0.203 | 95.7 | 0.00 |
| LARGE | pvac_style_prime | 0.410 | 0.770 | 1.164 | 0.200 | 0.058 | 0.202 | 97.6 | 0.00 |
| MEDIUM | genuine_prime | 0.772 | 1.293 | 2.082 | 0.344 | 0.104 | 0.347 | 49.9 | 0.00 |
| MEDIUM | frozen_epicurus | 0.694 | 1.157 | 1.977 | 0.344 | 0.099 | 0.299 | 50.1 | 0.00 |
| MEDIUM | pvac_style_prime | 0.656 | 1.129 | 1.959 | 0.333 | 0.098 | 0.286 | 51.2 | 0.00 |
| SMALL | genuine_prime | 1.143 | 1.874 | 3.241 | 0.533 | 0.162 | 0.455 | 27.2 | 0.02 |
| SMALL | frozen_epicurus | 1.031 | 1.791 | 3.214 | 0.537 | 0.161 | 0.415 | 27.3 | 0.02 |
| SMALL | pvac_style_prime | 1.009 | 1.778 | 2.951 | 0.493 | 0.148 | 0.397 | 27.9 | 0.02 |

**Epicurus − PRIME hits@20 (paired, mean [95% band over seeds]):**  
LARGE -0.131 [-0.1311, -0.1311]  ·  
MEDIUM -0.105 [-0.2459, -0.0295]  ·  
SMALL -0.026 [-0.1607, 0.0754]  ·  


**Variant B (label-blind EL-percentile gate; size-matched; positives may drop).**

| pool | pos. retention (mean / min) | all retained? | # patients losing a positive | epi−PRIME hits@20 |
|---|--:|:--:|--:|--:|
| LARGE | 1.000 / 1.000 | yes | 0 | -0.131 |
| MEDIUM | 0.598 / 0.000 | NO | 51 | -0.082 |
| SMALL | 0.357 / 0.000 | NO | 56 | -0.016 |

## multimer ⚠️ *Epicurus IN-SAMPLE (training cohort)*

CD8 pMHC-multimer; frozen Epicurus TRAINING cohort -> Epicurus IN-SAMPLE (flagged).  
Eligible patients **19** · positives **34** · tested-negatives **7091** · excluded 7 (no_positive).


**Variant A (oracle — reranker stress test; positives 100% retained by construction).** Mean over patients×seeds.

| pool | arm | hits@5 | hits@10 | hits@20 | recall@20 | prec@20 | MRR | med.pos.rank | sat≤20 |
|---|---|--:|--:|--:|--:|--:|--:|--:|--:|
| LARGE | genuine_prime | 0.632 | 0.895 | 1.105 | 0.730 | 0.056 | 0.319 | 65.8 | 0.05 |
| LARGE | frozen_epicurus | 0.737 | 1.211 | 1.263 | 0.733 | 0.064 | 0.337 | 52.9 | 0.05 |
| LARGE | pvac_style_prime | 0.684 | 0.842 | 1.053 | 0.635 | 0.053 | 0.348 | 62.7 | 0.05 |
| MEDIUM | genuine_prime | 0.912 | 1.143 | 1.299 | 0.782 | 0.068 | 0.445 | 33.7 | 0.05 |
| MEDIUM | frozen_epicurus | 1.078 | 1.265 | 1.347 | 0.754 | 0.071 | 0.459 | 27.6 | 0.05 |
| MEDIUM | pvac_style_prime | 0.821 | 1.017 | 1.312 | 0.742 | 0.069 | 0.484 | 32.1 | 0.05 |
| SMALL | genuine_prime | 1.095 | 1.257 | 1.368 | 0.807 | 0.077 | 0.597 | 17.6 | 0.16 |
| SMALL | frozen_epicurus | 1.225 | 1.320 | 1.457 | 0.850 | 0.081 | 0.608 | 14.6 | 0.16 |
| SMALL | pvac_style_prime | 0.989 | 1.282 | 1.404 | 0.797 | 0.078 | 0.628 | 16.7 | 0.16 |

**Epicurus − PRIME hits@20 (paired, mean [95% band over seeds]):**  
LARGE 0.158 [0.1579, 0.1579]  ·  
MEDIUM 0.048 [0.0, 0.1263]  ·  
SMALL 0.088 [0.0526, 0.1579]  ·  


**Variant B (label-blind EL-percentile gate; size-matched; positives may drop).**

| pool | pos. retention (mean / min) | all retained? | # patients losing a positive | epi−PRIME hits@20 |
|---|--:|:--:|--:|--:|
| LARGE | 1.000 / 1.000 | yes | 0 | 0.158 |
| MEDIUM | 0.814 / 0.000 | NO | 6 | 0.158 |
| SMALL | 0.718 / 0.000 | NO | 8 | 0.105 |

## Answers

**Q1 — Do all models improve just because there are fewer negatives?**  
recall@20 and precision@20 rise for *every* arm as the pool shrinks (mechanical: fewer competitors / smaller top-20 denominator). MRR and median positive rank are the denominator-robust checks — see the per-cohort tables. If MRR is roughly flat while recall climbs, the gain is the denominator, not better ranking.

**Q2 — Does Epicurus's advantage over PRIME increase as the pool shrinks?**  

- gartner: LARGE 0.115 → MEDIUM 0.175 → SMALL 0.123 (paired hits@20 delta)

- improve: LARGE -0.131 → MEDIUM -0.105 → SMALL -0.026 (paired hits@20 delta)

- multimer: LARGE 0.158 → MEDIUM 0.048 → SMALL 0.088 (paired hits@20 delta) — IN-SAMPLE

**Q3 — Does a real label-blind filter retain the positives, or only the oracle?**  

- gartner: MEDIUM gate retains 0.917 of positives (all retained: False); SMALL gate retains 0.667 (all retained: False, 12 patients lose ≥1). The oracle keeps 100% by construction; the deployable gate does not.

- improve: MEDIUM gate retains 0.598 of positives (all retained: False); SMALL gate retains 0.357 (all retained: False, 56 patients lose ≥1). The oracle keeps 100% by construction; the deployable gate does not.

- multimer: MEDIUM gate retains 0.814 of positives (all retained: False); SMALL gate retains 0.718 (all retained: False, 8 patients lose ≥1). The oracle keeps 100% by construction; the deployable gate does not.


> Oracle (variant A) retention is diagnostic, NOT validation: a real pipeline cannot pre-know positives.
> Multimer frozen-Epicurus scores are IN-SAMPLE (training cohort) -> optimistic; do not read as external.
> Cohorts are heterogeneous with fixed roles and are NEVER pooled into one headline number.
> Precision@20 and recall@20 rise mechanically as negatives thin; MRR / median positive rank are the shift-invariant checks for whether ranking QUALITY (not just the denominator) changed.
> No superiority is claimed; verdicts follow the project's challenger-vs-baseline convention.
