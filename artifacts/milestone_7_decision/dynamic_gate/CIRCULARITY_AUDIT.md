# Circularity / feature-overlap audit — dynamic gate v1

The v1 gate's **core veto axes are within-patient percentiles of EL and PRIME**. The downstream rankers are **genuine PRIME** (`-PRIME`) and **frozen Epicurus** (logistic on prime/el/expr percentiles; PRIME coef 1.04 dominates el 0.235 / expr 0.25). The gate's veto inputs are a SUBSET of — and the dominant terms in — the ranker inputs.


## Feature overlap: gate veto axes vs downstream ranker inputs

| downstream ranker | inputs | overlaps gate veto {EL, PRIME}? |
|---|---|---|
| genuine PRIME | PRIME %rank | **PRIME (100% overlap)** |
| frozen Epicurus v0.1 | PRIME, EL, expr percentiles | **EL + PRIME (2 of 3, both dominant)** |

## Empirical consequence (frozen t=0.25)

| cohort | removed | removed that were in ANY ranker top-20 | hard-decoy stratum n | stratum removed |
|---|--:|--:|--:|--:|
| gartner | 166 | **0** | 528 | **0.0** |
| improve | 362 | **0** | 2048 | **0.0** |
| multimer | 222 | **1** | 1030 | **0.0** |

**Spearman(gate keep-margin, ranker score)** — near +1 confirms the gate orders candidates the same way the rankers do:

- gartner: genuine_prime +0.84, frozen_epicurus +0.90
- improve: genuine_prime +0.71, frozen_epicurus +0.84
- multimer: genuine_prime +0.82, frozen_epicurus +0.90

## Conclusion
Essentially zero gate-removed candidates ever sat in a ranker's top-20 (0 / 0 / 1 across gartner / improve / multimer; the single multimer case is a saturated ≤20-candidate pool where every row is trivially 'top-20'), and the gate removes **0%** of the hard-decoy stratum (high-EL AND high-PRIME tested-negatives). This is a **structural tautology of feature overlap**, not evidence about label-blind gating in general: a same-feature monotone presentation gate CANNOT move a top-20 produced by those same features. The falsification is therefore scoped to **presentation/PRIME-derived gating with current peptide features** — the orthogonal-feature dynamic-gate hypothesis is untested (data-blocked), NOT falsified.


## Orthogonal probe (multimer, IN-SAMPLE, exploratory)
Within the hard-decoy stratum (EL%>0.75 AND PRIME%>0.75: 24 positives / 1181 tested-negatives), a cross-fitted leave-one-patient-out model on orthogonal features (Agretopicity, Foreignness score, Proteasomal processing score, Dissimilarity, RNA expression (TPM)) scores **AUROC 0.5919**; a purely-orthogonal veto would remove **0.0923** of stratum negatives at >=95% positive retention. A FAINT orthogonal signal in-sample (AUROC 0.59; ~9% of high-presentation decoys removable at >=95% stratum retention) — but this is multimer's own training cohort (IN-SAMPLE), tiny (24 stratum positives), and leave-one-patient-out is NOT leave-one-study-out. An INDEPENDENT sequence-only hard-decoy experiment (train-Gartner->test-IMPROVE) collapsed OOD (retained 1.5% of positives), so cross-study/assay shift is the real killer. => promising mechanism, NOT validated; requires real WES/RNA features + leave-one-study-out + OOD abstention.

> multimer is frozen Epicurus' training cohort (IN-SAMPLE) and tiny (34 positives) — this is a mechanism sanity check, never a headline. It is the ONLY orthogonal signal available now; a real test requires the WES/RNA features absent from every current eval cohort.
