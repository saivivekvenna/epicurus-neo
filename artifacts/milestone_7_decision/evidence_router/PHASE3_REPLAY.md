# Evidence router — Phase 3 locked replay + conditional no-regression

> Policy `epicurus-evidence-router-1.0.0`. **No fitting/tuning.** Frozen router + route-aware selection replayed over existing local artifacts. Mechanical reachability and learned no-regression are reported separately; neither is used to argue the other.

## 1. osteosarc.com / Sid — locked reachability replay

> Sid informed the router design and is **not** independent validation. No policy constant was tuned on it. Vaccine-selected peptides are not used to fill candidate generation; missing peptide/HLA is never charged to PRIME.

- Hudson IFNgamma/TCR recognized targets: **3** (ASPM-chr1-197102716, DYNC1H1-chr14-101980529, MAP2-chr2-209694772).
- **Multi-caller raw variant union recall: 3/3 = 1.00** (all three are called by DRAGEN/Sarek/oncoanalyser — the recoverable loss is candidate recall upstream of any ranker).

| Stage | Count |
|---|---:|
| generated (raw union) | 3 |
| peptide-generated (pVACtools 2025.01) | 1 |
| rankable (peptide+HLA) | 1 |
| selected (frozen Epicurus v0.1 top-20) | 1 |
| NEEDS_PEPTIDE_GENERATION | 2 |

- Peptide-generation gap (NEEDS_PEPTIDE_GENERATION, not a ranker miss): **ASPM-chr1-197102716, MAP2-chr2-209694772**.
- Reached ranking / shortlist: **DYNC1H1-chr14-101980529**.
- hits@20 conditional on rankability: **1/1**; end-to-end hits@20: **1/3**.

Router routing of the three peptide-free variants (offline): ASPM-chr1-197102716 -> LONGITUDINAL/NEEDS_PEPTIDE_GENERATION; MAP2-chr2-209694772 -> RESCUE/NEEDS_PEPTIDE_GENERATION; DYNC1H1-chr14-101980529 -> LONGITUDINAL/NEEDS_PEPTIDE_GENERATION. None is hard-removed (IMPOSSIBLE); none is charged to the ranker.

## 2. Reranker cohorts — conditional no-regression

> Requirement (§7): route-aware top-20 must not LOSE hits@20 vs the pure-score top-20. The router is recall-preserving.

| Cohort | Score | Patients | Positives | Rankable | HLA populated | Route-aware hits@20 | Pure-score hits@20 | Δ (CI) | Verdict |
|---|---|---:|---:|---:|---:|---:|---:|---|---|
| cd8_multimer | prime | 26 | 34 | 8095 | 1.00 | 21 | 21 | 0.000 [0.000, 0.000] | NO_REGRESSION_EXACT_PASS_THROUGH |
| gartner_nci | full_model | 26 | 46 | 0 | 0.00 | n/a | n/a | n/a | NOT_EVALUABLE_NO_RANKABLE_CANDIDATES |
| cd8_multimer_decision | full_model | 26 | 34 | 8095 | 1.00 | 13 | 13 | 0.000 [0.000, 0.000] | NO_REGRESSION_EXACT_PASS_THROUGH |

### Honest feature-availability caveats

- On the evaluable reranker artifacts (both CD8 multimer scorings) the router's recall-discriminating features (RNA/expression/HLA-LOH/multi-caller provenance) are **absent** and no diversity-cap key is populated, so all candidates collapse to a single route and route-aware selection is a **pure-score pass-through** (Δ=0 by construction). This satisfies no-regression but is **not** evidence of benefit; the router's value is realised at candidate **generation** (the Sid recall recovery), not at reranking a feature-poor list.
- The **Gartner NCI** stored artifact carries **no HLA allele** (hla populated = 0.00), so nothing is rankable and the top-20 no-regression is **NOT_EVALUABLE** there — a data limit of the file, not a ranker result.
- Router removals per cohort (route-verifiable only; positives lost must be 0): cd8_multimer: {'DUP_CANDIDATE': 8} removed, 0 positives lost; gartner_nci: {} removed, 0 positives lost; cd8_multimer_decision: {'DUP_CANDIDATE': 8} removed, 0 positives lost.

Route composition (valid candidates) per cohort: cd8_multimer={'UNCERTAIN': 8095}; gartner_nci={'UNCERTAIN': 3768}; cd8_multimer_decision={'UNCERTAIN': 8095}.

