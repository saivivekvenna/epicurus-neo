> **CORRECTION (supersedes the §2 IMPROVE row and §6 bottom line below).** The IMPROVE entry here audits the *reduced* export (prime/el/expr) and shows no orthogonal features. That is a scope artifact, not a fact about IMPROVE: the **raw 88-column** table (`data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt`) yields **36 deployable candidate-varying orthogonal features** across 11 families (agretopicity, clonality, expression, foreignness_selfsim, hla_expression, mutation_annotation, nn_align, physicochemical, rna_support, stability_processing, vaf_readsupport). See `IMPROVE_RAW_COLUMN_AUDIT.md` / `IMPROVE_DEPLOYABLE_WHITELIST.json`.

# Gate feature audit — orthogonal levers against high-PRIME false positives

**Why this exists.** The dynamic gate (`configs/frozen/dynamic_gate_v1.json`) was falsified: a label-blind *presentation* gate removes 0 pct of the high-presentation decoys that outrank positives, so downstream Δhits@20 = 0. A gate can only help via an **orthogonal** feature that separates the top-ranked TESTED_NEGATIVE decoys from POSITIVES *within the high-presentation stratum* (top-20 by presentation per patient). This audit measures exactly that, read-only, on all seven cohorts.

Presentation baseline on the stratum should sit near 0.5 by construction — that is the wall. Any orthogonal feature with |AUROC−0.5| meaningfully above 0 on the stratum is a candidate unlock. Study/patient identity is audited as a confound and is **never** a deployable feature.

## 1. Ranked feature-unlock matrix

| rank | family | best stratum AUROC | \|signal\| | best cohort | max coverage | availability | leakage risk | next experiment |
|---|---|---|---|---|---|---|---|---|
| 1 | expression | 0.8187 | 0.3187 | gartner | 1.0 | available_now (Gartner deciles); Miller needs RNA-seq; Sid sparse | medium (ascertainment/PU) | conditional expr AUROC on Miller after RNA quant; guard vs PU per m7-ascertainment-correction |
| 2 | predictor_disagreement | 0.7039 | 0.2039 | gartner | 1.0 | available_now (Gartner 5, Zhao 9); absent multimer/IMPROVE | low | test std-of-ranks on stratum with a pre-registered direction |
| 3 | vaf_readsupport | 0.5999 | 0.0999 | gartner | 1.0 | available_now (Gartner vaf_decile); Miller/Sid need WES | medium | join VAF/read-depth on Miller after somatic calling |
| 4 | physicochemical | 0.56 | 0.06 | zhao | 1.0 | available_now (Zhao); trivially computable everywhere | low | recompute gravy/charge for all cohorts (cheap) |
| 5 | stability_processing | 0.5369 | 0.0369 | zhao | 1.0 | available_now (Zhao NetMHCstab/NetCTLpan); add to Gartner/Miller | low | add NetMHCstab/NetCTLpan pass to Gartner+Miller universe |
| 6 | agretopicity | 0.5106 | 0.0106 | zhao | 1.0 | available_now (Zhao mut_kd_delta); Gartner/Miller need 1 WT predictor pass | low | score WT vs mutant EL for Gartner+Miller; test mutant/WT ratio on stratum |
| 7 | mutation_annotation | 0.4953 | 0.0047 | zhao | 1.0 | available_now (Gartner CGC/Cosmic/dbSNP) | medium (annotation-source leakage) | map driver/germline flags on Miller from public annotation |
| 8 | repeated_antigen | None | None | None | 0.0 | available_now (CEDAR n_subjects; cross-cohort recurrence) | HIGH (studied-because-immunogenic) | restrict to prospective recurrence only; never use assay-count of the same cohort |
| 9 | assay_context | None | None | None | 0.0 | available_now but STUDY-CONFOUNDED | audit-only (never deployable) | use as stratification/confound check, not a feature |
| 10 | llm_artifact_plausibility | None | None | None | 0.0 | computable_now (label-blind, annotation-only) | low (no labels sent) | score full universe; test whether artifact_risk_score down-ranks TESTED_NEGATIVE decoys on the stratum |

## 2. Per-cohort conditional signal (POSITIVE vs TESTED_NEGATIVE)

### gartner  (n=8777, patients=26)
- anchor: `prime` (lower=better) — genuine PRIME 2.1
- labels: {'UNTESTED': 5009, 'TESTED_NEGATIVE': 3722, 'POSITIVE': 46}
- **presentation baseline on stratum**: AUROC=0.5011 (pos=17, neg=418); marginal=0.7289  ← near 0.5 on stratum = the wall

| feature | family | cov | marginal AUROC | stratum AUROC | \|stratum signal\| | stratum n(pos/neg) |
|---|---|---|---|---|---|---|
| expr_decile | expression | 1.0 | 0.817 | 0.8187 | 0.3187 | 17/418 |
| seen_in_rna | expression | 1.0 | 0.751 | 0.772 | 0.272 | 17/418 |
| pred_rank_std | predictor_disagreement | 1.0 | 0.7196 | 0.7039 | 0.2039 | 17/418 |
| vaf_decile | vaf_readsupport | 1.0 | 0.6029 | 0.5999 | 0.0999 | 17/418 |
| is_cosmic | mutation_annotation | 1.0 | 0.5 | 0.5 | 0.0 | 17/418 |
| is_dbsnp | mutation_annotation | 1.0 | 0.5 | 0.5 | 0.0 | 17/418 |

- cross-fitted (patient-grouped OOF) on stratum: orthogonal-only=0.8551, presentation-only=0.344, combined=0.8377 (orthogonal-only full-cohort=0.8451)
- confound audit (AUDIT-ONLY, never deployable): best orthogonal `expr_decile` marginal=0.817 vs patient-generalising OOF=0.8069 (identity-parasitic gap=0.0101); tissue-label-only OOF=0.4924

### zhao  (n=2315, patients=352)
- anchor: `mixmhcpred3_score` (higher=better) — mixMHCpred 3.0 (not genuine PRIME)
- labels: {'TESTED_NEGATIVE': 2002, 'POSITIVE': 313}
- **presentation baseline on stratum**: AUROC=0.5447 (pos=313, neg=2002); marginal=0.5447  ← near 0.5 on stratum = the wall

| feature | family | cov | marginal AUROC | stratum AUROC | \|stratum signal\| | stratum n(pos/neg) |
|---|---|---|---|---|---|---|
| gravy | physicochemical | 1.0 | 0.56 | 0.56 | 0.06 | 313/2002 |
| aromatic_frac | physicochemical | 1.0 | 0.5506 | 0.5506 | 0.0506 | 313/2002 |
| pred_NetMHCstab | stability_processing | 1.0 | 0.5369 | 0.5369 | 0.0369 | 313/2002 |
| net_charge | physicochemical | 1.0 | 0.469 | 0.469 | 0.031 | 313/2002 |
| pred_NetCTLpan Cleavage | stability_processing | 1.0 | 0.5224 | 0.5224 | 0.0224 | 313/2002 |
| pred_MHCflurryProc | stability_processing | 1.0 | 0.5159 | 0.5159 | 0.0159 | 313/2002 |
| pred_z_std | predictor_disagreement | 1.0 | 0.5134 | 0.5134 | 0.0134 | 313/2002 |
| pred_NetCTLpan TAP | stability_processing | 1.0 | 0.512 | 0.512 | 0.012 | 313/2002 |
| mut_kd_delta | agretopicity | 1.0 | 0.5106 | 0.5106 | 0.0106 | 313/2002 |
| mut_charge_delta | physicochemical | 1.0 | 0.5078 | 0.5078 | 0.0078 | 313/2002 |
| mut_position_frac | mutation_annotation | 1.0 | 0.4953 | 0.4953 | 0.0047 | 313/2002 |
| mut_is_anchor | mutation_annotation | 1.0 | 0.4987 | 0.4987 | 0.0013 | 313/2002 |

- cross-fitted (patient-grouped OOF) on stratum: orthogonal-only=0.5402, presentation-only=0.5363, combined=0.5503 (orthogonal-only full-cohort=0.5402)
- confound audit (AUDIT-ONLY, never deployable): best orthogonal `gravy` marginal=0.56 vs patient-generalising OOF=0.549 (identity-parasitic gap=0.011); tissue-label-only OOF=None

### multimer  (n=8103, patients=26)
- anchor: `prime` (lower=better) — genuine PRIME 2.1 (rank-like: lower=better)
- labels: {'TESTED_NEGATIVE': 8069, 'POSITIVE': 34}
- **presentation baseline on stratum**: AUROC=0.6319 (pos=21, neg=494); marginal=0.8457  ← near 0.5 on stratum = the wall

- cross-fitted (patient-grouped OOF) on stratum: orthogonal-only=None, presentation-only=0.5379, combined=0.5379 (orthogonal-only full-cohort=None)

### improve  (n=17520, patients=70)
- anchor: `prime` (lower=better) — genuine PRIME 2.1 (rank-like: lower=better)
- labels: {'TESTED_NEGATIVE': 17053, 'POSITIVE': 467}
- **presentation baseline on stratum**: AUROC=0.5139 (pos=83, neg=1317); marginal=0.5916  ← near 0.5 on stratum = the wall

- cross-fitted (patient-grouped OOF) on stratum: orthogonal-only=None, presentation-only=0.4531, combined=0.4531 (orthogonal-only full-cohort=None)

### osteosarc  (n=7565, patients=1)
- anchor: `genuine_prime` (higher=better) — genuine PRIME 2.1

## 3. CEDAR (recognition prior — no gate stratum)

- recognition-prior asset: NO within-patient presentation anchor -> no gate stratum; repeated-antigen signal is reported but is HEAVILY leakage-prone (studied because immunogenic).
- labels: {'TESTED_NEGATIVE': 24974, 'POSITIVE': 7491}
- repeated-antigen `n_assays_for_peptide`: AUROC=0.7522 (HIGH leakage: peptides are re-assayed *because* immunogenic)

## 4. Miller (LOCKED_TEST — labels in hand, inputs not yet computed)

- LOCKED_TEST; labels ingested but NO presentation anchor and NO expression/VAF yet (requires the public WES/RNA download + HLA typing + RNA quant). No gate stratum computable now.
- labels: {'TESTED_NEGATIVE': 574, 'POSITIVE': 180}
- available now: {"agretopicity": {"raw": "ref_peptide present (coverage 0.995)", "status": "computable after ONE predictor pass scoring ref vs mutant"}, "mutation_annotation": {"gene/transcript/HGVS present": true}}
- requires Miller WES/RNA: ['expression', 'vaf_readsupport', 'predictor_disagreement', 'stability_processing', 'presentation_anchor']

## 5. LLM structured feature — feasibility (label-blind, artifact/transcript plausibility ONLY)

- status: **RAN (blind feasibility sample; annotation-only, no labels/identifiers)**
- sends only: ['gene_symbol', 'hgvs_c', 'hgvs_p', 'mutant_peptide', 'transcript_id', 'variant_type', 'wt_peptide'] (contains_labels=False, contains_patient_or_study_id=False)
- purpose: down-weight *annotation artifacts* (pseudogene / retained-intron / frame errors / NMD), an axis orthogonal to both presentation and recognition. **Not** an 'is this immunogenic?' guess.
- schema + prompt cached in `llm_feasibility_cache.json`.

| gene | variant | artifact_risk | coding_plausible | nmd_risk |
|---|---|---|---|---|
| ATRX | SNV | 0 | True | False |
| ATRX | SNV | 0 | True | False |
| CCAR2 | SNV | 0 | True | False |

## 6. Bottom line

- **Gartner and Zhao** carry both a presentation anchor and orthogonal features in their *loaded* frames, so they test the gate question directly here. **IMPROVE is NOT featureless** — the reduced export loaded in §2 kept only prime/el/expr, but its raw 88-column table is orthogonally rich (see the CORRECTION banner and `IMPROVE_RAW_COLUMN_AUDIT.md`); it simply needs the raw table wired into the gate frame. multimer's loaded frame is genuinely presentation-only. CEDAR has no anchor; Miller has no inputs yet; Sid has 3 positives and no clean negative denominator.
- The unlock matrix above ranks which family to invest in. Read the stratum AUROC, not the marginal: a feature strong marginally but ~0.5 on the stratum cannot remove high-PRIME decoys (that was expression's risk per `m7-ascertainment-correction`).
- The single highest-leverage, lowest-leakage NEW axis is **LLM artifact/transcript plausibility**, because it is label-blind, computable now on existing annotations, and orthogonal to the presentation wall.

### Tension & caveats (do not over-read the expression result)

- **Expression is the strongest orthogonal signal on the exact stratum where the gate failed** (Gartner stratum AUROC 0.82 while presentation is 0.50), it cross-fits across patients (OOF 0.855 vs presentation 0.34), and it is **not** identity-parasitic (marginal−OOF gap ≈ 0.01; tissue-label-only OOF ≈ 0.49). That is a real, generalising, orthogonal axis — exactly what a presentation gate lacks.
- **But this does NOT license an expression rank penalty.** The frozen expression policy (`configs/frozen/expression_policy_v1.json`, memory `expression-ranking-policy`) already showed that turning expression into a rank penalty REGRESSES conditional-ranking cohorts within a fixed top-20 budget. High discrimination on the stratum ≠ a safe deployable penalty; the unlock is a *conditional / confidence* use that spares strong positives, not a monotone demotion.
- **Underpowered:** only 17 POSITIVES survive onto the Gartner top-20 stratum (46 total). Per `m7-ascertainment-correction`, Gartner expression/VAF is 'promising but UNDERPOWERED'; this audit is CONSISTENT_WITH a real effect, it does not establish one. Negatives here are measured TESTED_NEGATIVE (not UNTESTED), which blunts but does not remove the ascertainment/PU risk.
- **Predictor disagreement** is the second orthogonal axis (Gartner stratum |signal| 0.20) and is cheap and low-leakage — but its direction must be pre-registered, not chosen from this table.
- **Zhao shows no orthogonal lever** (best stratum |signal| 0.06); its weak anchor (mixMHCpred) means its 'stratum' is barely selective. Only multimer's *loaded* frame is truly presentation-only; IMPROVE's apparent featurelessness was a loader-scope artifact, now corrected.