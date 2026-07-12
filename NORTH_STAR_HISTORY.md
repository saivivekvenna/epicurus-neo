# NORTH_STAR_HISTORY — canonical project ledger

Single source of truth for *where the project stands*: the north star, what has been tried, what was
falsified, what is frozen/preserved, which datasets are consumed, and what the next gate is. This file **links
to** the artifacts that hold the actual results; it does not restate or re-derive them. Keep it updated as the
last step of every milestone.

_Last updated: 2026-07-12 (Dynamic upstream gate — layered label-blind safe-rejection gate built, frozen & LOCKED-tested on CheckMate 153. It Pareto-dominates the incumbent EL gate on retention and is a real recall-preserving candidate-universe pruner, BUT the aggressive premise is falsified: at any safe retention the downstream Δhits@20 is 0 because a label-blind presentation gate removes 0% of the high-presentation decoys that outrank positives — the oracle gap is unclosable this way = the recognition wall re-derived. Unlock = WES/RNA rescue axes. Verdict B/C. Prior entry: RTTP SR24-58221 — FIRST end-to-end Epicurus deployment on a real clinical patient: WES+RNA+HLA+LOH → ranked neoantigen list. Deployment DEMO, not a benchmark — RTTP carries no experimental recognition label (its immunogenicity score is predicted). Epicurus raw-scores like PRIME (ρ 0.966) but its per-mutation shortlist diverges (2/20 shared w/ PRIME); no method agrees w/ the vendor's predicted immunogenicity → labels needed to adjudicate. Prior same day: CheckMate 153 first external PRIME-untouched test, frozen v0.1 TIES PRIME (Δ +0.071) = sixth straight parity; presentation is the ceiling.)_

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
| — | **CheckMate 153** — first external PRIME-untouched test | Frozen Epicurus v0.1 **TIES** genuine PRIME (Δ **+0.071**, CI[−0.357,+0.571]); 6th straight parity, now on EXTERNAL data. Among predicted binders recognition AUROC collapses to ~0.53–0.60 for ALL incl genuine PRIME → **presentation is the ceiling**. | TIE (14 pts, underpowered) | `artifacts/milestone_7_decision/checkmate153/`, memory `checkmate153-external-outcome` |
| — | **RTTP SR24-58221** — first real-patient tool deployment | First **end-to-end Epicurus run on a genuine clinical case** (WES+RNA+HLA+LOH → ranked neoantigen list). **DEMO, not a benchmark** (no experimental recognition label). Epicurus raw-scores like PRIME (per-candidate ρ **0.966**) but its per-mutation shortlist **diverges** (2/20 shared w/ PRIME, 7/20 w/ SHERPA); no method agrees w/ the vendor's *predicted* immunogenicity (~0.10). Handles HLA LOH. | deployment demo | `artifacts/…/rttp_sr24/` (gitignored, DUA), manifest `configs/source_manifests/rttp.yml`, memory `rttp-multiomics-asset` |
| M7 | **Evidence router + route-aware selection** (Phase 2 impl + Phase 3 locked replay) | Additive, non-destructive router (frozen `evidence_router_v1.json`; legacy gate byte-unchanged). Hard-removes ONLY route-verifiable impossibilities; keeps expression-N / TPM-zero / zero-read / atypical-class / single-caller candidates **eligible-but-flagged**; reports empty peptide / missing HLA as **`NEEDS_PEPTIDE_GENERATION`** (never a ranker miss). **Sid locked replay (informed the design → NOT independent validation): multi-caller RAW variant union recall of the 3 Hudson targets = 3/3 = 1.00** (built non-circularly from the real VAF table, 1,213 detected tumor rows → 184 variants), but only **1/3 peptide-generated → rankable → selected** (DYNC1H1; ASPM + MAP2 Gly868fs are NEEDS_PEPTIDE_GENERATION, absent from the single pVACtools 2025.01 set; vaccine peptides NOT used to fill, gap NOT charged to PRIME). **No-regression:** route-aware top-20 == pure-score top-20 **exactly (Δ=0)** on both CD8-multimer scorings; **Gartner NOT_EVALUABLE** (stored artifact has no HLA → nothing rankable); 0 positives removed by the router. Router value is at candidate **generation** (Sid recall), NOT reranking a feature-poor list. | recall-preserving; no superiority claimed | `artifacts/milestone_7_decision/evidence_router/{PHASE3_REPLAY.md,phase3_replay.json}`, `src/epicurus_neo/{evidence_router,variant_union}.py`, `configs/frozen/evidence_router_v1.json` |
| M7 | **MAP2 identity/ORF correction + peptide-recovery protocol** | Supersedes the earlier “different neo-frame” interpretation. The Leu867fs and Gly868fs deletions remain distinct genomic variants, but primary Ensembl VEP/CDS reconstruction shows a shared downstream mutant ORF containing `RVVPFTKAL`; `GYCVFNKYTV` is VEP reference context. Reading the Hudson label as reference-context naming is explicitly an inference. Lossless recovery is disclosed as post-hoc because Sid scores were already inspected before the protocol freeze. | exploratory protocol; not validation | `artifacts/milestone_7_decision/peptide_recovery/EXPLORATORY_PROTOCOL.md` |
| M7 | **Dynamic upstream gate** (layered safe-rejection) | Built a label-blind selective-prediction gate (`configs/frozen/dynamic_gate_v1.json`, `dynamic-gate` CLI, 16 safety tests): Layer0 deterministic impossibility + Layer1 AND-of-core-vetoes `{el,prime}` (missing→KEEP) + Layer2 expression rescue-only + Layer3 patient-adaptive rails; CP-lower-bound (Neyman–Pearson) LOCO calibration; frozen before LOCKED CheckMate 153. **Feasibility:** the AND-gate **Pareto-dominates** the incumbent pure-EL gate on positive retention at matched removal. **The aggressive premise is FALSIFIED but SCOPED** (post-review correction): (a) 50–75% removal at ≥95% retention is unreachable on peptide-only features (IMPROVE positives spread); (b) at any SAFE retention downstream Δhits@20 = **0**. **Circularity audit** (`CIRCULARITY_AUDIT.md`): the v1 veto axes {EL,PRIME} ARE the dominant downstream-ranker inputs, so 0 removed candidates ever sat in a ranker top-20 and 0% of the high-EL/high-PRIME decoy stratum is removed — the null is a **structural tautology of feature overlap**, falsifying **same-feature presentation gating ONLY**. The general **orthogonal-feature** gate is UNTESTED / data-blocked, NOT falsified. Two more falsified/inconclusive v2 probes: an **independent sequence-only** hard-decoy gate collapsed cross-study (train-Gartner→test-IMPROVE retained **1.5%** of positives → sequence motifs encode study shift → v2 forbids them, requires leave-one-study-out + OOD abstention); an **in-sample multimer** orthogonal probe shows only a faint signal (AUROC 0.59). CP-LB sample-capped (46 pos ⇒ 0.937). Honest v1 value = **recall-preserving candidate-universe pruning** (CheckMate LOCKED: 10.4% neg removed @ 95.7% retention, no regression), NOT top-20 lift. v2 reframed to **NET-utility reselection** (maximize top-20 hits; positives may be sacrificed if removing higher-ranked decoys backfills more): implemented fixed-budget + counterfactual-replacement policies, outer leave-one-study-out (`src/event_b/dynamic_gate_v2.py`, `scripts/dynamic_gate_v2.py`, `V2_REPORT.md`). **v2 ALSO fails, decisively:** the **random-matched-pool control matches/beats** the learned selection (Gartner Δ+0.038 vs random +0.04–0.07) ⇒ the only positive number is a **denominator × 25-mer-regime** artifact; IMPROVE (deployable 9mer) and multimer are NEGATIVE; the counterfactual policy failed to abstain because a bootstrap-logistic ensemble is over-confident (needs conformal uncertainty). Five independent angles agree (incl. two independent experiments: HistGBT direct-utility, sequence-only which collapsed OOD retaining 1.5% of positives) ⇒ **no source-invariant negative-risk signal among top ranks; recognition wall holds at top-20. NO v2 frozen.** True unlock (spec'd + pre-registered, `V2_CONTRACT.md`/`V2_PREREGISTRATION.md`) = orthogonal-residual gate on mutant-RNA-VAF / DNA-VAF+depth/CCF / processing / agretopicity in **minimal-peptide-regime (8-11mer)** cohorts, conditioned on but never vetoing with PRIME/EL, behind a feature/regime OOD router with calibrated abstention — blocked on Miller IPV `PRJNA980652` + Gartner class-I reconstruction. | B (safe pruner) / C SCOPED (same-feature gating) + v2 net-utility also data-limited | `artifacts/milestone_7_decision/dynamic_gate/{SPEC,FEASIBILITY,REPORT,VERDICT,CIRCULARITY_AUDIT,V2_CONTRACT,V2_PREREGISTRATION,FALSIFICATION_LEDGER}.md`, `src/event_b/dynamic_gate.py`, `scripts/{dynamic_gate,dynamic_gate_orthogonal_probe}.py`, `configs/frozen/dynamic_gate_v1.json` |
| M7 | **Rich-feature dynamic gate** (base-anchored residual; IMPROVE 88-col table) | The orthogonal WES/RNA features v2 called "missing" were **present locally all along** — the raw IMPROVE table has 88 cols (DNA VAF, mutant-RNA af/reads, Stability, DAI/RankEL_wt agretopicity, Foreigness/SelfSim, physchem); v1/v2 used only the 7-col pool export. **Rich-gate v1** (orthogonal-only q as REPLACEMENT utility) = **NULL/falsified** (HistGBT Δ−0.200, worse than random; discarded the base rank). **Rich-gate v2 = BASE-ANCHORED** (`U=base_pct+α·(feat_pct−0.5)`, α by inner-CV, Epicurus UNCHANGED anchor, exact no-op at α=0): **PropHydroAro Δ+0.371 CI[0.114,0.629] beats matched-random (+0.055)** on nested patient-disjoint IMPROVE Partitions = **FIRST real dynamic-gate signal**; source-invariant within IMPROVE (bladder/melanoma/Basket all +; positives systematically more hydrophobic than negatives). **BUT external transfer FALSIFIED** (frozen→multimer −0.42/−0.58, Gartner −0.15): hydrophobicity is IMPROVE-assay/regime-specific. ⇒ architecture VALIDATED; **regime-aware abstention mandatory** (activate only if local-significant ∧ source-invariant ∧ external-non-harmful); **DO NOT freeze hydrophobicity**. Development/local only; next = untouched external RICH 8-11mer cohort + transport-supported features. | DEV-positive local / FALSIFIED external (regime-local) | `artifacts/milestone_7_decision/rich_gate/{REPORT,BASE_ANCHORED_REPORT}.md`, `src/event_b/rich_gate.py`, `scripts/{rich_gate_experiment,rich_gate_base_anchored}.py`, ledger `dynamic_gate/FALSIFICATION_LEDGER.md` |
| M7 | **Lossless input-only peptide recovery** (Phase 2/3 impl + online/offline run) | Frozen input-only generator (`lossless-peptide-generation-1.0.0`): raw GRCh38 allele → Ensembl VEP MANE/canonical (URL+SHA offline cache, fail-closed) → CDS/protein → all standard-AA 8–14 windows spanning the mutant/novel-frame residue → HLA panel from pVAC → **genuine PRIME** → union with pVAC (+`_cache_prime`) on stable genomic identity → frozen router route-aware top-20 (`genuine_prime = −PRIME %rank`). Reads NO assay/vaccine/label input (test-enforced); Hudson labels joined only AFTER ranking (eval-only). **Reproduced the disclosed post-hoc feasibility exactly:** ASPM **77** windows (best %rank **0.088**), MAP2 Gly868fs **259** windows (**0.010**), DYNC1H1 control **77** windows reproducing pVAC `KRFHATISF` (0.002). Adding the recovered candidates moves the 3-target coverage from 1/3 → 3/3. **⛔ CORRECTED/WITHDRAWN (2026-07-12): the "3/3" was TARGET-LEAKED** — `osteosarc_peptide_recovery.py` hard-codes `TARGETS = {ASPM, MAP2, DYNC1H1}` (the exact positives) and generates only for those (covers 10.2% of the 147-variant label-blind universe; `assert_generation_label_blind` fails it). Superseded by the honest label-blind end-to-end run (next row). This entry stands only as a **target-conditioned reconstruction feasibility test**. | reachability feasibility; **L3 claim WITHDRAWN (leakage)** | `artifacts/milestone_7_decision/peptide_recovery/{…,NORTH_STAR_FINAL.md (CORRECTION banner)}`, superseded by `sid_benchmark/` |
| M7 | **Sid identical-input end-to-end benchmark** (label-blind, falsification-first) | Corrects the leakage above. Complete label-blind universe = **200 public → 147 class-I eligible** variants (`variant_vafs_long.tsv`; eligibility = coding consequence ∧ tumor alt-reads>0, declared before label join); exact eval-only positives = **3** (ASPM/DYNC1H1/MAP2 exact IDs), joined post-freeze; hard leakage guard (test-enforced) + generic `expected=None` generation. **Ran the full universe:** 137 supported (missense+fs) → **130 generated ok / 7 VEP-failed / 10 unsupported (stop_gained,inframe)** → 59,755 peptide×HLA scored by genuine MixMHCpred+PRIME (coverage **88.4%**, honestly <95%). **Honest mutation-level hits@20 = 2/3** (NOT 3/3): MixMHCpred DYNC1H1#2 + ASPM#20, misses MAP2#25; genuine PRIME DYNC1H1#3 + MAP2#10, misses ASPM#39. All 3 positives reach scoring ⇒ the miss is at **conditional RANKING within the full universe**, not generation; the leaked 3/3 was an artifact of pooling only the 3 targets. 2/3 is an optimistic upper bound (17 uncovered variants would only add competitors). Post-hoc n=1/3, descriptive. Competitors: pVAC 2025.01 = boundary-mismatch (not identical inputs); Vaxrank/nextNEOpi/NeoDisc = NOT_EVALUABLE (input boundary/FASTQ/license). | honest **2/3** end-to-end; no superiority claim | `artifacts/milestone_7_decision/sid_benchmark/{BENCHMARK_PROTOCOL,REPORT}.md,{generation,variant_universe}.json,per_variant.csv,top20.csv`, `src/event_b/sid_benchmark.py`, `scripts/sid_benchmark_generate.py`, tests `tests/test_sid_benchmark.py` |

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
- **Evidence router (2026-07-12) — the Sid "recall is the wall" finding is now CODE, and confirms the current
  reranker cohorts cannot test it.** The frozen router keeps weak-RNA/atypical/single-caller candidates
  eligible-but-flagged and reports missing peptide/HLA as `NEEDS_PEPTIDE_GENERATION` (never a ranker miss); the
  multi-caller variant union recovers **3/3** of Sid's Hudson targets at candidate *generation* vs the single
  pVACtools step's 1/3. But the no-regression replay shows the router is a **pure-score pass-through (Δ=0)** on
  the multimer reranker artifacts and **NOT_EVALUABLE on Gartner (no HLA in the stored file)** — i.e. the
  reranker cohorts lack the RNA/expression/LOH/multi-caller features the router acts on, so its value is
  realised at **candidate generation**, not reranking. **This sharpens lever #2:** the needed cohort must carry
  a **raw multi-caller variant callset + WES/RNA/HLA through to peptide generation** (not just a pre-generated
  peptide×HLA list), so recall recovery and route-aware selection can actually be measured.
- **CheckMate 153 (2026-07-12) — first independent PRIME-untouched external test executed; TIE.** Alban *Nat
  Med* 2024 combinatorial-tetramer NSCLC screen, downloaded openly (Nature static host), own model score
  discarded, genuine PRIME 2.1 + MHCflurry-EL + study RNA-seq computed here. 14 pts / 1,197 class-I HLA-resolved
  9-mers / 162 tetramer+ / 1,035 tested-neg. **Frozen Epicurus v0.1 vs genuine PRIME: Δhits@20 +0.071, CI spans
  0 → TIE** (both all-negatives and the binders-only hard cut). **Load-bearing:** restricted to predicted
  binders (TETRAMER+ vs TETRAMER−), recognition **AUROC collapses to ~0.53–0.60 for every method including
  genuine PRIME (0.597)** — presentation (MixMHCpred 0.735 w/ easy negs) is the ceiling, recognition-among-
  presented is near-chance for the whole field. Sixth straight parity, now external. Underpowered (14 pts, wide
  CIs) → directional confirmation of the wall, not a new gate. Runner `scripts/checkmate153_dev.py`; artifact
  `artifacts/milestone_7_decision/checkmate153/`; manifest `configs/source_manifests/checkmate153.yml`.
- **osteosarc.com / Sid Sijbrandij (2026-07-12) — first DEPLOYMENT-grade patient with a MEASURED recognition
  label; recall (not ranking) is the wall.** PUBLIC open-data osteosarcoma case (osteosarc.com / Research to the
  People; `b2://osteosarc-data`, no DUA), longitudinal T0–T3, full North-Star input. **Distinct patient from
  RTTP `SR24-58221`** (disjoint HLA A\*01:01/B\*08:01/B\*27:05/C\*01:02/C\*07:01; the on-disk
  `artifacts/osteosarc_audit/`+`osteosarc_product/` are that RTTP file MISLABELED as osteosarc → disregard).
  **Measured label = Hudson-Lab IFNγ peptide-expansion assay** (PBMC stim → IFNγ+/− sort → MiXCR TCR-seq →
  mutation-specific expansion): recognized = ASPM p.G2179R (May+Aug), DYNC1H1 p.V314I (May+Aug), MAP2 …868fs
  (May). Ran genuine PRIME + frozen Epicurus v0.1 (OUT-OF-SAMPLE, no fit) vs the full pVACtools ensemble on the
  21 curated candidates (14,780 pep×HLA). **Two results:** (1) **RECALL = 1/3** — only DYNC1H1 (TPM 357) reached
  the shortlist; MAP2 (TPM 5.2) was expression-filtered, ASPM was off the single 2025.01 callset → 2/3
  recognized neoantigens dropped *before* any ranking. (2) **Ranking is saturated** — presentation ranks
  DYNC1H1 #1 of 21, Epicurus matches it, but the two *immunogenicity* models do WORSE (BigMHC_IM #3, DeepImmuno
  #18). **Take:** the recoverable failure is the candidate FUNNEL (recall + single-timepoint), not the ranker;
  adding a learned recognition score on top of presentation hurts. Principled (unvalidated, NOT fit to n=3)
  levers: recall-first candidate policy, cross-timepoint mutanome union, high-recall + honest abstention.
  n=3 positives (1 in-universe) → DIAGNOSTIC, not a gate. Runner `scripts/osteosarc_rank.py`; artifact
  `artifacts/milestone_7_decision/osteosarc_sid/`; manifest `configs/source_manifests/osteosarc_sid.yml`.
  **Recall gap now shown mechanically recoverable (post-hoc):** the frozen input-only lossless generator
  regenerates ASPM (77) + MAP2 Gly868fs (259) windows from the raw allele + Ensembl alone and, scored by
  genuine PRIME, lifts 3-target coverage 1/3 → 3/3 at every stage (`peptide_recovery/`). This is a
  reachability fix feeding genuine PRIME, **not** a superiority result — it is post-hoc on the motivating
  patient. **Highest-value follow-up: ask Hudson Lab/RTTP for the stimulation-pool composition (true
  denominator) + more patients with the assay → turns n=3 descriptive into a real deployable-patient
  benchmark, and gives the lossless generator a genuinely untouched cohort to validate on.**
- **Indicated levers (model form now exhausted on these 3 cohorts):**
  1. ~~Assay/regime heads~~ — **DONE (v0.5); closes at a tie.** Do NOT iterate v0.6 model form against the
     current 3 cohorts expecting ACCEPT — v0.3/v0.4/v0.5 prove parity is data-bound, not a tuning miss.
  2. **DATA — acquire a second Gartner-like denominator cohort** (real per-patient candidate ranking + explicit
     immunogenicity negatives + WES/RNA). Now the single highest-value action; see the data-acquisition search
     brief (`artifacts/milestone_7_decision/external_validation/`). **Partial delivery: CheckMate 153 acquired +
     run (above) — independent + untouched but small (14 pts) and range-restricted; it confirms rather than
     clears. Still needed: a DENSE-denominator, well-powered cohort. Best open lead = Miller IPV `PRJNA980652`
     (label table behind the STM paywall → needs user download; raw seq 0.23TB). CheckMate raw WES/RNA (dbGaP)
     would complete its inputs. WashU TNBC = mostly CD4/long-peptide + no class-I HLA typing → weak for PRIME.**
  3. **External proof** — acquire Hu/Ott via dbGaP phs001451
     (`artifacts/milestone_7_decision/external_validation/ACQUISITION_PACKET.md`).
- **Not yet spent:** one pre-registered look at the SEMI-CONSUMED Gartner TEST holdout is still available, but
  only *after* a development candidate clears the gate — none of v0.3/v0.4/v0.5 did, so the holdout stays closed.
  **A superiority claim still requires a win on a genuinely untouched cohort.**

- **DECISION (2026-07-12) — the benchmark is a THREE-LEVEL hierarchy; never pool cohorts into one headline.**
  Product review: do NOT make heterogeneous pooled reranker performance the central metric. Three DISTINCT
  tasks, each reported/interpreted on its own: **(1) REACHABILITY** — raw variants/WES/RNA/HLA *through peptide
  generation*, with stage-loss attribution; **(2) CONDITIONAL RANKING** — only among generated/rankable
  candidates, strictly within each cohort's own assay/denominator; **(3) END-TO-END PATIENT UTILITY = the
  PRIMARY north star** — recognized mutations in the final top-20 from *common raw inputs* vs standard pVAC +
  genuine PRIME. **Fixed cohort roles (never merge):** CEDAR/Zhao = training/recognition-prior; cd8_multimer =
  presentation/T-cell compatibility (+ Epicurus v0.1 training ⇒ Epicurus arms leakage-invalid there); Gartner =
  conditional *broad-denominator* ranking; IMPROVE = *prefiltered-subset* ranking; CheckMate153 = external
  conditional (small); osteosarc/Sid = end-to-end diagnostic (all 3 levels, post-hoc n=3); RTTP SR24 =
  deployment only. The four-arm generation×scorer harness (`src/benchmark/four_arm.py`, tested) is **reusable
  infrastructure**; it yields an end-to-end HEADLINE only where a cohort is Level-3 eligible — **currently Sid
  only.** Eligibility audit encodes all three levels + roles + the no-pooling invariant
  (`src/benchmark/cohort_audit.py`; runner `scripts/four_arm_benchmark.py`; artifacts
  `artifacts/milestone_7_decision/four_arm/`). Sid L1 reachability = generation recall 1/3→3/3; L3 end-to-end:
  under genuine PRIME `lossless_prime` recovers all 3 (protected incumbent).
- **FAIR feature re-run (2026-07-12) — genuine MHCflurry presentation features on recovered candidates.**
  NetMHCpan (the frozen Epicurus `el`) is not locally runnable; recovered candidates were previously el=0.5
  imputed. Computed GENUINE MHCflurry presentation %rank/affinity/processing for every candidate
  (`src/benchmark/presentation_features.py`, tested; cached to `four_arm/mhcflurry_presentation_cache.csv`).
  Re-running the four-arm attribution with one consistent genuine predictor: frozen Epicurus scorer stage
  **−1** (was −2 under the 0.5 impute → the extra −1 was an imputation artifact), full stack **+1** vs
  pVAC+PRIME. Epicurus recovers ASPM once it has real presentation evidence but STILL drops low-expression
  MAP2 (a real expression-reweighting effect — a learned recognition score on presentation does not help).
  MHCflurry↔NetMHCpan-EL Spearman 0.52 (moderate) disclosed: the fair run approximates, not reproduces, the
  frozen NetMHCpan-EL feature. Descriptive, n=3, no constant tuned to Sid.
- **Expression ranking-policy frozen (2026-07-12) — RNA expression is CONFIDENCE-ONLY, not a rank penalty.**
  Label-blind analysis on the 3 Level-2 development cohorts (multimer/Gartner/IMPROVE), each within its own
  denominator, NEVER pooled, under strict no-regression vs the protected lossless+PRIME incumbent
  (`src/benchmark/expression_policy.py` tested; runner `scripts/expression_policy_analysis.py`; frozen
  `configs/frozen/expression_policy_v1.json`; artifact `.../expression_policy/`). Structure differs by
  denominator: Gartner recognition concentrates in the top expression quartile (penalty would help),
  IMPROVE is flat with 46% of positives low-expression, multimer weakly monotone with 35% low-expression.
  Result: **expr_rank_penalty** helps Gartner (recall .37→.57) but REGRESSES multimer (.62→.47) and IMPROVE
  (.18→.16); the fixed-budget **portfolio reserve** regresses Gartner (.37→.30). The only no-regression-
  everywhere forms are **prime_only (confidence-only)** and the presentation-protected **soft_saturating**,
  and the latter is IDENTICAL to prime_only on all 3 cohorts (it only clears low-presentation+low-expression
  junk that never reaches top-20). ⇒ Frozen decision: expression is **confidence-only** in the score;
  lossless+PRIME stays the protected ranking; soft-saturating kept as an equivalent route-dependent guard;
  portfolio reserve retained OPTIONAL/off-by-default. Constants (bottom-quartile stratum) fixed, NOT tuned to
  any eval. Sid DESCRIPTIVE (post-freeze, n=3): confidence-only keeps 3/3 recognized in top-20, the expr
  penalty drops low-expression MAP2 to 2/3 — consistent with the dev finding, nothing tuned to Sid. See
  `memory/expression-ranking-policy.md`.
- **Acquisition packet (2026-07-12) — canonical, primary-source-verified plan to obtain an untouched
  end-to-end patient.** `artifacts/milestone_7_decision/external_validation/`:
  `ACQUISITION_EXECUTION_PLAN.md` (cohorts ranked by benchmark level × leakage × denominator × negatives ×
  input-recoverability × effort × north-star value; split policy; go/no-go), `COHORT_ACQUISITION_TRACKER.csv`
  (13 cohorts, 19 columns), `MINIMUM_PATIENT_DATA_PACKAGE.md` (exact de-identified request schema),
  `AUTHOR_OUTREACH_DRAFTS.md` (7 tailored data-not-PHI requests), + verify-first manifests
  `configs/source_manifests/{miller_ipv,southampton_nsclc}.yml`. Web-verified corrections: **Miller "IPV" =
  identify-prioritize-validate PLATFORM** (raw WES+RNA OPEN at `PRJNA980652`, ELISpot negatives → the ONLY
  fully-open L3 build); **CheckMate raw = EGA `EGAD00001011302` controlled (BMS discretion), not dbGaP**;
  **EVX-01 (PMC11116868) = closed data**; GBM `GSE237936` = RNA-only n=4 (2024, not 2026); medRxiv aggregate
  UNLOCATABLE. Split: LOCK >=2 low-leakage cohorts (Miller IPV + CheckMate; Southampton secondary), DEV =
  phs001003/Gartner/IMPROVE/multimer, TRAIN = CEDAR/IEDB/Zhao; by patient & study, no random peptide split.
  **Go/no-go: no headline unless raw->generation->same universe->PRIME AND Epicurus->paired hits@20
  reproduces.** Single best user action: verify Miller IPV S1/S2 supplement is downloadable (science.org 403
  to bot), then build; parallel DARs for phs001003 + CheckMate EGA. See `memory/acquisition-packet.md`.
  Next = acquire/identify an **untouched end-to-end patient**, not more rerankers. See
  `memory/benchmark-three-level-hierarchy.md`.

## Maintenance

Update this ledger (row in the chronological table + any role/verdict change) as the final step of each
milestone. Link artifacts; never paste result tables here — they drift. Memory index: `memory/MEMORY.md`.
