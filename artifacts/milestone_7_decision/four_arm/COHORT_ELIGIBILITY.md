# Benchmark cohort eligibility audit — three-level hierarchy

> Policy `three-level-benchmark-1.0.0`. 9 cohorts; **2** end-to-end (Level-3) eligible.

> **No pooling.** These cohorts serve DIFFERENT tasks over DIFFERENT denominators (broad-somatic Gartner, prefiltered IMPROVE, presentation/T-cell multimer, external CheckMate, end-to-end Sid) and are NEVER pooled into one headline metric. A single cross-cohort reranker number would be a category error. Each level and each cohort is reported and interpreted on its own.

> **Four-arm harness** = reusable infrastructure; produces an end-to-end headline ONLY for Level-3-eligible cohorts (currently osteosarc_sid, post-hoc).

## The three levels

1. **reachability** — raw variants/WES/RNA/HLA through peptide generation, with stage-loss attribution (recognized mutations lost at generation vs rankable vs top-k)
2. **conditional_ranking** — ranking only among candidates actually generated/rankable, interpreted strictly within this cohort's own assay/denominator (never pooled)
3. **end_to_end_patient_utility** — PRIMARY north star: recognized mutations in the final top-20 from common raw inputs, vs a standard pVAC-style pipeline + genuine PRIME

## Per-cohort level eligibility

✅ = eligible · ❌ = not eligible (each cohort interpreted ONLY within its own denominator)

| cohort | role | L1 reachability | L2 conditional ranking | L3 end-to-end (north star) | denominator |
|---|---|---|---|---|---|
| **osteosarc_sid** | End-to-end diagnostic (all 3 levels; post-hoc, n=3, single patient) | ✅ | ✅ | ✅ | per-patient somatic pVAC universe (21 curated muts) + input-only lossless recovery |
| **miller_ipv** | End-to-end LOCKED_TEST (labels ingested; run blocked on WES download + HLA/expr compute) | ✅ | ✅ | ✅ | 13 patients; 754 tested 20-mers over 343 IPV-prefiltered mutations (re-enumerate full mutanome from open WES for a fair denominator) |
| **gartner_nci** | Conditional broad-denominator ranking (Level 2) | ❌ | ✅ | ❌ | Gartner TEST minimal peptide-HLA list (SEMI-CONSUMED holdout); no raw callset->gen |
| **improve_srhgroup** | Conditional prefiltered-subset ranking (Level 2) | ❌ | ✅ | ❌ | pre-screened ~200 candidates/patient — NOT the full somatic universe |
| **checkmate153** | External conditional ranking (Level 2; small, PRIME-untouched) | ❌ | ✅ | ❌ | HLA-resolved 9-mer candidate list (14 patients); raw WES/RNA behind dbGaP |
| **cd8_multimer** | Presentation / T-cell compatibility asset (and Epicurus v0.1 training set) | ❌ | ✅ | ❌ | pMHC-multimer peptide-HLA candidate list |
| **cedar_tcell** | Training / recognition-prior asset | ❌ | ❌ | ❌ | mutation-derived cancer T-cell recognition rows (not a per-patient ranking denom.) |
| **zhao_dc_2026** | Training / recognition-prior asset (genuine PRIME blocked) | ❌ | ❌ | ❌ | 352 pts / 2317 SNV peptide-HLA list |
| **rttp_sr24_58221** | Deployment only (no measured recognition label) | ❌ | ❌ | ❌ | complete North-Star INPUT (WES+RNA+HLA+candidate universe), 1 patient |

## Four-arm infrastructure evaluability (headline only where L3-eligible)

The generation×scorer matrix is reusable infra; read it as a HEADLINE only for the L3-eligible cohort. Elsewhere the evaluable arms are Level-2 conditional-ranking probes.

| cohort | pvac_prime | lossless_prime | lossless_epicurus | full_epicurus |
|---|---|---|---|---|
| **osteosarc_sid** | ✅ | ✅ | ✅ | ✅ |
| **miller_ipv** | ✅ | ✅ | ✅ | ✅ |
| **gartner_nci** | ✅ | ❌ no-lossless-gen | ❌ no-lossless-gen | ❌ no-lossless-gen |
| **improve_srhgroup** | ✅ | ❌ no-lossless-gen | ❌ no-lossless-gen | ❌ no-lossless-gen |
| **checkmate153** | ✅ | ❌ no-lossless-gen | ❌ no-lossless-gen | ❌ no-lossless-gen |
| **cd8_multimer** | ✅ | ❌ no-lossless-gen | ❌ no-lossless-gen, LEAKAGE:epicurus_v0_1_training_cohort | ❌ no-lossless-gen, LEAKAGE:epicurus_v0_1_training_cohort |
| **cedar_tcell** | ❌ no-PRIME | ❌ no-lossless-gen, no-PRIME | ❌ no-lossless-gen, no-Epicurus-feat | ❌ no-lossless-gen, no-Epicurus-feat, router_features |
| **zhao_dc_2026** | ❌ no-PRIME | ❌ no-lossless-gen, no-PRIME | ❌ no-lossless-gen, no-Epicurus-feat | ❌ no-lossless-gen, no-Epicurus-feat |
| **rttp_sr24_58221** | ❌ labels | ❌ labels | ❌ labels | ❌ labels |

## Notes

- **osteosarc_sid** (End-to-end diagnostic (all 3 levels; post-hoc, n=3, single patient)) — ONLY cohort with a raw callset carried to generation AND a measured label -> the only Level-3 (end-to-end) instance. POST-HOC, n=3, 1 patient: diagnostic, not a powered/blinded gate. Frozen Epicurus is OUT-OF-SAMPLE here (trained on cd8_multimer).
- **miller_ipv** (End-to-end LOCKED_TEST (labels ingested; run blocked on WES download + HLA/expr compute)) — INDEPENDENT, PRIME-untouched (2024, La Jolla). Raw WES/RNA OPEN (PRJNA980652); the S1/S2 label table is now INGESTED + validated (LOCKED_TEST, never used for development). ELIGIBLE but NOT-YET-EXECUTED: the four-arm run is blocked on the WES download + HLA typing + RNA expression + full-mutanome re-enumeration, not on any missing file. Labels are 20-mer/mutation-level (no HLA) -> benchmark runs at MUTATION granularity.
- **gartner_nci** (Conditional broad-denominator ranking (Level 2)) — Level-2 only. No lossless-generation arm: the stored artifact is a pre-generated peptide list, not a raw multi-caller callset. Holdout semi-consumed -> not a fresh test. Broadest somatic denominator of the ranking cohorts; interpret ONLY within itself.
- **improve_srhgroup** (Conditional prefiltered-subset ranking (Level 2)) — Level-2 only; denominator is range-restricted (prefiltered) so top-20 recall is optimistic and NOT comparable to Gartner's broad denominator. No raw callset -> no lossless arm.
- **checkmate153** (External conditional ranking (Level 2; small, PRIME-untouched)) — Level-2 only, underpowered (14 pts). Acquiring dbGaP WES/RNA would upgrade it toward Level 1/3. Independent + PRIME-untouched but range-restricted.
- **cd8_multimer** (Presentation / T-cell compatibility asset (and Epicurus v0.1 training set)) — Presentation/T-cell-compatibility role. It IS the frozen Epicurus v0.1 training cohort -> the Epicurus arms are leakage-INVALID here even with inputs present; only the genuine-PRIME scorer is honest. No raw callset -> no lossless arm. Level-2 for genuine PRIME only.
- **cedar_tcell** (Training / recognition-prior asset) — Recognition-PRIOR / training asset, not a ranking cohort: no per-patient somatic denominator, no genuine-PRIME incumbent, no raw callset. Feeds priors/representation, never a benchmark headline. Eligible for none of the three levels.
- **zhao_dc_2026** (Training / recognition-prior asset (genuine PRIME blocked)) — Recognition-prior/training asset. genuine PRIME is BLOCKED (only a MixMHCpred proxy; the incumbent guard forbids labeling a proxy PRIME) -> no Level-2/3 incumbent, no genuine-PRIME feature -> no Epicurus arm, no raw callset -> no Level 1. Eligible for none of the levels.
- **rttp_sr24_58221** (Deployment only (no measured recognition label)) — Deployment/DEMO asset, not a benchmark: no measured recognition label -> no denominator -> eligible for none of the three levels. Inputs would otherwise support all three.

## Interpretation

Exactly one cohort (osteosarc_sid) is END-TO-END (Level-3) eligible, and only post-hoc with n=3. Conditional-ranking (Level-2) cohorts each stand alone within their own denominator and are never pooled. Reachability (Level-1) needs a raw callset carried to generation, which only Sid provides among labelled cohorts. The binding constraint is a DENSE-denominator, PRIME-untouched, end-to-end patient — acquire/identify one; do not manufacture a pooled reranker headline.

