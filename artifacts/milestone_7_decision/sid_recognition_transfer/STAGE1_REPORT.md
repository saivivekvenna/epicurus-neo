# Recognition transfer — STAGE 1 (non-Sid freeze, leakage-clean)

_frozen SHA-256 `e44d3d1ac66a272e`; data: data/raw/gartner_nci/_cache_improve_prime.tsv, data/raw/improve/data.zip (NO Sid). Quarantined 731/17520 peptide-leaked held-out rows; selection on leakage-clean only. null hits clean=80.0._


## Nested CV per family (inner mask recomputed within outer-train; null included, conservative tie-break)

| family | nested total | null | Δ | per-cohort Δ | paired-boot Δ [CI95] p>0 | matched-random (beats) | eligible | chosen (fold:α,q) |
|---|--:|--:|--:|---|---|---|:--:|---|
| core_deployable | 83 | 80 | +3 | {'Basket': 1.0, 'bladder': 1.0, 'melanoma': 1.0} | +0.043 [-0.0429, 0.1286] 0.8004 | 77.2 (y) | yes | [(0, 0.1, 1), (1, 0.0, 3), (2, 0.5, 1), (3, 0.1, 1), (4, 0.2, 1)] |
| improve_rich_partial_bridge | 83 | 80 | +3 | {'Basket': 0.0, 'bladder': 2.0, 'melanoma': 1.0} | +0.043 [-0.0429, 0.1429] 0.7726 | 77.4 (y) | yes | [(0, 0.1, 1), (1, 0.0, 3), (2, 0.5, 3), (3, 0.3, 1), (4, 0.2, 1)] |

**Selected family (by nested evidence): core_deployable**. Deployment params chosen by full non-Sid CV within that family only:


| arm | C | α | q | clean hits (Δ) | per-cohort Δ | matched-random | beats rand |
|---|--:|--:|--:|--:|---|--:|:--:|
| core_deployable | 0.5 | 0.1 | 1 | 88 (+8) | {'Basket': 2.0, 'bladder': 3.0, 'melanoma': 3.0} | 82.55 | y |
| core_deployable | 1.0 | 0.1 | 1 | 88 (+8) | {'Basket': 2.0, 'bladder': 3.0, 'melanoma': 3.0} | 82.55 | y |
| core_deployable | 2.0 | 0.5 | 1 | 88 (+8) | {'Basket': 2.0, 'bladder': 3.0, 'melanoma': 3.0} | 83.25 | y |
| core_deployable | 0.5 | 0.25 | 1 | 87 (+7) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 4.0} | 82.25 | y |
| core_deployable | 1.0 | 0.25 | 1 | 87 (+7) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 4.0} | 82.2 | y |
| core_deployable | 1.0 | 0.5 | 1 | 87 (+7) | {'Basket': 2.0, 'bladder': 3.0, 'melanoma': 2.0} | 82.25 | y |
| core_deployable | 2.0 | 0.1 | 1 | 87 (+7) | {'Basket': 1.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.55 | y |
| core_deployable | 0.5 | 0.2 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.3 | y |
| core_deployable | 0.5 | 0.5 | 1 | 86 (+6) | {'Basket': 1.0, 'bladder': 3.0, 'melanoma': 2.0} | 81.25 | y |
| core_deployable | 1.0 | 0.2 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.3 | y |
| core_deployable | 2.0 | 0.2 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.4 | y |
| core_deployable | 2.0 | 0.25 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.2 | y |
| core_deployable | 0.5 | 0.0 | 1 | 85 (+5) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 2.0} | 79.55 | y |
| core_deployable | 0.5 | 0.3 | 1 | 85 (+5) | {'Basket': 1.0, 'bladder': 2.0, 'melanoma': 2.0} | 80.3 | y |

## FROZEN (§6, leakage-clean)
arm **core_deployable** α=0.1 q=1 C=0.5 clean_hits=88; multiplicity 180 configs; external transport {'gartner': {'base': 18, 'gated': 18, 'delta': 0, 'note': 'ANCHORED-COMPONENT-ONLY: q reserve DISABLED externally (no comparable mutant-RNA signal); VAF/RNA features neutral where absent. Not the full frozen policy.'}, 'multimer': {'base': 21, 'gated': 22, 'delta': 1, 'note': 'ANCHORED-COMPONENT-ONLY: q reserve DISABLED externally (no comparable mutant-RNA signal); VAF/RNA features neutral where absent. Not the full frozen policy.'}}. SHA-256 `e44d3d1ac66a272e89cb00aaf85dc856b905e76d8c8712665f6e3f5b13214890`.

_patient-level paired bootstrap CI for the nested per-family delta is COMPUTED and reported (see nested_cv_per_family[*].paired_bootstrap_delta_vs_null; fixed seed 12345), NOT gated. Stage 1 was NARROWED to the exact runner vs the original contract (see CONTRACT.md 'STAGE-1 v2 PROTOCOL CORRECTION'): PRIME base only (MixMHCpred base deferred), 2 anchored full-feature arms (not the physchem/expression/absence families), 5 IMPROVE official partitions (Gartner/multimer are an ANCHORED-COMPONENT-ONLY external transport check, q disabled there). Sid previously inspected -> Stage 2 is exploratory confirmation, not pristine validation._

> PRE-SID: no Sid file/label accessed; STOP for audit before Stage 2.
