# Lossless peptide recovery — osteosarc.com / Sid (EXPLORATORY, post-hoc)

> **⛔ CORRECTION (2026-07-12): TARGET LEAKAGE — not end-to-end.** The generator hard-codes the 3 recognized
> positives as `TARGETS` and generates only for them (covers 10.2% of the 147-variant label-blind universe;
> `assert_generation_label_blind` fails it). The recall/coverage figures here are a **target-conditioned
> sensitivity** test, NOT proof of end-to-end recovery. L3 claim withdrawn. See
> `../sid_benchmark/BENCHMARK_PROTOCOL.md`.

> Generator policy `lossless-peptide-generation-1.0.0` composed with router policy `epicurus-evidence-router-1.0.0`. Mode: **online**. Genuine PRIME commit `7b18d4e110`.

> **Status:** post-hoc reachability diagnostic on the patient that motivated it — NOT preregistered / blind / independent. This does **not** show Epicurus beats PRIME: the selection score IS genuine PRIME (`genuine_prime = -PRIME %rank`); better candidate GENERATION lets genuine PRIME score targets it previously never received.

## Generation (input-only; no assay/vaccine/label input)

| Variant | Windows | Unique peptides | Peptide×HLA pairs | Best genuine-PRIME %rank | Role |
|---|---:|---:|---:|---:|---|
| ASPM | 77 | 77 | 385 | 0.0880 | recovered |
| MAP2 | 259 | 259 | 1295 | 0.0100 | recovered |
| DYNC1H1 | 77 | 77 | 385 | 0.0020 | positive control |

- pVAC candidate rows: **14780**; union rows: **7565**; recovered pairs scored by genuine PRIME: **2065/2065**.

## Coverage of the 3 Hudson-expanded targets (labels joined AFTER ranking, evaluation only)

| Stage | pVAC-only | Augmented (pVAC + lossless recovery) |
|---|---|---|
| candidate generation recall | 1/3 (DYNC1H1-chr14-101980529) | 3/3 (ASPM-chr1-197102716, DYNC1H1-chr14-101980529, MAP2-chr2-209694772) |
| rankable recall (peptide+HLA) | 1/3 (DYNC1H1-chr14-101980529) | 3/3 (ASPM-chr1-197102716, DYNC1H1-chr14-101980529, MAP2-chr2-209694772) |
| pure genuine-PRIME top-20 | 1/3 (DYNC1H1-chr14-101980529) | 3/3 (ASPM-chr1-197102716, DYNC1H1-chr14-101980529, MAP2-chr2-209694772) |
| route-aware top-20 | 1/3 (DYNC1H1-chr14-101980529) | 3/3 (ASPM-chr1-197102716, DYNC1H1-chr14-101980529, MAP2-chr2-209694772) |

- Content hash (mode-invariant): `74bbd8fd6947275d61c51b0ab3f8e88d5f6b133026c667719871e21ec6e183e1`.

## Interpretation guardrail

Post-hoc reachability fix on the motivating patient. NOT a claim that Epicurus beats PRIME: the selection score IS genuine PRIME (genuine_prime = -PRIME %rank); better candidate GENERATION lets genuine PRIME score targets it previously never received. A benefit claim needs a future untouched cohort.
