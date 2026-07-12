# Sid full Epicurus filter-stack benchmark

> Post-hoc n=1 patient / 3 recognized mutations. All selections were serialized before the exact Hudson labels were imported. Upstream generation covers 130/147 eligible mutations (88.4%), so this is a diagnostic, not a general superiority claim.

## Filter effects (label-blind)

- Deterministic validity: removed **0** candidate rows and **0** mutations. Sid HLA-LOH is unavailable.
- Dynamic gate v1: removed **1,561** rows (2.6%) but **0 mutations**.
- Product RNA policy: excluded **27,875** rows; 62 mutations remained.

## Recognized mutations in the selected top 20

| Arm | Hits / 3 | Selected routes | Unique mutations | Hit mutations |
|---|---:|---:|---:|---|
| `baseline_genuine_prime` | **2/3** | 20 | 20 | DYNC1H1-chr14-101980529, MAP2-chr2-209694772 |
| `baseline_frozen_epicurus_v0_1` | **1/3** | 20 | 20 | DYNC1H1-chr14-101980529 |
| `deterministic_then_prime` | **2/3** | 20 | 20 | DYNC1H1-chr14-101980529, MAP2-chr2-209694772 |
| `dynamic_then_prime` | **2/3** | 20 | 20 | DYNC1H1-chr14-101980529, MAP2-chr2-209694772 |
| `deterministic_dynamic_then_prime` | **2/3** | 20 | 20 | DYNC1H1-chr14-101980529, MAP2-chr2-209694772 |
| `deterministic_dynamic_then_frozen_epicurus` | **1/3** | 20 | 20 | DYNC1H1-chr14-101980529 |
| `mutation_level_dynamic_sensitivity_then_prime` | **2/3** | 20 | 20 | DYNC1H1-chr14-101980529, MAP2-chr2-209694772 |
| `deterministic_dynamic_router_prime` | **2/3** | 20 | 17 | DYNC1H1-chr14-101980529, MAP2-chr2-209694772 |
| `deterministic_dynamic_router_frozen_epicurus` | **1/3** | 20 | 13 | DYNC1H1-chr14-101980529 |
| `product_v1_as_shipped` | **1/3** | 20 | 10 | DYNC1H1-chr14-101980529 |
| `full_stack_mutation_level_fair` | **1/3** | 20 | 20 | DYNC1H1-chr14-101980529 |
| `full_current_epicurus_stack` | **1/3** | 20 | 10 | DYNC1H1-chr14-101980529 |

## Verdict

- Genuine PRIME baseline: **2/3**.
- Complete currently runnable Epicurus stack: **1/3**.
- Filters improve PRIME on this patient: **False**.
- Full stack beats PRIME on this patient: **False**.

The deterministic gate has nothing impossible to remove in the generated pool. The dynamic gate removes low-scoring peptide×HLA routes but no whole mutations and does not change PRIME's mutation-level top-20. With the newly recovered matched T2 RNA counts, the legacy product policy hard-excludes zero-mutant-read mutations and falls to 1/3. That policy is over-aggressive here: MAP2 is experimentally recognized despite zero observed mutant-allele RNA reads.

Matched T2 RNA depth, mutant reads, and RNA VAF are now evaluated for every generated mutation. Sid-specific HLA-LOH remains unavailable and is not silently imputed. A separate absence-safe RNA gate is reported in sid_benchmark/rna_support_arm.json.
