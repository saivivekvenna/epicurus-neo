# Canonical Epicurus product-path end-to-end audit

> Actual production normalization, deterministic gate, evidence score, eligibility policy, and capped portfolio. Both patients were previously inspected; this is an integration audit, not independent validation.

| Patient | Generated positive mutations | Deterministic-valid | Product-eligible | Product top 20 | PRIME plain | PRIME cap-2 |
|---|---:|---:|---:|---:|---:|---:|
| Hu_287 | 3/3 | 3/3 | 3/3 | **3/3** | 0/3 | 2/3 |
| Sid | 3/3 | 3/3 | 3/3 | **1/3** | 2/3 | 2/3 |

## Patient funnels

### Hu_287

- Candidate rows: 3,192
- Generated mutations: 14
- Deterministic-valid rows: 3,192
- Product-eligible rows: 2,280
- Selected routes / unique mutations: 20 / 10
- Exclusions: {'NO_MUTANT_RNA_SUPPORT': 912}
- Positive last-reached stages: {'16:10907146:C:T': 'selected', '17:7673535:C:G': 'selected', '20:18511052:C:T': 'selected'}

### Sid

- Candidate rows: 62,540
- Generated mutations: 137
- Deterministic-valid rows: 62,540
- Product-eligible rows: 34,665
- Selected routes / unique mutations: 20 / 10
- Exclusions: {'NO_MUTANT_RNA_SUPPORT': 20960, 'NO_RNA_EXPRESSION': 6915}
- Positive last-reached stages: {'ASPM-chr1-197102716': 'product_eligible', 'DYNC1H1-chr14-101980529': 'selected', 'MAP2-chr2-209694772': 'product_eligible'}

## Honest interpretation

This report is the deliverable-level check: every headline number comes from the same shipped product logic, not from a mix-and-match research arm. Hu_287 tests a complete local raw-data reconstruction; Sid accounts for all 147 eligible variants (137 generated plus 10 documented non-enumerable) and exposes downstream product losses. Neither patient is blind, so this establishes runnable behavior and patient-specific outcomes—not general superiority.
