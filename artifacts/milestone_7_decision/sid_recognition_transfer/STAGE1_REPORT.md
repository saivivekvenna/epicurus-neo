# Recognition transfer — STAGE 1 (non-Sid freeze, leakage-clean)

_frozen SHA-256 `d69c498d04abe9ef`; data: data/raw/gartner_nci/_cache_improve_prime.tsv, data/raw/improve/data.zip (NO Sid). Quarantined 731/17520 peptide-leaked held-out rows; selection on leakage-clean only. null hits clean=80.0._


| arm | C | α | q | clean hits (Δ) | per-cohort Δ | matched-random | beats rand |
|---|--:|--:|--:|--:|---|--:|:--:|
| core_deployable | 0.1 | 0.1 | 1 | 88 (+8) | {'Basket': 2.0, 'bladder': 3.0, 'melanoma': 3.0} | 82.3 | y |
| core_deployable | 1.0 | 0.1 | 1 | 88 (+8) | {'Basket': 2.0, 'bladder': 3.0, 'melanoma': 3.0} | 82.3 | y |
| core_deployable | 0.1 | 0.25 | 1 | 87 (+7) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 4.0} | 82.35 | y |
| core_deployable | 1.0 | 0.25 | 1 | 87 (+7) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 4.0} | 82.35 | y |
| improve_rich_partial_bridge | 1.0 | 0.3 | 1 | 87 (+7) | {'Basket': 1.0, 'bladder': 3.0, 'melanoma': 3.0} | 82.3 | y |
| core_deployable | 0.1 | 0.2 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.2 | y |
| core_deployable | 1.0 | 0.2 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.2 | y |
| improve_rich_partial_bridge | 0.1 | 0.1 | 1 | 86 (+6) | {'Basket': 1.0, 'bladder': 3.0, 'melanoma': 2.0} | 80.25 | y |
| improve_rich_partial_bridge | 0.1 | 0.2 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.25 | y |
| improve_rich_partial_bridge | 0.1 | 0.25 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.25 | y |
| improve_rich_partial_bridge | 0.1 | 0.3 | 1 | 86 (+6) | {'Basket': 1.0, 'bladder': 3.0, 'melanoma': 2.0} | 81.2 | y |
| improve_rich_partial_bridge | 1.0 | 0.2 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.25 | y |
| improve_rich_partial_bridge | 1.0 | 0.25 | 1 | 86 (+6) | {'Basket': 0.0, 'bladder': 3.0, 'melanoma': 3.0} | 81.25 | y |
| improve_rich_partial_bridge | 1.0 | 0.3 | 2 | 86 (+6) | {'Basket': 1.0, 'bladder': 3.0, 'melanoma': 2.0} | 81.55 | y |

## FROZEN (§6, leakage-clean)
arm **core_deployable** α=0.1 q=1 C=0.1 clean_hits=88; multiplicity 77 configs; external transport {'gartner': {'base': 18, 'gated': 18, 'delta': 0, 'note': 'VAF/RNA neutral where absent'}, 'multimer': {'base': 21, 'gated': 22, 'delta': 1, 'note': 'VAF/RNA neutral where absent'}}. SHA-256 `d69c498d04abe9ef4d44cf8731545bb2a1448a8943d1508b6860a6dabb96afd1`.

_bootstrap patient CI on +hits spans 0 (underpowered, n=70) -> reported not gated; Sid previously inspected -> Stage-2 is exploratory confirmation, not pristine validation._

> PRE-SID: no Sid file/label accessed; STOP for audit before Stage 2.
