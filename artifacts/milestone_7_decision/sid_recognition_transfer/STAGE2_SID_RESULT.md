# Stage 2 — frozen gate applied ONCE to Sid

_code commit `19f505b284f6`; frozen SHA `e44d3d1ac66a272e`; payload SHA `580339d7b6bcc8da`; PASS (config SHA + payload SHA + stage2_must_not_refit); no refit._

_EXPLORATORY confirmation (Sid previously inspected), n=1 patient / 3 positives; not pristine external validation. One shot; no tuning/rerun after results._


**Accounting:** {'generated': 137, 'unrepresentable_documented': 10, 'total_accounted': 147, 'mutations_scored': 137} (147 = generated + documented-unrepresentable).

**Input hashes:** {'frozen_config': '19a7cc38d0566c7dd080ac411adb8c0731c266701eb1f813d9c0972e004c3d61', 'scored_candidates': 'ee6d69635327d45d17faa04b965d8ddd51cbae5a92632dc081a7a7709355cda2', 'variant_vafs_long': '265224cae8bb7c43c85a65a7d4512195edbd7bbd9741d1cbe935f6027acbe93c', 'per_variant': 'ec1173f96836a34da6f9f9713658f41989bf82bacd54fd545d5f09ea37de5762'}


## Tie-aware mutation-level hits@20 (labels joined post-freeze)

| arm | nominal | guaranteed | hit IDs |
|---|--:|--:|---|
| genuine PRIME | 2/3 | 2/3 | ['DYNC1H1-chr14-101980529', 'MAP2-chr2-209694772'] |
| frozen Epicurus v0.1 | 1/3 | 1/3 | ['DYNC1H1-chr14-101980529'] |
| **frozen gate** | **2/3** | **2/3** | ['DYNC1H1-chr14-101980529', 'MAP2-chr2-209694772'] |

## Per-positive (frozen gate)

- `ASPM-chr1-197102716`: {'score_rank_interval': [41, 41], 'rna_af': 0.5, 'protect_guaranteed': False, 'reserve_guaranteed': False, 'nominally_selected': False}
- `DYNC1H1-chr14-101980529`: {'score_rank_interval': [3, 3], 'rna_af': 0.4275, 'protect_guaranteed': True, 'reserve_guaranteed': False, 'nominally_selected': True}
- `MAP2-chr2-209694772`: {'score_rank_interval': [9, 9], 'rna_af': 0.0, 'protect_guaranteed': True, 'reserve_guaranteed': False, 'nominally_selected': True}

PRIME rank intervals of positives: {'ASPM-chr1-197102716': [41, 41], 'DYNC1H1-chr14-101980529': [3, 3], 'MAP2-chr2-209694772': [9, 10]}
