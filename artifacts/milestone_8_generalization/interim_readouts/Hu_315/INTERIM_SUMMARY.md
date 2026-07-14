# Hu_315 immutable interim result

Status: **development-only calibration readout**. The frozen portfolios were
created before this outcome join. This patient may not be used to tune or select
policy, and the final-held-out cohort remains sealed.

## Top-20 result

| Frozen arm | Unique recognized mutations in top 20 | Unique mutations represented | Duplicate slots | Delta vs PRIME |
|---|---:|---:|---:|---:|
| Epicurus mutation-cap-1 | **8** | 20 | 0 | **+7** |
| PRIME + Epicurus rank fusion, cap-1 | **8** | 20 | 0 | **+7** |
| Shipped Epicurus product | **5** | 11 | 9 | **+4** |
| Epicurus evidence-lane portfolio | **5** | 20 | 0 | **+4** |
| Genuine PRIME | **1** | 19 | 1 | baseline |
| Genuine PRIME, mutation-cap-1 | **1** | 20 | 0 | 0 |
| Epicurus route-level plain | **0** | 1 | 19 | -1 |

The strongest Epicurus-only preregistered portfolio therefore placed **8 unique
recognized mutations in 20 unique mutation slots (40% precision by the locked
mutation-level endpoint)**. Genuine PRIME placed 1.

## End-to-end reachability

- Experimentally recognized mutations recorded for Hu_315: **33**.
- Reconstructed by matched tumor/normal WES: **18 / 33 (54.5%)**.
- Passed the frozen evidence gate: **18 / 18 reachable**.
- Losslessly generated into the peptide/HLA universe: **18 / 18 reachable**.
- Epicurus mutation-cap-1 selected: **8 / 18 reachable (44.4%)**.
- Shipped Epicurus selected: **5 / 18 reachable (27.8%)**.
- Genuine PRIME selected: **1 / 18 reachable (5.6%)**.

For this patient, the dominant remaining loss is therefore **upstream somatic
reachability** (15 recognized mutations absent from the reconstructed callset),
not the validity/presentation gate. Among reachable mutations, mutation-level
portfolio diversification is the major reason Epicurus outperforms PRIME.

## Interpretation boundary

This is one calibration patient and is not a generalization claim. It is a valid
frozen-pipeline result showing that the Epicurus-only mutation-cap portfolio
substantially outperformed genuine PRIME on the identical reachable universe.
The headline full-pipeline comparison against nextNEOpi still requires sealed
multi-patient evaluation.

