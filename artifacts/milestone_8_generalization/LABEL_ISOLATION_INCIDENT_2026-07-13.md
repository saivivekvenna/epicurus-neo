# Calibration label-isolation incident — 2026-07-13

## What happened

At approximately 13:09 MDT, a metadata-acquisition search used a broad path:

```text
rg -n 'Hu_182|Hu_268|Hu_254|Hu_343' artifacts/milestone_8_generalization/INPUT_CROSSWALK.csv data/raw/miller_ipv ...
```

Because `data/raw/miller_ipv` also contains `miller_recognition_labels.csv`, the
command opened that file before the six calibration reconstructions were frozen.
Terminal output visibly included recognition-label rows for `Hu_182` and
`Hu_254`. The query also named `Hu_268` and `Hu_343`; even though their rows were
not visible in the captured/truncated output, all four are conservatively treated
as early-unsealed.

No final-held-out patient identifier was part of the query. The final cohort
(`Hu_333`, `Hu_159`, `Hu_344`, `Hu_048`, `Hu_293`, `Hu_250`) remains sealed.
`Hu_277` and `Hu_315` calibration labels were also not queried.

## Impact

The original statement that all six calibration labels would remain unopened
until every calibration portfolio was frozen is no longer true. Results on the
calibration cohort must therefore be described as development evidence with an
early-unseal protocol deviation, never as label-isolated validation.

The decisive six-patient final evaluation remains potentially valid because its
IDs and outcomes were not accessed. That claim is conditional on keeping the
pre-incident semantic pipeline fixed and opening final labels only through the
committed fail-closed evaluator after all final portfolios are frozen.

## Containment

1. No exposed gene, mutation, peptide, outcome, or patient-specific information
   may be used in reconstruction, gating, scoring, arm design, thresholds, or
   portfolio selection.
2. The registered arms, simplicity order, comparator, metrics, and policy
   selection objective remain those committed before the incident.
3. The semantic pipeline is pinned to pre-incident `HEAD`
   `f23672314924c31b3e87c026d0a3ac801672e4e1`. No semantic changes are permitted
   during calibration reconstruction. Operational recovery, if unavoidable,
   must preserve byte-identical scientific commands and be documented separately.
4. All six calibration patients will still be reconstructed and frozen without
   consulting labels. The evaluator will select the universal arm mechanically
   from the pre-registered set; the calibration result will carry this incident
   disclosure.
5. No final-held-out ID may be searched in any outcome-bearing path. Final
   evaluation remains a single once-only unseal after the universal lock and all
   six final freezes exist.

## Pre-incident semantic file pins

| File | SHA-256 |
|---|---|
| `src/benchmark/miller_product_freeze.py` | `526c17d1927a13ed0954e8e988913ae72797fb96b317b7d546068980e45da124` |
| `src/benchmark/miller_universe_core.py` | `2e02bf1dab53a8166df873b88529d19a70c0c9217539953e37e4b17cf56bd517` |
| `src/epicurus_neo/product.py` | `d17400dd2c90ce4ab779cacc3e0f8d8099077fee2d533331b5b635cf4d4c15f2` |
| `src/epicurus_neo/gates.py` | `75d7ddb2aeff1368ec2e56ed3b81338e26b3934f80548785560417d082ff0842` |
| `src/benchmark/universal_portfolio.py` | `4be47d22b1767619197f9386cdf6fd3ecc91e954bdbab233451e14e3baf7e577` |
| `scripts/miller_patient_reconstruct.py` | `1b7f54bb8225a03282d326b554ffad93f77f5afb1b86886fd9294f34b95174bb` |
| `scripts/miller_hu287_somatic.sh` | `4afc9635e9b2cc9ba3d639e5de82108fd67b8573c73ff78e08c2b5d97364d609` |
| `scripts/miller_hu287_hla.sh` | `564b61a567d3b3d8f8a360a8f2004ac17ee2034108fd6a3b7255ae8884d669ec` |
| `scripts/miller_hu287_rna.sh` | `20d7ef6f249031e5ad6b39faebf5992b51b7c74b08e241633198348b018076a5` |

This document records the failure rather than retroactively redefining the
original freeze as successful.
