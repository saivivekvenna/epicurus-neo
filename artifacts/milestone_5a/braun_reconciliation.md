# Braun RCC 2025 — Event-B ingestion reconciliation

Source: Braun et al., Nature 2025 (PMC11903305, NCT02950766), CC BY-NC-ND. Immunogenicity
is recomputed from raw IFN-gamma ELISpot replicates using the paper's stated rule
(P<0.05 two-sided t-test AND mean spot count >=3x no-stim DMSO control). No counts are
hard-coded; no negative is inferred from omission.

## Expected (research report) vs accepted

| Quantity | Expected | Accepted | Match |
|---|---|---|---|
| Patients | 9 | 9 | True |
| Positives | 61 | 61 | True |
| Tested-negatives | 68 | 68 | True |

## Independent recomputation vs the paper's summary (sheet 2e)

- Source In Vitro peptides: 130 across 9 patients (unscorable rows: 0).
- Recomputed by rule, split by mutation type: `{'driver': {'immunogenic': 11, 'non_immunogenic': 6}, 'passenger': {'immunogenic': 50, 'non_immunogenic': 62}, 'unclassified': {'immunogenic': 1, 'non_immunogenic': 0}}`
- Paper sheet-2e targets: `{'driver': {'immunogenic': 11, 'non_immunogenic': 6}, 'passenger': {'immunogenic': 50, 'non_immunogenic': 62}}`
- Driver/passenger splits reconcile exactly: **True**

## Why accepted positives are 61, not 62

The rule scores 62 peptides immunogenic across all 130 In Vitro rows. One immunogenic
peptide (AMACR|p.Y41N, patient 104) has a blank `Mutation_type` and is excluded from the
paper's driver/passenger summary (Fig. 2e / Supplementary Table 2). Rather than silently
inflate the accepted count to 62, it is routed to the review queue; the 129 accepted assays
(61 positive + 68 tested-negative) match the paper's reported figures.

## Review queue, event typing, and completeness

- Review queue: 1 record(s) by code `{'UNCLASSIFIED_MUTATION': 1}`.
- Accepted event counts: `{'EVENT_B_VACCINE_INDUCED_RESPONSE': 129}` (Event-B only; no Event-A relabelled).
- De-novo basis: paper states 'no pre-existing immune responses were detected'; ex-vivo week-0 pool baselines are background (below the ELISpot positivity floor).
- MHC class breakdown (candidates): `{'UNKNOWN': 130}` (long SLPs left UNKNOWN; class not resolved per peptide by the assay).
- HLA-resolved fraction (candidates): 0.00 (per-peptide HLA is a prediction, not an assay restriction; not stored as such).
- Candidate identity completeness: `{'gene': 1.0, 'protein_change': 1.0, 'genomic_variant': 1.0, 'transcript': 0.0, 'wildtype_peptide': 0.0}`.

## Model-readiness (this slice alone)

```json
{'sufficient_for_recognition_model_development': False, 'event_b_patient_n': 9, 'event_b_study_n': 1, 'event_b_positive_patient_n': 9, 'registered_minimums': {'event_b_patients': 100, 'event_b_studies': 2, 'positive_patients': 30}, 'decision': 'EVENT_B_VERTICAL_SLICE_VALIDATED_NOT_YET_SUFFICIENT_FOR_GENERAL_MODEL'}
```
