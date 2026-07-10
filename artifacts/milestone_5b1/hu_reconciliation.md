# Hu 2021 melanoma NeoVax — Event-B ingestion reconciliation

Source: Hu et al., Nat Med 2021 (PMC8273876, NCT01970358), author manuscript. Per-peptide
CD8 (class-I minimal epitopes) and CD4 (class-II assay peptides) IFN-gamma ELISpot calls are
ingested as reported (Hu scored a response positive at >=2.5x DMSO). They are reconciled
against Ott 2017's independently published totals (PMC5577644 (Ott 2017; per-peptide screen for patients 1-6)).

## Vaccine-peptide recognition vs Ott 2017 (patients 1-6, week 16)

| Channel | Positive neoantigens (observed) | Ott reported | Total neoantigens | Reconciles |
|---|---|---|---|---|
| CD8 | 15 | 15 | 97 (Ott 97) | True |
| CD4 | 56 | 58 | 96 (Ott 97) | within tol (delta 2): True |

The CD8 positive-neoantigen count matches Ott exactly. CD4 differs by a small margin that is
not closed by any other assay condition in the table (minigene/tumor add nothing), consistent
with Ott counting at the immunizing-peptide rather than neoantigen granularity; it is reported,
not tuned away. Per-channel totals differ slightly from Ott's 97 because CD8 and CD4 cover
different neoantigen subsets; the positive count is the anchor.

## Accepted Event-B corpus

- Patients: 8 (`['1', '2', '3', '4', '5', '6', '11', '12']` — Ott's 1-6 plus new 11-12; one consolidated
  source, so no cross-study patient double-counting).
- Accepted assays: 623 — 128 POSITIVE / 495 TESTED_NEGATIVE / 0 UNTESTED.
- Event counts: `{'EPITOPE_SPREADING': 82, 'EVENT_B_VACCINE_INDUCED_RESPONSE': 541}`.
- MHC class (candidates): `{'CLASS_I': 304, 'CLASS_II': 334}` (CD8 class I, CD4 class II).

## Epitope spreading kept separate

- 82 epitope-spreading assays; all non-vaccine: **True**; any labelled a vaccine-candidate positive: **False**.
- week-16 epitope-spreading responses in Datasets 11a-c are uniformly non-reactive; the paper's spreading positives arise at later / post-checkpoint timepoints not ingested here. Kept as EPITOPE_SPREADING, never vaccine-candidate recognition.

## Reliability is not flattened

```json
{
  "FUNCTIONAL_T_CELL_ASSAY": {
    "assay_directness": [
      0.9
    ],
    "candidate_specificity": [
      1.0
    ],
    "n": 82,
    "temporal_clarity": [
      0.6
    ],
    "vaccine_relevance": [
      0.0
    ]
  },
  "LONGITUDINAL_PERSISTENCE": {
    "assay_directness": [
      0.7
    ],
    "candidate_specificity": [
      0.8
    ],
    "n": 128,
    "temporal_clarity": [
      0.8
    ],
    "vaccine_relevance": [
      1.0
    ]
  },
  "VACCINE_EVENT_B": {
    "assay_directness": [
      0.6,
      0.9
    ],
    "candidate_specificity": [
      0.6,
      1.0
    ],
    "n": 541,
    "temporal_clarity": [
      0.7
    ],
    "vaccine_relevance": [
      1.0
    ]
  }
}
```

De-novo basis: Ott 2017 and Hu 2021 both report no pre-vaccination neoantigen reactivity; this table carries no week-0 column, so de-novo is author-asserted (not baseline-verified) and recognition_evidence temporal_clarity is set lower than a baseline-verified de-novo claim.

## Model-readiness (this slice alone)

```json
{'sufficient_for_recognition_model_development': False, 'event_b_patient_n': 8, 'event_b_study_n': 1, 'event_b_positive_patient_n': 8, 'registered_minimums': {'event_b_patients': 100, 'event_b_studies': 2, 'positive_patients': 30}, 'decision': 'EVENT_B_VERTICAL_SLICE_VALIDATED_NOT_YET_SUFFICIENT_FOR_GENERAL_MODEL'}
```
