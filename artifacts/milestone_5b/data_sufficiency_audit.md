# Public Event-B data sufficiency audit

**Verdict: INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA**

No recognition model was fitted. The registered gate is evaluated on candidate-resolved (peptide-level) patients only; patient-level-only cohorts are reported but never counted toward peptide-ranking readiness. CD4/CD8 counts below use source-resolved candidate class, not inferred cellular phenotype.

## Evidence tiers (peptide-ranking sample vs. total)

- `total_event_b_patient_n`: 82
- `candidate_resolved_patient_n`: 45
- `patient_level_only_patient_n`: 37
- `candidate_resolved_positive_patient_n`: 37
- `candidate_resolved_study_n`: 4
- `candidate_level_primary_label_n`: 974

## Global counts

- `total_studies`: 5
- `event_b_studies`: 5
- `publications`: 12
- `unique_patients`: 82
- `event_b_patients`: 82
- `event_b_positive_patients`: 74
- `candidate_resolved_patient_n`: 45
- `candidate_resolved_positive_patient_n`: 37
- `candidate_resolved_study_n`: 4
- `patient_level_only_patient_n`: 37
- `patients_with_explicit_tested_negatives`: 38
- `vaccine_components`: 1072
- `unique_patient_candidate_pairs`: 1072
- `primary_candidate_labels`: 974
- `assay_observations`: 1398
- `event_b_assay_observations`: 1012
- `event_b_positives`: 272
- `event_b_tested_negatives`: 693
- `event_b_untested_candidates`: 9
- `event_a_observations`: 304
- `epitope_spreading_observations`: 82
- `class_i_observations`: 207
- `class_ii_observations`: 334
- `unknown_class_observations`: 433
- `cd4_responses`: 110
- `cd8_responses`: 18
- `shared_antigen_observations`: 72
- `personalized_antigen_observations`: 902
- `review_queue_size`: 38

## Registered minimum

- `thresholds`: {'candidate_resolved_patients': 100, 'candidate_resolved_studies': 2, 'candidate_resolved_positive_patients': 30}
- `met`: False
- `study_holdout_feasible`: True
- `no_overwhelming_primary_label_dominance`: True
- `explicit_tested_negatives_available`: True
- `class_i_and_class_ii_coverage`: True
- `label_comparability_for_primary_candidate_set`: True

## Dominance

- `largest_study_patient_fraction`: 0.4512
- `largest_study_observation_fraction`: 0.4456
- `largest_study_primary_label_fraction`: 0.5554
- `largest_study_positive_fraction`: 0.4706
- `largest_study_negative_fraction`: 0.5960

## Per-study status

- `braun_rcc_2025`: ACCEPTED; patients=9; primary=129 (+61/-68/untested=0); patient-level-only=0
- `fukuoka_dc`: BLOCKED_SOURCE_UNAVAILABLE; patients=0; primary=0 (+0/-0/untested=0); patient-level-only=0
- `hu_neovax_2021`: ACCEPTED; patients=8; primary=541 (+128/-413/untested=0); patient-level-only=0
- `mkras_vax_2026`: ACCEPTED; patients=12; primary=72 (+60/-12/untested=0); patient-level-only=0
- `nous_209_2025`: ACCEPTED; patients=37; primary=0 (+0/-0/untested=0); patient-level-only=37
- `pdac_neovax_2023`: ACCEPTED; patients=16; primary=232 (+23/-200/untested=9); patient-level-only=1

## Split feasibility

- `PATIENT_HOLDOUT`: viable
- `STUDY_HOLDOUT`: viable
- `HLA_HOLDOUT`: viable
- `PEPTIDE_CLUSTER_HOLDOUT`: viable
- `CANCER_TYPE_HOLDOUT`: viable
- `TEMPORAL_HOLDOUT`: not viable
- `SHARED_ANTIGEN_GROUP_HOLDOUT`: viable
