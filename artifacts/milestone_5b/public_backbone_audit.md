# Public Event-B backbone audit

## Outcome

**`INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA`**

Five independent accepted Event-B cohorts reproduce in one deterministic corpus. Fukuoka is
blocked and contributes no records. No recognition model was trained.

## Corpus

- Event-B studies: 5
- Unique Event-B patients: 82 (total tier)
- Candidate-resolved patients: 45 across 4 studies (peptide-ranking tier)
- Patient-level-only patients: 37 (Nous-209; eligibility/abstention only, no peptide labels)
- Event-B-positive patients: 74 total; 37 candidate-resolved
- Patients with explicit tested negatives: 38
- Unique patient-candidate pairs: 1,072
- Primary candidate labels: 974
- Assay observations: 1,398
- Event-B assay observations: 1,012
- Primary positives: 272
- Primary tested negatives: 693
- Primary untested: 9
- Event-A observations: 304
- Epitope-spreading observations: 82
- Class-I primary observations: 207
- Class-II primary observations: 334
- Unknown/both-class primary observations: 433
- Review queue: 38 deterministic Hu source-ambiguity records; contradictions: 0

## Independence and dominance

Nous-209 is the largest patient contributor (37/82, 45.12%) but contributes no candidate labels.
Hu/NeoVax is the largest primary-label contributor (541/974, 55.54%), positive contributor
(128/272, 47.06%) and negative contributor (413/693, 59.60%). No study exceeds the registered 70%
primary-label dominance threshold. Repeated peptides account for 9.96% of primary labels.

## Label quality

All 693 accepted candidate negatives are explicit source calls; inferred negatives are zero. The
nine untested targets are the seven non-decomposable PDAC combined-pool members and two no-data
targets. Nous-209 contributes 37 patient-level responses and no peptide labels. Fukuoka contributes
none. Missing accepted-record provenance and accepted-label contradictions are both zero.

## Split feasibility

Patient, study, HLA, peptide-cluster, cancer-type and shared-antigen-group holdouts can each retain
positive patients and explicit negatives on both sides. Temporal holdout is not viable because the
public candidate tables do not expose enough calendar dates. A full study holdout is viable, but the
registered minimum still fails at 45 candidate-resolved versus 100 patients (the 82-patient headline
includes 37 patient-level-only Nous-209 participants that cannot train peptide ranking).
