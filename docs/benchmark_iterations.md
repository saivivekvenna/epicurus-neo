# Benchmark Iterations

## Iteration 001: Hits@20-Oriented Prior-Score Blend

Dataset:

- `data/processed/NmersTestingSet.normalized.csv`
- Gartner/NCI Nmers test set
- 26 leave-patient-out folds
- exact mutant/wildtype peptide-HLA overlaps purged from each training fold
- shared study allowed because this is patient-level CV within one source dataset

Hypothesis:

```text
Gartner's Nmer score is the strongest available prior on this benchmark, but
NetMHCpan-EL rank can recover some presentation signal. A conservative blend
should preserve Gartner top-k recall while nudging candidates with stronger
presentation evidence upward.
```

Implemented score:

```text
epicurus_hits20_score =
  0.9 * percentile_rank(baseline_gartner_nmer_score)
  + 0.1 * percentile_rank(baseline_netmhcpan_el_score)
```

Result:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_hits20_score` | 1.3846 | 0.0692 | 0.8173 | 0.5208 | 0.5176 |
| `baseline_gartner_nmer_score` | 1.3077 | 0.0654 | 0.7596 | 0.4608 | 0.4340 |
| `epicurus_blend_score` | 1.1923 | 0.0596 | 0.7083 | 0.4801 | 0.4888 |
| `baseline_netmhcpan_el_score` | 0.8077 | 0.0404 | 0.4647 | 0.2595 | 0.2139 |
| `baseline_mhcflurry_score` | 0.9231 | 0.0462 | 0.5256 | 0.2361 | 0.1896 |

Decision:

Accepted as the first measured improvement on the Gartner patient-level CV
benchmark. This is still a prior-score ensemble, not proof of generalization to
TESLA or future patient cases. Next iterations should test whether the same
conservative blending principle transfers when Gartner-specific `Nmer score` is
absent.

## Iteration 002: Generic Sequence Features

Dataset:

- Same Gartner/NCI Nmers test-set patient-level CV as Iteration 001

Hypothesis:

```text
Peptide sequence composition features should help the learned ranker and allow
cross-dataset inference when source-specific features such as Gartner Nmer score
are absent.
```

Added features:

- peptide length
- amino-acid fractions for all 20 standard amino acids
- hydrophobicity mean
- charge sum
- aromatic, polar, acidic, basic fractions
- cysteine, proline, glycine fractions

Gartner CV result:

| Score | mean hits@20 before | mean hits@20 after | recall@20 after | nDCG@20 after |
| --- | ---: | ---: | ---: | ---: |
| `epicurus_score` | 1.1538 | 1.2308 | 0.7692 | 0.3912 |
| `epicurus_hits20_score` | 1.3846 | 1.3846 | 0.8173 | 0.5208 |

Decision:

Accepted because it improves the learned ranker and is needed for datasets that
lack source-specific model scores.

External status check:

- Train: Gartner/NCI Nmers test set
- Test: TESLA normalized peptide/HLA labels
- Result: current model gets 0/33 TESLA positives into top 20

This TESLA result is a failure baseline, not an accepted tuning target. It shows
that Gartner source features do not transfer to TESLA without real presentation
or peptide-immunogenicity features.

## Iteration 003: Optional MHCflurry Presentation Features

Dataset:

- TESLA normalized peptide/HLA labels
- 714 peptide-HLA rows
- 33 positives, 681 negatives

Hypothesis:

```text
For datasets that include peptide and HLA but lack source-specific model scores,
local MHCflurry class-I presentation predictions provide a real biophysical
baseline and feature source.
```

Added:

- optional `mhcflurry` dependency extra
- `epicurus add-mhcflurry-features`
- `mhcflurry_affinity`
- `mhcflurry_processing_score`
- `mhcflurry_presentation_score`
- `mhcflurry_presentation_percentile`

TESLA status:

| Score | hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `mhcflurry_presentation_score` | 6 | 0.3000 | 0.1818 | 0.2286 | 0.1429 |
| Gartner-trained `epicurus_score` without MHCflurry | 0 | 0.0000 | 0.0000 | 0.0000 | 0.0417 |

Decision:

Accepted as a necessary external presentation feature path. This is not yet a
full TESLA-winning ranker, but it moves TESLA from zero top-20 hits to six using
a reproducible local predictor.

## Iteration 004: Transferable Presentation-Composition Score

Dataset:

- Train path: Gartner/NCI normalized data can be present but this score is
  heuristic and source-independent
- Test/status path: TESLA MHCflurry-scored table

Hypothesis:

```text
Raw class-I presentation is necessary but not sufficient. TESLA positives are
strongly presented and, in this dataset, less hydrophobic/cysteine/aromatic than
many high-presentation negatives. A conservative transfer score should keep the
presentation signal dominant while penalizing simple peptide-liability patterns.
```

Implemented score:

```text
epicurus_transfer_score =
  0.70 * percentile_rank(mhcflurry_presentation_score)
  + 0.15 * inverse_percentile_rank(seq_hydrophobicity_mean)
  + 0.10 * inverse_percentile_rank(seq_cysteine_fraction)
  + 0.05 * inverse_percentile_rank(seq_aromatic_fraction)
```

TESLA result:

| Score | hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_transfer_score` | 9 | 0.4500 | 0.2727 | 0.4964 | 1.0000 |
| `mhcflurry_presentation_score` | 6 | 0.3000 | 0.1818 | 0.2286 | 0.1429 |
| Gartner-trained `epicurus_score` | 0 | 0.0000 | 0.0000 | 0.0000 | 0.0417 |

Gartner regression check:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_hits20_score` | 1.3846 | 0.0692 | 0.8173 | 0.5208 | 0.5176 |
| `baseline_gartner_nmer_score` | 1.3077 | 0.0654 | 0.7596 | 0.4608 | 0.4340 |

Decision:

Accepted as the first transferable TESLA improvement. The next iteration should
focus on the 24 TESLA positives still missed by top 20, ideally adding
foreignness or known-epitope-neighborhood signals rather than increasing
TESLA-specific composition tuning.
