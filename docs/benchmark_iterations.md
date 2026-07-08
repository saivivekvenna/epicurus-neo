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
