# Benchmark Contract

This document defines what must be true before we claim Epicurus Neo beats an
existing neoantigen ranking method.

## Task

For each patient or case, rank candidate neoantigens so experimentally reactive
peptides appear as high as possible, especially within the final submitted
top-k list.

The hackathon deliverable is a peptide ranking. Therefore the benchmark is a
grouped ranking benchmark, where each group is a patient/case/study unit.

## Labels

Rows must be assigned one of three label states:

| State | Meaning | Training use |
| --- | --- | --- |
| `positive` | Experimentally tested and T-cell reactive | positive label |
| `negative` | Experimentally tested and non-reactive | negative label |
| `unknown` | Not experimentally tested or ambiguous | excluded from supervised labels |

Unknown rows must not be silently converted to negatives.

Optional label confidence weights should reflect assay quality, for example:

- high: patient autologous T-cell response, ELISpot/multimer/functional assay
- medium: compatible ex vivo response with partial metadata
- low: weak, ambiguous, pooled, or non-patient-specific assay

## Splits

Headline claims require at least one of:

- leave-patient-out
- leave-study-out
- leave-peptide-cluster-out
- official benchmark train/test split

Random peptide-row splits are allowed only for smoke tests and debugging.

## Leakage Rules

A split is invalid if train and test share any of the following exact keys:

- normalized mutant peptide + HLA allele
- normalized wildtype peptide + HLA allele
- source patient/case ID
- source study ID when doing study-level holdout

For strict external benchmarks, also flag near-duplicate peptide clusters. The
first implementation uses exact keys; later versions should add sequence-cluster
and embedding-neighbor leakage checks.

## Metrics

Primary:

- `hits@20`: number of positives in top 20 per group
- `precision@20`: positive fraction in top 20 per group
- `recall@20`: fraction of group positives recovered in top 20
- `ndcg@20`: ranking quality with positive labels as relevance
- `mrr`: reciprocal rank of first positive

Secondary:

- AUROC
- average precision
- expected calibration error
- Brier score

## Baselines

Minimum baselines before claiming improvement:

- binding-affinity-only rank
- presentation-score-only rank
- pVAC-style heuristic rank when features are available
- existing model scores available in the dataset, such as BigMHC, PRIME,
  DeepImmuno, or NeoRanking outputs

## Fair Goal

The first fair quantitative goal is:

```text
On at least one locked patient/study-level benchmark, Epicurus Neo improves
mean hits@20 and mean precision@20 over the strongest available baseline without
degrading recall@20 by more than 10% relative.
```

The stronger goal is:

```text
Improve hits@20 on TESLA and Gartner/NCI official held-out sets after all model
choices are frozen.
```

