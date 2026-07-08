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

## Iteration 005: BigMHC Immunogenicity External Benchmark

Dataset:

- `data/raw/bigmhc/datasets.zip`
- `im_train.csv`: 6,185 rows, 1,407 positives, 4,778 negatives
- `im_val.csv`: 688 rows, 173 positives, 515 negatives
- `im_test.csv`: 937 rows, 198 positives, 739 negatives
- Evaluation group: `hla_allele`
- Metric: mean hits@20 / precision@20 / recall@20 across 54 HLA groups

Purpose:

```text
Add an independent hard-part benchmark where the test set includes published
scores from BigMHC and other presentation/immunogenicity tools. This gives us a
fair target beyond Gartner and TESLA: beat the best existing score column on the
same rows without training or tuning on im_test labels.
```

Published-score baselines on `im_test.csv`:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `TransPHLA` | 2.4259 | 0.2700 | 0.5006 | 0.3906 | 0.3732 |
| `MHCnuggets-2.4.0` | 2.4074 | 0.2691 | 0.5259 | 0.3927 | 0.3785 |
| `HLAthena_Scores` | 2.3889 | 0.2682 | 0.4946 | 0.3818 | 0.3699 |
| `netmhcpan_41_score` | 2.3519 | 0.2663 | 0.5100 | 0.4001 | 0.4015 |
| `mhcflurry_20_score` | 2.3519 | 0.2663 | 0.4981 | 0.4018 | 0.3881 |
| `bigmhc_el_score` | 2.3519 | 0.2663 | 0.4917 | 0.3916 | 0.3800 |
| `prime_20_score` | 2.3333 | 0.2654 | 0.5221 | 0.4019 | 0.3787 |
| `bigmhc_elim_score` | 2.3333 | 0.2654 | 0.5133 | 0.3915 | 0.3819 |
| `bigmhc_im_score` | 2.3148 | 0.2645 | 0.5041 | 0.3856 | 0.3778 |

Epicurus sequence-only check:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_score`, trained on `im_train + im_val` only | 2.2963 | 0.2635 | 0.5080 | 0.3679 | 0.3426 |

Decision:

Accepted as a fair quantitative target, not as an Epicurus win. The target for
the next modeling iteration is to exceed `TransPHLA` on mean hits@20 (>2.4259)
and exceed the strongest published recall@20 (`MHCnuggets-2.4.0`, 0.5259)
without using `im_test` labels for model selection. Sequence-only learning is
not enough; next work should add train-side presentation predictions, epitope
retrieval density, and foreignness/self-similarity features.

## Iteration 006: Positive-Neighborhood Retrieval Score

Dataset:

- Reference for validation: BigMHC `im_train`
- Validation queries: BigMHC `im_val`
- Reference for locked test: BigMHC `im_train + im_val`
- Locked test queries: BigMHC `im_test`
- Evaluation group: `hla_allele`

Hypothesis:

```text
T-cell recognition may transfer locally in peptide/HLA space. A candidate whose
mutant peptide is close to a known immunogenic peptide for the same HLA allele
should rank higher, even when presentation scores are similar.
```

Implemented features:

- `retrieval_max_positive_similarity`
- `retrieval_max_negative_similarity`
- `retrieval_positive_minus_negative_similarity`
- `retrieval_topk_positive_similarity_mean`
- `retrieval_topk_negative_similarity_mean`
- `retrieval_topk_positive_fraction`

Accepted score:

```text
epicurus_retrieval_score = retrieval_max_positive_similarity
```

BigMHC `im_test` result:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_retrieval_score` | 2.5000 | 0.2737 | 0.5135 | 0.3963 | 0.4098 |
| `TransPHLA` | 2.4259 | 0.2700 | 0.5006 | 0.3906 | 0.3732 |
| `MHCnuggets-2.4.0` | 2.4074 | 0.2691 | 0.5259 | 0.3927 | 0.3785 |
| `epicurus_transfer_score` | 2.4259 | 0.2700 | 0.5032 | 0.3997 | 0.3692 |
| `epicurus_score` with retrieval features | 2.4074 | 0.2691 | 0.5101 | 0.3633 | 0.3405 |

Decision:

Accepted as the first BigMHC hits@20 win over published score columns. This is
not yet a complete win over every metric: `MHCnuggets-2.4.0` still has better
recall@20. The next iteration should preserve the retrieval hits gain while
recovering recall, likely through a score that balances positive-neighborhood
similarity with a recall-oriented negative-neighborhood or presentation term.

## Iteration 007: Biochemical Retrieval Similarity

Dataset:

- Same BigMHC validation/test setup as Iteration 006
- Validation queries: `im_val` against `im_train`
- Locked test queries: `im_test` against `im_train + im_val`

Hypothesis:

```text
Exact residue identity may be too brittle for T-cell recognition. Conservative
substitutions should count partially, so biochemical peptide-neighborhood
features may improve recall without losing the retrieval hits@20 gain.
```

Added:

- `retrieval_biochemical_max_positive_similarity`
- `retrieval_biochemical_max_negative_similarity`
- `retrieval_biochemical_positive_minus_negative_similarity`
- `retrieval_biochemical_topk_positive_similarity_mean`
- `retrieval_biochemical_topk_negative_similarity_mean`
- `retrieval_biochemical_topk_positive_fraction`

Validation result:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `retrieval_positive_minus_negative_similarity` | 1.7455 | 0.3043 | 0.5307 | 0.4420 | 0.4286 |
| `retrieval_max_positive_similarity` | 1.7273 | 0.3034 | 0.5257 | 0.4444 | 0.4394 |
| `retrieval_biochemical_topk_positive_similarity_mean` | 1.6545 | 0.2997 | 0.5228 | 0.4499 | 0.4658 |
| `retrieval_biochemical_max_positive_similarity` | 1.6182 | 0.2979 | 0.5170 | 0.4261 | 0.4271 |

Locked BigMHC `im_test` check:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `retrieval_max_positive_similarity` | 2.5000 | 0.2737 | 0.5135 | 0.3963 | 0.4098 |
| `retrieval_biochemical_positive_minus_negative_similarity` | 2.5000 | 0.2737 | 0.5066 | 0.3899 | 0.3557 |
| `retrieval_biochemical_topk_positive_similarity_mean` | 2.3704 | 0.2672 | 0.5147 | 0.4039 | 0.4015 |
| `retrieval_biochemical_max_positive_similarity` | 2.3148 | 0.2645 | 0.5052 | 0.4058 | 0.3963 |

Decision:

Rejected as the default ranking signal. The biochemical features are retained as
model inputs because they expose a different neighborhood view, but validation
does not support replacing exact positive-neighborhood similarity. The next
iteration should focus on score selection/calibration rather than making peptide
similarity softer.

## Iteration 008: Validation-Selected Score Family

Dataset:

- Validation: BigMHC `im_val`
- Locked test: BigMHC `im_test`
- Selection unit: `hla_allele`
- Candidate score families: exact retrieval, biochemical retrieval, MHCflurry
  presentation/processing, and `epicurus_transfer_score`

Hypothesis:

```text
Different HLA alleles may favor different evidence types. Choose the best score
family per HLA allele on validation labels, with a global validation winner as
fallback, then apply that selection to the locked test set.
```

Validation-selected default:

```text
retrieval_positive_minus_negative_similarity
```

BigMHC `im_test` result:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_selected_score` | 2.5370 | 0.2756 | 0.5159 | 0.3870 | 0.3710 |
| `epicurus_retrieval_score` | 2.5000 | 0.2737 | 0.5135 | 0.3963 | 0.4098 |
| `TransPHLA` | 2.4259 | 0.2700 | 0.5006 | 0.3906 | 0.3732 |
| `MHCnuggets-2.4.0` | 2.4074 | 0.2691 | 0.5259 | 0.3927 | 0.3785 |

Decision:

Accepted as the current best BigMHC hits@20/precision@20 result. It improves
the hits target further using validation-only selection, but still does not beat
the best recall@20. The next hard problem remains recall recovery without
sacrificing top-20 precision.

## Iteration 009: nDCG-Selected Score Family

Dataset:

- Same validation/test setup as Iteration 008
- Selection objective: validation `nDCG@20`
- Candidate score families: exact retrieval, biochemical retrieval, MHCflurry
  presentation/processing, and `epicurus_transfer_score`

Hypothesis:

```text
Hits-first validation selection improves precision but under-recovers recall.
nDCG-first selection should reward putting positives early while still allowing
more recall-oriented score families to win on HLA groups where retrieval alone
is too narrow.
```

Validation-selected default:

```text
mhcflurry_presentation_score
```

BigMHC `im_test` result:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_selected_score` (`objective=ndcg`) | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 |
| `epicurus_selected_score` (`objective=mrr`) | 2.5185 | 0.2746 | 0.5308 | 0.4088 | 0.3880 |
| `epicurus_selected_score` (`objective=hits`) | 2.5370 | 0.2756 | 0.5159 | 0.3870 | 0.3710 |
| `epicurus_retrieval_score` | 2.5000 | 0.2737 | 0.5135 | 0.3963 | 0.4098 |
| `TransPHLA` | 2.4259 | 0.2700 | 0.5006 | 0.3906 | 0.3732 |
| `MHCnuggets-2.4.0` | 2.4074 | 0.2691 | 0.5259 | 0.3927 | 0.3785 |
| `mhcflurry_20_score` | 2.3519 | 0.2663 | 0.4981 | 0.4018 | 0.3881 |
| `netmhcpan_41_score` | 2.3519 | 0.2663 | 0.5100 | 0.4001 | 0.4015 |

Decision:

Accepted as the new BigMHC headline result. The nDCG objective gives the best
hits@20/precision@20 operating point, while the MRR objective gives the best
observed Epicurus recall@n and nDCG@n operating point. Both beat the published
BigMHC comparison columns on mean hits@20, precision@20, recall@20, and nDCG@20
using validation-only score-family selection. Neither beats the best MRR
observed among all available scores, so the next iteration should focus on
early-first ranking without giving up the top-20 gains.

## Iteration 010: Balanced Validation Objective

Dataset:

- Same BigMHC validation/test setup as Iterations 008-009
- Candidate score families: exact retrieval, biochemical retrieval, MHCflurry
  presentation/processing, and `epicurus_transfer_score`
- Selection objective: validation
  `precision@20 + recall@20 + nDCG@20 + MRR`

Hypothesis:

```text
A balanced validation objective may recover the strong early-ranking behavior
of MRR-oriented selectors without giving up the top-20 hit and recall gains
from nDCG-oriented selection.
```

Validation-selected default:

```text
mhcflurry_processing_score
```

BigMHC `im_test` result:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_selected_score` (`objective=balanced`) | 2.5185 | 0.2746 | 0.5308 | 0.4088 | 0.3880 |
| `epicurus_selected_score` (`objective=ndcg`) | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 |
| `epicurus_retrieval_score` | 2.5000 | 0.2737 | 0.5135 | 0.3963 | 0.4098 |
| `TransPHLA` | 2.4259 | 0.2700 | 0.5006 | 0.3906 | 0.3732 |
| `MHCnuggets-2.4.0` | 2.4074 | 0.2691 | 0.5259 | 0.3927 | 0.3785 |

Decision:

Rejected as the headline operating point. The balanced objective is useful as a
CLI-supported selector because it exposes a validation-only multi-metric policy,
but on BigMHC it converges to the same held-out operating point as the MRR
objective: strong recall/nDCG, not enough hits/precision to replace the
`objective=ndcg` headline.

## Iteration 011: Validation-Selected Rank Blends

Dataset:

- Same BigMHC validation/test setup as Iterations 008-010
- Candidate score families: exact retrieval, biochemical retrieval, MHCflurry
  presentation/processing, and `epicurus_transfer_score`
- Blend search: HLA-local percentile-rank fusion over single score families and
  pairwise blends with weights 0.25/0.75, 0.50/0.50, and 0.75/0.25
- Selection still uses validation only; test labels are used only for final
  reporting

Hypothesis:

```text
Single-family selection leaves complementary signal on the table. Rank-normalized
pairwise blends may recover early positives from presentation-like signals while
keeping retrieval's immunogenic-neighborhood signal in the top 20.
```

Best BigMHC `im_test` blend results:

| Score | objective | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_blend_score` | nDCG | 2.3889 | 0.2682 | 0.5134 | 0.4082 | 0.4160 |
| `epicurus_blend_score` | MRR | 2.4259 | 0.2700 | 0.5159 | 0.4098 | 0.4006 |
| `epicurus_blend_score` | balanced | 2.4444 | 0.2709 | 0.5171 | 0.4006 | 0.3830 |
| `epicurus_blend_score` | hits | 2.4259 | 0.2700 | 0.4999 | 0.3804 | 0.3775 |
| `epicurus_selected_score` | nDCG | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 |
| `epicurus_retrieval_score` | single score | 2.5000 | 0.2737 | 0.5135 | 0.3963 | 0.4098 |
| `netmhcpan_41_score` | published | 2.3519 | 0.2663 | 0.5100 | 0.4001 | 0.4015 |

Decision:

Accepted as a secondary operating point, not the top-20 headline. The
validation-selected nDCG blend is the first Epicurus score to beat the strongest
published MRR baseline and the retrieval-only MRR while also improving nDCG.
It does not replace the nDCG score-selector headline because it sacrifices
hits@20 and precision@20, which remain the primary hackathon submission
constraint.

## Iteration 012: Dense Rank-Blend Grid

Dataset:

- Same BigMHC validation/test setup as Iteration 011
- Same candidate score families
- Blend search expanded from 0.25/0.50/0.75 pairwise weights to a denser
  validation grid: 0.10, 0.20, ..., 0.90

Hypothesis:

```text
The coarse pairwise blend grid may be missing a better validation-transfer
tradeoff. A denser grid could recover the headline selector's top-20 strength
while preserving some of the rank-blend MRR gain.
```

Best BigMHC `im_test` dense-grid results:

| Score | objective | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_blend_score` | nDCG dense | 2.5000 | 0.2737 | 0.5189 | 0.4004 | 0.3924 |
| `epicurus_blend_score` | balanced dense | 2.4815 | 0.2728 | 0.5177 | 0.4050 | 0.3891 |
| `epicurus_blend_score` | MRR dense | 2.4815 | 0.2728 | 0.5177 | 0.4050 | 0.3891 |
| `epicurus_blend_score` | hits dense | 2.4630 | 0.2719 | 0.5004 | 0.3848 | 0.3846 |
| `epicurus_selected_score` | nDCG selector | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 |
| `epicurus_blend_score` | coarse nDCG blend | 2.3889 | 0.2682 | 0.5134 | 0.4082 | 0.4160 |

Decision:

Rejected. The denser grid improves over several published top-20 baselines, but
it does not beat the nDCG selector on the primary top-20 metrics and it loses
the coarse nDCG blend's MRR advantage. This is a useful overfitting warning:
more validation search capacity is not automatically better transfer.

## Iteration 013: Guarded Rank-Blend Selector

Dataset:

- Same BigMHC validation/test setup as Iterations 011-012
- Same candidate score families
- Guarded policy:
  - choose the best single score family on validation as the baseline
  - search single-score and pairwise rank blends
  - accept a blend only if it does not fall below the baseline on a validation
    guard metric, then optimize a secondary objective

Hypothesis:

```text
The previous rank blends improved early-hit ranking but sometimes gave up too
many top-20 hits. A validation guard on hits or recall should preserve the
current selector's primary top-20 behavior while allowing MRR-oriented blends
only where validation says they are not costly.
```

BigMHC `im_test` guarded results:

| Score | policy | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_blend_score` | baseline nDCG, optimize MRR, guard hits | 2.4259 | 0.2700 | 0.5159 | 0.4098 | 0.4006 |
| `epicurus_blend_score` | baseline nDCG, optimize MRR, guard recall | 2.4259 | 0.2700 | 0.5159 | 0.4098 | 0.4006 |
| `epicurus_blend_score` | baseline nDCG, optimize balanced, guard hits | 2.4444 | 0.2709 | 0.5171 | 0.4006 | 0.3830 |
| `epicurus_blend_score` | dense grid, baseline nDCG, optimize MRR, guard hits | 2.4815 | 0.2728 | 0.5177 | 0.4050 | 0.3891 |
| `epicurus_selected_score` | nDCG selector | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 |
| `epicurus_blend_score` | coarse nDCG blend | 2.3889 | 0.2682 | 0.5134 | 0.4082 | 0.4160 |

Decision:

Rejected as a headline replacement. The guard makes the search policy more
auditable and prevents intentionally selecting validation hit-losing blends, but
it still does not transfer into a better held-out top-20 operating point than
the validation-selected nDCG score selector. Keep the guarded selector as
infrastructure for future datasets, not as the current BigMHC headline.

## Iteration 014: Crossfit Retrieval Stacker

Dataset:

- BigMHC `im_train`, `im_val`, and locked `im_test`
- Training retrieval features:
  - `im_train`: 5-fold out-of-fold retrieval features for validation tuning
  - `im_train + im_val`: 5-fold out-of-fold retrieval features for final test
    training
- Validation/test retrieval features:
  - `im_val` retrieves only from `im_train`
  - `im_test` retrieves only from `im_train + im_val`
- Stacker: existing `HistGradientBoostingClassifier` ranker over presentation,
  sequence, and retrieval features
- Artifact columns (`rand`, `retrieval_fold`) are excluded from model features

Hypothesis:

```text
The selector wins are driven by retrieval-neighborhood signal. A learned stacker
trained on leakage-safe out-of-fold retrieval features may learn when retrieval,
presentation, and biochemical-neighborhood features should dominate instead of
choosing one score family per HLA.
```

Validation read:

On `im_val`, the corrected stacker matched retrieval on hits@20/precision@20 and
improved nDCG over retrieval, but still trailed presentation-style scores on MRR.
That was enough signal to run the locked test check, but not enough to promote it
without test evidence.

BigMHC `im_test` result:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_score` (crossfit retrieval stacker) | 2.3519 | 0.2663 | 0.5102 | 0.3712 | 0.3364 |
| `epicurus_retrieval_score` | 2.5000 | 0.2737 | 0.5135 | 0.3963 | 0.4098 |
| `epicurus_selected_score` | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 |
| `TransPHLA` | 2.4259 | 0.2700 | 0.5006 | 0.3906 | 0.3732 |
| `MHCnuggets-2.4.0` | 2.4074 | 0.2691 | 0.5259 | 0.3927 | 0.3785 |
| `netmhcpan_41_score` | 2.3519 | 0.2663 | 0.5100 | 0.4001 | 0.4015 |

Decision:

Rejected as a ranking method, accepted as infrastructure. The crossfit retrieval
feature path gives the learned stacker a leakage-safer training view of the
retrieval signal, but the gradient-boosted probability model does not transfer
to the locked BigMHC test set. The current evidence continues to favor
validation-selected score policies over direct supervised stacking.

## Iteration 015: Minimum-Positive Evidence Selector

Dataset:

- Same BigMHC validation/test setup as Iterations 008-014
- Same candidate score families as the nDCG selector:
  - exact retrieval
  - biochemical retrieval
  - MHCflurry presentation/processing
  - `epicurus_transfer_score`
- Selection objective: validation `nDCG@20`
- Sweep: require at least 2, 3, 5, or 10 validation positives before trusting an
  HLA-specific score-family choice; otherwise use the global validation winner

Hypothesis:

```text
Per-HLA validation selection may overfit low-evidence HLA groups. Raising the
minimum-positive threshold should reduce noisy HLA-specific choices and may
improve held-out transfer, especially nDCG/MRR, even if it gives up some
score-family specialization.
```

BigMHC `im_test` result:

| Score | min positives | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `epicurus_selected_score` | 1 | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 |
| `epicurus_selected_score` | 2 | 2.5000 | 0.2737 | 0.5278 | 0.4104 | 0.3947 |
| `epicurus_selected_score` | 3 | 2.4630 | 0.2719 | 0.5031 | 0.4057 | 0.3986 |
| `epicurus_selected_score` | 5 | 2.4630 | 0.2719 | 0.5031 | 0.4057 | 0.3986 |
| `epicurus_selected_score` | 10 | 2.4074 | 0.2691 | 0.5003 | 0.4069 | 0.4109 |
| `epicurus_retrieval_score` | single score | 2.5000 | 0.2737 | 0.5135 | 0.3963 | 0.4098 |
| `TransPHLA` | published | 2.4259 | 0.2700 | 0.5006 | 0.3906 | 0.3732 |
| `MHCnuggets-2.4.0` | published | 2.4074 | 0.2691 | 0.5259 | 0.3927 | 0.3785 |
| `netmhcpan_41_score` | published | 2.3519 | 0.2663 | 0.5100 | 0.4001 | 0.4015 |

Decision:

Accepted as a secondary operating point, not the headline. `min_positive=2`
improves held-out nDCG over the current headline and keeps recall above the
published comparison columns, but it gives up hits@20 and precision@20. The
current `min_positive=1` nDCG selector remains the top-20 headline because
hits/precision/recall are the primary hackathon submission constraints. The
sweep does show that evidence-gating HLA-specific choices is a useful control
against validation overfitting.

## Iteration 016: Precision-Target Abstention

Dataset:

- Same BigMHC validation/test setup as Iterations 008-015
- Validation threshold calibrated to target at least 50% precision
- `min_selected=20` on validation, so the threshold must support at least a
  vaccine-sized candidate set before being accepted
- Applied threshold to the locked test split without using test labels during
  calibration

Hypothesis:

```text
If the model cannot honestly make 20 high-confidence calls, it should abstain
instead of filling all slots. A validation-calibrated threshold may identify a
smaller candidate core with precision closer to the desired 50% hit-rate target.
```

BigMHC `im_test` threshold results:

| Score used for threshold | validation selected | validation precision | test selected | test hits | test precision | test recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `retrieval_topk_positive_fraction` | 59 | 0.5932 | 102 | 36 | 0.3529 | 0.1818 |
| `mhcflurry_processing_score` | 76 | 0.5000 | 70 | 18 | 0.2571 | 0.0909 |
| `mhcflurry_presentation_score` | 142 | 0.5000 | 95 | 23 | 0.2421 | 0.1162 |
| `retrieval_positive_minus_negative_similarity` | 90 | 0.5889 | 100 | 18 | 0.1800 | 0.0909 |
| `retrieval_max_positive_similarity` | 122 | 0.5410 | 60 | 8 | 0.1333 | 0.0404 |
| `retrieval_biochemical_max_positive_similarity` | 107 | 0.5047 | 57 | 7 | 0.1228 | 0.0354 |
| `retrieval_biochemical_topk_positive_similarity_mean` | 108 | 0.5000 | 57 | 6 | 0.1053 | 0.0303 |

Decision:

Rejected as a route to a 50% held-out hit rate. The abstaining threshold can
raise precision above the current top-20 average (`retrieval_topk_positive_fraction`
reaches 35.3% on test), but validation-calibrated 50% thresholds do not transfer
to the locked BigMHC split. This is strong evidence that a 50% true-positive
rate cannot be claimed from score thresholding alone; reaching that target will
require either materially better patient-specific signal, stronger external
training data, or a much smaller/stricter nomination budget.
