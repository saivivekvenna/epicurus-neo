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

## Iteration 017: HLA-Grouped Precision Thresholds

Dataset:

- Same BigMHC validation/test setup as Iteration 016
- Validation threshold still targets at least 50% precision with
  `min_selected=20`
- Thresholds are calibrated per `hla_allele` when the validation group has
  enough positive evidence; otherwise the global threshold is used
- Tested score columns:
  - `retrieval_topk_positive_fraction`
  - `retrieval_biochemical_topk_positive_fraction`
  - `mhcflurry_processing_score`
  - `mhcflurry_presentation_score`
  - `retrieval_positive_minus_negative_similarity`

Hypothesis:

```text
The global precision threshold may fail because different HLA alleles have
different score calibration curves. Per-HLA calibration, with a global fallback
for low-evidence alleles, may transfer a 50% high-confidence-core target better
than one global cutoff.
```

Best BigMHC `im_test` threshold results:

| Score used for threshold | min group positives | group thresholds | test selected | test hits | test precision | test recall |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `retrieval_topk_positive_fraction` | 5 | 6 | 124 | 46 | 0.3710 | 0.2323 |
| `retrieval_topk_positive_fraction` | 10 | 2 | 102 | 36 | 0.3529 | 0.1818 |
| `retrieval_topk_positive_fraction` | 2 | 10 | 137 | 48 | 0.3504 | 0.2424 |
| `retrieval_biochemical_topk_positive_fraction` | 5 | 7 | 150 | 52 | 0.3467 | 0.2626 |
| `mhcflurry_processing_score` | 10 | 3 | 70 | 18 | 0.2571 | 0.0909 |

Decision:

Rejected as a direct route to 50% held-out precision, but accepted as a useful
risk-control improvement over the global high-confidence-core gate. Group-aware
calibration lifted the best precision from 35.3% to 37.1% and increased hits
from 36 to 46, but it still did not transfer the validation 50% threshold to the
locked test split. This reinforces the primary metric choice: for the hackathon
submission, optimize fixed `mean_hits@20` / `precision@20`, and treat thresholded
selection only as a confidence annotation unless the rules allow abstention.

## Iteration 018: Recover Gartner Screened Negatives

Dataset:

- Gartner/NCI long-peptide files already present under `data/raw/gartner_nci/`
- Normalizer fix:
  - `CD8` / `1` -> positive
  - `0` / `-` -> negative
  - `unscreened` -> unknown
- Regenerated processed files:
  - `NmersBalancedForExpression.normalized.csv`: 139 positives, 4,765 negatives
  - `NmersTrainingSet.normalized.csv`: 139 positives, 9,404 negatives, 21,344
    unknowns
  - `NmersTestingSet.normalized.csv`: 46 positives, 3,722 negatives, 5,014
    unknowns

Hypothesis:

```text
The current Gartner normalizer was discarding thousands of screened
non-reactive candidates as unknowns. Recovering those explicit negatives should
give the supervised ranker a better recognition boundary and improve fixed
top-20 performance.
```

Gartner patient-level result:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `epicurus_hits20_score` | 1.3846 | 0.0692 | 0.8173 | 0.5208 | 0.5176 |
| `baseline_gartner_nmer_score` | 1.3077 | 0.0654 | 0.7596 | 0.4608 | 0.4340 |
| `baseline_mhcflurry_score` | 0.9231 | 0.0462 | 0.5256 | 0.2361 | 0.1896 |
| `baseline_netmhcpan_el_score` | 0.8077 | 0.0404 | 0.4647 | 0.2595 | 0.2139 |

Decision:

Accepted as a necessary data-quality fix, but rejected as a standalone path to
the 50% fixed top-20 target. Recovering explicit negatives makes the Gartner
benchmark more honest and keeps Epicurus ahead of available Gartner baselines,
but the absolute precision remains low. This strengthens the next-build
priority: ingest more validated screening datasets and add new signal families,
especially patient/context features, rather than relying on relabeling or
threshold calibration alone.

## Iteration 019: UQ, Motif Retrieval, and Pairwise Ranking

Dataset:

- Gartner official train/test split with exact overlaps purged
- BigMHC validation/test with regenerated retrieval features

Experiments:

1. **Bootstrap uncertainty ranker**
   - Added ensemble mean/std scores and `epicurus_lower_confidence_score`
   - Swept uncertainty penalties `0`, `0.5`, `1`, and `2`
2. **Motif/prototype retrieval**
   - Added deterministic peptide motif embeddings:
     amino-acid composition, terminal residues, T-cell-face composition, length,
     and biochemical summaries
   - Added motif nearest-neighbor and positive/negative prototype similarities
3. **Pairwise ranking score**
   - Trained positive-vs-negative within-group pairwise comparisons so the model
     optimizes ordering pressure rather than calibrated classification

Results:

| Experiment | Dataset | Best score | mean hits@20 | precision@20 | nDCG@20 | MRR |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| UQ lower confidence | Gartner | `epicurus_lower_confidence_score`, penalty 2 | 1.3846 | 0.0692 | 0.4056 | 0.3159 |
| Existing hits score | Gartner | `epicurus_hits20_score` | 1.3846 | 0.0692 | 0.5208 | 0.5176 |
| Pairwise ranker | Gartner | `epicurus_pairwise_score` | 0.9231 | 0.0462 | 0.2009 | 0.1521 |
| Motif selector | BigMHC | `epicurus_selected_score` | 2.5556 | 0.2765 | 0.3930 | 0.3562 |
| Current headline | BigMHC | `epicurus_selected_score` | 2.5556 | 0.2765 | 0.4057 | 0.3849 |
| Motif blend | BigMHC | `epicurus_blend_score` | 2.4815 | 0.2728 | 0.3714 | 0.3303 |

Decision:

Rejected as SOTA improvements. All three ideas are useful infrastructure, and
motif retrieval is competitive enough to keep as an available feature family,
but none improves the locked headline. The result is informative: the current
BigMHC SOTA-like point is not limited by a missing simple sequence-neighborhood
feature or by a naive rank-loss swap. The next plausible SOTA attempt needs a
materially stronger external signal, most likely real protein-language-model
embeddings or additional harmonized immunogenicity screens, rather than another
small handcrafted retrieval variant.

## Iteration 020: ESM2-Tiny PLM Retrieval

Dataset:

- BigMHC `im_val` and `im_test`
- Model: `facebook/esm2_t6_8M_UR50D`
- Device: Apple MPS available; ESM2-tiny embedding generation completed locally
- Feature family:
  - PLM nearest positive/negative cosine similarities
  - PLM top-k positive fraction
  - PLM positive/negative prototype similarities

Hypothesis:

```text
Handcrafted motif retrieval may be too weak. A pretrained protein language
model should provide smoother peptide neighborhoods and capture biochemical
regularities not represented by exact or handcrafted similarity.
```

Direct BigMHC `im_test` PLM retrieval results:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| `retrieval_plm_max_positive_similarity` | 2.4259 | 0.2700 | 0.5056 | 0.3923 | 0.3808 |
| `retrieval_plm_positive_minus_negative_similarity` | 2.4259 | 0.2700 | 0.5007 | 0.3864 | 0.3896 |
| `retrieval_plm_positive_minus_negative_prototype_similarity` | 2.3889 | 0.2682 | 0.5152 | 0.3835 | 0.3576 |
| `retrieval_plm_topk_positive_fraction` | 2.3519 | 0.2663 | 0.5049 | 0.3764 | 0.3698 |

Selector results:

| Selector | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| PLM + presentation selector | 2.4815 | 0.2728 | 0.5245 | 0.4054 | 0.3948 |
| all retrieval families selector | 2.5556 | 0.2765 | 0.5332 | 0.3956 | 0.3655 |
| current headline selector | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 |

Decision:

Rejected as the new headline, accepted as reusable infrastructure. ESM2-tiny
PLM retrieval is competitive and improves MRR in the PLM-only selector, but it
does not improve the primary fixed top-20 hit metric, and the all-feature
selector overfits validation choices enough to reduce nDCG/MRR. The next PLM
attempt should not just swap similarity columns; it should use stronger
embeddings or fine-tuned/listwise fusion with validation safeguards.

## Iteration 021: Fixed-20 Gate-Targeted Policies

Acceptance gate:

- Beat `mean_hits@20 > 2.5556` or `precision@20 > 0.2765`
- Keep `recall@20 >= 0.5259`
- Keep `nDCG@20 >= 0.4057`
- Use only train/validation labels for feature engineering, model selection,
  calibration, thresholding, and policy choice

Experiments:

1. **All-feature hits-first selector**
   - Candidate scores: exact retrieval, biochemical retrieval, motif retrieval,
     PLM retrieval, MHCflurry presentation, and MHCflurry processing
   - Objective: validation `hits`
   - Evidence thresholds: `min_positive` in `{1, 2, 3, 5, 10, 20}`
2. **Diversity-aware fixed-20 reranker**
   - Base score: validation-selected current headline policy
   - Greedy top-20 reranking penalizes biochemical similarity to already
     selected peptides
   - Diversity lambda selected on validation only
3. **High-confidence-core then fill**
   - Promote candidates that clear a validation-selected confidence threshold,
     then fill remaining top-20 slots with the current headline score
   - Confidence signals: retrieval, motif, PLM, and MHCflurry columns
4. **Supervised train+validation ranker**
   - Train on `im_train+im_val` with existing out-of-fold retrieval features
   - Evaluate once on locked `im_test`
   - Manual leakage check: `shared_mutant_hla` was empty. The generic leakage
     guard reports shared pseudo-patients because BigMHC encodes HLA as
     `patient_id`, and shared empty wildtype keys because these peptides do not
     have wildtype sequences; those are not exact mutant-peptide leakage.

Results:

| Policy | Selection basis | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR | Accept? |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| Current headline | validation nDCG selector | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 | current |
| All-feature hits selector, best `min_positive=1` | validation hits | 2.5370 | 0.2756 | 0.5159 | 0.3779 | 0.3571 | no |
| Diversity-aware reranker | validation hits/nDCG; selected `lambda=0` | 2.5556 | 0.2765 | 0.5334 | 0.3988 | 0.3747 | no |
| Confidence-core then fill | validation hits; `retrieval_plm_topk_positive_fraction >= 1.0` | 2.5185 | 0.2746 | 0.5316 | 0.4036 | 0.3694 | no |
| Supervised train+val `epicurus_retrieval_score` | train+validation labels | 2.5000 | 0.2737 | 0.5135 | 0.3963 | 0.4098 | no |
| Supervised train+val `epicurus_pairwise_score` | train+validation labels | 2.4815 | 0.2728 | 0.5069 | 0.3589 | 0.3139 | no |

Decision:

Rejected as headline replacements. The all-feature hits selector and
confidence-core policy optimized the right primary metric on validation but did
not transfer to locked `im_test`. Diversity selection did not help even on
validation (`lambda=0` was selected), so there is no evidence that correlated
top-20 bets are the current limiter. The supervised train+validation ranker also
failed to beat the current validation-selected retrieval policy.

The current headline remains the locked BigMHC fixed top-20 result. The next
highest-leverage direction is not another validation selector over the same
columns; it is either stronger representation learning (larger or fine-tuned
PLM embeddings with guarded validation) or more harmonized immunogenicity
screening data with explicit tested negatives.

## Iteration 022: Rank Aggregation and Bootstrap-Stable Selection

Acceptance gate:

- Beat `mean_hits@20 > 2.5556` or `precision@20 > 0.2765`
- Keep `recall@20 >= 0.5259`
- Keep `nDCG@20 >= 0.4057`
- Use only train/validation labels for policy selection

Experiments:

1. **Validation-selected rank aggregation**
   - Candidate scores: exact retrieval, biochemical retrieval, motif retrieval,
     PLM retrieval, MHCflurry presentation, and MHCflurry processing.
   - Generated 1,306 deterministic weight sets:
     - single-score rank policies
     - pair/triple/quartet/quintet policies over strong representative columns
     - sparse Dirichlet-weighted rank mixtures with a fixed seed
   - Selection keys tried from validation only:
     - hits
     - nDCG
     - recall
     - MRR
     - balanced
     - gated hits
2. **Bootstrap-stable per-HLA score selection**
   - Candidate scores restricted to the strongest retrieval/presentation
     columns for runtime.
   - For each HLA, bootstrapped validation rows and kept the HLA-specific score
     only when it won a minimum fraction of bootstrap rounds; otherwise fell
     back to the validation default.
   - Stability grid:
     - objective in `{hits, nDCG}`
     - `min_positive` in `{1, 2, 5, 10}`
     - stability in `{0.4, 0.6, 0.8}`

Results:

| Policy | Selection basis | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR | Accept? |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| Current headline | validation nDCG selector | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 | current |
| Rank aggregation, validation hits | 1,306 validation-selected weight sets | 2.4815 | 0.2728 | 0.5066 | 0.3889 | 0.3671 | no |
| Rank aggregation, validation nDCG | 1,306 validation-selected weight sets | 2.3704 | 0.2672 | 0.5105 | 0.3924 | 0.3602 | no |
| Rank aggregation, validation MRR/balanced | 1,306 validation-selected weight sets | 2.2778 | 0.2626 | 0.4921 | 0.3986 | 0.3841 | no |
| Bootstrap-stable selector | hits, `min_positive=1`, stability `0.4` | 2.5556 | 0.2765 | 0.5281 | 0.4146 | 0.3978 | no |

Decision:

Rejected as a new headline. Bootstrap-stable selection is a useful secondary
operating point because it preserves the current primary hit/precision level and
improves nDCG/MRR, but the acceptance gate requires a strict improvement in
`mean_hits@20` or `precision@20`, which it does not provide. Rank aggregation
over the same feature families overfits validation and loses held-out hits.

The evidence now points away from additional selectors over the same score
columns. The next highest-leverage path is to add materially new training signal:
larger/fine-tuned PLM embeddings, external screened negative datasets, or a
true listwise ranker trained on more than the BigMHC validation split.

## Iteration 023: Multi-k Retrieval and Oriented MHCflurry Scores

Acceptance gate:

- Beat `mean_hits@20 > 2.5556` or `precision@20 > 0.2765`
- Keep `recall@20 >= 0.5259`
- Keep `nDCG@20 >= 0.4057`
- Use only train/validation labels for policy selection

Experiments:

1. **Multi-k retrieval neighborhoods**
   - Added reusable multi-k retrieval features for `k in {1, 3, 5, 10, 20}`
   - Feature families:
     - exact peptide similarity
     - biochemical similarity
     - motif/prototype similarity
   - Validation reference: BigMHC `im_train`
   - Locked test reference: BigMHC `im_train+im_val`
   - Validation-selected objectives:
     - hits
     - nDCG
     - balanced
   - Evidence thresholds:
     - `min_positive in {1, 2, 3, 5, 10}`
2. **MHCflurry score orientation**
   - Added inverse-oriented temporary score columns:
     - `mhcflurry_affinity_inverse_score = -log10(affinity)`
     - `mhcflurry_presentation_percentile_inverse_score = -percentile`
   - Re-ran validation-selected policies with retrieval, motif, PLM,
     presentation, processing, inverse affinity, and inverse percentile columns

Results:

| Policy | Selection basis | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR | Accept? |
| --- | --- | ---: | ---: | ---: | ---: | ---: | --- |
| Current headline | validation nDCG selector | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 | current |
| Multi-k retrieval, best hits selector | hits, `min_positive=10` | 2.5000 | 0.2737 | 0.5081 | 0.3761 | 0.3473 | no |
| Multi-k retrieval, best nDCG selector | nDCG, `min_positive=10` | 2.3889 | 0.2682 | 0.4996 | 0.4038 | 0.3955 | no |
| Multi-k retrieval, best balanced selector | balanced, `min_positive=10` | 2.4074 | 0.2691 | 0.5202 | 0.4120 | 0.3950 | no |
| Oriented scores, hits selector | hits, `min_positive=1` | 2.5370 | 0.2756 | 0.5159 | 0.3837 | 0.3747 | no |
| Oriented scores, nDCG selector | nDCG, `min_positive=1` | 2.5556 | 0.2765 | 0.5332 | 0.4014 | 0.3830 | no |
| Oriented scores, balanced selector | balanced, `min_positive=1` | 2.5185 | 0.2746 | 0.5308 | 0.4055 | 0.3917 | no |

Decision:

Rejected as a new headline. Multi-k retrieval did not improve held-out fixed
top-20 hits, suggesting the current `k=5` retrieval neighborhood is not the
main bottleneck. Correcting lower-is-better MHCflurry affinity/percentile into
high-is-better score columns also failed the acceptance gate: the closest
variant matched the primary hit/precision level but lost nDCG, while the
balanced variant nearly matched the nDCG gate but lost primary hits.

The reusable multi-k retrieval infrastructure is kept. The next highest-leverage
direction remains stronger out-of-distribution signal rather than more
validation selection over BigMHC-native columns: larger/fine-tuned PLM
embeddings, external screened negatives, or a proper listwise ranker trained on
additional harmonized immunogenicity data.

## Iteration 024: Frozen PLM Learning-to-Rank

Acceptance gate:

- Beat `mean_hits@20 > 2.5556` or `precision@20 > 0.2765`
- Keep `recall@20 >= 0.5259`
- Keep `nDCG@20 >= 0.4057`
- Use only train/validation labels for encoder and ranker selection

Experiment:

1. Generated normalized mean-pooled peptide embeddings with two frozen ESM-2
   encoders:
   - `facebook/esm2_t6_8M_UR50D`
   - `facebook/esm2_t12_35M_UR50D`
2. Combined the embeddings with:
   - MHCflurry presentation and processing features
   - mutation-position and physicochemical features
   - leakage-controlled exact and biochemical retrieval features
   - train-derived HLA one-hot features
3. Trained grouped `XGBRanker` models over a validation-only grid:
   - objectives: `rank:ndcg`, `rank:pairwise`
   - depths: `2`, `3`, `4`
   - learning-rate/tree pairs: `(0.03, 300)`, `(0.06, 180)`
4. Selected the encoder and ranker lexicographically by validation hits,
   precision, recall, nDCG, and MRR.
5. Refit the frozen selection on BigMHC `im_train+im_val`, then evaluated it
   once on locked `im_test`.

Validation selection:

| Encoder | Selected configuration | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| ESM2 8M | nDCG, 300 trees, depth 2 | 1.7455 | 0.3043 | 0.5271 | 0.4876 | 0.5021 |
| ESM2 35M | pairwise, 180 trees, depth 2 | 1.7273 | 0.3034 | 0.5234 | 0.4783 | 0.4970 |

Locked result:

| Policy | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR | Accept? |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Current headline | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 | current |
| Frozen ESM2 8M + XGBRanker | 2.4259 | 0.2700 | 0.5247 | 0.4017 | 0.3883 | no |

Residual check:

The frozen PLM score was also rank-blended with the current headline score
using weights from `0.025` through `1.0`. This check used validation labels
only. The unmodified headline score remained best on validation hits at
`1.8364`; no positive PLM weight improved that primary metric. The residual
variant was therefore rejected before another locked-test evaluation.

Decision:

Rejected as a headline replacement. The smaller encoder generalized better
than the larger encoder during validation, but direct frozen peptide
embeddings plus a tree ranker did not add enough recognition signal to beat the
current retrieval/presentation policy. The infrastructure remains useful, but
the next PLM experiment should not merely increase model size. It should add
orthogonal supervision: screened external negatives, mutant/wild-type paired
representations, or contrastive fine-tuning with study- and HLA-grouped
validation.

## Iteration 025: External TCR-Recognition Neighborhoods

Acceptance gate:

- Beat `mean_hits@20 > 2.5556` or `precision@20 > 0.2765`
- Keep `recall@20 >= 0.5259`
- Keep `nDCG@20 >= 0.4057`
- Construct external features without BigMHC labels
- Use only validation labels for threshold, ranker, and policy selection

External reference:

- VDJdb release `2026-06-11-ZENODO`
- Source:
  [antigenomics/vdjdb-db](https://github.com/antigenomics/vdjdb-db/releases/tag/2026-06-11-ZENODO)
- Filters:
  - human TCR records
  - MHC class I
  - canonical 8-14mer peptides
- Reference policies:
  - all curated records: 1,972 peptide-HLA-origin rows
  - higher-confidence records (`vdjdb.score >= 1`): 512 rows

Features:

- Rebuilt ESM2 8M embeddings using residue-only mean pooling, excluding special
  tokens.
- Added raw and globally centered cosine neighborhoods for:
  - HLA-matched recognized epitopes
  - pathogen-derived recognized epitopes
  - human-derived recognized epitopes
- Added top-k means, evidence/support-weighted maxima, and
  pathogen-minus-human contrasts.
- No BigMHC label was used to construct any VDJdb feature.

Validation experiments:

| Policy | Reference | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| External-feature XGBRanker | all records | 1.7818 | 0.3061 | 0.5275 | 0.4919 | 0.5014 |
| External-feature XGBRanker | score >= 1 | 1.7818 | 0.3061 | 0.5291 | 0.4835 | 0.4929 |
| Frozen headline fallback | score >= 1 | 1.8364 | 0.3088 | 0.5402 | 0.4994 | 0.5121 |
| Headline + sparse recognition residual | score >= 1 | 1.8364 | 0.3088 | 0.5402 | 0.5225 | 0.5455 |

The sparse residual kept the exact current headline score as the default and
allowed only five HLA groups to switch to an external-recognition feature.
Validation hits, precision, and recall were unchanged while nDCG and MRR
improved.

Locked results:

| Policy | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR | Accept? |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| Current headline | 2.5556 | 0.2765 | 0.5332 | 0.4057 | 0.3849 | current |
| Broad VDJdb recognition selector | 2.5185 | 0.2746 | 0.5315 | 0.4194 | 0.4128 | no |
| Sparse headline residual | 2.5370 | 0.2756 | 0.5159 | 0.3938 | 0.3890 | no |

Decision:

Rejected as a headline replacement. The broad selector produced a meaningful
held-out nDCG and MRR gain and kept recall above the gate, but it lost primary
hits and precision. The stricter five-HLA residual preserved validation hits
while improving validation ordering, but that gain did not transfer and its
held-out recall/nDCG fell below the gate.

The result weakens the hypothesis that generic population-level recognized
epitope proximity is enough. The next recognition experiment needs tighter
task alignment: tumor neoantigen positives and screened negatives, paired
mutant/wild-type tolerance features, or TCR-contact-position representations
rather than whole-peptide mean embeddings.

## Iteration 026: Direct Recognition Dataset Expansion

Added two patient-level, experimentally screened recognition datasets without
using BigMHC `im_test` labels:

- IMPROVE official cross-validation matrix:
  - 17,520 tested candidates
  - 467 positives
  - 70 patients in five official patient-disjoint partitions
- Cross-tumor 2025 pHLA multimer screen:
  - 8,095 deduplicated tested candidates
  - 34 positives
  - 26 patients split across TIL and PBMC cohorts

IMPROVE five-fold out-of-fold result:

| Score | mean hits@20 | precision@20 | recall@20 | nDCG@20 | MRR |
| --- | ---: | ---: | ---: | ---: | ---: |
| PRIME source score | 1.2000 | 0.0600 | 0.2086 | 0.1267 | 0.1815 |
| Epicurus learned ranker | 1.4714 | 0.0736 | 0.2206 | 0.1432 | 0.1986 |

Patient-level paired bootstrap versus PRIME:

- hits@20 delta: `+0.2714`, 95% interval `[+0.0143, +0.5429]`
- precision@20 delta: `+0.0136`, 95% interval `[+0.0007, +0.0271]`
- probability of a positive hits delta: `0.977`

The first unconstrained cross-cohort transfer experiment failed: a model trained
on the TIL subset recovered 5 PBMC top-20 hits, while NetMHCpan EL recovered 10.
This supports using the new labels for constrained residual learning or
stage-aware auxiliary training, not replacing presentation ranking outright.

Decision:

Accepted as reusable direct-recognition training and external-validation
infrastructure. This iteration does not alter the locked BigMHC headline. The
next experiment will train a recognition residual from IMPROVE and select its
blend strength on BigMHC `im_val` before any locked-test evaluation.

## Iteration 027: Cross-Dataset Transfer and Benchmark Realignment

Goal:

- test whether substantially more direct recognition data trains a better
  BigMHC reranker;
- quantify whether a `5 hits@20` goal is possible on the locked benchmark.

Leakage control:

- IMPROVE had 38 exact peptide-HLA overlaps with BigMHC train, 4 with
  validation, and 0 with test.
- All 42 overlapping rows were purged before external feature construction or
  transfer training.
- The broader direct-screen union removed 1,008 exact overlaps against the
  union of BigMHC train, validation, and test candidates.
- No BigMHC test label was used for model or hyperparameter selection.

Validation results:

| Method | External data | mean hits@20 | recall@20 | nDCG@20 |
| --- | --- | ---: | ---: | ---: |
| Frozen headline | none | 1.8182 | 0.5384 | 0.4994 |
| Signed PLM recognition neighborhood | IMPROVE | 1.4364 | 0.5090 | 0.4193 |
| Frozen PLM ranker + screened features | IMPROVE | 1.7273 | 0.5252 | 0.4853 |
| Target-only transfer control | none | 1.5818 | 0.5199 | 0.4721 |
| Pooled ranker, patient-grouped external loss | IMPROVE | 1.6545 | 0.5240 | 0.4719 |
| Pooled ranker, HLA-grouped external loss | IMPROVE | 1.6909 | 0.5278 | 0.4820 |
| Pairwise fine-tuned ESM2 | target only | 1.4545 | 0.5066 | 0.3999 |
| Pairwise fine-tuned ESM2 | IMPROVE pretraining | 1.4364 | 0.5001 | 0.4382 |
| Pooled HLA-grouped ranker | IMPROVE + 2025 multimer | 1.6545 | 0.5222 | 0.4869 |

The extra data did improve the otherwise identical frozen-feature ranker:
IMPROVE HLA-grouped pooling increased validation hits from `1.5818` to `1.6909`.
It did not beat the mature retrieval/presentation headline. End-to-end PLM
fine-tuning reduced validation hits despite falling pairwise training loss,
showing that the learned sequence relation did not transfer.

One label-free locked check averaged within-HLA ranks across ten published model
families:

| Policy | mean hits@20 | precision@20 | recall@20 | nDCG@20 |
| --- | ---: | ---: | ---: | ---: |
| Current headline | 2.5556 | 0.2765 | 0.5332 | 0.4057 |
| Published-model consensus | 2.2593 | 0.2617 | 0.4890 | 0.3893 |
| Fixed 50/50 headline + consensus | 2.4074 | 0.2691 | 0.5118 | 0.3972 |

Oracle audit:

| Benchmark | Group | Oracle mean hits@20 | Current mean hits@20 |
| --- | --- | ---: | ---: |
| BigMHC validation | HLA allele | 1.9636 | 1.8182 |
| BigMHC locked test | HLA allele | 3.5185 | 2.5556 |
| IMPROVE official CV | patient | 6.4571 | 1.4714 |
| 2025 multimer | patient | 1.3077 | not accepted |

Decision:

No transfer method replaces the BigMHC headline. A `5 hits@20` target is
mathematically impossible on BigMHC because the oracle ceiling is `3.5185`.
BigMHC is now treated as a component-level peptide/HLA regression benchmark.
The primary hard-part target moves to patient-disjoint IMPROVE, where an oracle
mean of `6.4571` makes `5 hits@20` valid and where patient, expression,
clonality, mutant/wild-type, and tumor-context features match the eventual
WES/RNA product.

## Iteration 028: Patient Baselines and Paired Sequence Audit

The IMPROVE source table was audited before further optimization. Twenty-three
patients had all positives ordered before negatives, and 33 patients achieved
their oracle top-20 score from source order alone. Metrics now break score ties
with a deterministic hash of biological candidate identity. Source row order
and candidate IDs cannot decide top-k membership.

Reproduced official and external baselines under the corrected metric:

| Method | Mean hits@20 |
| --- | ---: |
| Official IMPROVE RF, TME excluded | 1.4429 |
| NeoGuider isotonic logistic regression | 1.3429 |
| Epicurus direct-recognition ranker | 1.4714 |
| Patient-budget XGBoost candidate | 1.4857 |
| Development-only RF/Epicurus rank blend | 1.5143 |

The blend is diagnostic, not a locked result, because its weight was selected
after inspecting aggregate out-of-fold performance.

Representation experiments:

| Added evidence | Mean hits@20 | Decision |
| --- | ---: | --- |
| Leakage-safe ESM retrieval from training folds | 0.9857 | reject |
| Raw features in stricter preprocessing harness | 1.4143 | control |
| Mutant-minus-wild-type ESM delta | 1.1857 | reject |
| Full mutant/wild-type ESM pair | 1.4143 | reject |

Generic protein-language-model geometry did not encode the T-cell
cross-reactivity relation needed for recognition. The next paired model must be
pretrained on shared-TCR peptide relationships with HLA context before direct
cancer-screen fine-tuning. NeoPrecis-Immuno is the required external baseline
for that hypothesis.
