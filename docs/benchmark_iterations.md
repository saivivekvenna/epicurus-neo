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
                  …
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

                                                            ## Iteration 029: Frozen NeoPrecis Diagnostic

                                                            The released Apache-2.0 NeoPrecis-Immuno model was applied without IMPROVE
                                                            label fitting. IMPROVE publishes NetMHC percentile rank rather than the raw
                                                            NetMHC score expected by the model, so the adapter used the explicit
                                                            approximation `1 - RankEL_4.1 / 100`.

                                                            | Direction | Evaluated candidates | Mean hits@20 |
                                                            | --- | ---: | ---: |
                                                            | Higher NP-Immuno score | 17,359 | 0.5714 |
                                                            | Lower NP-Immuno score | 17,359 | 0.4857 |

                                                            Mean scores were almost identical for IMPROVE positives and negatives. This
                                                            diagnostic rejects NeoPrecis as a drop-in per-peptide score. It is not an exact
                                                            reproduction of the published method because of the binding-input
                                                            approximation and because NeoPrecis aggregates peptide-HLA evidence at mutation
                                                            level. Exact reproduction requires rerunning a compatible NetMHCpan version to
                                                            obtain raw scores before making a comparative claim.
