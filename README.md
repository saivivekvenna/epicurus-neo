# Epicurus Neo

Epicurus Neo is a competition-grade neoantigen ranking system for the prescribed
hackathon task: given patient tumor data, rank candidate neoantigens so the final
submitted peptide set maximizes experimentally validated T-cell responses.

This repository starts with the hard part: a fair benchmark and training harness.
The first target is not a UI or report generator. It is a leakage-controlled,
top-k-optimized model development loop that can prove whether a method beats
existing neoantigen ranking baselines.

## Quantitative Goal

Primary goal:

```text
Improve validated hits@20 and precision@20 over strong public baselines on
locked held-out neoantigen immunogenicity benchmarks.
```

Primary locked benchmarks:

- TESLA / Wells et al. 2020
- Gartner/NCI official held-out sets
- External held-out datasets such as BigMHC `im_test`, ITSNdb-style sets, or
  study-level CEDAR/NEPdb/dbPepNeo holdouts after strict deduplication

Primary metrics:

- `hits@20`
- `precision@20`
- `recall@20`
- `ndcg@20`
- `mrr`
- calibration error for reported probabilities

AUC is tracked only as a secondary diagnostic because the deliverable is a
ranked peptide list, not a generic binary classifier.

## Method Thesis

The core bet is not another single immunogenicity predictor. Epicurus Neo trains
a top-k ranker that learns when existing predictors fail.

Initial model stack:

1. Presentation gate
2. Immunogenicity ranker
3. False-positive / dud detector
4. Probability calibrator
5. Final top-20 selector

The highest-priority feature families are:

- existing predictor scores: NetMHCpan/MHCflurry/BigMHC/PRIME/DeepImmuno/pVAC-style
- expression and clonality
- mutant-vs-wildtype binding and embedding deltas
- self-similarity / foreignness
- known reactive and known failed epitope neighborhoods
- anchor-vs-TCR-facing mutation position
- model disagreement and out-of-domain indicators

## Non-Negotiables

- No random row splits for headline claims.
- Unknown/unassayed candidates are not negatives.
- Locked test sets are not used for feature search, threshold tuning, or model
  selection.
- Report top-k patient/group-level metrics, not just pooled AUC.
- Every benchmark run records data versions, split definitions, feature set, and
  model config.

