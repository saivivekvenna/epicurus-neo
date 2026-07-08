# Big-Move Method Scan

The local BigMHC iterations show selector-style methods beating published
comparison columns on top-20 metrics, but not reaching a 50% held-out hit rate.
That changes the strategy: incremental reweighting is unlikely to be enough.

## What the Literature Suggests

- Current immunogenicity tools still have limited ability to predict cytotoxic
  activation against neoantigens. The VACCIMEL evaluation reports 94
  experimentally validated melanoma neopeptides and concludes that reviewed
  methods have limited activation-prediction ability:
  <https://www.explorationpub.com/Journals/ei/Article/100391>
- A high-value data direction is the reprocessed multi-assay dataset described
  by Ludwig Lausanne: three large immunogenicity screens from 131 cancer
  patients, with reported cross-dataset ranking gains:
  <https://www.ludwigcancerresearch.org/ludwig-link/april-2024/machine-learning-method-improves-neoantigen-selection-for-immunotherapy/>
- Uncertainty quantification is directly relevant because immunogenicity data is
  scarce and distribution shift is severe. ImmUQBench argues that UQ methods can
  improve reliability and robustness for immunogenicity prediction:
  <https://academic.oup.com/ooim/article/7/1/iqag003/8505599>

## Big Changes Worth Trying

1. **Risk-controlled nomination instead of forced top-20 filling**
   - Current precision-threshold experiment reached 35.3% held-out precision,
     not 50%.
   - Next version should be group-aware: calibrate thresholds per HLA family or
     allele cluster, not globally.
   - Output should explicitly split candidates into:
     - high-confidence core
     - lower-confidence filler
     - abstained cases where 50% is unsupported

2. **External validated-screening data ingestion**
   - Prioritize datasets with tested negatives, not only positive epitopes.
   - Target sources:
     - Ludwig/Immunity reprocessed screens if accessible
     - VACCIMEL supplementary table
     - CEDAR/IEDB cancer T-cell assay exports
     - NEPdb validated neoepitopes
   - The immediate engineering need is a normalizer that keeps assay context:
     pre-vaccine/post-vaccine, TIL/PBMC, cancer type, HLA, peptide length, and
     tested negative status.

3. **Conformal or UQ-gated ranker**
   - Treat “50% precision” as a risk-control target, not just a threshold.
   - Train multiple simple rankers over different feature families, estimate
     disagreement/uncertainty on validation, and nominate only when confidence
     clears a calibrated risk bound.
   - This matches the ImmUQBench direction and avoids overclaiming from a single
     overconfident model.

4. **TCR-aware mode when patient repertoire data exists**
   - The core missing variable is still patient-specific TCR availability.
   - If hackathon patient data includes TCR-seq, build a second-mode scorer:
     peptide/HLA presentation score plus TCR-pMHC retrieval/similarity to known
     reactive motifs.
   - If no TCR data exists, do not pretend this signal is available.

5. **Foreignness via embedding neighborhoods**
   - Replace crude sequence similarity with protein-language-model embedding
     neighborhoods against:
     - human self peptides
     - known immunogenic viral peptides
     - known cancer-reactive peptides
   - Use this as a retrieval feature, not as a standalone classifier.

## Current Verdict

The next most credible path to 50% is not a larger supervised stacker. It is:

```text
more validated screening data
+ explicit negative assays
+ uncertainty/risk-gated nomination
+ optional TCR-aware scoring if patient repertoire data exists
```

Until those signals are available, the honest product behavior is to produce a
ranked top-20 for the hackathon while separately reporting which subset clears a
calibrated high-confidence threshold.
