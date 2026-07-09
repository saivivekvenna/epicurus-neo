# Ranking Engine Plan

## Exact Goal

Epicurus Neo is the ranking component inside a larger pipeline:

```text
tumor WES + matched-normal WES + tumor RNA
    -> somatic variants, HLA, expression, clonality, peptide candidates
    -> Epicurus Neo
    -> ranked top-20 neoantigen portfolio
```

The primary research metric is patient-level mean `hits@20`. The stretch target
is `5.0` on the IMPROVE official patient-disjoint folds. Five hits in 20 is 25%
precision, not 50%. A 50% precision target would require 10 hits per patient and
is impossible on IMPROVE, whose oracle mean is 6.46.

Current credible development results:

| Policy | Mean hits@20 | Oracle capture |
| --- | ---: | ---: |
| Published IMPROVE RF | 1.44 | 22.3% |
| Epicurus patient model | 1.47 | 22.8% |
| Patient-budget XGBoost candidate | 1.49 | 23.0% |
| Development-only RF/Epicurus blend | 1.51 | 23.5% |
| Oracle | 6.46 | 100% |

The `5.0` target requires 77.4% oracle capture and a 3.3x improvement over the
current credible result. It remains the target, but it must not be presented as
an expected result until a locked external evaluation supports it.

## Component Contract

Epicurus does not consume FASTQ directly. An upstream adapter must emit one row
per patient, mutation, transcript, mutant peptide, wild-type peptide, and HLA
allele. Required evidence families are:

- patient, sample, mutation, transcript, gene, peptide, and HLA identifiers
- DNA depth, tumor VAF, caller confidence, purity, copy number, and clonality
- RNA depth, mutant RNA reads, RNA VAF, transcript expression, and expression
  uncertainty
- mutant and wild-type sequences, mutation position, consequence, and phase
- binding, processing, presentation, and stability predictions
- HLA expression, HLA loss, and antigen-processing machinery state
- provenance and missingness for every derived value

Optional RNA-derived patient evidence:

- reconstructed TCR beta CDR3 sequences, clone abundance, and reconstruction
  confidence
- repertoire depth, clonality, diversity, and tumor infiltration
- immune activation, exhaustion, and suppression signatures

The output contains one calibrated score per candidate, stage-specific
probabilities, uncertainty, out-of-distribution flags, evidence provenance, and
the selected top-20 portfolio.

## Biological Factorization

The model should preserve the causal stages:

```text
P(useful target)
  = P(translated)
  * P(presented | peptide, HLA)
  * P(recognized | peptide, HLA, patient)
  * P(tumor coverage | clonality, HLA viability)
```

A learned residual may correct this product, but no single assay label should be
used as if it supervised every stage. For example, an eluted-ligand negative is
not a T-cell non-response, and an untested candidate is not a negative.

## Two Recognition Routes

### Universal Route

This route is always available from WES/RNA-derived candidate evidence:

1. HLA-aware mutant/wild-type cross-reactivity representation.
2. Direct cancer-screen recognition score.
3. Self-proteome and validated-epitope neighborhood evidence.
4. Presentation, abundance, clonality, and tumor viability gates.
5. Patient-level ranking model focused on the top-20 boundary.

Generic ESM embeddings are not sufficient. The paired representation must first
be trained on TCR cross-reactivity triplets and HLA context, then fine-tuned on
direct T-cell assay outcomes.

### TCR-Aware Route

When raw tumor RNA contains enough immune-receptor reads:

1. Reconstruct a partial TCR repertoire with TRUST4 or an equivalent adapter.
2. Score candidate peptide-HLA pairs against each confident TCR clonotype.
3. Aggregate interaction scores using clone abundance, reconstruction
   confidence, and a multiple-comparisons correction.
4. Combine the TCR evidence as a constrained residual over the universal score.
5. Fall back to the universal route when repertoire coverage is insufficient.

TCR integration itself is not new. The specific product opportunity is a
confidence-gated route derived from the same bulk tumor RNA already required by
the end-to-end pipeline, rather than requiring dedicated TCR sequencing.

## Training Program

### Phase 1: Stage-Specific Pretraining

Use large datasets only for the biological edge they actually label:

| Data | Training task |
| --- | --- |
| binding and eluted ligands | peptide-HLA presentation |
| IEDB/VDJdb shared-TCR peptide groups | HLA-aware cross-reactivity metric |
| paired TCR-peptide-HLA records | peptide-HLA-TCR interaction |
| healthy immunopeptidome/proteome | self-likeness and tolerance |

Use peptide-, HLA-, and study-cluster holdouts. Generated negatives are
auxiliary contrastive examples, never equivalent to assay-confirmed negatives.

### Phase 2: Direct Recognition Fine-Tuning

Fine-tune on experimentally tested cancer candidates from CEDAR, NCI,
IMPROVE, TESLA, and compatible direct screens. Preserve assay type and study as
explicit variables. Use:

- positive-unlabeled or selection-aware loss when negatives were not uniformly
  screened
- hard negatives that were expressed and strongly presented but non-reactive
- study-balanced sampling
- exact and near-peptide leakage purges
- mutant/wild-type pairs and HLA context

### Phase 3: Patient Ranking

Train only on patient-grouped candidate lists:

- optimize a differentiable top-k/listwise objective around ranks 1-30
- weight each patient's negatives by the 20 available portfolio slots
- use nested patient-disjoint validation for model and blend selection
- calibrate probabilities out of fold
- estimate epistemic uncertainty across folds or model seeds

The final score is learned from out-of-fold stage predictions. Raw labels from a
patient must never train a feature later evaluated on that patient.

### Phase 4: Portfolio Selection

Select 20 candidates by expected utility, not score alone. Penalize redundant
peptides from the same mutation and account for:

- uncertainty
- clonal tumor coverage
- HLA loss/expression
- gene and HLA diversity
- synthesis and assay feasibility

Run both unconstrained top-20 and portfolio selection. Diversity constraints are
accepted only if they increase held-out hits or improve coverage without losing
hits.

## Benchmark Program

1. **IMPROVE official CV:** primary development benchmark; all choices made
   inside nested patient-disjoint folds.
2. **Leave-cohort-out IMPROVE:** domain-shift check across melanoma, bladder,
   and basket cohorts.
3. **NCI/CEDAR/TESLA:** frozen external recognition tests, with overlap purged.
4. **TCR-aware cohort:** separate benchmark requiring patient TCRs and direct
   peptide outcomes.
5. **End to end:** candidate-generation recall plus final wet-lab hits@20 from
   WES/RNA inputs.

Every report includes mean hits@20, precision@20, recall@20, nDCG@20, MRR,
oracle hits@20, oracle capture, confidence intervals, and per-patient results.

## Research Order

1. Reproduce NeoPrecis-Immuno as an external Apache-2.0 baseline on IMPROVE.
2. Add its HLA-aware cross-reactivity score to the patient ranker without
   retraining on IMPROVE labels.
3. Rebuild cross-reactivity pretraining with stricter peptide/TCR/study
   holdouts and direct top-20 downstream selection.
4. Add selection-aware fine-tuning on direct cancer screens.
5. Build the TRUST4 adapter and a TCR-aware benchmark before using TCR evidence
   in a headline claim.
6. Freeze the universal model, then evaluate the conditional TCR residual.
7. Integrate the frozen ranker with the upstream pVACtools-compatible adapter.

An auto-research loop may propose configurations and run these experiments.
Only predeclared metrics on nested validation can promote a hypothesis. An LLM
does not assign production candidate scores.
