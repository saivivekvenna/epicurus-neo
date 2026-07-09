# Epicurus Neo System Architecture

## Product Boundary

Epicurus Neo accepts tumor/normal whole-exome sequencing, tumor RNA sequencing,
and patient metadata, then emits a ranked, evidence-backed top-20 neoantigen
portfolio.

The system should not train one opaque model directly from raw sequencing reads
to peptide rank. Direct T-cell labels are far too scarce for that. The correct
architecture is a mechanistic candidate-generation pipeline followed by a
patient-aware learned ranker.

```text
tumor WES + normal WES + tumor RNA
    -> QC, alignment, somatic variants, HLA typing
    -> annotation, phasing, expression, VAF, clonality, HLA viability
    -> mutant/wild-type peptide candidates
    -> processing and presentation predictors
    -> patient-aware recognition ranker
    -> uncertainty-aware top-20 portfolio
```

## Upstream Genomics

Use established tools as replaceable adapters rather than rebuilding their core
logic:

1. QC and alignment for tumor WES, matched-normal WES, and tumor RNA.
2. Somatic SNV/indel calling with caller agreement and artifact filters.
3. HLA class-I typing from normal WES, with RNA support when available.
4. Variant annotation and local phasing so nearby variants produce the correct
   mutant protein sequence.
5. RNA evidence:
   - gene and transcript expression
   - mutant RNA read count and RNA VAF
   - allele-specific expression when available
6. Copy number, purity, cancer cell fraction, and mutation clonality.
7. HLA loss of heterozygosity, HLA expression, and antigen-processing machinery
   expression.
8. Candidate generation for missense variants, indels, frameshifts, fusions,
   and splice alterations.

pVACtools is the initial candidate-generation adapter because pVACseq already
integrates annotated somatic variants, patient HLA alleles, DNA/RNA evidence,
mutant/wild-type sequences, and multiple binding predictors:
<https://pvactools.readthedocs.io/en/stable/>.

## Epicurus-Owned Hard Part

Epicurus starts from a canonical candidate table. Each row must preserve:

```text
patient -> sample -> mutation -> transcript -> mutant/wild-type peptide
-> HLA allele -> evidence source -> assay -> outcome
```

The ranker should model separate biological stages rather than flattening every
proxy into an immunogenicity label:

```text
P(useful target)
  = P(translated)
  * P(presented | peptide, HLA, tumor)
  * P(recognized | peptide-HLA, patient)
  * P(tumor coverage | clonality, HLA viability)
```

The final learned model may correct this factorization, but each component must
remain observable for debugging and explanation.

### Required Feature Families

- variant confidence, tumor VAF, RNA VAF, expression, and clonality
- mutant and wild-type peptide representations
- mutation position: HLA anchor versus TCR-facing residue
- binding, processing, presentation, and stability ensembles
- mutant-versus-wild-type presentation delta
- normal-proteome and healthy-ligand similarity
- validated reactive and screened-negative recognition neighborhoods
- HLA expression/loss and antigen-processing state
- patient immune context when available from RNA
- model disagreement, missingness, and out-of-distribution indicators

### Ranking Objective

Train and select models on patient-level `hits@20`. AUROC is diagnostic only.
A useful training recipe is:

1. Auxiliary pretraining on stage-specific data.
2. Direct recognition training on experimentally screened cancer candidates.
3. Patient- and study-disjoint validation.
4. Listwise or hard-negative ranking loss focused on the top-20 boundary.
5. Calibration and uncertainty estimation.
6. Portfolio selection with optional gene, mutation, and HLA diversity
   constraints.

An LLM may propose experiments, normalize evidence, and explain predictions. It
must not directly assign the peptide score used for a benchmark claim.

## Benchmark Hierarchy

Do not force one dataset to validate the whole system.

| Level | Question | Primary evidence |
| --- | --- | --- |
| Candidate generation | Did the true mutation/peptide enter the candidate set? | truth sets and known validated patient targets |
| Presentation | Is the peptide displayed by the patient HLA? | ligand and binding benchmarks |
| Recognition ranker | Are reactive candidates high in each patient's list? | patient-level direct screens |
| External generalization | Does the frozen policy transfer to new studies and assays? | study-held-out cohorts |
| End to end | Do WES/RNA inputs produce validated top-20 hits? | hackathon patient cases and wet-lab results |

### Current Quantitative Reality

| Benchmark | Unit | Current hits@20 | Oracle hits@20 | Role |
| --- | --- | ---: | ---: | --- |
| BigMHC `im_test` | HLA allele | 2.5556 | 3.5185 | component regression |
| IMPROVE official CV | patient | 1.4714 | 6.4571 | primary hard-part target |
| 2025 multimer cohort | patient | not accepted | 1.3077 | external robustness |

The `5 hits@20` goal is mathematically impossible on BigMHC but valid on
IMPROVE. Therefore:

- freeze the current BigMHC headline as a component baseline;
- optimize patient-disjoint IMPROVE mean hits@20 toward 5;
- require gains to survive leave-study/cohort-out validation;
- use the multimer cohort as a domain-shift test, not as a 5-hit benchmark.

Always report:

- raw mean hits@20
- oracle mean hits@20
- oracle capture rate
- precision, recall, nDCG, and MRR
- per-patient positive availability
- candidate-generation recall for end-to-end runs

## Build Order

1. Freeze and document the canonical WES/RNA-to-candidate contract.
2. Reproduce the strongest source baselines on IMPROVE official folds.
3. Build a patient-aware learning-to-rank baseline using every feature that can
   be derived at inference time from WES/RNA.
4. Add mutant/wild-type paired representations and patient context.
5. Add assay-aware auxiliary pretraining only when it improves patient-disjoint
   hits@20.
6. Run leave-cohort-out and independent multimer validation.
7. Integrate the frozen ranker behind the pVACtools candidate adapter.
8. Produce the final top-20 table with uncertainty, provenance, and rationale.
