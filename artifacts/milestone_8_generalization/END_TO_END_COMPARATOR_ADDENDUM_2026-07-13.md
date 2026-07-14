# End-to-end comparator addendum — 2026-07-13

## Objective refinement

The north-star comparison is a patient-level, raw-input, end-to-end benchmark.
Epicurus Neo and the closest runnable competing workflow receive equivalent raw
tumor WES, matched-normal WES, and tumor RNA evidence and must each produce a
ranked vaccine portfolio without access to recognition outcomes.

The primary question is:

> From equivalent raw patient evidence, which complete pipeline places more
> experimentally recognized neoantigens in its final top-20 portfolio on sealed
> held-out patients?

This addendum changes the headline comparator, not the frozen Epicurus biology,
patient split, label-isolation rules, or final-held-out seal.

## Comparator hierarchy

1. **Primary end-to-end comparator: nextNEOpi.** It is the closest open-source
   raw WES/WGS plus optional RNA neoantigen workflow and therefore tests the
   product-level claim.
2. **Secondary workflow comparators: pVACtools/pVACseq and OpenVax/Vaxrank.**
   These begin from later-stage artifacts and are evaluated only in an
   input-controlled track where their required VCF, expression, RNA, and HLA
   inputs are supplied without recognition labels.
3. **Secondary component diagnostic: genuine PRIME.** PRIME remains useful for
   isolating reranker behavior on a common reachable peptide universe, but it is
   not the headline end-to-end opponent.

## Two-track design

### Track A — full-pipeline autonomy (primary)

- Input: the same patient tumor WES, matched-normal WES, and tumor RNA FASTQs.
- Each complete workflow performs its own supported alignment, somatic calling,
  HLA inference, expression integration, peptide construction, filtering, and
  ranking.
- No artifact from Epicurus may silently replace a failed or missing competitor
  stage.
- Before local lifecycle cleanup, every FASTQ consumed by Epicurus is pinned by
  byte size, read count, and SHA-256 in `CONVERT_PROVENANCE.json`. The preserved
  public SRA may be converted again on the nextNEOpi host, but the comparator
  bundle must reject every regenerated FASTQ whose basename, size, or SHA-256
  differs. Matching read counts alone do not establish identical raw input.
- A pipeline failure, unreachable measured target, or absent top-20 portfolio is
  reported as such rather than repaired after labels are opened.

### Track B — controlled downstream comparison (diagnostic)

- Input: one frozen common somatic-variant, RNA-expression, and HLA evidence
  bundle produced before recognition labels are loaded.
- Epicurus, pVACtools/pVACseq, Vaxrank, and genuine PRIME operate only on inputs
  each tool officially supports.
- This track diagnoses whether differences arise in candidate generation,
  presentation filtering, recognition ranking, or portfolio construction. It
  cannot replace Track A as evidence for the end-to-end claim.

## Locked evaluation unit and metrics

- Unit: patient, not pooled peptide row.
- Primary endpoint: recognized positives in the final top 20, reported per
  patient and as a paired mean across the held-out cohort.
- Required supporting endpoints: recall@20, measured-positive reachability
  before ranking, tested-negative count in the top 20, candidate-universe size,
  portfolio size, and pipeline failure/abstention rate.
- Diagnostic endpoints: rank of each reachable measured positive, loss stage
  for each unreachable measured positive, runtime, and peak storage.
- No pooled AUROC can substitute for the patient-level primary endpoint.
- Ties, missing outputs, and fewer-than-20 portfolios are reported explicitly.

## Fairness and anti-cherry-picking rules

Before any final-held-out recognition label is opened, the benchmark must freeze:

1. exact container/tool versions and reference assets;
2. raw-input manifests and checksums;
3. allowed configuration for every workflow;
4. deterministic conversion from each workflow's native output to one ranked
   top-20 portfolio;
5. peptide/variant identity normalization and matching rules;
6. timeout, retry, failure, and abstention policy;
7. output hashes for every label-blind candidate universe and portfolio.

Default configurations are preferred. Any non-default configuration must be
motivated without recognition outcomes and applied uniformly to every patient.
No comparator-specific rescue is permitted after label access.

## Interpretation

Epicurus is not novel merely because it processes raw sequencing. The intended
contribution is demonstrated only if its complete, patient-agnostic decision
policy improves experimentally recognized top-20 recovery under the sealed
end-to-end comparison, or if the diagnostic track identifies a reproducible and
clinically meaningful advantage at a specific stage.

The prior genuine-PRIME protocol remains valid as a component-level analysis.
Where its verdict vocabulary conflicts with this addendum, the final report must
label it `COMPONENT_LEVEL` and reserve the headline end-to-end verdict for Track
A against nextNEOpi.
