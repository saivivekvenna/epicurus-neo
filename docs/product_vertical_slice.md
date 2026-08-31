# Product vertical slice

## Supported boundary

The first production path starts at a pVACseq `all_epitopes` table, not raw FASTQ. It accepts an
optional RNA-evidence table and emits a deterministic, evidence-backed portfolio. Alignment,
somatic calling, HLA typing, VEP annotation, and peptide generation remain replaceable upstream
adapters.

This boundary is deliberate: Epicurus Neo owns evidence normalization, transparent stage scoring,
portfolio selection, abstention, and reporting. It does not duplicate mature genomics tools.

## Contracts

Production inference and research supervision are separate:

- `CandidateEvidence` contains only values available before a vaccine-response assay.
- `AssayOutcome`, represented by the existing Event-B corpus, contains labels and provenance.
- `RankedCandidate` adds component scores, completeness, uncertainty, tier, rank, and rationale.

The machine-readable contracts live in `configs/schemas/`. A patient run never requires `label`,
`study_id`, `label_weight`, or `assay_type`.

## Run the demonstration

```bash
epicurus-neo validate-patient-input \
  --input examples/demo_patient/pvacseq_all_epitopes.tsv \
  --patient-id DEMO-001

epicurus-neo run-patient \
  --input examples/demo_patient/pvacseq_all_epitopes.tsv \
  --patient-id DEMO-001 \
  --output-dir outputs/demo_patient
```

The output directory contains:

- `ranked_candidates.csv`: every candidate, including exclusions and unselected rows.
- `report.json`: schema/policy metadata and patient-level abstention state.
- `report.md`: a human-readable selected portfolio with evidence rationale.

An RNA table may be added with `--rna-evidence`. It is joined using the most specific safe shared
key, in this order: candidate, mutation plus transcript, mutation, transcript, then gene. Ambiguous
RNA rows are rejected rather than silently aggregated.

## Evidence policy

Before scoring, the default deterministic validity gate removes only rule-verifiable invalid
candidate routes while retaining every rejected row and its reason in the output:

- `HLA_LOH_LOST_ALLELE`: the source explicitly says the presenting HLA allele is lost;
- `GENE_NOT_EXPRESSED`: the source explicitly calls the gene unexpressed;
- malformed peptide, duplicate candidate, invalid class-I length, or mutation/peptide mismatch only
  when the metadata required to establish that condition is present.

These are eligibility rules, not recognition predictions. Vendor binary calls are preserved as
provenance; Epicurus Neo does not silently recreate their thresholds from TPM. The gate is default-on
and may be disabled with `--disable-validity-gate` solely for before/after audits.

The v1 policy exposes four mechanistic components:

1. translated: expression, RNA VAF, and mutant RNA reads;
2. presented: a supplied presentation score, binding percentile, or binding affinity;
3. recognized: a supplied frozen recognition score or mutant-versus-wild-type binding delta;
4. coverage: clonality or tumor DNA VAF.

The weighted geometric score is a prioritization score, **not a probability of vaccine response**.
Unavailable components receive a neutral value and increase reported uncertainty. This avoids
turning missing evidence into an artificial negative.

Hard exclusions are intentionally narrow:

- measured zero RNA expression;
- zero mutant RNA reads at RNA depth of at least 10.

Selected candidates are labelled `CORE`, `SUPPORTING`, or `FILLER`. A patient is marked as
abstained when no candidate clears the core evidence policy, even if lower-confidence candidates are
returned to fill the portfolio.

The core and supporting thresholds are policy defaults, not clinically validated cutoffs. They are
CLI-configurable and must not be described as response-risk bounds until locked external evidence
supports that interpretation.

## Next adapters

The intended upstream deployment profile is:

```text
nf-core/sarek (tumor/normal WES)
  + nf-core/rnaseq (tumor RNA)
  + HLA caller
  -> VEP-annotated somatic VCF + quantified RNA evidence
  -> pVACseq
  -> epicurus-neo run-patient
```

Copy number, purity, HLA loss, and reconstructed bulk-RNA TCR evidence should be added as optional
adapters after the base path has been exercised on real patient exports. Recognition-model work
remains gated on materially better candidate-resolved Event-B data.
