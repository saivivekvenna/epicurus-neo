# Dataset Plan

The repository will normalize public neoantigen immunogenicity datasets into a
single canonical table. Data files are not committed unless licensing permits it.

## Candidate Sources

| Dataset | Role | Notes |
| --- | --- | --- |
| NeoRanking / Immunity 2023 | main training and validation backbone | harmonized patient-level neo-peptide features and labels |
| TESLA / Wells et al. 2020 | locked test | widely recognized held-out exam; do not tune on it |
| Gartner/NCI datasets | locked or secondary test | official train/test sets for neoantigen ranking models |
| CEDAR | auxiliary train/validation and retrieval features | requires careful assay-label normalization and dedupe |
| NEPdb | auxiliary train/validation and negative mining | useful for non-reactive tested neoepitopes |
| dbPepNeo2.0 | positive/retrieval features and external checks | strong overlap risk; dedupe strictly |
| BigMHC data | baseline scores and optional training | keep `im_test` locked when used for comparison |

## Canonical Table

Minimum columns:

- `candidate_id`
- `source_dataset`
- `study_id`
- `patient_id`
- `hla_allele`
- `mutant_peptide`
- `wildtype_peptide`
- `label`
- `label_weight`
- `assay_type`

## Dataset-Specific Label Notes

- Gartner/NCI `Screening Status` is normalized as:
  - `CD8` / `1` -> `positive`
  - `0` / `-` -> `negative`
  - `unscreened` -> `unknown`
- The `-` value is handled only in the Gartner normalizer because other public
  files may use dash-like placeholders for missing values rather than screened
  non-reactivity.

Optional feature columns are allowed and should be prefixed by feature family:

- `binding_*`
- `presentation_*`
- `expression_*`
- `clonality_*`
- `foreignness_*`
- `epitope_neighbor_*`
- `mutation_*`
- `model_*`
