# PDAC NeoVax 2023 source and Event-B semantics

## Identity and overlap

- Primary publication: Rojas et al., *Nature* (2023),
  `doi:10.1038/s41586-023-06063-y`, `PMCID:PMC10171177`.
- Trial/cohort: `NCT04161755`, 16 vaccinated biomarker-evaluable participants.
- The 2025 Nature report (`doi:10.1038/s41586-024-08508-4`) is long-term follow-up of the same
  cohort. It is provenance for persistence, not a second independent study or 16 new patients.

## Frozen source

`41586_2023_6063_MOESM4_ESM.xlsx` (Supplementary Table 5), SHA-256
`c9942caa0de461e87c3725ae42377fea2a283cd170074fb2e36820303b429595`.
The compact public workbook remains under the ignored `data/raw/pdac_neovax_2023/` tree.

## Locked semantics

The table lists 232 manufactured personalized targets for 16 vaccinated patients. Each target
preserves its encoded mutant sequence, wild-type sequence, mutation, transcript, best predicted
class-I epitope/HLA and best predicted class-II epitope/HLA. Predictions are not recorded as
experimental restrictions.

The source values map as follows:

- `De novo response`: Event B positive at the encoded-target assay level (23 targets).
- `No response`: explicit Event B tested negative (200 targets).
- `De novo response in pool`: Event B untested at candidate level (7 patient-25 targets).
- `No data`: Event B untested (2 targets).

Patient 25 had two positive combined pools spanning seven targets. The public table does not reveal
the 2-versus-5 membership partition. One pool-level positive observation is retained with a value of
two positive pools, while all seven component labels remain `UNTESTED`. This reconciles the paper's
25 positive response entities as 23 single targets plus 2 combined pools without label inflation.

The encoded mRNA antigen, its overlapping 15-mer stimulation pool, and predicted minimal class-I
and class-II epitopes are separate entities linked explicitly. Clinical recurrence and the response
to concurrent atezolizumab are not Event-B labels.
