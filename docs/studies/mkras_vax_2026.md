# mKRAS-VAX 2026 source and Event-B semantics

## Identity and overlap

- Primary publication: Huff et al., *Nature Communications* (2026),
  `doi:10.1038/s41467-026-68324-4`.
- Trial/cohort: `NCT04117087`, 12 immunogenicity-evaluable participants with resected PDAC.
- This adapter covers the resected-PDAC cohort only. A later metastatic-CRC publication using
  the same trial identifier is a distinct cohort and is not included or counted as independent.

## Frozen sources

- `41467_2026_68324_MOESM1_ESM.pdf`, SHA-256
  `8bc2ad373da7dabcd6efb2b25926f2eacc04456333e9474f06a4fb9ebd6ff61b`.
- `41467_2026_68324_MOESM4_ESM.xlsx`, SHA-256
  `b7ac94e8c1427b63eb35b503c60ea12bea863165402b7a11cc2d83d927e09a7a`.

Both files are public Springer supplementary objects. Raw copies remain under the ignored
`data/raw/mkras_vax_2026/` directory.

## Locked semantics

The vaccine contains the same six mutant-KRAS 21-mers for every patient. Figure 2c is an explicit
12-by-6 matrix and defines `1=positive response, 0=negative response`. These 72 cells—not the
patient-level immune-responder endpoint—are the primary candidate labels. The resulting source
calls are 60 positive and 12 tested negative.

All six peptides were individually tested by ex-vivo IFN-gamma ELISpot at baseline and after
vaccination. The paper reports no pre-vaccination response. Baseline observations are retained as
Event A tested negatives; the Figure 2c within-17-week calls are Event B. Clinical recurrence is
not ingested as an immune label.

Six global antigen entities hold the shared sequences. Seventy-two patient-candidate entities hold
patient-specific vaccination and assay observations. Only the matching component receives that
patient's directly confirmed tumor-KRAS mutation; the other five are shared vaccine components,
not inferred tumor mutations. No HLA restriction is inferred from long-peptide response.
