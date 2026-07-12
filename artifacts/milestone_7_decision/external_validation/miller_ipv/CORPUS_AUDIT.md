# Miller IPV — recognition-label corpus audit

> Status **LABELS_INGESTED_AND_VALIDATED** · LOCKED_TEST · DOI 10.1126/scitranslmed.abj9905 (PMID 38416845).
> Canonical label source: S2 (Data_file_S2.xlsx, 'Sheet 1'): 754 tested 20-mers x 13 patients.

## Counts
- peptides tested: **754** · patients: **13** · mutations: **343**
- IFN-g label (primary): {'TESTED_NEGATIVE': 574, 'POSITIVE': 180}
- positive peptides — IFN-g 180 · IL-5 143 · any 199 (paper's 199)
- recognized mutations — IFN-g 134 · any 145
- peptide lengths: {20: 743, 21: 8, 22: 3} · variant types: {'SNV': 711, 'DEL': 14, 'COMPLEX': 14, 'MNP': 8, 'INS': 7}

## Contract validation
- ok=**True** · rows=754 · invalid_labels=[] · invalid_peptides=0 · conflicting_keys=0

## SRA crosswalk
- all 13 label patients have public inputs: **True** · input trios complete: True

## Per-patient (IFN-g)

| patient | peptides | IFN-g+ peptides | mutations | recognized mutations |
|---|---:|---:|---:|---:|
| Hu_048 | 89 | 15 | 43 | 11 |
| Hu_159 | 33 | 9 | 15 | 6 |
| Hu_182 | 74 | 40 | 35 | 24 |
| Hu_250 | 27 | 5 | 12 | 4 |
| Hu_254 | 51 | 11 | 23 | 9 |
| Hu_268 | 51 | 13 | 21 | 7 |
| Hu_277 | 11 | 2 | 6 | 2 |
| Hu_287 | 19 | 3 | 5 | 3 |
| Hu_293 | 17 | 6 | 9 | 4 |
| Hu_315 | 135 | 40 | 65 | 33 |
| Hu_333 | 75 | 14 | 35 | 11 |
| Hu_343 | 111 | 14 | 50 | 12 |
| Hu_344 | 61 | 8 | 30 | 8 |

## Reconciliation notes

- Paper reports 349 tested variants; S2 has 343 distinct (gene,chr,pos) mutations — a 6-variant gap (likely a few multi-transcript/indel rows collapsing on coordinates); documented, not resolved.
- Paper's '199 (26 pct) induced T cell responses' = the 'any' (IFN-g OR IL-5) set (verified 199); the PRIMARY class-I label here is IFN-g (180 pos peptides).

## Remaining gaps for the four-arm run

- HLA per patient: NOT in S2 (20-mer ELISpot, no MHC restriction) -> type from WES (OptiType).
- RNA expression/TPM: NOT in S2 -> quantify from the RNA-seq runs.
- Full candidate universe: S2 has 343 IPV-PREFILTERED tested mutations -> re-enumerate the full class-I mutanome from WES for a fair denominator (both arms share it).
- 20-mer -> class-I minimal epitope: recover mutant residue position (present: derivable from mut vs ref peptide diff) and enumerate 8-11mers spanning it for PRIME/Epicurus scoring.
- => the LABEL half is DONE; the ranking arms are BLOCKED on the SRA download + bioinformatics (HLA + expression + mutanome), not on any missing file.
