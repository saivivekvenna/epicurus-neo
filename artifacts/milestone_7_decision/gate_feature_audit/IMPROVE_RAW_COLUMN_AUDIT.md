# IMPROVE raw 88-column feature audit (CORRECTION to commit 45af791)

**Correction.** supersedes commit 45af791 'IMPROVE has no orthogonal features' — that audited the reduced prime/el/expr export; the raw 88-col table is feature-rich.

Source: `data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt` — 17520 rows × 88 cols, 70 patients / 81 samples.

**Selection policy.** whitelist chosen from PRE-DECLARED biology + within-patient variation + coverage; held-out response used ONLY as a leakage screen on suspicious/forbidden (to reject).

## Bucket counts

- **deployable**: 39
- **forbidden**: 22
- **context_only**: 13
- **presentation**: 7
- **suspicious_derived**: 4
- **split_only**: 3

## Deployable candidate-varying whitelist (by orthogonal family)

- **agretopicity**: DAI_4.1, DAI
- **clonality**: CelPrev
- **expression**: Expression
- **foreignness_selfsim**: SelfSim, Foreigness
- **hla_expression**: HLAexp
- **mutation_annotation**: Mutation_Consequence, Cancer_Driver_Gene, Misense_mutation, Frameshift_mutation, Inframe_deletion_mutation
- **nn_align**: Of, Gp, Gl
- **physicochemical**: PeptLen, HydroAll, HydroCore, PropSmall, PropAro, PropBasic, PropAcidic, mw, Aro, Inst, PropHydroAro, CysRed, pI
- **rna_support**: rna_confirm, rna_var, rna_total, rna_af, ValMutRNACoef, rna_bin
- **stability_processing**: Stability
- **vaf_readsupport**: VarAlFreq

## Full column classification

| column | bucket | family | cov | within-pt var | n_uniq | response-leak AUROC | note |
|---|---|---|---|---|---|---|---|
| DAI | deployable | agretopicity | 1.0 | 1.0 | 16823 |  | differential agretopicity index (mut vs WT binding) |
| DAI_4.1 | deployable | agretopicity | 0.9908 | 1.0 | 15756 |  | DAI on NetMHCpan-4.1 |
| CelPrev | deployable | clonality | 1.0 | 0.9714 | 5965 |  | cancer-cell prevalence / clonality of the mutation |
| Expression | deployable | expression | 1.0 | 1.0 | 5974 |  | gene expression (TPM-like) |
| Foreigness | deployable | foreignness_selfsim | 1.0 | 1.0 | 496 |  | foreignness / dissimilarity-to-self score |
| SelfSim | deployable | foreignness_selfsim | 1.0 | 1.0 | 14788 |  | self-similarity to self-proteome (lower=more foreign) |
| HLAexp | deployable | hla_expression | 1.0 | 0.9571 | 202 |  | expression of the restricting HLA allele |
| Cancer_Driver_Gene | deployable | mutation_annotation | 1.0 | 0.8857 | 2 |  | driver-gene flag |
| Frameshift_mutation | deployable | mutation_annotation | 1.0 | 0.8429 | 2 |  | frameshift flag |
| Inframe_deletion_mutation | deployable | mutation_annotation | 1.0 | 0.5714 | 2 |  | inframe-del flag |
| Inframe_insertion | deployable | mutation_annotation | 1.0 | 0.2714 | 2 |  | inframe-ins flag |
| Misense_mutation | deployable | mutation_annotation | 1.0 | 0.8857 | 2 |  | missense flag |
| Mutation_Consequence | deployable | mutation_annotation | 1.0 | 0.8857 | 4 |  | categorical consequence class |
| Gl | deployable | nn_align | 1.0 | 1.0 | 3 |  | gap length (low value) |
| Gp | deployable | nn_align | 1.0 | 1.0 | 9 |  | gap position (low value) |
| Il | deployable | nn_align | 1.0 | 0.3714 | 2 |  | insertion length (low value) |
| Ip | deployable | nn_align | 1.0 | 0.3714 | 9 |  | insertion position (low value) |
| Of | deployable | nn_align | 1.0 | 0.5857 | 3 |  | NetMHC alignment offset (low value) |
| Aro | deployable | physicochemical | 1.0 | 1.0 | 20 |  | aromaticity |
| CysRed | deployable | physicochemical | 1.0 | 1.0 | 18 |  | reduced-cysteine count |
| HydroAll | deployable | physicochemical | 1.0 | 1.0 | 1770 |  | hydrophobicity (whole peptide) |
| HydroCore | deployable | physicochemical | 1.0 | 1.0 | 1303 |  | hydrophobicity (core) |
| Inst | deployable | physicochemical | 1.0 | 1.0 | 3334 |  | instability index |
| PeptLen | deployable | physicochemical | 1.0 | 1.0 | 4 |  | peptide length |
| PropAcidic | deployable | physicochemical | 1.0 | 1.0 | 14 |  | acidic fraction |
| PropAro | deployable | physicochemical | 1.0 | 1.0 | 13 |  | aromatic fraction |
| PropBasic | deployable | physicochemical | 1.0 | 1.0 | 16 |  | basic fraction |
| PropHydroAro | deployable | physicochemical | 1.0 | 1.0 | 35 |  | hydrophobic+aromatic fraction |
| PropSmall | deployable | physicochemical | 1.0 | 1.0 | 17 |  | small-residue fraction |
| mw | deployable | physicochemical | 1.0 | 1.0 | 13634 |  | molecular weight |
| pI | deployable | physicochemical | 1.0 | 1.0 | 894 |  | isoelectric point |
| ValMutRNACoef | deployable | rna_support | 1.0 | 0.9857 | 1895 |  | validated-mutation RNA coefficient |
| rna_af | deployable | rna_support | 1.0 | 0.9857 | 1643 |  | RNA allele frequency |
| rna_bin | deployable | rna_support | 1.0 | 0.9857 | 3 |  | RNA-confirmed (binarised) |
| rna_confirm | deployable | rna_support | 0.8466 | 0.9857 | 3067 |  | ref/alt RNA read-count string (needs parse) |
| rna_total | deployable | rna_support | 1.0 | 0.9857 | 819 |  | RNA total read depth at locus |
| rna_var | deployable | rna_support | 1.0 | 0.9857 | 306 |  | RNA reads supporting the variant |
| Stability | deployable | stability_processing | 1.0 | 1.0 | 97 |  | pMHC stability (NetMHCstab) |
| VarAlFreq | deployable | vaf_readsupport | 1.0 | 1.0 | 706 |  | DNA variant allele frequency |
| IB_CB | suspicious_derived |  | 1.0 | 1.0 | 16823 | 0.5007 | unknown composite magnitude (CB mean~17.6 vs IB~0.78); provenance needed |
| IB_CB_cat | suspicious_derived |  | 1.0 | 1.0 | 2 |  | IB vs CB category of IB_CB; provenance needed |
| NetMHCExp | suspicious_derived |  | 1.0 | 1.0 | 16987 | 0.5195 | NetMHC x Expression composite — double-counts presentation+expression primitives |
| PrioScore | suspicious_derived |  | 1.0 | 1.0 | 95 | 0.5038 | IMPROVE's OWN 0-100 priority bucket — model output, circular; provenance needed |
| Blinage | context_only |  | 1.0 | 0.1571 | 81 |  | B lineage (per sample) |
| CYT | context_only |  | 1.0 | 0.1571 | 81 |  | cytolytic activity (per sample) |
| CytoxLympho | context_only |  | 1.0 | 0.1571 | 81 |  | cytotoxic lymphocytes (per sample) |
| Endothelial | context_only |  | 1.0 | 0.1571 | 81 |  | endothelial (per sample) |
| Fibroblasts | context_only |  | 1.0 | 0.1571 | 81 |  | fibroblasts (per sample) |
| MCPmean | context_only |  | 1.0 | 0.1571 | 81 |  | MCP-counter mean (per sample) |
| Monocytes | context_only |  | 1.0 | 0.1571 | 81 |  | monocytes (per sample) |
| MyeloidDC | context_only |  | 1.0 | 0.1571 | 81 |  | myeloid DC (per sample) |
| NKcells | context_only |  | 1.0 | 0.1571 | 81 |  | NK cells (per sample) |
| Neutrophils | context_only |  | 1.0 | 0.1571 | 81 |  | neutrophils (per sample) |
| Sample_TME | context_only |  | 1.0 | 0.1571 | 81 |  | TME class (per sample) |
| Tcells | context_only |  | 1.0 | 0.1571 | 81 |  | MCP-counter T cells (per sample) |
| TcellsCD8 | context_only |  | 1.0 | 0.1571 | 78 |  | CD8 T cells (per sample) |
| Prime | presentation |  | 1.0 | 1.0 | 16200 |  | PRIME immunogenicity+presentation score (incumbent scorer) |
| RankBA | presentation |  | 1.0 | 1.0 | 14990 |  | NetMHCpan BA %rank (mut) |
| RankBA_4.1 | presentation |  | 1.0 | 1.0 | 4722 |  | NetMHCpan-4.1 BA %rank |
| RankEL | presentation |  | 1.0 | 1.0 | 10886 |  | NetMHCpan-4.0 EL %rank (mut) |
| RankEL_4.1 | presentation |  | 1.0 | 1.0 | 3673 |  | NetMHCpan-4.1 EL %rank |
| RankEL_wt | presentation |  | 1.0 | 1.0 | 13415 |  | EL %rank (WT) |
| RankEL_wt_4.1 | presentation |  | 0.9908 | 1.0 | 6332 |  | EL %rank WT 4.1 |
| Core | forbidden |  | 1.0 | 1.0 | 15063 |  | binding-core residues (sequence-derived) |
| CoreNonAnchor | forbidden |  | 1.0 | 1.0 | 15074 |  | non-anchor core residues (sequence-derived) |
| Gene_ID | forbidden |  | 1.0 | 1.0 | 4423 |  | gene identifier |
| Gene_Symbol | forbidden |  | 1.0 | 1.0 | 4409 |  | gene symbol (identity; driver captured via Cancer_Driver_Gene) |
| Genomic_Position | forbidden |  | 1.0 | 1.0 | 6016 | 0.5105 | genomic coordinate (identity) |
| HLA_allele | forbidden |  | 1.0 | 0.9857 | 36 |  | restricting allele (identity) |
| HLA_num | forbidden |  | 1.0 | 0.9857 | 35 |  | allele-count/identity string |
| HLA_type | forbidden |  | 1.0 | 0.9571 | 3 |  | hetero/homo HLA status string |
| Mut_peptide | forbidden |  | 1.0 | 1.0 | 15293 |  | mutant peptide sequence (identity) |
| Norm_peptide | forbidden |  | 0.9908 | 1.0 | 15000 |  | WT peptide sequence (identity) |
| Patient | forbidden |  | 1.0 | 0.0 | 70 |  | patient id |
| PeptNorm | forbidden |  | 1.0 | 1.0 | 15143 |  | normalised peptide sequence (identity) |
| Protein_position | forbidden |  | 1.0 | 1.0 | 4993 | 0.5221 | protein coordinate (identity) |
| Sample | forbidden |  | 1.0 | 0.1571 | 81 | 0.5296 | sample id |
| Transcript_ID | forbidden |  | 1.0 | 1.0 | 5246 |  | transcript identifier |
| identi_pep_patient | forbidden |  | 1.0 | 1.0 | 17520 |  | peptide-patient row id |
| identity | forbidden |  | 1.0 | 1.0 | 17520 |  | row identity string |
| norm_pMHC | forbidden |  | 1.0 | 1.0 | 17087 |  | allele_peptide IDENTITY string (WT) — NOT a score |
| pMHC | forbidden |  | 1.0 | 1.0 | 17353 |  | allele_peptide IDENTITY string (mutant) — NOT a score |
| response | forbidden |  | 1.0 | 0.8714 | 2 | 1.0 | THE LABEL (immunogenic 0/1) |
| sample_hla | forbidden |  | 1.0 | 0.9571 | 202 |  | sample-allele identity |
| validation | forbidden |  | 1.0 | 0.0 | 1 |  | QC flag (constant 'Sufficient'); outcome-adjacent, unusable |
| Loci | split_only |  | 1.0 | 0.0 | 1 |  | constant (single value) |
| Partition | split_only |  | 1.0 | 0.0 | 5 |  | CV fold assignment |
| cohort | split_only |  | 1.0 | 0.0 | 3 |  | tumour cohort (bladder/melanoma/Basket) — context/split, not deployable |

## Notes on the columns the previous session flagged

- **PrioScore** — IMPROVE's own 0–100 priority bucket (95 distinct values); marginal response AUROC ≈ 0.50 so it is not leaking the outcome, but it is the pipeline's OWN composite output (circular / double-counts primitives) → suspicious_derived, excluded from the clean whitelist.
- **CelPrev** — cancer-cell prevalence / clonality of the mutation, candidate-varying (within-patient var ≈ 0.97); legitimate clonality axis → deployable.
- **IB_CB / IB_CB_cat** — a magnitude with a two-level category (CB≈17.6 vs IB≈0.78); semantics undocumented in the reduced export → suspicious_derived pending provenance, not deployed.
- **NetMHCExp** — NetMHC × Expression composite: double-counts presentation and expression primitives that are separately available → suspicious_derived (use the primitives instead).
- **pMHC / norm_pMHC** — these are `allele_peptide` IDENTITY strings (mutant / WT), NOT composite scores → forbidden/identity.
- **validation** — constant 'Sufficient' (an RNA-QC flag), not an experimental outcome and useless as a feature → forbidden/constant.
- **Immune-deconvolution block** (Tcells…MCPmean, CYT, Sample_TME) is per-sample constant (within-patient var ≈ 0.16): it can move a whole patient's prior but cannot re-order that patient's candidates → context_only, never a within-patient gate feature.