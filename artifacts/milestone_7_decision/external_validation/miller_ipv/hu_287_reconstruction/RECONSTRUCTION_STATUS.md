# Miller Hu_287 — T1 reconstruction status (LOCKED_TEST)

_LOCKED_TEST: labels never consulted for download/convert/HLA/expression/calling/generation/ranking_


**Download complete:** True


| run | role | bytes | sha256 |
|---|---|--:|---|
| SRR24836184 | normal_exome | 1632487104 | `be2f87043b534858` |
| SRR24836169 | tumor_exome | 3404956280 | `b1c69a3b3184a8bf` |
| SRR24836183 | tumor_rna | 2469815696 | `7665834788b643aa` |

## Reachable now
- download: **True**
- sra_to_fastq: **True**
- expression_tpm: **True**

## Expression (salmon, label-blind)
60883 genes / 251955 transcripts, ΣTPM=1000000.0. Top genes: ENSG00000283907=14666.9, NPIPB3=14370.67, RNA5-8SN2=9070.18, COL1A1=6690.77, UBC=5488.6, ACTB=5450.25.


## NOT_EVALUABLE stages (machine-actionable — see DEPENDENCY_MANIFEST.md)

- **hla_typing_classI**: missing tool(s): need [OptiType OR arcasHLA OR hla-la] on PATH
- **wes_alignment**: tools present (bwa+samtools) but missing reference/index: ['data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa', 'data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bwt']
- **somatic_calling**: missing tool(s): need [gatk OR strelka] on PATH
- **mutanome_enumeration**: missing tool(s): need [pvacseq OR pvactools] on PATH
- **scoring_prime_epicurus**: local PRIME/MixMHCpred present but UPSTREAM-BLOCKED: needs the re-enumerated candidate universe + class-I HLA before any peptide can be scored

## North-star loop
NOT_EVALUABLE — hard-blocked at HLA typing (no OptiType/arcasHLA) and somatic calling (no Mutect2/Strelka); without somatic variants the shared candidate universe cannot be re-enumerated, so genuine-PRIME-vs-Epicurus hits@20 is not computable. Expression (gene TPM) is the only north-star input reconstructed. Install commands in DEPENDENCY_MANIFEST.md.
