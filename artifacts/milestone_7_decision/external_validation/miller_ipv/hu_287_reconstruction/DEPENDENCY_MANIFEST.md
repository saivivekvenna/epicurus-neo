# Miller Hu_287 — reconstruction dependency/reference manifest

_LOCKED_TEST. PRIME/MixMHCpred binaries are present on disk (gitignored) but scoring is upstream-blocked until the candidate universe is re-enumerated. HLA typing, somatic calling, and mutanome enumeration are the hard blockers on this toolchain; expression via salmon is reachable once the GENCODE index is built. No released processed Hu_287 VCF/TPM/HLA exists (web 2026-07-12)._


## Installed tools

- **fasterq-dump**: fasterq-dump : 3.4.1
- **bwa**: Program: bwa (alignment via Burrows-Wheeler transformation)
- **samtools**: samtools 1.24
- **bcftools**: bcftools 1.24
- **salmon**: salmon 2.3.3
- **PRIME**: /Users/saivivekvenna/conductor/repos/epicurus-neo/data/raw/tools/PRIME/PRIME
- **MixMHCpred**: /Users/saivivekvenna/conductor/repos/epicurus-neo/data/raw/tools/MixMHCpred/MixMHCpred

## Per-stage status (machine-actionable)

| stage | status | produces | missing / reason | install hint |
|---|---|---|---|---|
| sra_to_fastq | **RUNNABLE** | paired FASTQ per run |  |  |
| hla_typing_classI | **NOT_EVALUABLE** | 4-digit class-I HLA (A/B/C) from normal exome (needed for PRIME/EL ranking) | OptiType: missing tools ['OptiTypePipeline.py', 'razers3'], missing refs ['data/raw/refs/hla/hla_reference_dna.fasta']; arcasHLA: missing tools ['arcasHLA'], missing refs ['data/raw/refs/hla/IMGTHLA']; T1K: missing tools ['run-t1k'], missing refs ['data/raw/refs/hla/t1k_hlaidx'] | conda install -c bioconda optitype razers3   # or: conda install -c bioconda arcas-hla |
| wes_alignment | **NOT_EVALUABLE** | coordinate-sorted, dup-marked tumor & normal BAM | bwa-mem+samtools: missing refs ['data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.bwt', 'data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.sa'] |  |
| somatic_calling | **NOT_EVALUABLE** | PASS somatic SNV/indel VCF with tumor VAF + depth | Mutect2: missing refs ['data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.dict']; Strelka2: missing tools ['configureStrelkaSomaticWorkflow.py'] | conda install -c bioconda gatk4   # Mutect2 tumor-vs-normal; or bioconda strelka |
| rna_quant | **RUNNABLE** | per-gene TPM (+ mutant-allele RNA evidence via genome align) |  |  |
| mutanome_enumeration | **NOT_EVALUABLE** | full class-I 8-11mer lossless peptide universe (shared by the lossless arms; the pvac arm generates its own set from the same base variants) | VEP+lossless: missing tools ['vep'], missing refs ['data/raw/refs/vep/homo_sapiens', 'data/raw/refs/gencode/gencode.v44.annotation.gtf']; pvacseq: missing tools ['pvacseq', 'vep'], missing refs ['data/raw/refs/vep/homo_sapiens'] | pip install pvactools && conda install -c bioconda ensembl-vep   # + VEP cache |
| scoring_prime_epicurus | **NOT_EVALUABLE** | genuine PRIME AND frozen Epicurus over the IDENTICAL universe | local PRIME/MixMHCpred present but UPSTREAM-BLOCKED: needs the re-enumerated candidate universe + class-I HLA before any peptide can be scored |  |

## References required (not auto-fetched)

| name | purpose | present | dest |
|---|---|:--:|---|
| GRCh38_primary_assembly | WES alignment + somatic calling | y | `data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa` |
| GENCODE_v44_transcripts | RNA quantification (salmon index) | y | `data/raw/refs/gencode/gencode.v44.transcripts.fa` |
| VEP_cache_GRCh38 | variant annotation for mutanome enumeration (pVACtools) | n | `data/raw/refs/vep/homo_sapiens` |
| IMGT_HLA_dna | class-I HLA typing (OptiType/arcasHLA reference) | n | `data/raw/refs/hla/hla_reference_dna.fasta` |

### Acquire commands

- **GRCh38_primary_assembly**: `curl -L https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | gunzip > <dest>; bwa index <dest>; samtools faidx <dest>`
- **GENCODE_v44_transcripts**: `curl -L https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.transcripts.fa.gz | gunzip > <dest>; salmon index -t <dest> -i data/raw/refs/gencode/salmon_index -k 31`
- **VEP_cache_GRCh38**: `curl -L https://ftp.ensembl.org/pub/release-110/variation/indexed_vep_cache/homo_sapiens_vep_110_GRCh38.tar.gz | tar xz -C data/raw/refs/vep`
- **IMGT_HLA_dna**: `OptiType ships hla_reference_dna.fasta; arcasHLA: 'arcasHLA reference --version latest'`
