# Epicurus Neo

**Prioritize personalized cancer-vaccine neoantigens from raw sequencing — end to end.**

Epicurus takes a patient's tumor/normal whole-exome sequencing and tumor RNA-seq and produces a
ranked, at-most-20-candidate neoantigen portfolio for a personalized cancer vaccine. It runs as a
single command on your own machine:

```bash
epicurus run-pipeline --config patient.yaml --output-dir out/PATIENT-001
```

```
tumor WES + normal WES + tumor RNA  ──▶  ranked top-20 vaccine neoantigen portfolio
```

## What Epicurus is (and is not)

Epicurus **orchestrates a complete pipeline** and owns the final prioritization. It does **not**
reimplement variant callers or peptide generators — those are established, validated tools that
Epicurus drives:

| Stage | Tool | Does |
|-------|------|------|
| align | BWA-MEM2 + MarkDuplicates | FASTQ → tumor/normal BAM |
| call | GATK Mutect2 + FilterMutectCalls | BAMs → filtered somatic VCF |
| annotate | Ensembl VEP | annotate variants |
| express | Salmon | RNA transcript TPM |
| hla | OptiType / arcasHLA | class-I HLA typing |
| generate | pVACtools + MHCflurry | candidate neoantigen table |
| **prioritize** | **Epicurus** | validity gate → calibrated ranking → ≤20 portfolio |
| report | Epicurus | portfolio CSV + JSON summary + provenance |

**What Epicurus itself contributes** is the `prioritize` stage: a deterministic biological-validity
gate (lost-HLA routes and unexpressed genes cannot take a top-20 slot), a transparent evidence
score combining translation / presentation / recognition / expression evidence, and a
diversity-constrained portfolio selection with patient-level abstention.

**Honest positioning.** In head-to-head evaluation Epicurus's *ranking* is at parity with strong
published rerankers (e.g. PRIME); its measured advantage comes from **portfolio diversification and
full-evidence routing**, not from a novel immunogenicity model. The scores are transparent
evidence-prioritization scores — **not** validated response probabilities. Neoantigen recognition
remains an open scientific problem, and this tool does not claim to have solved it.

## Requirements

Neoantigen calling from raw reads is heavyweight; these requirements are inherent, not incidental:

- **OS:** Linux, or any host via the provided container.
- **Reference data:** GRCh38 + GATK resource bundle + VEP cache + Salmon index — a one-time
  download (tens–hundreds of GB) via `epicurus fetch-references`.
- **Compute:** a 30× tumor/normal WES pair takes hours of CPU.

Check your machine is ready before a run:

```bash
epicurus doctor --bundle-dir ~/.epicurus/references/GRCh38
```

## Install

```bash
# Container (recommended — brings every external tool):
docker run --rm -v "$PWD":/work epicurus run-pipeline --config /work/patient.yaml --output-dir /work/out

# Or a local bioconda environment on Linux:
conda env create -f environment.yml
pip install -e .
```

## Run a patient (full pipeline)

`patient.yaml`:

```yaml
patient_id: PATIENT-001
inputs:
  tumor_wes:  [tumor_R1.fastq.gz, tumor_R2.fastq.gz]
  normal_wes: [normal_R1.fastq.gz, normal_R2.fastq.gz]
  tumor_rna:  [rna_R1.fastq.gz, rna_R2.fastq.gz]
references:
  bundle_dir: ~/.epicurus/references/GRCh38
prioritize:
  k: 20
  max_per_mutation: 1
  max_per_gene: 4
```

```bash
epicurus run-pipeline --config patient.yaml --output-dir out/PATIENT-001
```

Every stage is resumable: re-running reuses valid cached artifacts. Use `--start`/`--stop` to run a
subrange and `--force` to recompute.

## Already have candidates? Start at the ranking stage

If you have already run pVACseq, skip the front-end and prioritize directly:

```bash
epicurus run-patient \
  --input examples/demo_patient/pvacseq_all_epitopes.tsv \
  --patient-id DEMO-001 \
  --output-dir out/demo_patient
```

This writes a full candidate CSV, a machine-readable JSON summary, and a human-readable portfolio
report. See [`docs/product_vertical_slice.md`](docs/product_vertical_slice.md) for the input
contract, RNA merging, abstention, and the upstream integration boundary.

## Status

- The **prioritize stage** and **pipeline orchestration** (config, stage contracts, provenance,
  resume, `doctor`) are implemented and unit-tested; `pytest` runs green on any machine.
- The reads-level stages wrap external tools and are validated end-to-end on a Linux host with the
  reference bundle installed — that run is the release-acceptance gate, not a laptop claim.

## Design

The full architecture and scope decisions are documented in
[`docs/superpowers/specs/2026-07-17-epicurus-v0-full-pipeline-design.md`](docs/superpowers/specs/2026-07-17-epicurus-v0-full-pipeline-design.md).

## License

MIT — see [`LICENSE`](LICENSE).
