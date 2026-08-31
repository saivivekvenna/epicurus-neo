# Epicurus Neo

**Turn tumor/normal WES and tumor RNA-seq into a defensible, patient-specific cancer-vaccine
portfolio—not just another peptide score.**

Epicurus Neo is an open-source, end-to-end neoantigen prioritization pipeline. It orchestrates the
established genomics stack, then applies its own biological-validity gates, evidence ranking, and
diversity-aware portfolio selection to produce at most 20 vaccine candidates. One command takes a
patient from raw reads to an auditable shortlist:

```bash
epicurus-neo run-pipeline --config patient.yaml --output-dir out/PATIENT-001
```

```
tumor WES + normal WES + tumor RNA  ──▶  ranked top-20 vaccine neoantigen portfolio
```

## Why Epicurus Neo

Most neoantigen tools stop at a long ranked list. Epicurus Neo treats the actual problem as a
constrained portfolio decision: every slot is scarce, duplicate routes to the same mutation crowd out
coverage, and candidates unsupported by the patient's biology should not survive merely because they
scored well in isolation.

That design has produced concrete wins:

- **8 recognized mutations in 20 slots versus PRIME's 1.** In a frozen, label-blind calibration
  readout on Hu_315, Epicurus Neo's mutation-diversified portfolio captured 8 of 18 reachable
  recognized mutations; genuine PRIME captured 1. That is **40% top-20 precision versus 5%**, with no
  duplicate mutation slots. This is a one-patient development result, not a population-level claim.
- **75 → 85 validated top-20 hits on IMPROVE.** A predeclared hierarchical reranking policy added
  10 hits across 70 patients (**+13.3%**) over global Epicurus Neo ranking. The lift survived **100/100
  randomized bracket assignments**.
- **A reranker that knows when not to rerank.** The same tournament mechanism reversed direction on
  Gartner (20 → 17 hits) and the multimer cohort (24 → 23). Epicurus Neo therefore fails closed: the
  guarded policy activates only for its explicit screened-candidate regime and preserves the global
  ranking elsewhere (20/20 Gartner hits and 24/24 multimer hits retained).

The takeaway is deliberately narrower—and more useful—than “we solved immunogenicity.” Epicurus Neo's
edge is **patient-level prioritization**: validity-aware routing, mutation-level diversification, and
regime-aware abstention that converts strong component scores into a better 20-slot decision.

## How it works

Epicurus Neo **orchestrates a complete pipeline** and owns the final prioritization. It does **not**
reimplement variant callers or peptide generators — those are established, validated tools that
Epicurus Neo drives:

| Stage | Tool | Does |
|-------|------|------|
| align | BWA-MEM2 + MarkDuplicates | FASTQ → tumor/normal BAM |
| call | GATK Mutect2 + FilterMutectCalls | BAMs → filtered somatic VCF |
| annotate | Ensembl VEP | annotate variants |
| express | Salmon | RNA transcript TPM |
| hla | OptiType / arcasHLA | class-I HLA typing |
| generate | pVACtools + MHCflurry | candidate neoantigen table |
| **prioritize** | **Epicurus Neo** | validity gate → calibrated ranking → ≤20 portfolio |
| report | Epicurus Neo | portfolio CSV + JSON summary + provenance |

Epicurus Neo owns the part shown in bold: the `prioritize` stage. It combines a deterministic validity
gate (lost-HLA routes and unexpressed genes cannot take a slot), a transparent score spanning
translation, presentation, recognition, and expression evidence, mutation- and gene-level diversity
constraints, and patient-level abstention. Every run emits both the portfolio and its provenance.

The scores are evidence-prioritization scores, **not validated probabilities of vaccine response**.
Ranking performance is cohort-dependent, and the benchmark results above include development
evidence. Untouched, multi-patient external validation is still required before making a clinical or
general superiority claim.

## Requirements

Neoantigen calling from raw reads is heavyweight; these requirements are inherent, not incidental:

- **OS:** Linux, or any host via the provided container.
- **Reference data:** GRCh38 + GATK resource bundle + VEP cache + Salmon index — a one-time
  download (tens–hundreds of GB). Scaffold the bundle and get exact sources with
  `epicurus-neo references --dest ~/.epicurus-neo/references/GRCh38`.
- **Compute:** a 30× tumor/normal WES pair takes hours of CPU.

Check your machine is ready before a run:

```bash
epicurus-neo doctor --bundle-dir ~/.epicurus-neo/references/GRCh38
```

## Install

```bash
# Container (recommended — brings every external tool):
docker build -t epicurus-neo .
docker run --rm -v "$PWD":/work epicurus-neo run-pipeline --config /work/patient.yaml --output-dir /work/out

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
  bundle_dir: ~/.epicurus-neo/references/GRCh38
prioritize:
  k: 20
  max_per_mutation: 1
  max_per_gene: 4
```

```bash
epicurus-neo run-pipeline --config patient.yaml --output-dir out/PATIENT-001
```

Every stage is resumable: re-running reuses valid cached artifacts. Use `--start`/`--stop` to run a
subrange and `--force` to recompute.

## Already have candidates? Start at the ranking stage

If you have already run pVACseq, skip the front-end and prioritize directly:

```bash
epicurus-neo run-patient \
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

## License

MIT — see [`LICENSE`](LICENSE).
