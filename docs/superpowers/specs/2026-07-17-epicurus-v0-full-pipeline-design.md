# Epicurus V0 — full WES/RNA → neoantigen-portfolio pipeline

_Design spec. Date: 2026-07-17. Status: DRAFT for review._

## 1. Goal

Ship an open-source tool a researcher installs and runs **on their own computer** that
takes a patient's raw sequencing data and produces a ranked, ≤20-candidate personalized
cancer-vaccine neoantigen portfolio — end to end, A to Z:

```
tumor WES FASTQ + normal WES FASTQ + tumor RNA FASTQ
  → ranked top-20 vaccine neoantigen portfolio (+ full provenance)
```

One command:

```bash
epicurus run-pipeline --config patient.yaml --output-dir out/PATIENT-001
```

## 2. Scope and honest non-goals

**In scope (what V0 owns and ships):**
- A single orchestrator that drives every stage from raw reads to final portfolio.
- Contracted, resumable stages with explicit artifact hand-offs and provenance.
- The **Epicurus prioritize stage** (validity gate → calibrated evidence ranking →
  diversity-constrained portfolio → patient-level abstention) — our own, unit-tested code.
- A reproducible environment (container recipe + pinned conda env) so a user obtains all
  external tools, and a reference-data bootstrap step.
- Honest docs: hardware/reference requirements, per-stage tool provenance, and what the
  scores do and do not mean.

**Explicit non-goals (and why):**
- **We do not reimplement variant callers, aligners, or peptide generators.** Somatic
  calling and neoantigen generation are delegated to established, validated tools
  (BWA-MEM2, GATK/Mutect2, VEP, Salmon, an HLA typer, pVACtools + MHCflurry/NetMHCpan).
  Rebuilding these would be reckless for a clinical-adjacent tool and adds no value.
- **We do not claim a recognition breakthrough.** The research record (M6–M8) is explicit:
  Epicurus's *ranking* sits at parity with genuine PRIME; its demonstrated edge is
  portfolio diversification + full-evidence routing, not a superior immunogenicity model.
  Docs state this plainly. The pipeline's value is being *complete, reproducible, honest,
  and locally runnable*, with a principled final prioritizer.
- **Not required to run on the developer Mac.** It targets the end user's machine
  (Linux workstation/server, or Docker on any host). Orchestration + the Epicurus stage
  are unit-tested on the dev Mac with fixtures/mocked tool calls; the first true
  FASTQ→portfolio run is validated on a real Linux host.

## 3. Runtime requirements (documented, inherent)

- **OS:** Linux, or any host via the provided container (Docker/Apptainer).
- **Reference data:** GRCh38 genome + GATK resource bundle + VEP cache + Salmon index —
  a one-time download of tens–hundreds of GB via a `epicurus fetch-references` helper.
- **Compute:** a 30× tumor/normal WES pair is hours of CPU; documented, not hidden.
- These are intrinsic to somatic neoantigen calling, not artifacts of our design.

## 4. Architecture

### 4.1 Stage graph

Each stage is a pure function `run(inputs, config, workdir) -> StageResult` that (a) checks
its external tool is available and fails loudly if not, (b) shells out to that tool with a
frozen, logged command, (c) writes a declared artifact, (d) records provenance
(tool version, command, input hashes). Stages are resumable: an existing valid artifact is
reused unless `--force`.

| # | Stage | Tool (external) | Input → Output |
|---|-------|-----------------|----------------|
| 1 | `align` | BWA-MEM2 + samtools/GATK MarkDuplicates | FASTQ → tumor/normal BAM |
| 2 | `call` | GATK Mutect2 + FilterMutectCalls | BAMs → filtered somatic VCF |
| 3 | `annotate` | Ensembl VEP | VCF → annotated VCF |
| 4 | `express` | Salmon | RNA FASTQ → transcript TPM |
| 5 | `hla` | OptiType (class I) / arcasHLA | normal BAM/FASTQ → HLA alleles |
| 6 | `generate` | pVACtools (pVACseq) + MHCflurry (NetMHCpan optional) | annotated VCF + HLA + TPM → candidate table |
| 7 | `prioritize` | **Epicurus (in-repo)** | candidate table → scored + gated + ≤20 portfolio |
| 8 | `report` | in-repo | portfolio CSV + JSON summary + human report + provenance manifest |

Stages 3/4/5 are independent and may run in any order (or in parallel) once stage 2 (and
raw RNA) exist; stage 6 joins them.

### 4.2 The Epicurus prioritize stage (existing, tested)

`epicurus_neo.product` already implements this and stays the differentiator:
- `normalize_product_candidates` — canonicalize a pVACseq-style table to the product schema.
- `merge_rna_evidence` — attach RNA expression as **confidence-only** evidence
  (frozen policy: expression never acts as a rank penalty — `configs/frozen/expression_policy_v1.json`).
- `score_product_candidates` — combine translated / presentation / recognition / coverage
  evidence into a transparent evidence-prioritization score (not a validated response probability).
- deterministic validity gate (`apply_deterministic_gate`) — lost-HLA routes and unexpressed
  genes cannot consume top-20 slots; rejected rows keep exact reason codes.
- portfolio selection — `k`, `max_per_mutation`, `max_per_gene`, `max_per_hla` diversity
  constraints; deterministic tie-break on `md5(mutant_peptide|hla_allele)`.

### 4.3 Configuration

`patient.yaml`:
```yaml
patient_id: PATIENT-001
inputs:
  tumor_wes:  [tumor_R1.fastq.gz, tumor_R2.fastq.gz]
  normal_wes: [normal_R1.fastq.gz, normal_R2.fastq.gz]
  tumor_rna:  [rna_R1.fastq.gz, rna_R2.fastq.gz]
  hla_alleles: null          # optional override; else typed in stage 5
references:
  bundle_dir: ~/.epicurus/references/GRCh38
generate:
  predictors: [MHCflurry]     # NetMHCpan added if licensed+installed
  epitope_lengths: [8,9,10,11]
prioritize:
  k: 20
  max_per_mutation: 1
  max_per_gene: 4
resources:
  threads: 8
```

### 4.4 CLI surface (product)

- `epicurus run-pipeline --config patient.yaml --output-dir DIR [--from STAGE] [--to STAGE] [--force]`
- `epicurus run-patient ...` — **retained** shortcut: start at stage 7 from an existing
  pVACseq candidate table (the honest A3→Z entry; useful when a user already has candidates).
- `epicurus validate-patient-input`, `validate-schema`, `select-portfolio` — retained utilities.
- `epicurus fetch-references --bundle GRCh38 --dest DIR` — reference bootstrap.
- `epicurus doctor` — checks each external tool + reference presence, prints a readiness table.

### 4.5 Distribution

- `Dockerfile` (or Apptainer def) that installs all external tools at pinned versions +
  the `epicurus` package. `docker run epicurus run-pipeline ...` runs anywhere.
- `environment.yml` (bioconda) as the non-container path for Linux users.
- `epicurus doctor` gives a clear, actionable readiness report before a long run.

## 5. Testing strategy

- **Unit (dev Mac, CI):** stage contract logic, config parsing/validation, artifact
  hand-off + resume, provenance manifest, tool-missing failure paths (external tools
  **mocked** — assert on the constructed command, not real execution), and the full
  Epicurus prioritize stage on small fixtures. This is the bulk of the safety net and runs
  green on the Mac.
- **Stage smoke (dev Mac):** any stage whose tool is pip/bioconda-installable and CPU-cheap
  (MHCflurry, Salmon on a tiny index) exercised on a miniature fixture where feasible.
- **End-to-end (Linux host, gated):** one real small WES/RNA pair from raw reads to
  portfolio, validated on a provisioned Linux host. Not part of Mac CI; documented as the
  release-acceptance gate.

## 6. Repository shape after V0 (hyper-clean)

Master ships only the product; the M5–M8 research corpus lives in the local `research-archive`
branch + `v0-archive` tag (already created; not pushed).

**Keep on master:**
- `src/epicurus_neo/` — prioritize stage closure (`product`, `contracts`, `gates`,
  `normalize`, `schema`, `portfolio`, `features`) + new `pipeline/` orchestration package
  (stage runners, config, provenance, references) + slim product `cli.py`.
- Generation-adjacent code productized from the archive **only as needed** (e.g. a cleaned
  MHCflurry feature helper), each with tests. Research pipelines (Miller/Sid/Zhao, v0.2–v0.5
  rankers, decision benchmarks, four-arm harness) stay archived.
- `tests/` — product + pipeline unit tests only.
- `examples/demo_patient/` — one tiny worked example (candidate-table entry; and a
  reads-level example pointer/fixture where size allows).
- `docs/` — a rewritten product README, a pipeline guide, and this spec.
- `LICENSE` (MIT), `pyproject.toml` (deps + `[pipeline]`/`[mhc]` extras), `Dockerfile`,
  `environment.yml`.

**Remove from master** (recoverable at `v0-archive`): `src/benchmark/`, `src/event_b/`,
the research modules in `src/epicurus_neo/` (m6/, plm_*, transfer_*, retrieval_*,
auto_research, decision/benchmark scaffolding), research `scripts/` + `tests/`,
`artifacts/`, `experiments/`, `notes/`, `comparators/`, `outputs/`.

## 7. Build sequence (incremental, each step committed)

1. Repo carve-out: MIT `LICENSE`; rewrite `README.md` to the honest full-pipeline framing;
   remove research surface from master (keep-set per §6); slim `cli.py`; green product tests.
2. `pipeline/` skeleton: config schema + loader, `StageResult`/stage protocol, provenance
   manifest, workdir/resume logic, `doctor`. Unit-tested (no external tools yet).
3. Wire stage 7 (`prioritize`) + stage 8 (`report`) into the orchestrator end-of-chain;
   `run-pipeline --from prioritize` works against a candidate table on the Mac.
4. Stages 6 (`generate`) then 1–5, each as a mocked-command stage runner + provenance +
   tests. Real execution deferred to Linux.
5. Distribution: `Dockerfile`, `environment.yml`, `fetch-references`; docs pass.
6. Push clean master to public origin. (Archive branch/tag stay local.)
7. Linux-host acceptance run (separate, gated) → first real end-to-end portfolio.

## 8. Risks

- **Not testable end-to-end on the dev Mac** — mitigated by heavy unit coverage of
  orchestration + a Linux acceptance gate before any "validated end-to-end" claim.
- **External-tool/reference drift** — pinned versions in the container; `doctor` verifies.
- **Overclaiming** — README states generation is delegated and ranking is at parity with
  PRIME; scores are evidence-prioritization, not validated response probabilities.
- **Large surface** — sequenced so every committed step leaves master shippable (the A3→Z
  prioritizer is always runnable even before the reads-level stages land).
