# Miller Hu_287 — LOCKED reconstruction-method preregistration (outcome-isolated)

**Frozen 2026-07-12.** This fixes every tool, version, threshold, and deviation so no reconstruction choice
can be tuned to the recognition outcomes.

**Honest isolation statement (NOT a pristine-blinding claim).** The Miller label table has ALREADY been
ingested in this project, so recognized-mutation identities may be known to us. This is therefore
**outcome-isolated / label-VALUE-unused execution**, not prospective blinding: the upstream code, tools, and
every threshold in this document are frozen and committed BEFORE the run, and the run does not READ any
label value (POSITIVE/NEGATIVE, cytokine, or which mutation was recognized) for HLA typing, alignment,
somatic calling, annotation, candidate generation, filtering, expression, or ranking. File access is
audited mechanically where possible (the reconstruction code path does not open
`miller_recognition_labels.csv`; only the final scoring join does). A resulting 3/3 is **locked one-shot
external reconstruction evidence**, explicitly NOT pristine prospective blinding. No parameter is tuned on
Hu_287.

## 1. Original Miller method (from the paper, for equivalence — Tempus xE/xO)
- Tumor/normal **exome** aligned to **GRCh38** with **SpeedSeq Align 0.1.0**; somatic variants via
  **SpeedSeq Somatic**; annotation **SnpEff 4.3i**.
- **RNA** aligned with **STAR 2.4.1c**; RNA alternate-allele observations via **FreeBayes 0.9.21**.
- Candidate variant filters: **exonic VAF ≥ 2%**, **tumor/normal alt-frequency ratio ≥ 1**, **RNA alt
  observations ≥ 1**, **tumor VAF ≥ 5%**, **normal VAF ≤ 5%**, **depth ≥ 10 in tumor AND normal**.
- Ranking: **RNA alt-frequency descending, gene TPM descending, tumor genotype rank ascending**.
- HLA typing **OptiType 1.3.1**; **NetMHCpan 4.0** used for 9–11mer comparisons.

## 2. Modern reproducible equivalent (FROZEN here; deviations enumerated in §4)
- **Reference:** Ensembl GRCh38 primary assembly, release 110 (`.fa` + `.fai` + BWA index sentinels
  `.amb/.ann/.bwt/.pac/.sa` + GATK `.dict`).
- **Alignment:** BWA-MEM (brew `bwa`) for tumor + normal exome with **explicit read groups carrying
  distinct sample names** `SM=Hu_287_T` (tumor) and `SM=Hu_287_N` (normal) → `samtools sort` (coordinate) →
  `samtools index` → GATK4 `MarkDuplicates`. (Replaces SpeedSeq Align.)
- **Somatic calling:** GATK4 **Mutect2** in matched tumor-vs-normal mode, passing `--normal Hu_287_N`, with
  the GRCh38 FASTA + `.fai` + sequence `.dict` → `FilterMutectCalls` → PASS variants only.
  **GATK 4.5 + Java compatibility:** GATK 4.5 targets Java 17; this machine has Java 22 (+ Java 8). Java 22
  will be tested first and, if Mutect2 misbehaves, **JDK 17 (`brew install openjdk@17`) will be installed and
  used** (recorded as a runtime note, not a method change). The `gatk` wrapper needs `python` on PATH — a
  `python`→venv shim is provided for the run.
  **Germline resource / PoN:** the correct-provenance concern is **specificity** — omitting a germline
  resource/PoN mainly **increases false positives** (does not merely reduce sensitivity). Broad/gnomAD
  resources use a **different (Broad `hg38`) contig naming** than the Ensembl primary assembly used here, so
  they are **NOT mixed in**. Either a reference-matched af-only resource is acquired, or — **pre-registered
  here as the default** — Mutect2 runs **matched-normal-only** (no germline resource/PoN), with the resulting
  higher-FP rate explicitly disclosed and partially mitigated by FilterMutectCalls + the tumor/normal VAF
  and depth filters in §3. (Replaces SpeedSeq Somatic.)
- **Exome calling-region policy (frozen BEFORE execution):** the Tempus xE/xO capture BED is not published
  in the SRA metadata. Policy, in order: (1) if a capture BED is recoverable label-blind from study/SRA
  metadata, use it; (2) else restrict calling to a **label-blind coding-exon interval list** derived from
  the GENCODE v44 CDS (a `.interval_list`), which is independent of Hu_287's variants; (3) else call over
  the whole reference. The choice actually used is recorded in the run provenance. This is fixed now so the
  region cannot be chosen after seeing calls.
- **Annotation:** Ensembl VEP (release 110, GRCh38) picking the canonical transcript per variant; if VEP
  cannot be installed on this platform, SnpEff 5.x GRCh38 is the pre-registered fallback. (Replaces
  SnpEff 4.3i.)
- **Shared BASE somatic-variant universe (NOT a shared peptide universe).** The four arms do **not** all
  share one peptide set — that would be wrong. What is frozen as identical is the **base PASS somatic-variant
  set** and its per-variant eligibility evidence (VAF/depth/consequence/RNA). From that base set:
  - `pvac_prime`: **standard pVAC-style generation** (its own peptide set).
  - `lossless_prime` and `lossless_epicurus`: the **exact same lossless peptide universe**
    (`event_b.lossless_peptide_generation`, `expected=None`, full class-I 8–11mers spanning each mutation;
    NOT the paper's 349 IPV-prefiltered set). These two arms MUST receive byte-identical candidate rows.
  - `full_epicurus`: the lossless universe with Epicurus's **frozen gates/portfolio** applied.
  Never hand pVAC candidates to the lossless arms, and never call two different peptide universes
  "identical". Differences are **attributed separately** to (i) variant calling, (ii) generation
  reachability, (iii) gating survival, and (iv) ranking — via the additive decomposition
  `full_epicurus − pvac_prime = (generation) + (scorer) + (selection)` already enforced by
  `benchmark.four_arm`.
- **RNA mutant evidence:** align tumor RNA with a splice-aware aligner (STAR if installable, else HISAT2)
  and count mutant-allele reads at each somatic site (`bcftools mpileup`/pileup). Gene TPM already
  reconstructed via salmon (label-blind). (Replaces STAR + FreeBayes for alt observations.)
- **HLA typing:** OptiType (class-I A/B/C, 4-digit) from the normal exome; if OptiType/razers3 has no
  working osx-arm64 build, the pre-registered fallbacks in order are **arcasHLA**, then **T1K**. If none
  runs on this platform, HLA is `NOT_EVALUABLE` and the class-I-restricted analysis (§5 secondary) is
  reported as NOT_EVALUABLE — the primary mutation-level endpoint (§5) does not require patient HLA for
  reachability but does for PRIME/EL ranking; that dependency is stated in §5.
- **Scoring:** genuine **PRIME 2.1 + MixMHCpred 3.0** (local binaries) and **frozen Epicurus v0.1**
  (`configs/frozen/epicurus_v0_1.json`), expression **confidence-only** per
  `configs/frozen/expression_policy_v1.json`. (Replaces NetMHCpan 4.0 as the presentation scorer; no model
  retuning.)

## 3. Candidate-universe filters (FROZEN, applied identically to ALL FOUR arms)
Applied to the shared base somatic variant set that seeds enumeration for every arm (mirrors the paper
where sensible):
`tumor VAF ≥ 0.05`, `normal VAF ≤ 0.05`, `tumor depth ≥ 10`, `normal depth ≥ 10`,
`tumor/normal alt-frequency ratio ≥ 1`, coding/exonic consequence. Expression/RNA-alt is recorded as an
**evidence annotation** (RNA alt obs ≥ 1 flags "expressed"), NOT a hard universe filter, so unexpressed
recognized mutations are not silently dropped before the reachability endpoint. Class-I peptide lengths
**8–11** (superset of the paper's 9–11). Any peptide with a non-standard residue or invalid length is
`NOT_EVALUABLE`, never a negative.

## 4. Frozen deviations from the paper (each with rationale, decided now)
| paper | ours | rationale |
|---|---|---|
| SpeedSeq Align | BWA-MEM + samtools + GATK4 MarkDuplicates | SpeedSeq is unmaintained; BWA-MEM is the reproducible standard |
| SpeedSeq Somatic | GATK4 Mutect2 + FilterMutectCalls | maintained, matched-normal, widely validated |
| SnpEff 4.3i | Ensembl VEP r110 (SnpEff 5.x fallback) | canonical-transcript consequences; VEP is the pVAC-standard |
| STAR + FreeBayes RNA alt | STAR/HISAT2 align + bcftools mpileup pileup | equivalent mutant-allele counting |
| NetMHCpan 4.0 | genuine PRIME 2.1 + MixMHCpred 3.0 | the benchmark's frozen presentation scorer for both PRIME arms |
| SpeedSeq germline handling | Mutect2 matched-normal-only (default) or reference-matched af-only if acquired | no germline resource ⇒ **more false positives (specificity)**, not merely less sensitivity; Broad/gnomAD contig naming is incompatible with the Ensembl primary assembly, so it is not mixed in |
Any FURTHER deviation forced at runtime (e.g. a fallback aligner/caller/HLA method, or JDK 17) will be
appended here with its reason BEFORE it is used, and committed, never chosen after seeing results.

## 5. Endpoints (FROZEN)
**hits@20 definition (frozen).** Each arm selects **up to 20** peptide-vaccine slots (the final ranked
list; an arm may yield **fewer than 20** rankable/selected candidates after generation/caps/gating — always
report `n_selected` and whether the list **saturated** at 20). `hits@20` = the number of **UNIQUE
lab-recognized mutations represented among the selected slots** — a recognized mutation counts **once** no
matter how many of its peptides occupy slots, and 3 peptides from a single mutation is **1** hit, not 3.
**3/3** means all three distinct recognized mutations each have ≥1 peptide among the selected slots.
**Peptide-level** results (which peptide/HLA filled each slot) are reported alongside, but the headline
count is unique mutations. Never assert "all 20 slots" when caps/generation yield fewer. Duplicate/overlapping peptides: the frozen portfolio/selection keeps its
deterministic tie-break (`md5(mutant_peptide|hla_allele)`); if the frozen selection dedupes or diversifies
across mutations that behavior is used as-is (not re-tuned here). This matches `benchmark.four_arm`, whose
metric granularity is already the mutation (a mutation is a hit iff ≥1 of its candidate rows is selected).

**Frameshift / indel enumeration (honesty, correction).** SNV missense → single mutant residue, fully
enumerable. **Frameshift and inframe indels** require transcript + reading-phase-aware translation of the
altered ORF; the lossless generator enumerates these only where transcript/phase context is available.
Where a consequence cannot be faithfully enumerated (ambiguous phase, missing transcript), those variants
are marked **`NOT_ENUMERABLE`** and reported — never silently treated as if every consequence were
enumerable, and never dropped from the reachability denominator without a stated reason.

## 6. HLA-restriction reporting (Miller 20-mers are NOT HLA-I restricted)
Miller's IFN-γ 20-mer ELISpot labels can be recognized by **CD4/class-II** as well as CD8/class-I, so a
class-I top-20 is **not** the assay's full biological recall. Three reports, clearly separated:
- **(a) HLA-agnostic mutation reachability** — of the lab-recognized mutations, how many are present as
  called somatic variants and survive the frozen §3 base filters (no HLA needed). This is the honest
  upper-level endpoint and is computable even if HLA typing fails.
- **(b) class-I top-20 mechanistic benchmark** — the four-arm hits@20 using patient class-I HLA + PRIME/
  Epicurus, presented explicitly as a CD8/class-I mechanistic view. `NOT_EVALUABLE` if HLA typing fails.
- **(c) limitation** — (b) undercounts CD4/class-II recognition by construction; it is never reported as
  the assay's full recall. IL-5 / `any`-cytokine readouts are sensitivity only, never the headline.
- **Four arms** (`pvac_prime`, `lossless_prime`, `lossless_epicurus`, `full_epicurus`) over the shared
  **base variant** universe (§2). Here n=1 (Hu_287) → a single-patient exploratory point with full
  per-mutation/per-peptide detail, not a cohort claim.

## 7. Isolation & GO/NO-GO
Labels remain sealed (`data/raw/miller_ipv/miller_recognition_labels.csv` is read ONLY at the final scoring
join, after the universe + ranks are frozen). No Hu_287 label value influences any upstream choice; the
reconstruction code path is audited to not open the label file. If a required stage cannot run, its output
is `NOT_EVALUABLE` with a machine-actionable reason; a partial loop yields at most the §6(a) HLA-agnostic
reachability diagnostic, never a class-I headline. This document is frozen; deviations require an appended,
separately-committed entry BEFORE use.
