# Acquisition packet — a genuinely NEW denominator-level external cohort

Goal: obtain an independent cohort with a **full/reconstructable candidate denominator + patient HLA +
expression + VAF + genuine POSITIVE and TESTED_NEGATIVE labels**, to test frozen Epicurus v0.1 vs
MHCflurry/NetMHCpan/genuine PRIME at within-patient top-20 on patients **no part of the pipeline has
touched**. Ranked by (independence × denominator-defensibility × executability now).

---

## PRIMARY (controlled access, slow) — Hu/Ott NeoVax melanoma via dbGaP `phs001451`
Full machine-readable manifest already exists: **`configs/source_manifests/neovax_reconstruction.yml`**
(checksum-pinned local sources + per-element status). This packet is the actionable request.

**Exact accession / who to request from**
- dbGaP study **phs001451**; versions **`phs001451.v1.p1`** (Ott 2017, patients 1–6) and
  **`phs001451.v3.p1`** (Hu 2021, patients 1–6, 11, 12). Data types: tumor WES, matched-normal WES, tumor
  RNA-seq. Access level: **CONTROLLED**.
- Mechanism: eRA Commons account → **Data Access Request (DAR)** to the study's Data Access Committee →
  signed **Data Use Certification (DUC)**. Individual-level data may NOT be redistributed; only aggregate
  summary stats may be published (confirm against the DUC).
- Papers: Ott PA et al. Nature 2017;547:217 (PMC5577644); Hu Z et al. Nat Med 2021;27:515 (PMC8273876).

**What is already local (no download):** somatic SNV/indel identities for Hu patients 11–12
(`data/raw/hu_melanoma_2021/manual/NIHMS1707651-supplement-Suppl_DataSet.xlsx`), vaccine epitopes +
post-vaccine ELISpot calls (patients 1–6, 11–12). **Not local / must be regenerated:** per-mutation VAF,
tumor RNA TPM, complete HLA genotype, and the full candidate universe (never published).

**Reconstruction pipeline once WES/RNA land (destination `data/raw/hu_melanoma_2021/regenerated/`,
`data/raw/ott_melanoma_2017/regenerated/`):**
1. HLA-type class I+II from WES (e.g. OptiType / HLA-HD) → `hla_genotype.tsv` (`patient\talleles`).
2. Somatic calling on tumor vs matched-normal WES → per-mutation VAF/read-depth.
3. Quantify tumor RNA-seq → per-gene/transcript TPM.
4. Candidate universe: mutations → all class-I 8–11mer neoepitopes (pVACtools/NetMHCpan) restricted to the
   patient's alleles = the **denominator**.
5. Join labels: vaccine-included + ELISpot-reactive → POSITIVE; screened non-reactive → TESTED_NEGATIVE;
   everything else in the universe → UNTESTED (never coerce to negative).
6. Emit the unified schema and run `event_b.prime_transfer.external_validate(frame)`.

**Unified schema expected by `external_validate`** (one row per candidate):
`patient_id, mutant_peptide, hla_allele, label∈{POSITIVE,TESTED_NEGATIVE,UNTESTED}, prime(%rank via tool),
el(NetMHCpan-EL %rank), expr(TPM or decile)` (+ optional `cmp_*` comparator columns).

**Blocker status:** BLOCKED_CONTROLLED_ACCESS. Timeline is DAR-approval-bound (weeks). This is the only path
to a *fully independent, denominator-complete, VAF+RNA* external cohort — worth starting now in parallel.

---

## OPEN EXECUTABLE ALTERNATIVES (ranked by independence × denominator-defensibility × executability)

**Honest bottom line: no open, executable, denominator-defensible NEW external cohort exists locally today.**
Every candidate fails on at least one of {ground-truth labels, real denominator, open access}. Details:

| Rank | Cohort | What's local / accessible | Fatal gap for a frozen top-20 validation | Verdict |
|---|---|---|---|---|
| 1 | **TNBC 2024 vaccine** (dbGaP `phs002787` — *accession named by requester; NOT verifiable from the repo, confirm it maps to the intended TNBC vaccine study before committing*) | **Nothing local** (repo-wide grep for `phs002787`/`tnbc` returns only the handoff note). Published supplement *may* list 18 pts / 198 candidates / 45 immunogenic (peptide, HLA, immunogenic label). WES+RNA controlled. | The 198 are the **vaccinated-selected** subset (~11/pt), NOT a full denominator; per-candidate expression/VAF need the controlled WES/RNA → frozen residual (needs `expr`) cannot run; only PRIME/EL peptide-HLA arms could. | Recognition-scale, not a denominator. **Not defensible now.** Acquire supplement to at least run PRIME/EL top-k; full run needs phs002787 (controlled). |
| 2 | **"Osteosarc" n=1** (`artifacts/osteosarc_audit/class_I.tsv` — **byte-identical (md5 751d16…) to the Sijbrandij report** `data/raw/sijbrandij/DNA_SR24-58221_C1_…tsv`; + `somatic.vcf.gz`, 1,213 MuTect variants) | FULL candidate universe for 1 patient with HLA, SHERPA/NetMHCpan ranks, DNA/RNA allelic fraction, TPM, *predicted* immunogenicity score. | **Same single patient as Sijbrandij** (not a separate cohort); **no T-cell assay label column** (prediction-only) and **n=1**. | Denominator present, **labels absent, n=1. Not usable.** Product/demo input only. |
| 3 | **Müller/neoranking NCI** (`neoranking_corpus`) | 56 pts full universe, HLA + expression, `Score_EL`. | **No genuine TESTED_NEGATIVE** (VALIDATED=0 is UNTESTED) → PU-only; and **same NCI lineage as Gartner** (not independent). | PU top-k only, not independent. **Development lineage, not external.** |
| 4 | Zhao DC / vaccine Event-B cohorts (Braun, mKRAS, PDAC-NeoVax, NOUS-209) | vaccinated-subset recognition rows. | No per-patient candidate denominator (vaccinated peptides only). | Recognition-scale. **Not a denominator.** |

**Conclusion for the north star:** the binding constraint is confirmed empirically — a genuinely new,
denominator-complete, independently-labeled cohort is **not** obtainable from open/local files. It requires a
controlled-access acquisition. The two realistic controlled paths are **Hu/Ott `phs001451`** (primary packet
above; full WES+RNA, denominator reconstructable) and **TNBC `phs002787`** (18 TNBC pts). Both are DAR-gated.
Recommended parallel action: submit DARs for both; meanwhile pull the TNBC published supplement to check
whether it exposes enough (peptide, HLA, immunogenic label) for at least a PRIME/EL-only external top-k while
the WES/RNA remains controlled.
