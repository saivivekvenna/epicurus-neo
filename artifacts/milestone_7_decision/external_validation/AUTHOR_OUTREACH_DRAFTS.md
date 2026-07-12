# Author outreach drafts

Short, tailored requests for **data, not PHI**. Each: states the exact scientific purpose, accepts
controlled access / DUAs, asks only for the [MINIMUM_PATIENT_DATA_PACKAGE](MINIMUM_PATIENT_DATA_PACKAGE.md),
offers to share aggregate results, and **promises no authorship** (offer acknowledgement only). Fill the
`[BRACKETS]` (your name, affiliation, the specific supplementary table numbers once verified) before sending.
Send from an institutional address. These are drafts for the user to send — nothing here is auto-sent.

**Shared purpose paragraph (reuse verbatim):**
> We are building an open, reproducible benchmark that asks a single question: starting from the *same*
> WES/RNA/HLA inputs a clinician would have, does a candidate neoantigen-ranking method place more
> *experimentally recognized* neoantigens in the final top-20 than a standard pVACtools-style pipeline plus
> genuine PRIME? The primary metric is a patient-level paired difference in recognized hits@20. To evaluate
> this without bias we need, per patient, the full pre-selection candidate universe plus explicit
> POSITIVE / TESTED-NEGATIVE / UNTESTED labels — not only the validated hits. We publish only aggregate,
> patient-level statistics; individual-level data stays within its controlled environment.

---

## 1 — Hudson Lab / osteosarc.com (Sid Sijbrandij case)

Subject: Request: tumor-peptide stimulation manifest + pool deconvolution for the osteosarc.com case

> Dear [Hudson Lab / Dr. Hudson],
>
> Thank you for openly publishing the osteosarc.com longitudinal multi-omics and the IFNγ
> peptide-expansion assay. Using only the public inputs, an input-only pipeline recovers all three
> IFNγ-expanded targets (ASPM, DYNC1H1, MAP2) into a top-20 shortlist — but we cannot yet compute a clean
> denominator because the **tested-peptide stimulation manifest** is not public.
>
> [Shared purpose paragraph.]
>
> Could you share, de-identified: (1) the exact tumor-peptide **stimulation pool composition / sample
> sheet** (which peptides were in which stimulation pool at each timepoint — the true tested denominator),
> (2) the **per-peptide deconvolution** from pool-level IFNγ readout to individual peptide outcomes, and
> (3) the positivity threshold used. The public TCR/expander files alone are not the tested denominator.
> This would convert a single-patient recovery demonstration into a properly-denominated benchmark case.
> Happy to work under any data-use terms you prefer; we would acknowledge (not co-author) the contribution.
>
> [Name, affiliation]

---

## 2 — Miller IPV (Sci Transl Med 2024, scitranslmed.abj9905)

Subject: Question on the PRJNA980652 deposit + per-peptide immunogenicity table

> Dear [Dr. Miller / corresponding author],
>
> We admire the individualized-peptide-vaccine immunogenicity screen and note the WES/RNA appear openly
> deposited under BioProject PRJNA980652. [Shared purpose paragraph.]
>
> Two requests: (1) confirmation that the **full per-peptide immunogenicity table** (every 20-mer assayed,
> its HLA restriction where known, and its POSITIVE/NEGATIVE ELISpot call) is the supplementary table
> [S#], and (2) if available, the **full pre-selection candidate list** (the 6,237 detected neoantigens
> before down-selection to the 349 tested), so we can reconstruct the complete denominator. We plan to
> re-enumerate the mutanome from the deposited WES ourselves; a pointer to the exact HLA genotypes per
> patient would save error. Aggregate results only; acknowledgement, not authorship.
>
> [Name, affiliation]

---

## 3 — GBM neoantigen vaccine cohort (GEO GSE237936 / BioProject PRJNA997267)

Subject: Request: matched WES/HLA + full tested-peptide panel for the GSE237936 GBM cohort

> Dear [corresponding author],
>
> Your GBM cohort deposits multi-region tumor RNA-seq (GSE237936 / PRJNA997267) and reports T-cell
> reactivity to roughly half the vaccine peptides. [Shared purpose paragraph.]
>
> To use it end-to-end we would need, de-identified and per subject: matched **tumor/normal WES** (or the
> somatic VCF), **4-digit HLA** genotypes, and the **full tested peptide panel with explicit negative
> calls** (peptides assayed and scored non-reactive), plus the pre-vaccine candidate universe if it exists.
> Controlled access / a DUA is fine. We would publish only aggregate patient-level statistics and
> acknowledge your contribution.
>
> [Name, affiliation]

---

## 4 — Zhao 2026 dendritic-cell recognition dataset

Subject: Availability of the full pre-selection candidate universe + raw inputs

> Dear [Dr. Zhao / corresponding author],
>
> Your dataset is, to our knowledge, the largest mutation-derived neoantigen recognition corpus with
> explicit negatives, and we already use its tested peptides as a recognition-scale training/prior asset.
> [Shared purpose paragraph.]
>
> To move it from a training corpus to an end-to-end benchmark we would need two additions: (1) the **full
> per-patient pre-selection candidate universe** (all enumerated neoantigens, not only the tested subset),
> and (2) the **raw WES/RNA inputs + 4-digit HLA** (controlled access acceptable). With those, the cohort
> could support the generation stage, not only ranking among pre-chosen peptides. Aggregate reporting only;
> acknowledgement, not authorship.
>
> [Name, affiliation]

---

## 5 — EVX-01 barcoded-multimer study (PMC11116868)

Subject: Per-patient candidate universe + explicit negatives from the EVX-01 barcoded-multimer screen

> Dear [corresponding author / Evaxion],
>
> Your barcoded-pMHC-multimer screen of ~145–231 predicted binders per patient is an unusually dense
> per-patient immunogenicity readout. [Shared purpose paragraph.]
>
> We would like to know whether the following are (or can be made) available, de-identified: the **full
> per-patient list of screened binders** with their HLA restriction, the **barcode/pool → peptide
> deconvolution**, the **explicit per-peptide POSITIVE/NEGATIVE table**, and the underlying **WES/RNA/HLA**
> inputs (any controlled-access mechanism is fine). If a patient-level crosswalk (inputs → candidates →
> outcomes) can be shared under a DUA, this could be one of the most informative external test cohorts in
> the field. We publish aggregate statistics only and would acknowledge the contribution.
>
> [Name, affiliation]

---

## 6 — PGV001 + atezolizumab (Nature Cancer 2025; OpenVax / Mount Sinai)

Subject: Explicit per-peptide ELISpot negatives + consent group for phs003922

> Dear [Dr. Bhardwaj / Dr. Gnjatic / corresponding author],
>
> Your PGV001 personalized-vaccine trial pairs OpenVax candidate selection with per-peptide ELISpot, which
> is exactly the structure our benchmark needs. [Shared purpose paragraph.]
>
> Could you confirm: (1) the **explicit per-peptide POSITIVE/NEGATIVE ELISpot counts** and where they are
> tabulated, and (2) the **consent/data-use group** for the deposited inputs (e.g. phs003922) — general
> research vs disease-specific vs non-commercial — since our tool is open-source. We would submit any
> required Data Access Request and report only aggregate results. Acknowledgement, not authorship.
>
> [Name, affiliation]

---

## 7 — 2026 medRxiv "automating neoantigen selection" aggregate authors

Subject: Data release + patient crosswalk for the 21-patient neoantigen-selection dataset

> Dear [authors],
>
> Your preprint describes ~927 training + ~310 development peptides across 21 patients — a valuable
> structured resource. [Shared purpose paragraph.]
>
> We would like to know your **data-release plan**: whether the **per-patient raw inputs (WES/RNA/HLA)**,
> the **full candidate universe**, and a **patient crosswalk** (inputs → candidates → tested outcomes with
> explicit negatives) will be deposited, and under what access model. Until a crosswalk and raw inputs
> exist, we would treat the peptide tables as a development/prior resource rather than an independent
> end-to-end test; a release would let us use it properly. Aggregate reporting only; acknowledgement, not
> authorship.
>
> [Name, affiliation]
