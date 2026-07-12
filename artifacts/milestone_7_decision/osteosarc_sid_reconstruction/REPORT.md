# osteosarc.com (Sid) — public reconstruction REPORT

> Reconstructed from the public osteosarc.com site (182 variant pages + VAF TSVs) plus the local public pVACtools/RSEM/Hudson inputs, per the frozen preregistration (`docs/superpowers/specs/2026-07-12-osteosarc-sid-reconstruction-preregistration.md`). No model is fit, tuned, or compared here. Supersedes the `dd3efd1` diagnostic (assumed-negatives + single-positive AUROC — descriptive only).

## 1. Evidence-graded counts

- Unique site variants: **182**  (vaccine-targeted **44**, site-ELISPOT-positive **14**).
- Peptide blocks (long vaccine peptides): **94**;  assay-ledger rows: **128** (68 real experiment rows + 60 UNTESTED vaccine peptides).
- Site ELISPOT tests by resolution: individual-peptide **0**, long-peptide **29**, pool **39**.
- Site positives: strong **7**, weak **15**, unqualified **10**;  negatives total **36** (defensible individual/long-peptide **12**, pool-only **24**);  ambiguous **0**.
- Hudson IFNγ/TCR stream (SEPARATE modality): **15** clonotype rows across **5** (timepoint,mutation) tests; mutation-specific recognized genes **['ASPM', 'DYNC1H1', 'MAP2']**.
- Contradictions (same peptide positive *and* negative across protocol/timepoint): **10** — reported, never collapsed (see AUDIT.json).

## 2. Do the 14 site-ELISPOT-positive variants overlap the Hudson positives?

Overlap (by site coordinate) = **['DYNC1H1-chr14-101980529']**. The two streams are different assays (peptide ELISPOT vs IFNγ/TCR expansion); overlap is reported without asserting equivalence.

## 3. True evidence-supported denominator for ranking today

The defensible within-patient evidence set is NOT the 21 curated mutations with the other 20 called negative. Excluding pool-only and untested rows leaves **29 assay rows** (17 positive, 12 negative) across **15 unique (variant, peptide) units**. At the unit level, 13 are ever-positive and 10 are ever-negative, with **8 in both groups** across protocol/timepoint. Therefore this is an evidence ledger, not a clean binary reranker denominator; ordinary AUROC requires a separately frozen longitudinal label policy.

## 4. Where recognized targets are lost (reachability)

| recognized gene | first failure stage in the automated funnel |
|---|---|
| ASPM | `pvactools_2025_candidate` |
| DYNC1H1 | `reached_shortlist` |
| MAP2 | `pvactools_2025_candidate` |

Adjudications:

- **ASPM G2179R** — ASPM p.G2179R (Hudson-recognized May+Aug) = site variant ASPM-chr1-197102716 (p.Gly2179Arg), called by DRAGEN/Sarek/oncoanalyser, vaccine-included and site-ELISPOT-listed, but ABSENT from the pVACtools 2025.01 candidate universe (pVACtools has no ASPM). Lost at candidate generation, NOT variant calling — corrects dd3efd1's 'off-callset' claim.
- **DYNC1H1 V314I** — DYNC1H1 p.V314I (recognized May+Aug) = site variant DYNC1H1-chr14-101980529 (p.Val314Ile); present in pVACtools 2025 AND curated-21 -> reaches the shortlist (the only recognized target Epicurus could rank).
- **MAP2 GYCVFNKYTV868FS** — MAP2 frameshift (Hudson label p.GYCVFNKYTV868fs, recognized May): the site carries two overlapping frameshift annotations ['MAP2-chr2-209694768(p.Leu867fs)', 'MAP2-chr2-209694772(p.Gly868fs)'] 4bp apart. They have different genomic alleles and different neo-frame sequences, so they remain distinct records. The vaccine (JLF V1/V2/V3) + ELISPOT-strong record is on MAP2-chr2-209694768 (p.Leu867fs); the Hudson '868fs' label position-matches MAP2-chr2-209694772 (p.Gly868fs, no vaccine/experiments). The Hudson neo-frame GYCVFNKYTV differs from the Leu867fs vaccine neo-frame (RVVPFTKAL), supporting the Gly868fs mapping. BOTH are called by DRAGEN/Sarek/oncoanalyser and absent from pVACtools 2025 candidate universe -> lost at candidate generation. The prior 'low-TPM expression-filter drop' attribution is not confirmable from these files (MAP2 is simply absent from the pVACtools output).

## 5. Changes justified for Epicurus (proposed, NOT fit to this patient)

1. **Multi-caller / longitudinal union at candidate generation.** 2 of 3 Hudson-recognized neoantigens (ASPM, MAP2 Gly868fs) were *called by DRAGEN/Sarek/oncoanalyser* yet dropped by the single pVACtools 2025.01 candidate step. The recoverable loss is candidate **recall**, upstream of any ranker.
2. **No hard TPM/tier drop; carry low-evidence candidates with a flag** rather than filtering them out before ranking.
3. **Evidence tiers + honest abstention + diversity/uncertainty in the top-20**, since the recognition axis is unobserved for most candidates.

Each must be validated on the independent labeled cohorts (multimer/Gartner/IMPROVE/CheckMate), never on this single patient.

## Provenance

- 185 fetched/hashed URLs (see PROVENANCE.json); rerun is network-free from `data/raw/osteosarc/site_cache/`.

