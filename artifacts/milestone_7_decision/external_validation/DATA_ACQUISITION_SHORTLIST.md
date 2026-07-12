# Neoantigen decision-problem cohort acquisition shortlist

> **SUPERSEDED FOR EXECUTION (2026-07-12) →** the canonical, primary-source-re-verified execution ranking,
> tracker, data-request schema, and outreach drafts now live in
> [`ACQUISITION_EXECUTION_PLAN.md`](ACQUISITION_EXECUTION_PLAN.md) +
> [`COHORT_ACQUISITION_TRACKER.csv`](COHORT_ACQUISITION_TRACKER.csv). This file is retained for its
> **lineage analysis** (the field = 3 lineages). Corrections since it was written: CheckMate 153 raw seq is
> **EGA-controlled `EGAD00001011302`, not dbGaP**; **Miller "IPV" = identify-prioritize-validate platform**
> (not a vaccine); **EVX-01 (PMC11116868) is closed data**; GBM `GSE237936` is RNA-only n=4 (2024). See the
> plan's "Verification corrections" section.

_Compiled 2026-07-11 from a 5-agent parallel web scout (vaccine trials · TIL/TCR screens · immunogenicity
databases · controlled-access repos · 2024–2026 newest/HT), merged + deduped + leakage-adjudicated against
a prior deep-research pass. Two linchpin facts independently re-verified (PRJNA980652 open; CheckMate 153 real)._

## The bar (a qualifying cohort needs ALL four)
1. **Denominator** — many candidates screened per patient (not a handful of validated hits).
2. **Explicit tested-negatives** — peptides assayed for T-cell recognition and scored NEGATIVE are recorded (not "untested"). ← the hard filter.
3. **Recognition readout** — ELISpot / ICS / multimer-tetramer / TIL reactivity (functional), not predicted or MS-eluted presentation alone.
4. **Recoverable inputs** — mutation calls (WES), RNA/expression, 4-digit+ HLA (controlled access OK if documented).

## The load-bearing finding — the field is 3 lineages, and we hold or are beaten by all three
- **NCI / Rosenberg Surgery Branch** (TMG/minigene whole-mutanome screens: Parkhurst, Tran, Gros, Cafri, Lowery, Zacharakis, Leidner, Gustafson 2025) = **our Gartner lineage**, and partly inside PRIME's training set.
- **DTU / Hadrup** (barcoded pMHC-multimer baskets: Kristensen, 26-pt melanoma) = **our IMPROVE cohort**.
- **Lausanne / NeoDisc / HiTIDE** (Bassani-Sternberg/Gfeller) = **inside PRIME 2.0's training set** (the Müller/Arnaud *Immunity* 2023 harmonized set IS PRIME's training data).
- IEDB / CEDAR / NEPdb / dbPepNeo2 / TSNAdb / Neodb / TumorAgDB = training-expansion or positives-only; **IEDB `tcell_full_v3` (downloaded 2021-03-27) literally IS PRIME's negative set**.
- HTAN / ICGC-ARGO / PCAWG / iAtlas / GENIE = genomics without a per-candidate recognition screen → fail #2/#3.

**Consequence:** genuinely-independent, negatives-bearing, inputs-recoverable cohorts that PRIME has NOT
touched are **few and small**. They are the Tier-A list below.

---

## TIER A — new, independent, low-leakage PROOF candidates (rank order)

| # | Cohort | Cancer | n (functional) | Denominator | Explicit negatives | Assay | Inputs | Access | Leakage | Note |
|---|--------|--------|----:|---|---|---|---|---|---|---|
| **A1** | **CheckMate 153 / Alban, Riaz** (Nat Med 2024, `s41591-024-03240-y`, PMC12066197) | NSCLC (nivolumab) | ~14 (of 80 seq'd) | **1,453 candidates** — Gartner-density | ~196 recognized / **~1,257 negatives** (Supp Table 2) | combinatorial tetramer | WES+RNA+HLA (Polysolver) | **supp OPEN now**; raw seq likely dbGaP/BMS-controlled (accession UNVERIFIED) | **LOW** (Oct 2024, post-cutoff); independent | Largest independent denominator found. Tetramer = binding+recognition. Verify raw-seq accession. |
| **A2** | **Miller IPV** (Sci Transl Med 2024, `abj9905`) | 8 solid types | **13** | 6,237 detected → 349 tested (~27/pt, IPV-prefiltered) | ~199 pos / **~555 negatives** (754 20-mers, supp) | IFN-γ ELISpot | WES+RNA+HLA | **FULLY OPEN — BioProject `PRJNA980652` (SRA), no DAR** ✓verified | **LOW**; independent (La Jolla/UCSD) | **Fastest path — days, free.** Small n; denominator IPV-prefiltered (re-enumerate full mutanome from WES to widen). |
| **A3** | **WashU TNBC neoantigen DNA vaccine** (Genome Med 2024, `s13073-024-01388-3`, PMC11562513) | Triple-neg breast | 18 | 198 tested (~11/pt, vaccine-slate); full candidate list in Table S1 | 45 pos / **153 negatives** | ELISpot+ICS+tetramer+TCR-seq | WES+RNA+HLA+TCR (deposited) | **dbGaP `phs002787`** (controlled DAR); supp labels open now | **LOW** (Nov 2024); independent; **pVACtools** (our open lineage) | New cancer type (TNBC). Vaccine range-restricted denominator (see caveat). |
| **A4** | **TESLA** (Wells, Cell 2020; Synapse `syn21048999`) | 3 melanoma + 3 NSCLC | 6 | 608 peptides (~97/pt) | 37 pos / **571 negatives** | pMHC multimer | WES+RNA+6-digit HLA | Synapse governed (DUC) — downloadable | **MED** — field's most-recycled benchmark; partial ingestion by 2021+ predictors | Downloadable now; pre-register leakage caveat; tiny n. |
| **A5** | **PGV001 + atezolizumab** (Nat Cancer 2025, `s43018-025-00966-7`; PMC handled) | Urothelial | 10 | ~10/pt (vaccine-slate) | Y (per-peptide ELISpot); **exact count UNVERIFIED** | ELISpot | WES+RNA+HLA | **dbGaP/EGA `phs003922`** (controlled) | **LOW** (2025); independent; **OpenVax** | Verify negative counts in supplement before committing. |
| **A6** | **Southampton NSCLC immunopeptidomics** (npj Prec Onc 2026, `s41698-026-01539-2`; bioRxiv 2024.05.30.596609) | NSCLC | ~6 | ~70 MS-presented peptides tested | ~9 pos / **~61 negatives** (recognition-neg among presented) | IFN-γ | WES+RNA (**EGA `EGAS00001005499`**, controlled); MS at PRIDE `PXD028990` | supp open (Zenodo 10.5281/zenodo.12820423); inputs gated | **LOW** (2024/26); independent | Underpowered; independent external *test* only, not training. The prior pass's lead — confirmed real. |
| **A7** *(probe)* | **HANSolo** (Sci Adv 2024, `sciadv.ado6491`; Kula/Elledge) | solid tumors | UNVERIFIED | full-HLA genetic minigene library | inherent (screened vs not) | genetic functional screen | recoverable UNVERIFIED | data-availability UNVERIFIED | **LOW** (2024); independent | High-upside emerging platform — verify cohort size + open per-antigen table. |

## TIER B — development / input-completion (overlaps our labels; PRIME did NOT train on the genomics)

| # | Cohort | What it is | Value |
|---|--------|-----------|-------|
| **B1** | **NCI Prospective Procurement — dbGaP `phs001003.v2.p1`** (145 pts WES+RNA) | The **deposited input half of the Gartner/NCI functional labels we already hold** (Parkhurst/Lowery/Zacharakis all draw on it) | **Highest-value development action.** Directly fixes our documented "0 WES/RNA/HLA" gap → converts our *label corpus* into a real *decision-problem corpus*. Not independent proof (overlaps labels) but PRIME never saw its genomics. Related: `phs002735` (Lowery), `phs002748` (scRNA/TCR), `phs002928` (p53). HLA recoverable via OptiType/arcasHLA. |

## TIER C — skip / already-have / contaminated (with reason)
- **HiTIDE `EGAS00001007101`** — clean 4-point bar but **same lab as PRIME** → leakage control only.
- **NCI Gustafson 2025** — huge screen, but **no per-neoantigen negatives released** + Gartner overlap.
- **Autogene cevumeran / BNT122** (Rojas 2023: 230 peptides, 205 neg, LOW leakage) — **inputs NOT deposited** (BioNTech proprietary) → can't recover WES/RNA.
- **NEO-PV-01** (Ott 2020 / Awad 2022: 570 vaccine peptides) — proprietary inputs; peptide tables usable as training-expansion only.
- **Ott 2017 NeoVax `phs001451`** (our known target), **Keskin GBM `phs001519`**, Sahin MUTANOME, Carreno, GAPVAC — pre-2021 → **HIGH leakage** into PRIME/IEDB.
- **mRNA-4157/V940 (KEYNOTE-942)** — efficacy-only; no public per-peptide immunogenicity or inputs.
- **Nous-209 / FixVac / UV1** — shared/off-the-shelf antigens, not per-patient WES-derived → out of scope.
- **NEPdb / dbPepNeo2 / TSNAdb / Neodb / TumorAgDB / ITSNdb / IEDB / CEDAR** — pooled peptide-level, mostly PRIME/IEDB training data.

---

## Recommended two-track plan
**Track 1 — DEVELOPMENT (break the ceiling): acquire `phs001003`.** Completes the Gartner labels we already
trust into a full per-patient decision problem with real inputs. Single highest-value action for actually
moving the gate; dbGaP DAR ~2–4 wk.

**Track 2 — INDEPENDENT PROOF (new, PRIME-untouched cohorts), in order of speed-to-signal:**
1. **Miller IPV (`PRJNA980652`) — start now, free.** Days to ingest; the only fully-open independent cohort.
   Widen its denominator by re-enumerating the full mutanome from the deposited WES.
2. **CheckMate 153 — pull open Supp Table 2 now; file dbGaP DAR for raw WES/RNA.** Biggest independent
   denominator (1,453), lowest leakage. Highest-ceiling proof cohort if raw inputs are obtainable.
3. **WashU TNBC (`phs002787`) — start on open supp labels; file DAR.** New cancer type; open pipeline.
4. **TESLA (`syn21048999`) — download as a secondary benchmark**, pre-registering the MED-leakage caveat.

## Verify before committing effort
- CheckMate 153 raw-sequencing accession + license (supp is open; raw = ?).
- `phs002787` / `phs003922` consent group (GRU vs non-commercial/DS) — matters for an open, possibly-commercial tool.
- PGV001 (`phs003922`) explicit negative counts (supplement).
- HANSolo data-availability + cohort size.
- Whether Miller IPV's 754 peptides are already in IEDB (would contaminate NetMHCpan-EL/BigMHC comparators, not PRIME 2.1 directly).

## Structural caveat (governs Track 2)
Vaccine cohorts (A3, A5) only T-cell-test the ~10–20 **deployed top-ranked** peptides → a **range-restricted**
denominator; they add independent negatives *within the deployable slate* but cannot reproduce Gartner's
full-candidate-universe screen. The genuinely Gartner-dense independent denominators are the **whole-screen /
tetramer** cohorts — **CheckMate 153 (1,453)** and **Miller IPV** (widen via WES). Even so, functional n stays
small (14 / 13 / 18 / 6), so power remains a real constraint — consistent with the v0.2→v0.5 conclusion that
the wall is data, and that the scarce resource is *independent dense-denominator patients*, not model form.
