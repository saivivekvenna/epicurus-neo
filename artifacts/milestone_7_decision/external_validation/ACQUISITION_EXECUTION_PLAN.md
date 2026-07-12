# Acquisition execution plan — north-star external-cohort acquisition

**North star (do not dilute):** across *untouched* patients, using *identical* WES/RNA/HLA inputs, does
Epicurus place **more experimentally recognized neoantigens in the final top-20** than standard
pVAC-style generation + genuine PRIME? **Primary metric = patient-level paired Δ in recognized hits@20.**
Reachability (raw→generation) and conditional ranking are **diagnostic only**. The Sid 3/3 is **recovery of
3 known positives in one patient with no clean negative denominator — never call it "accuracy."**

This plan ranks cohorts by which benchmark level they can support and by expected north-star value. It is
the canonical acquisition doc; it corrects and supersedes the fact rows in
[DATA_ACQUISITION_SHORTLIST.md](DATA_ACQUISITION_SHORTLIST.md) (kept for the lineage analysis) and pairs
with [COHORT_ACQUISITION_TRACKER.csv](COHORT_ACQUISITION_TRACKER.csv),
[MINIMUM_PATIENT_DATA_PACKAGE.md](MINIMUM_PATIENT_DATA_PACKAGE.md),
[AUTHOR_OUTREACH_DRAFTS.md](AUTHOR_OUTREACH_DRAFTS.md), and the four-arm
[SYNTHESIS.md](../four_arm/SYNTHESIS.md). All accession claims below were primary-source re-verified
2026-07-12 (see the "Verification corrections" list at the end).

## Benchmark-level definitions (which a cohort can serve)
- **L3 end-to-end (the north-star gate)** — needs raw WES/RNA/HLA → candidate GENERATION → the same
  candidate universe scored by genuine PRIME *and* Epicurus → paired hits@20, with explicit
  POSITIVE/TESTED-NEGATIVE/UNTESTED labels. Requires open-or-obtainable raw inputs **and** a
  full/reconstructable denominator.
- **L2 conditional ranking (diagnostic)** — rank among the peptides the authors already tested; needs the
  tested peptide-HLA set + labels, interpreted within that cohort's own denominator.
- **L1 reachability (diagnostic)** — raw → generation stage-loss attribution; needs raw inputs + labels.

## Ranked cohorts (by expected north-star value)

| rank | cohort | best level attainable | leakage | denominator quality | explicit negatives | input recoverability | effort | expected north-star value |
|---|---|---|---|---|---|---|---|---|
| **1** | **Miller IPV** (`PRJNA980652`) | **L3 (only fully-OPEN build)** | LOW (2024, La Jolla) | 349 tested = IPV-prefiltered; **full 6,237 detected + open WES → re-enumerable to a real universe** | **YES** (per-peptide ELISpot, S1/S2) | **OPEN** raw WES+RNA; HLA from WES | **MED** (0.23 TB DL + re-enumerate + HLA-type) | **Highest immediate.** Only cohort where the full L3 loop can be built today with no DAR. n=13. |
| **2** | **CheckMate 153** (EGA `EGAD00001011302`) | L2 done (TIE); **L3 if EGA DAR clears** | LOW (Sep-2024, post-PRIME) | **1,453 tetramer (densest independent)** | 1,257 (derived) | raw = **CONTROLLED EGA, BMS discretion** | HIGH (BMS DAC, discretionary) | **Best independent LOCKED test** if raw obtainable; L2 already a TIE. n=14. |
| **3** | **NCI Procurement** (`phs001003`) | DEV completion (L1/L3 for our Gartner labels) | DEV (Gartner-overlap; PRIME didn't train genomics) | Gartner full-mutanome (broad) | Gartner tested-negatives (held) | CONTROLLED dbGaP | MED (dbGaP DAR) | **Highest DEV value** — converts our trusted Gartner *labels* into a decision-problem corpus. Not independent proof. |
| **4** | **WashU TNBC** (`phs002787`) | L2 now (open labels); L3 partial | LOW (Nov-2024) | 198 vaccine-slate = **range-restricted** | 45 pos / 153 neg (open Table S1) | raw CONTROLLED dbGaP | MED | New cancer type; open labels enable an L2 arm immediately; not a full-universe proof alone. n=18. |
| **5** | **Southampton NSCLC** (Zenodo `12820423` + EGA `EGAS00001005499`) | L2 now (open labels); L3 if EGA DAR | LOW | 70 MS-presented tested | **YES** (S9 None/Weak/Strong) | raw CONTROLLED EGA | MED | Small independent LOCKED test; **open negative labels** already minable. n=6. |
| **6** | **PGV001** (`phs003922`) | L2/L3 (controlled) | LOW (2025) | ~103/pt OpenVax | **UNVERIFIED** (paywalled supp) | CONTROLLED dbGaP | MED | Promising but **verify the explicit-negative table exists before spending a DAR**. n=10. |
| **7** | **Hudson/Sid** (osteosarc.com) | descriptive only (n=3) | already used post-hoc | 21 curated; no clean negatives | **NO** (needs pool manifest) | OPEN | LOW (one request) | Upgrades the descriptive case if the **stimulation-pool manifest + deconvolution** is obtained; still n=3, never a headline. |
| **8** | **Hu/Ott NeoVax** (`phs001451`) | DEV reconstruction | MED/HIGH (pre-2021 vaccine peptides leaky) | reconstructable | screened→negatives | CONTROLLED dbGaP | HIGH | Full VAF+RNA reconstruction (existing [ACQUISITION_PACKET.md](ACQUISITION_PACKET.md)); labels leak into IEDB/PRIME → DEV, not clean proof. |
| — | **EVX-01** (PMC11116868) | none (data closed) | moot | 145–231/pt dense | not public | **NOT deposited** | — | **BLOCKED — closed data** ("not publicly available"). Dense but not minable; author-request only, low probability. |
| — | **Zhao 2026** (Front Immunol) | TRAIN only | TRAIN | tested-only (no universe) | 2,004 | not public | — | Recognition-scale **TRAINING** corpus; no raw inputs / no universe → cannot go end-to-end. |
| — | **GBM NeoVax** (`GSE237936`) | incomplete | LOW | n=4, RNA-only | UNVERIFIED | RNA open; **WES/HLA missing** | — | Incomplete; author-request for WES/HLA/panel. n=4 anecdotal. |
| — | **medRxiv aggregate** | unknown | unknown | unknown | unknown | unknown | — | **Unlocatable** (not found after 8 searches); get a DOI before treating as anything. |
| — | **TESLA** (`syn21048999`) | L2 secondary | **MED (recycled)** | 608 multimer | 571 | Synapse DUC | LOW | Secondary benchmark only; pre-register the leakage caveat. |

## Split policy (assign BEFORE any label touches the model)
- **LOCKED_TEST (independent, low-leakage, never used for development):** **Miller IPV** and **CheckMate 153**
  are the two primary locked cohorts (keep ≥2 locked). **Southampton** is a small secondary locked test.
- **DEV (development / ceiling-breaking):** `phs001003` (Gartner completion), plus existing
  multimer/Gartner/IMPROVE; TNBC and Hu/Ott as reconstruction DEV. PGV001 pending negative-table check.
- **TRAIN (priors only):** CEDAR / IEDB / Zhao 2026. Shared-antigen vaccines are recognition priors, not
  personalized end-to-end proof.
- **Leakage controls:** split and hold out **by patient and preferably by study** — never a random peptide
  80/20. Apply exact + near-peptide (k-mer) de-duplication across TRAIN/DEV/LOCKED boundaries. Record each
  cohort's leakage status in the tracker. A cohort's split is fixed the moment it is acquired.

## Go / no-go rule (headline gate)
**No benchmark headline is published unless, for a LOCKED_TEST cohort, the full loop is reproducible:**
raw input → candidate generation → the **same** candidate universe → genuine PRIME **and** Epicurus →
patient-level paired hits@20. If only tested peptides exist (no reconstructable universe), the cohort is
**L2 diagnostic only** and cannot produce a north-star headline. NOT_EVALUABLE must be explicit wherever a
requirement is missing.

## Blockers (current)
1. **Miller IPV S1/S2 supplement** returned HTTP 403 to automated fetch (science.org/Atypon); the per-peptide
   label table's open-download status needs a human browser check. HLA table not confirmed (recoverable from
   WES). This is the only thing between us and a fully-open L3 build.
2. **CheckMate raw seq** is EGA-controlled at **BMS discretion** (harder than a standard dbGaP DAR).
3. **PGV001** explicit-negative table is behind the Nature Cancer paywall — unverified whether a labeled
   per-peptide table exists at all.
4. **EVX-01 / Zhao / GBM / medRxiv** cannot go end-to-end from public data (closed / train-only / incomplete /
   unlocatable).

## Single best action for the user
**Manually open the Miller IPV supplement — `https://www.science.org/doi/suppl/10.1126/scitranslmed.abj9905`
— and confirm data files S1/S2 (the per-peptide ELISpot label table) are downloadable.** Miller IPV is the
only fully-open end-to-end (L3) build prospect (raw WES+RNA already open at `PRJNA980652`, explicit
negatives recorded); that one check unblocks the fastest path to a real north-star run. **In parallel** file
two DARs — `phs001003` (breaks the DEV ceiling on our Gartner labels) and the CheckMate EGA
`EGAD00001011302` (best independent locked test) — since both are DAR-timeline-bound.

---

## Verification corrections (primary-source, 2026-07-12) vs the prior shortlist
- **Miller "IPV" = "identify-prioritize-validate" PLATFORM, not "individualized peptide vaccine."** It is a
  spontaneous-response T-cell identification study; `PRJNA980652` confirmed OPEN (WES+RNA FASTQ, ~0.23 TB).
- **CheckMate 153 raw sequencing is EGA-controlled, NOT dbGaP** — `EGAD00001011302` (studies
  `EGAS00001007508` WXS + `EGAS00001007509` RNA), DAC `EGAC00001003387`, BMS discretion. Supp Table 2 (labels)
  is open on PMC.
- **EVX-01 (PMC11116868 = JITC 2024, 10.1136/jitc-2024-008817) is CLOSED data** — "not publicly available,"
  no accessions, no per-candidate negative table. Downgraded from "urgent dense asset" to author-request-only.
- **GBM `GSE237936` is RNA-seq only, n=4, and the paper is 2024 (Clin Cancer Res, PMID 38639919), not 2026.**
- **medRxiv "Automating neoantigen selection" could not be located** — treat as unverifiable until a DOI is
  supplied.
- **Southampton negatives are OPEN** in Zenodo `10.5281/zenodo.12820423` table S9; TNBC `phs002787` string is
  printed "phs0002787" in-paper (verify on dbGaP); Zhao 2026 = Front Immunol `10.3389/fimmu.2026.1829509`
  (352 pts / 2,317 tested / 313 pos / 2,004 neg), TRAIN-only.
