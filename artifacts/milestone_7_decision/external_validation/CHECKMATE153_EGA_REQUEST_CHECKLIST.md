# CheckMate 153 raw WES/RNA — controlled-access request checklist (EGA)

**No application has been submitted and no one has been contacted.** This is the exact checklist for the
user to file the EGA Data Access Request when ready. CheckMate 153 is a **LOCKED_TEST** cohort; its L2
conditional-ranking arm is already ingested (a TIE) — this raw acquisition upgrades it to an L1/L3
end-to-end test.

## Verified metadata (public EGA API, 2026-07-12)
- **Dataset:** `EGAD00001011302` — "Neoantigen Immunogenicity Landscapes and Evolution of Tumor Ecosystems
  During Immunotherapy with Nivolumab." Whole-Exome (WXS) + RNA-Seq, matched tumor/normal, pre- and
  on-therapy. Verified via `https://metadata.ega-archive.org/datasets/EGAD00001011302` (HTTP 200).
- **Parent studies:** `EGAS00001007508` (WXS) + `EGAS00001007509` (RNA-Seq).
- **DAC:** `EGAC00001003387` ("Senior Director, Grants and Contracts"). Access is at **Bristol Myers Squibb
  discretion** — the dataset policy directs requesters to Dr. William Geese (`william.geese@bms.com`) per
  the EGA dataset landing page. This is a **sponsor-discretion** gate, harder/slower than a standard NIH
  dbGaP DAR.
- **Sample counts (from the dataset page):** WXS 58 pre + 42 on-therapy; RNA-seq 24 pre + 12 on-therapy.
  The 14 functionally-screened patients are a subset.
- **License:** article is subscription (Springer Nature exclusive licence); the per-candidate label
  supplement (Supp Table 2) is already OPEN on PMC (`PMC12066197`) and ingested.

## Request checklist (do in order)
1. **EGA account** at `ega-archive.org` (institutional email).
2. **Identify the DAC/contact:** DAC `EGAC00001003387`; email the sponsor contact
   (`william.geese@bms.com`) to ask the current mechanism (BMS may route via a data-sharing portal or a
   direct DUA rather than the standard EGA form).
3. **Draft the Data Access Agreement / DUA** with:
   - Exact scientific purpose (the north-star benchmark; see `MINIMUM_PATIENT_DATA_PACKAGE.md`).
   - Commitment to publish **aggregate patient-level statistics only** (paired Δ hits@20, CIs); no
     individual-level redistribution; controlled data stays in a controlled environment.
   - Confirmation of the **LOCKED_TEST** use (never for development/tuning) and patient/study-level holdout.
   - Note the tool is **open-source** — confirm the consent group permits non-commercial open-source
     research use (BMS discretion may restrict commercial use).
4. **Request scope:** dataset `EGAD00001011302` (WXS `EGAS00001007508` + RNA `EGAS00001007509`). Ask
   whether processed VCF/BAM and the Polysolver HLA calls can be shared to avoid re-alignment.
5. **On approval:** download to a controlled environment; reconstruct the mutanome; join the already-open
   Supp Table 2 labels **after** ranking; run the frozen four-arm protocol.

## Reproducibility already in hand (no DAR needed)
The L2 conditional-ranking arm is done: `configs/source_manifests/checkmate153.yml` (1,197 candidates, 14
patients, 196 tetramer+ / 1,001 tested-neg after binder filtering; genuine PRIME vs frozen Epicurus = TIE).
The DAR only adds the raw-input L1/L3 (generation) half.

## Status
`EXTERNAL ACTION REQUIRED` — user files the DAR/DUA. Metadata verification is **DONE**; no automated route
can obtain sponsor-discretion controlled data, and none was attempted beyond public metadata.
