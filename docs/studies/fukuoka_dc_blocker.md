# Fukuoka dendritic-cell corpus blocker

## Canonical study decision

The registry anchors the requested Fukuoka dataset to Morisaki et al., *Anticancer Research* 2021,
`doi:10.21873/anticanres.15213`: a retrospective analysis of 17 advanced-cancer patients treated
with intranodal personalized neoantigen peptide-pulsed dendritic-cell vaccine monotherapy.

The original task's `10.21873/anticanres.15215` identifier was verified to belong to an unrelated
hepatitis-C frailty paper and has been corrected. The 2021 clinic paper is the only publication in
the Fukuoka series that could plausibly serve as the requested multi-patient backbone dataset.

## Why it is blocked

An open abstract verifies 17 patients and reports that some nonresponding patients had no ELISpot
reaction, but it does not expose the patient-by-peptide denominator, sequences, baseline status,
post-vaccine assay inclusion or missing-follow-up mechanism. The publisher endpoint returns an HTML
access challenge rather than the article PDF. No restricted access was bypassed.

Later open clinic publications cannot substitute safely:

- The 2023 hybrid class-I/class-II paper contains four stage-IV patients and cites the earlier clinic
  cohort; cross-publication patient reuse is unresolved.
- The 2024 breast study contains five postoperative patients and explicitly states that no dataset
  can be made publicly available.
- Ovarian, salivary-duct, gastric and other reports are single cases that may overlap the 2021 series.
- The clinic protocol measures ELISpot after three doses and may use ELISpot-positive peptides for
  later doses, creating response-adaptive observation and selection that must be reconstructed.

Without the full source, `No response` could be patient-level, positive-only follow-up could bias the
peptide denominator, and repeated publications could duplicate patients. The study is therefore
`BLOCKED_SOURCE_UNAVAILABLE`, contributes zero patients and zero labels, and remains outside both the
primary corpus and dominance calculations.

## Unblocking contract

Follow `configs/source_manifests/fukuoka_dc.yml`. Required inputs are the lawful full 2021 article
and supplements plus an author-confirmed patient map across publications. Acceptance requires an
explicit per-peptide tested denominator at the locked primary timepoint and a documented account of
response-adaptive vaccination and follow-up.
