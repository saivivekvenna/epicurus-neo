# Milestone 5a — Braun RCC Event-B vertical slice

## Why this milestone

Milestone 4 built the Event-A/B/C corpus substrate and correctly found **zero Event-B
patients in IMPROVE**, because IMPROVE measures pre-existing reactivity (Event A), not
vaccine-induced response (Event B). Milestone 5a proves the substrate can ingest a *real*
public personalized-vaccine study and produce **accepted, provenance-backed Event-B
observations without weakening label semantics** — the first nonzero Event-B in the corpus.

It fits no recognition model. It ingests exactly one study, cleanly, then stops for review.

## Source (open access, frozen)

- Braun DA, Moranzoni G, Chea V, et al. *A neoantigen vaccine generates antitumour immunity
  in renal cell carcinoma.* **Nature 2025;639:474-482.** doi:10.1038/s41586-024-08507-5.
- **PMC11903305**, license **CC BY-NC-ND**; trial **NCT02950766** (phase I, high-risk fully
  resected clear cell RCC; 9 vaccinated patients).
- Vaccine: personalized **synthetic long peptides** (4 pools, ≤5 SLPs each, 300 µg/peptide) +
  **poly-ICLC (Hiltonol)**; cohorts "Vaccine + ipilimumab" (5) and "Vaccine alone" (4).
- Supplements are fetched reproducibly from the Europe PMC open-access endpoint and pinned by
  sha256 (`EXPECTED_SHA256` in `event_b.adapters.braun_rcc`). If the network is unavailable,
  ingestion fails with an actionable manual-placement message rather than fabricating data.

## What the adapter does (and refuses to do)

`BraunRCCAdapter` recomputes the per-peptide immunogenicity call from the raw IFN-γ ELISpot
replicates (`MOESM4`, sheet *In Vitro*) using the paper's **own stated rule**, verbatim from
the Methods / Extended Data Fig. 3 legend:

> **P < 0.05 by two-sided t-test AND mean spot count at least three-fold higher than the DMSO
> (no-stim) control.**

- **POSITIVE** = immunogenic by that rule; **TESTED_NEGATIVE** = assayed but not immunogenic
  (every peptide is explicitly in the deconvolution denominator, so negatives are real, never
  inferred from omission).
- **Event typing is de-novo Event-B**, justified by the paper's explicit statement *"For all
  peptides, no pre-existing immune responses were detected"* and corroborated by ex-vivo
  week-0 (pre-vaccine) pool baselines that sit below the ELISpot positivity floor. If a week-0
  pool ever crossed that floor, those peptides would be demoted to `UNKNOWN_EVENT`; none do.
- The **long vaccine peptide is the tested entity** (`mutant_peptide`), with `mhc_class`
  left `UNKNOWN` (a >20mer SLP is not a class-I minimal epitope; the assay does not resolve
  CD4/CD8 per peptide). The predicted best short epitope + HLA are *not* stored as the assay
  restriction.
- Pools are never decomposed to the peptide level; per-peptide resolution comes only from the
  paper's own in-vitro deconvolution table.

## Reconciliation (expected vs. reproduced)

Applying the rule to the raw replicates reproduces the paper's summary sheet `2e` **exactly**:

| Split | Immunogenic | Non-immunogenic |
|---|---|---|
| Driver | 11 | 6 |
| Passenger | 50 | 62 |
| **Accepted total** | **61** | **68** |

The In Vitro table has 130 assayed peptides; the rule scores 62 immunogenic. The 130th
(`AMACR|p.Y41N`, patient 104) is immunogenic but has a **blank `Mutation_type`** and is
excluded from the paper's driver/passenger summary — so it is **routed to the review queue**
(`UNCLASSIFIED_MUTATION`) rather than silently inflating the accepted count to 62. The 129
accepted assays (61 positive + 68 tested-negative across **9 patients**) match the paper.

## Effect on the global audit

Combining the frozen IMPROVE (Event-A) corpus with Braun (Event-B):

| | Before (Milestone 4) | After (this slice) |
|---|---|---|
| Event-B patients | 0 | **9** |
| Event-B positive-patients | 0 | **9** |
| Event-B studies | 0 | 1 |
| Event-B observations | 0 | **129** (61 POSITIVE / 68 TESTED_NEGATIVE) |
| Verdict | `INSUFFICIENT_DATA_DO_NOT_FIT_RECOGNITION_MODEL` | `EVENT_B_VERTICAL_SLICE_VALIDATED_NOT_YET_SUFFICIENT_FOR_GENERAL_MODEL` |

Event-A (17,082, IMPROVE) and Event-B (129, Braun) remain **separate events**; nothing is
relabelled. The verdict stays conservative and honest: nonzero and validated, still far below
the registered minimum (100 Event-B patients / 2 studies / 30 positive-patients).

## Integrity gates (CI-enforced)

`event_b.braun_pipeline.assert_quality_gates` fails on any accepted Event-B positive without a
post-vaccine timepoint, any tested-negative without an explicit assay denominator, any missing
patient id or label provenance, any clinical-outcome-as-assay, or any contradictory accepted
labels. A non-empty review queue is **not** a failure. Determinism, manifest reproducibility,
and the 61/68 reconciliation are covered by `tests/test_braun_rcc_event_b.py`; a dedicated
`milestone-5a` CI job fetches the open-access supplements and runs them.

## Reproduce

```bash
python scripts/event_b_corpus.py import-braun-rcc
pytest -q tests/test_braun_rcc_event_b.py
```

Outputs: `outputs/event_b_braun/` (normalized tables, event labels, review queue, funnel
links, source manifest, split manifest, reconciliation, Braun audit, model-ready parquet),
`outputs/event_b_corpus_combined/` (IMPROVE+Braun), and `artifacts/milestone_5a/`
(reconciliation + combined audit, tracked in git). The IMPROVE and combined corpora live under
gitignored `outputs/`; regenerate IMPROVE with `import-improve-event-a` first if absent.

## Not in scope / next

No recognition model is fitted; IMPROVE labels are untouched; no negative is inferred from
omission. The next step (**M5B**) scales adapters to the remaining open Event-B studies
(NeoVax melanoma/PDAC, mKRAS-VAX, Nous-209, Fukuoka DC) toward the registered minimum — each
with its own definitions and traps, added one at a time and reconciled the same way.
