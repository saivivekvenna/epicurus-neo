# Milestone 5b.1 — Hu 2021 melanoma NeoVax Event-B vertical slice

## Why this milestone

Milestone 5a ingested the first real Event-B study (Braun RCC) and moved the corpus from
zero to nine vaccine-induced-response patients, with the single-study verdict
`EVENT_B_VERTICAL_SLICE_VALIDATED_NOT_YET_SUFFICIENT_FOR_GENERAL_MODEL`. Milestone 5b.1
adds the **second independent Event-B study** — the Hu 2021 melanoma NeoVax follow-up —
one study at a time, reconciled the same way, and trips the study-diversity gate of the
registered win condition (two Event-B studies) while remaining honestly below the
patient-count minimum.

It fits no recognition model. It ingests exactly one new study, cleanly, then stops.

## Source (manual, checksum-pinned)

- Hu Z, Leet DE, Allesøe RL, et al. *Personal neoantigen vaccines induce persistent memory
  T cell responses and epitope spreading in patients with melanoma.* **Nat Med
  2021;27:515-525.** doi:10.1038/s41591-020-01206-4.
- **PMC8273876**; trial **NCT01970358** (the same NeoVax melanoma trial reported in Ott et
  al. 2017). Long-term follow-up of the original patients plus two new ones.
- Hu 2021 and Ott 2017 (**PMC5577644**) are PMC *author manuscripts* behind a JavaScript
  proof-of-work download gate, so they are not fetched programmatically. The one workbook the
  slice needs (`NIHMS1707651-supplement-Suppl_DataSet.xlsx`, ~2.2 GB uncompressed) is placed
  manually and pinned by sha256 (`EXPECTED_SHA256` in `event_b.adapters.hu_neovax`). Ingestion
  refuses to run — with an actionable download/placement message — rather than fabricating data
  if the file is absent or its checksum differs. See `data/raw/MANUAL_SOURCES.md`.

## What the adapter does (and refuses to do)

`HuNeoVaxAdapter` reads only the small recognition sheets from the giant workbook via a
**streaming reader** that never materialises the multi-hundred-MB shared-string table or the
bulk TCR-repertoire sheets. It ingests the source's own consolidated per-peptide ELISpot calls
(Hu scored a response positive at ≥2.5× the DMSO control); the table carries no raw replicates,
so calls are ingested **as reported** and reconciled against Ott 2017's independently published
totals rather than recomputed from replicates (contrast Braun, where the rule was re-derived).

- **CD8 (Dataset 4a)** is class I: the tested entity is the predicted minimal epitope
  (`mutant_peptide`) with its restricting HLA (`mhc_class = CLASS_I`); the parent long
  immunizing peptide is kept in `candidate_source`, preserving the long-peptide/minimal-epitope
  relationship.
- **CD4 (Dataset 4b)** is class II: the tested entity is the overlapping assay peptide
  (`mhc_class = CLASS_II`, no per-peptide class-II restriction resolved by the source); the
  week-16 label is positive if reactive ex-vivo or after pre-stimulation.
- **Vaccine inclusion vs assay target** is explicit: Datasets 4a/4b peptides are vaccine-included
  Event-B; the assay target (minimal epitope / assay peptide) is distinct from the vaccine
  candidate (the long immunizing peptide).
- **Epitope spreading (Datasets 11a–c)** — responses to neoantigens that were **not** in the
  vaccine — is typed `EPITOPE_SPREADING` with `vaccine_inclusion = NOT_INCLUDED`, and is **never**
  counted as vaccine-candidate recognition, even when a spreading response is positive.
- **Patient de-duplication is structural**: patients 1–6 are the Ott 2017 cohort and 11–12 are
  new; a single consolidated source is used for all eight, so there is no cross-source double
  count.
- **Week 16** is the primary recognition timepoint; a scored 3–4.5-year readout is recorded as
  `LONGITUDINAL_PERSISTENCE` recognition evidence rather than a second label.

## Reconciliation (observed vs. Ott 2017)

Vaccine-peptide recognition for the original cohort (patients 1–6, week 16), aggregated to the
neoantigen level, against Ott's published totals:

| Channel | Positive neoantigens | Ott reported | Total neoantigens |
|---|---|---|---|
| CD8 | **15** | 15 | 97 (Ott 97) |
| CD4 | 56 | 58 | 96 (Ott 97) |

The **CD8 positive count and denominator match Ott exactly (15/97)**. CD4 differs by two
neoantigens; the gap is not closed by any other assay condition in the table (minigene/tumor add
nothing), consistent with Ott counting at the immunizing-peptide rather than the neoantigen
granularity. It is reported, not tuned away. Ott's `Supp_5/6` PDFs are recorded as the deeper
per-patient audit reference.

## Evidence strength is not flattened

Per the 5a feedback, distinct recognition channels carry **distinct** `recognition_evidence`
reliability vectors instead of one uniform strength:

| Channel (evidence family) | candidate_specificity | assay_directness | vaccine_relevance | temporal_clarity |
|---|---|---|---|---|
| CD8 minimal epitope (VACCINE_EVENT_B) | 1.0 | 0.9 | 1.0 | 0.7 |
| CD4 assay peptide, ex-vivo (VACCINE_EVENT_B) | 0.6 | 0.9 | 1.0 | 0.7 |
| CD4 assay peptide, pre-stim (VACCINE_EVENT_B) | 0.6 | 0.6 | 1.0 | 0.7 |
| Epitope spreading (FUNCTIONAL_T_CELL_ASSAY) | 1.0 | 0.9 | **0.0** | 0.6 |
| Persistence 3–4.5 yr (LONGITUDINAL_PERSISTENCE) | 0.8 | 0.7 | 1.0 | 0.8 |

`temporal_clarity` is 0.7 for vaccine recognition (not 1.0 as a baseline-verified de-novo claim
would earn): this table has **no pre-vaccine week-0 column**, so de-novo status is author-asserted
(Ott/Hu both report no pre-vaccination reactivity), not re-verified from a baseline as it was for
Braun. Epitope spreading carries zero vaccine relevance by construction.

## Effect on the global audit

Combining frozen IMPROVE (Event-A), Braun (Event-B), and Hu (Event-B):

| | After Milestone 5a | After this slice |
|---|---|---|
| Event-B patients | 9 | **17** |
| Event-B studies | 1 | **2** |
| Event-B observations | 129 | **670** (Braun 129 + Hu 541) |
| Epitope-spreading observations | 0 | **82** (kept separate) |
| Verdict | `EVENT_B_VERTICAL_SLICE_VALIDATED_NOT_YET_SUFFICIENT_FOR_GENERAL_MODEL` | `EVENT_B_MULTI_STUDY_CORPUS_VALIDATED_INSUFFICIENT_PATIENTS_FOR_GENERAL_MODEL` |

Event-A (17,082, IMPROVE), Event-B (670), and epitope spreading (82) remain **separate**;
nothing is relabelled. The study-diversity gate (≥2 Event-B studies) is now met; the sample-size
gate (≥100 patients / ≥30 positive-patients) is not, and the verdict says so.

## Integrity gates (CI-enforced)

`event_b.hu_pipeline.assert_quality_gates` shares the Braun Event-B checks and adds the
melanoma-specific invariants that **no epitope-spreading response is tied to a vaccine-included
candidate** and **every Event-B recognition is tied to a vaccine-included candidate**. Three-state
labels, class-I/II handling, the differentiated reliability vectors, epitope-spreading separation,
the multi-study verdict, streaming-reader correctness, determinism, and the CD8 15/97
reconciliation are covered by `tests/test_hu_neovax_event_b.py`. The synthetic gates run in the
standard suite; the real-source reconciliation, determinism, and combined-audit tests run when the
manually-placed supplement is present.

## Reproduce

```bash
# Place the source first (see data/raw/MANUAL_SOURCES.md):
#   data/raw/hu_melanoma_2021/manual/NIHMS1707651-supplement-Suppl_DataSet.xlsx
python scripts/event_b_corpus.py import-hu-neovax
pytest -q tests/test_hu_neovax_event_b.py
```

Outputs: `outputs/event_b_hu/` (normalized tables, event labels, recognition evidence, review
queue, source manifest, split manifest, reconciliation, Hu audit, model-ready parquet),
`outputs/event_b_corpus_combined/` (IMPROVE + Braun + Hu), and `artifacts/milestone_5b1/`
(reconciliation + combined audit, tracked in git). Prior corpora (`import-improve-event-a`,
`import-braun-rcc`) are combined automatically when materialised.

## Not in scope / next

No recognition model is fitted; IMPROVE labels are untouched; no negative is inferred from
omission. Epitope-spreading Datasets 11d/11e (MSEC-nominated, wild-type-ligand only; 11e has
`#REF!`-corrupted cells) and Dataset 12 (class-II spreading predictions, no measured reactivity)
are deferred, not mis-parsed. The next studies (**mKRAS-VAX, PDAC NeoVax, Nous-209, Fukuoka DC**)
are added one at a time and reconciled the same way, toward the registered patient-count minimum.
