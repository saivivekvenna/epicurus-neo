# Milestone 4 — Event-B corpus and recognition-evidence substrate

## Why Event B is the label

Epicurus designs a vaccine. The relevant target is therefore a response newly induced or meaningfully
expanded after vaccination, with candidate-resolved evidence that the source attributes to the
vaccinated target. Pre-existing reactivity (Event A), vaccine-induced response (Event B), clinical
outcome (Event C), and presentation-only evidence are stored as separate biological events. They are
never collapsed into a generic `immunogenic` label.

The corpus is a structured evidence substrate. **It is not itself proof of clinical benefit.**

## Three-state response semantics

Every candidate-assay observation is `POSITIVE`, `TESTED_NEGATIVE`, or `UNTESTED`.
`TESTED_NEGATIVE` requires explicit inclusion in the relevant assay denominator. A generated,
vaccine-included, or unreported candidate remains `UNTESTED` when no candidate-resolved assay result
exists. Ambiguity enters the review queue; it never becomes a negative by convenience.

## Canonical entities and provenance

`event_b.models` defines versioned Study, Patient, Vaccine, Candidate, Assay, Clinical Outcome,
Recognition Evidence, Funnel Link, and field-level Provenance tables. Candidate IDs hash study,
patient, sample/timepoint, variant, transcript, mutant peptide, and HLA context; peptide sequence alone
is not an identity.

Every reported or derived field can point to a provenance record containing document, page, table,
figure, supplement, row, column, source fragment, extraction method, confidence, value origin, and
review status. Source manifests checksum the raw documents.

## LLM extraction and review

LLMs may emit schema-constrained extraction records, but the schema forces LLM-derived records to
remain `NEEDS_REVIEW`. Raw outputs are content-addressed and cached. When no provider is configured,
the pipeline writes `extraction_tasks.jsonl` and the expected JSON Schema and reports `PENDING`; it
does not fabricate records. Deterministic validators and human review are the only routes to an
accepted label.

## Validation and contradictions

Validation checks post-vaccine evidence for Event B, negative denominators, chronology, cross-study
links, HLA compatibility, peptide length/class consistency, vaccine-inclusion provenance,
mutant/wild-type roles, clinical/assay separation, contradictory accepted labels, and field origins.
Contradictory or unsafe records are removed from the accepted assay table and written to
`review_queue.jsonl`; source records remain preserved in the normalized corpus.

## Funnel linkage and evidence channels

Event-B candidates link to the existing reachability ledger from mutation calling through top-k.
Vaccine inclusion and functional assay status are appended without inferring missing upstream stages.
Recognition evidence remains separated into presentation, immunopeptidomics, tolerance, primability,
TCR, functional, clonality, spatial, and longitudinal channels. Reliability dimensions are stored
independently; no universal evidence score is imposed.

## Leakage-safe partitions and timepoint safety

Deterministic manifests support patient, study, HLA, peptide-cluster, cancer-type, and temporal
holdouts. Patient rows are atomic, peptide clusters use normalized edit similarity, and study holdout
keeps trials disjoint. Pre-selection evidence is explicitly separated from outcome-only evidence;
post-vaccine ELISPOT, TCR expansion, imaging, and clinical outcomes cannot enter the pre-vaccine
feature matrix.

## Audits and exports

Exports preserve normalized tables as Parquet, repeated assays and timepoints as separate rows,
review issues as JSONL, source manifests as JSON, deterministic split manifests, and both JSON and
Markdown audits. Audits report peptide, patient, and study counts separately, plus event/label/assay
coverage, linkage, missingness, contradictions, and a registered model-readiness decision.

## Data currently missing

No candidate-resolved vaccine-trial corpus is currently stored in this repository. IMPROVE can be
imported reproducibly through `ImproveEventAAdapter`, but it remains Event A and contributes **zero
Event-B patients**. Osteosarc has a declaration-only case-study adapter and refuses to fabricate data
when source files are absent. It remains an n=1 integration/longitudinal case, never a population
training cohort.

The registered minimum for even transparent recognition diagnostics is 100 evaluable Event-B
patients across at least two studies, including 30 patients with an Event-B positive. Until curated
trial sources satisfy those thresholds, the audit decision is:

```text
INSUFFICIENT_DATA_DO_NOT_FIT_RECOGNITION_MODEL
```

This milestone feeds the later recognition-model milestone by supplying accepted labels, independent
patient/study groups, evidence availability, provenance, and pre-selection-safe features. It does not
begin that modeling work.

## Commands

```bash
python scripts/event_b_corpus.py import-improve-event-a /tmp/IMPROVE_paper outputs/event_b_corpus

python scripts/event_b_corpus.py emit-extraction-tasks \
  --study-id trial-id \
  --source extracted-paper-text.txt \
  --output outputs/event_b_extraction_tasks/trial-id
```
