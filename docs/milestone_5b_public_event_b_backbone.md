# Milestone 5b: public Event-B backbone

This milestone builds a provenance-first, resumable ingestion path for public evidence of
vaccine-induced or meaningfully expanded target-specific T-cell response. It does not train a
recognition model.

## Factory commands

```bash
.venv/bin/python scripts/event_b_corpus.py ingest-study <study_id>
.venv/bin/python scripts/event_b_corpus.py ingest-backbone --resume
.venv/bin/python scripts/event_b_corpus.py ingest-study <study_id> --rebuild
.venv/bin/python scripts/event_b_corpus.py rebuild-combined
.venv/bin/python scripts/event_b_corpus.py audit-sufficiency
```

The versioned registry is `configs/event_b_studies.yml`. A job checkpoint records the registry,
schema, adapter, and source-manifest fingerprint. A completed output is reused only when all of
those values still match and its success marker exists. Changed inputs are marked `STALE` and are
not silently reused. An interrupted adapter run is marked `INTERRUPTED`; rerunning repeats the
deterministic transform from frozen sources.

Source validation checks the detected file type rather than trusting the suffix, rejects HTML
access/error pages disguised as supplements, streams checksums in 1 MiB chunks, records file size,
and validates JSON/JSONL and OOXML containers. Raw sources remain under the ignored `data/raw/`
tree.

Braun RCC and Hu melanoma reproduce through this factory without changing their accepted labels.
The mandatory backbone studies are accepted only if public sources support the required patient
and assay granularity; otherwise the registry and study note retain a reproducible blocker.

## Final disposition

mKRAS-VAX and PDAC NeoVax have candidate-resolved accepted adapters. Nous-209 is accepted only at
patient level because patient-by-pool identities are absent. Fukuoka is blocked by source access and
patient-overlap uncertainty. Together with Braun and Hu, the combined corpus contains 82 total
Event-B patients, but only 45 are candidate-resolved (the peptide-ranking sample) across 4
candidate-resolved studies; the other 37 are Nous-209 patient-level-only participants. All 974
primary candidate labels are nested inside those 45 patients. The corpus supports a viable full-study
holdout, but because the registered gate counts candidate-resolved patients only, 45 < 100 yields
`INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA`.
