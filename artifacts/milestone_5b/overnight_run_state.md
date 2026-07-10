# Public Event-B backbone run state

- Branch: `milestone-5b-public-event-b-backbone`
- Starting merged commit: `8731e1d`
- Completed phases: baseline integrity gate; registry and resumable ingestion factory
- Completed studies: existing Braun RCC and Hu/NeoVax remain accepted and unchanged
- Backbone studies: source and label-semantics review in progress
- Source files pinned: existing Braun and Hu sources only
- Files currently being edited: registry, source notes, backbone audit implementation
- Latest focused verification: 27 tests passed; Ruff passed; Braun factory run reused cleanly
- Exact next task: freeze primary-source semantics and resolve each mandatory backbone study

## Reproduction

```bash
.venv/bin/python scripts/event_b_corpus.py ingest-study braun_rcc_2025 --resume
.venv/bin/python scripts/event_b_corpus.py ingest-backbone --resume
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
```
