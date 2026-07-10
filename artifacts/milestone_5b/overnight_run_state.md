# Public Event-B backbone run state

- Branch: `milestone-5b-public-event-b-backbone`
- Starting merged commit: `8731e1d`
- Completed phases: baseline integrity gate; registry and resumable ingestion factory; mKRAS-VAX
- Completed studies: Braun RCC, Hu/NeoVax, and mKRAS-VAX
- Backbone studies: PDAC NeoVax, Nous-209, and Fukuoka source/semantics review in progress
- Source files pinned: existing Braun and Hu sources; mKRAS PDF and XLSX supplements
- Files currently being edited: registry, source notes, backbone audit implementation
- Latest focused verification: mKRAS factory export reconciled 12 patients and 72 primary labels
- Exact next task: freeze PDAC NeoVax sources and resolve candidate-level denominators

## Reproduction

```bash
.venv/bin/python scripts/event_b_corpus.py ingest-study braun_rcc_2025 --resume
.venv/bin/python scripts/event_b_corpus.py ingest-backbone --resume
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
```
