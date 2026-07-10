# Public Event-B backbone run state

- Branch: `milestone-5b-public-event-b-backbone`
- Starting merged commit: `8731e1d`
- Completed phases: baseline integrity gate; ingestion factory; mKRAS-VAX; PDAC NeoVax
- Completed studies: Braun RCC, Hu/NeoVax, mKRAS-VAX, and PDAC NeoVax
- Backbone studies: Nous-209 and Fukuoka source/semantics review in progress
- Source files pinned: Braun, Hu, mKRAS, and PDAC NeoVax source sets
- Files currently being edited: registry, source notes, backbone audit implementation
- Latest focused verification: PDAC factory export reconciled 16 patients and 232 target labels
- Exact next task: audit Nous-209 pool resolution and public patient mapping

## Reproduction

```bash
.venv/bin/python scripts/event_b_corpus.py ingest-study braun_rcc_2025 --resume
.venv/bin/python scripts/event_b_corpus.py ingest-backbone --resume
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
```
