# Public Event-B backbone run state

- Branch: `milestone-5b-public-event-b-backbone`
- Starting merged commit: `8731e1d`
- Completed phases: baseline gate; factory; mKRAS-VAX; PDAC NeoVax; Nous-209
- Completed studies: Braun RCC, Hu/NeoVax, mKRAS-VAX, PDAC NeoVax, and Nous-209
- Backbone studies: Fukuoka source/semantics review in progress
- Source files pinned: Braun, Hu, mKRAS, PDAC NeoVax, and Nous-209 source sets
- Files currently being edited: registry, source notes, backbone audit implementation
- Latest focused verification: Nous-209 export reconciled 37 patient-level observations, 0 candidate labels
- Exact next task: resolve the canonical Fukuoka dataset and label comparability

## Reproduction

```bash
.venv/bin/python scripts/event_b_corpus.py ingest-study braun_rcc_2025 --resume
.venv/bin/python scripts/event_b_corpus.py ingest-backbone --resume
.venv/bin/python -m pytest -q
.venv/bin/python -m ruff check .
```
