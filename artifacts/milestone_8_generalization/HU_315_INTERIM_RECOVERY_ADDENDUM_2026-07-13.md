# Hu_315 interim evaluator recovery addendum — 2026-07-13

The first outcome-mediated attempt created its exclusive claim and then stopped
with `LABELS_INVALID` before evaluation or metric persistence. The reason was a
generic evaluator defect: it treated positive and negative assays for the same
genomic mutation as contradictory. The locked endpoint is mutation-level
recognition, and the downstream evaluator already defined recognition as the
presence of any positive assay. Therefore the correct cohort-wide normalization
is fixed **without inspecting Hu_315 values**:

- `POSITIVE` if any recorded assay for patient + genomic mutation is positive;
- otherwise `TESTED_NEGATIVE`.

No portfolio, gate, score, candidate, or patient rule changes. A second outcome
read is permitted only through a new durable recovery claim that pins the first
claim, the first failure record, this addendum, and the corrected evaluator. If
that recovery attempt fails after its claim, the interim readout remains failed
closed. This is development-only evidence and cannot select policy.

