# Gate feature audit

Isolated, read-only-first audit that asks the one question the falsified dynamic
gate could not answer: **which currently-available ORTHOGONAL feature can remove
the high-presentation TESTED_NEGATIVE decoys that outrank positives** (the top-20
presentation stratum where a label-blind presentation gate removes 0% and
Δhits@20 = 0).

## Files
- `GATE_FEATURE_AUDIT.md` — full narrative: ranked unlock matrix, per-cohort
  conditional/cross-fitted signal, confound audit, LLM feasibility, caveats.
- `FEATURE_UNLOCK_MATRIX.json` — the ranked matrix (family → best stratum AUROC,
  availability now vs requires-Miller-WES/RNA vs unavailable, leakage, next experiment).
- `FEATURE_AUDIT.json` — full machine-readable results across all 7 cohorts.
- `llm_feasibility_cache.json` — cached blind LLM schema + prompt + a 3-row
  feasibility sample (annotation-only; NO labels, NO patient/study IDs; raw
  locked-test sequences redacted).

## Reproduce
`.venv/bin/python scripts/gate_feature_audit.py` (pure core + tests:
`src/benchmark/gate_feature_audit.py`, `tests/test_gate_feature_audit.py`).

## One-line finding
Only **Gartner** and **Zhao** carry both a presentation anchor and orthogonal
features. On Gartner the presentation baseline sits exactly on the wall at the
top-20 stratum (AUROC 0.50) while **expression still separates positives from
decoys there (0.82), cross-fits across patients (OOF 0.855 vs 0.34), and is not
identity-parasitic** — the orthogonal lever the gate lacked — but it is
underpowered (17 stratum positives) and must NOT be used as a rank penalty (that
form is already frozen-falsified). multimer/IMPROVE have no orthogonal column at
all; CEDAR has no anchor; Miller has labels but no inputs yet; Sid has 3
positives and no clean negative denominator. Highest-leverage NEW axis =
label-blind LLM artifact/transcript plausibility.

Does NOT touch any `dynamic_gate` file. Never uses study identity as a deployable feature.
