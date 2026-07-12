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
- `IMPROVE_RAW_COLUMN_AUDIT.md` / `.json` — **correction** to the first pass:
  full classification of every one of the 88 columns in the raw IMPROVE table
  (`data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt`)
  into deployable / suspicious-derived / context-only / presentation / forbidden
  / split, with per-column coverage + within-patient variation + a leakage screen.
- `IMPROVE_DEPLOYABLE_WHITELIST.json` — the deployable candidate-varying feature
  whitelist (36 features / 11 families) for the rich-feature gate session, plus
  the explicit excluded buckets.

## Correction (supersedes the first pass on IMPROVE)
The first run concluded "IMPROVE has no orthogonal features." That was a
**loader-scope artifact**: the loaded IMPROVE frame is the *reduced* export
(prime/el/expr only). The **raw 88-column** table is orthogonally rich — DNA VAF,
RNA support (reads/AF/confirmation), agretopicity (DAI), foreignness/self-similarity,
clonality (CelPrev), stability, HLA expression, driver/mutation-class flags, and a
full physicochemical block: **36 deployable candidate-varying features across 11
families**. Suspicious/derived (`PrioScore`, `IB_CB`/`IB_CB_cat`, `NetMHCExp`) are
withheld pending provenance (not leaky, but circular/composite); `pMHC`/`norm_pMHC`
are `allele_peptide` identity strings not scores; `validation` is a constant QC
flag; the immune-deconvolution block is per-sample constant (context-only, cannot
re-order a patient's candidates). The whitelist is chosen from pre-declared
biology + within-patient variation + coverage; the held-out `response` is used
only as a leakage screen on the suspicious/forbidden buckets, never to select.

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
