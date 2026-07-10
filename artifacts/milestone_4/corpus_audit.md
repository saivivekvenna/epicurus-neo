# Event-B corpus audit — available repository data

Generated from the official IMPROVE archive via `ImproveEventAAdapter`.

This corpus is a recognition-evidence substrate. It is not proof of clinical benefit.

## Independent sample sizes

- Source studies: 1
- Patients: 70
- Unique peptide sequences: 15,293
- Canonical candidates: 17,520
- Accepted candidate-assay observations: 17,082
- Vaccines: 0

The peptide count is not an independent sample size. The only available compatible source has 70
patients and contains Event A, not Event B.

## Events and labels

| Biological event | Positive | Tested negative | Untested |
|---|---:|---:|---:|
| Event A — pre-existing reactivity | 458 | 16,624 | 0 |
| Event B — vaccine-induced response | 0 | 0 | 0 |
| Event C — clinical outcome | 0 | 0 | 0 |
| Presentation only | 0 | 0 | 0 |

The source contains 17,520 assay rows before canonical review (467 positive, 17,053 tested negative).
Validation queued 438 candidates whose stored mutant and wild-type peptide sequences are identical;
their linked assays are preserved in the normalized source layer but excluded from the accepted/model-
ready assay export. Nine of those quarantined rows were source positives.

## Model-readiness decision

```text
INSUFFICIENT_DATA_DO_NOT_FIT_RECOGNITION_MODEL
```

- Evaluable Event-B patients: 0
- Event-B studies: 0
- Patients with an Event-B positive: 0
- Registered diagnostic minimum: 100 Event-B patients, 2 studies, 30 positive patients

No recognition baseline or model was fit.
