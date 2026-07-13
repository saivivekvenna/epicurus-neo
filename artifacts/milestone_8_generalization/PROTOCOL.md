# Epicurus Neo multi-patient generalization protocol

**Locked before any additional Miller patient reconstruction or outcome join.**

> **Protocol deviation (2026-07-13):** a broad metadata search accidentally
> opened the recognition-label file and exposed calibration rows before all six
> calibration freezes existed. See
> `LABEL_ISOLATION_INCIDENT_2026-07-13.md`. Calibration is consequently
> development evidence with an early-unseal deviation. The six final-held-out
> IDs were not queried and remain the only clean generalization evaluation.

## Objective

Build one patient-agnostic end-to-end Epicurus policy that generalizes across
patients. The same code, gates, score weights, portfolio constraints, and
abstention behavior must run for every patient.

## Roles

- **Known development stress cases:** Hu_287 and Sid. Their outcomes have already
  been inspected. Hu_287's 3/3 is preserved as development evidence; Sid exposes
  abundance/coverage over-penalization and is not a validation patient.
- **Miller calibration:** six patients selected by the deterministic ID-only
  split below. Their labels may be opened only after label-free reconstruction
  and portfolio freeze. They may inform one subsequent universal-policy update.
- **Miller final held-out:** six patients selected by the same ID-only split.
  Their label values remain unused until the final universal policy and all six
  portfolios are frozen. No policy changes after final unseal.

## Label-independent split

Starting population is the 12 Miller patients other than Hu_287 in
`INPUT_CROSSWALK.csv`. Sort by
`sha256("epicurus-generalization-v1|" + patient_id)`; first six are calibration,
last six are final held-out. No assay counts or outcomes participate.

## Universal product boundary

```text
tumor/normal WES + tumor RNA
  -> somatic calling + patient HLA + RNA evidence
  -> lossless class-I mutation-spanning peptide generation
  -> genuine PRIME/MixMHCpred + Epicurus evidence
  -> deterministic validity and absence-safe evidence handling
  -> one universal 20-slot portfolio policy
```

All patient-specific values are evidence, never policy branches. No gene,
mutation, peptide, patient ID, or known-positive whitelist may affect generation,
gating, scoring, or selection.

## Metrics

Primary:

- macro mean unique recognized mutations in the final top 20 per patient.

Co-primary:

- fraction of patients with at least one recognized mutation in the top 20;
- macro recognized-mutation recall@20.

Required safeguards and diagnostics:

- worst-patient hits@20;
- generated -> valid -> eligible -> selected positive reachability;
- paired patient-level delta against genuine PRIME on the identical universe;
- mean and worst-patient duplicate-slot burden;
- number of patients helped, tied, and harmed;
- per-patient abstention and missing-evidence status.

Patients with zero measured positives are retained for hits@20 and `P(>=1)` as
zero but are excluded from recall denominators.

## Policy development rule

Candidate policies may use Hu_287, Sid, existing non-Miller development cohorts,
and the six Miller calibration patients. Selection uses a lexicographic objective:

1. maximize worst-patient hits@20;
2. maximize macro mean hits@20;
3. maximize `P(>=1 hit)`;
4. minimize patients harmed relative to genuine PRIME;
5. prefer the simpler policy.

No candidate is accepted merely for producing 3/3 on Hu_287 or Sid. A universal
policy must not contain patient- or target-specific constants.

## Final freeze and verdict

Before final held-out labels are loaded, commit:

- exact policy/config hashes;
- reconstruction and scoring code hashes;
- raw/derived input hashes per patient;
- stage funnels and ordered selected candidate IDs;
- genuine-PRIME comparator selections;
- a `FROZEN_NO_FINAL_LABELS` manifest.

Final verdict vocabulary: `GENERALIZES`, `TIES_PRIME`, `DOES_NOT_GENERALIZE`, or
`NOT_EVALUABLE`. With only six final patients, all claims remain retrospective
research evidence rather than clinical validation.

## Executable label-isolation gate

The protocol is enforced by `scripts/miller_generalization_eval.py`. Its two
read-only preflight commands re-hash every frozen input, output, historical code
blob, product selection, and ordered candidate list without reading outcomes.

```text
PYTHONPATH=src .venv/bin/python -m scripts.miller_generalization_eval preflight-calibration
PYTHONPATH=src .venv/bin/python -m scripts.miller_generalization_eval calibrate
PYTHONPATH=src .venv/bin/python -m scripts.miller_generalization_eval preflight-final
PYTHONPATH=src .venv/bin/python -m scripts.miller_generalization_eval finalize
```

`calibrate` and `finalize` each create a durable exclusive unseal claim before
their sole label-file read. If a process crashes after that boundary, rerunning
fails closed instead of reopening labels. A selectable policy and genuine PRIME
must both be evaluable on all six relevant patients; partial-cohort averages do
not qualify. The calibration lock pins the exact split, evaluator/CLI bytes,
arm registry, selection objective, per-patient freeze manifests, calibration
result, and recomputed winning policy.
