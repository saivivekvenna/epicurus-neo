# Milestone 3 — Candidate reachability and label headroom

This is the next focused milestone after the frozen Milestone 1 instrumentation commit. It implements
the useful part of the proposed headroom work without reopening recognition-model iteration.

## Question

For every experimentally validated positive, where did it first become unreachable?

```text
validated positive
→ mutation called
→ correct transcript represented
→ peptide generated
→ survives expression/clonality gates
→ restricting HLA included
→ presentation candidate
→ ranking stage
→ top k
```

The unit is a validated peptide/target, not a patient. Every stage reports recall with a 95% Wilson
interval. Patient counts remain visible but are not misrepresented as the sampling unit for this
candidate-generation measurement.

## Evidence rules

- `reached`: the target identity is present in a complete stage artifact.
- `lost`: the stage artifact is complete and the target identity is absent.
- `not_assessed`: the required artifact is missing or incomplete.
- Missing evidence is never converted to `lost`.
- Once a target is lost, it cannot reappear downstream.
- Identity keys are declared independently per stage. Mutation, transcript, peptide, and peptide-HLA
  membership are never treated as interchangeable.

The implementation is in `benchmark.funnel`. It accepts either an explicit status ledger or stage
tables with declared identity keys.

## Label audit

Every evaluated row must also carry:

```yaml
event_type: pre_existing_reactivity | vaccine_induced_response | post_treatment_reactivity | unknown
assay: elispot | ics | tetramer | manafest | cytokine_release | other | unknown
label: POSITIVE | TESTED_NEGATIVE | UNTESTED
timepoint: pre_vaccine | post_prime | post_boost | post_vaccine | post_treatment | unknown
provenance: non-empty source location
```

`benchmark.label_metadata.validate_label_metadata` rejects tested labels without an assay and rejects
pre-vaccine observations labelled as vaccine-induced responses.

## Current dependency boundary

This repository does not currently contain pVACtools or complete raw WES/RNA-to-stage artifacts.
Therefore no empirical funnel result is claimed yet. The first run requires:

1. a frozen validated-positive registry with mutation, transcript, peptide and restricting-HLA IDs;
2. complete stage exports from variant calling through ranking;
3. a per-stage completeness declaration so absence can honestly mean loss;
4. NetMHCpan licensing or a pre-registered MHCflurry fallback for presentation generation.

Until those inputs exist, the framework reports `not_assessed`; it does not manufacture headroom.

## Out of scope

- training or selecting a recognition model;
- adding contact-residue, anchor-creation, PLM, or hand-reasoned features;
- opening TESLA, the 2025 multimer set, or Sijbrandij early;
- treating Event A and Event B as interchangeable;
- reporting funnel recall without a confidence interval.
