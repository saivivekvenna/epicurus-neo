# Hu_315 immutable interim readout addendum — 2026-07-13

## Authorization and scope

The user explicitly requested an immediate Hu_315 performance result rather
than waiting for all six Miller calibration reconstructions. Hu_315 belongs to
the **calibration/development** split, not the final-held-out split. Calibration
was already classified as development evidence after the label-isolation
incident recorded in `LABEL_ISOLATION_INCIDENT_2026-07-13.md`.

## Rules registered before this readout

1. The existing Hu_315 `FROZEN_NO_LABELS` universe and all registered ordered
   portfolios must independently verify before any outcome access.
2. All frozen evaluation inputs are byte-snapshotted before the label read.
3. One durable exclusive claim is created before exactly one evaluator-mediated
   label-file read. The outcome table is never searched or manually inspected.
4. Every registered arm is evaluated. No arm is selected from this one patient.
5. The result is immutable descriptive development evidence and **must not**
   change gates, weights, reconstruction, portfolio policy, or patient handling.
6. Hu_315 is thereafter explicitly observed development data. It cannot be
   called held out, and it cannot contribute to the final generalization claim.
7. The six patients in `SPLIT.json.final_held_out` remain sealed and unchanged.

The result answers only: *given the already-frozen Hu_315 reconstruction and
portfolios, how many experimentally recognized mutations did each registered
top-20 contain, and at which upstream stage were recognized mutations lost?*

