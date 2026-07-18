# Demo patient

A tiny, synthetic pVACseq-style candidate table for trying the **prioritize** stage without
running the full read-level pipeline.

```bash
epicurus run-patient \
  --input examples/demo_patient/pvacseq_all_epitopes.tsv \
  --patient-id DEMO-001 \
  --output-dir out/demo_patient
```

This writes `ranked_candidates.csv`, `report.json`, and `report.md`. The demo deliberately includes
an unexpressed gene (CTNNB1, TPM 0) so you can see the validity gate keep it out of the top-20.

`patient.yaml` is a template for the **full** pipeline (`epicurus run-pipeline`); point the FASTQ and
reference paths at real files on a machine with the external tools and reference bundle installed
(see the top-level README for requirements).
