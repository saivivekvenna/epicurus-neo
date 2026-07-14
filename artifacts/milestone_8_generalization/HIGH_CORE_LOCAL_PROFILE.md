# High-core local reconstruction profile

Registered before the remaining Miller reconstructions. This is an operational
acceleration only and cannot change variants, peptide semantics, scores, gates,
or portfolios.

- The reusable patient driver defaults to 12 workers on hosts with at least 12
  logical CPUs, capped at the available CPU count.
- `EPICURUS_THREADS` or `--threads` may lower the count for memory pressure, but
  the chosen value is passed uniformly to FASTQ conversion, WES alignment, RNA
  alignment, and Salmon quantification.
- `fasterq-dump` now receives the same worker count via `-e`; previously it ran
  at its implicit default even when the driver requested more cores.
- On the current 15-core, 24-GiB Mac, 12 workers leave three logical CPUs for
  the OS and orchestration. The shell stages retain their existing integrity,
  resume, and provenance checks.
- Ensembl annotation remains network-bound; adding CPU workers to that stage
  would not improve throughput and could increase API pressure.
