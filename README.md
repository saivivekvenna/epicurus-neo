# Epicurus Neo

Epicurus prioritizes at most 20 neoantigen candidates for a personalized cancer vaccine. The product
target is vaccine-inducible response (Event B), not pre-existing T-cell reactivity. Candidate
generation remains the responsibility of pVACtools; Epicurus owns gating, calibrated ranking,
portfolio selection, and patient-level abstention.

Milestone 1 instrumentation is frozen in git. The next focused work is the candidate-reachability
funnel: identify whether validated positives are lost at mutation calling, transcript selection,
peptide generation, gating, HLA inclusion, presentation, ranking, or top-k selection. It trains no
recognition model and does not touch the external TESLA, 2025 multimer, or Sijbrandij sets early.

## Registered evaluation contract

- Primary: patient-level mean `hits@20`, including zero-positive patients as zero.
- Co-primary: capture fraction, excluding zero-positive patients as unevaluable.
- Clinical gate and eventual headline: `P(≥1 hit in top 20)`.
- Diagnostics: `precision@20`, MRR, unreachable patients, candidate-list random expectation, and
  reranking headroom.
- Every result carries a 20,000-resample paired bootstrap interval against a named baseline.
- Ties always break on `md5(mutant_peptide|hla_allele)`; source row order never decides membership.
- Labels have three states: `POSITIVE`, `TESTED_NEGATIVE`, and `UNTESTED`.
- A hand-reasoned or LLM-derived rule cannot ship unless it beats PRIME on the same rows with a paired
  confidence interval excluding zero.

`benchmark.scorecard.scorecard()` is the reporting path. It emits all five metrics, paired deltas,
retained sample sizes, the candidate-universe fingerprint and random baseline, unreachable-patient
count, current-n MDE, and a computed `ACCEPT`, `CONSISTENT_WITH_NO_EFFECT`, or `REJECT` verdict.

## Reproduce Milestone 1

The official IMPROVE repository contains both required archives. A normal clone is sufficient:

```bash
git clone https://github.com/SRHgroup/IMPROVE_paper.git /tmp/IMPROVE_paper
python scripts/milestone_1.py /tmp/IMPROVE_paper verify
pytest -q
```

Generate the ten blind masking-ablation question sets with:

```bash
python scripts/milestone_1.py /tmp/IMPROVE_paper generate-ablation \
  artifacts/milestone_1/masking_ablation
```

The frozen-score results are in
[`docs/milestone_1_reaudit.md`](docs/milestone_1_reaudit.md). Historical research iterations remain
in [`docs/benchmark_iterations.md`](docs/benchmark_iterations.md); they are a record, not the current
evaluation contract.

The next-stage evidence contract is in
[`docs/milestone_3_funnel_spec.md`](docs/milestone_3_funnel_spec.md). Once complete stage exports are
available, `epicurus funnel-report ledger.csv` reports per-stage candidate recall with confidence
intervals and explicit bounds for missing evidence.

## Benchmark roles

- IMPROVE official five-fold patient CV: primary Event-A component/ranking regression.
- BigMHC `im_test`: HLA-grouped component regression only.
- Gartner/NCI Nmers: frozen TIL-reactivity component regression.
- TESLA and the 2025 multimer screen: external domain-shift sets, opened once per later milestone.
- Curated vaccine trials: Event-B validation target.
- Sijbrandij: end-to-end acceptance test only; never fitted and never treated as a patient benchmark.

See [`docs/system_architecture.md`](docs/system_architecture.md) for the product boundary and
[`docs/data_workflow.md`](docs/data_workflow.md) for the legacy ingestion commands.
