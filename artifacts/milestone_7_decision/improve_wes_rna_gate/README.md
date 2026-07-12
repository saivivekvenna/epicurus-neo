# Biology-first WES/RNA gate on IMPROVE

Isolated investigation: can an audited WES/RNA biology **gate** (label-blind
demotion on the frozen Epicurus base order, freed slots backfilled by base order,
reranker untouched) raise experimentally recognized hits@20 on IMPROVE?

## Files
- `WES_RNA_GATE_MECHANISM.md` — mechanism report: reachability ceiling, within-patient
  partial effects (expression / RNA / clonality / HLAexp interactions), pre-declared
  gate rules with matched-random control + leave-one-cancer-cohort-out + per-partition
  stability, ascertainment cross-check, verdict.
- `MECHANISM_RESULTS.json` — full machine-readable results.
- Code: `scripts/improve_wes_rna_gate.py` (runner), `src/benchmark/improve_wes_rna_gate.py`
  (pure tested core), `tests/test_improve_wes_rna_gate.py` (8 tests).

## Reproduce
`.venv/bin/python scripts/improve_wes_rna_gate.py`

## Bottom line
**No deployable gate.** There is real headroom (111 promotable positives at ranks
21–60), but every pre-declared biology demotion gate is net-**negative** and none
beats matched-random removal — within the frozen top-20 boundary, positives fail
the RNA/expression prerequisites as often as decoys. The marginal HLA-expression
signal dissolves within-patient (a confound). Mechanism: (1) **ascertainment,
verified** — IMPROVE's expression slope is ~0 at the boundary while the same axis
is clearly positive in Gartner (0.039); IMPROVE's ~200-candidate denominator was
pre-screened on expression, so the residual boundary is recognition-limited, not
presentation-limited; (2) **clonality > expression** — the within-patient VAF×expr
2×2 puts clonal/low-expr highest and subclonal/low-expr lowest, and DNA VAF is the
only axis whose weak positive direction is consistent across IMPROVE (+0.016) and
Gartner (+0.011). One candidate direction for **external** testing (not established):
a clonality/truncality PROMOTE-side prior on a NON-pre-screened full-mutanome
denominator (Miller re-enumeration / Gartner), with the same control discipline.

Uses only audited deployable pre-outcome primitives (commit 47b2064). Excludes
outcome, identity, PrioScore, IB_CB*, NetMHCExp, circular composites. Observational
— no causal claims. Touches no `dynamic_gate` / other base-anchored gate files.
