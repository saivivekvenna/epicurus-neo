# PREREGISTERED PROTOCOL — osteosarc.com (Sid) public reconstruction

**Frozen pre-fit.** Canonical design: `docs/superpowers/specs/2026-07-12-osteosarc-sid-reconstruction-preregistration.md`
(committed together with this file). This is the machine-checkable contract Phase B must satisfy **before** any
model is fit, tuned, or compared. Corrects commit `dd3efd1` (superseded — see spec §0/§11).

## Frozen site invariants (build aborts on any mismatch)
- `site_index` (`https://osteosarc.com/variants/`) has **182** variant rows (182 `/variant/` links).
- rows with `data-vaccines > 0` = **44**.
- rows with `data-elispot > 0` (0/1 flag) = **14**.
- variant pages resolve only at trailing-slash URL `…/variant/<GENE-chrN-POS>/` (bare path → HTTP 308).

## Frozen label states
`POSITIVE_STRONG | POSITIVE_WEAK | POSITIVE | NEGATIVE | AMBIGUOUS | UNTESTED`
Mapping on raw result **text** (class preserved separately): `positive (strong[**])`→STRONG,
`positive (weak)`→WEAK, bare `positive`→POSITIVE, `negative`→NEGATIVE, else→AMBIGUOUS; vaccine-included peptide
with zero experiments→UNTESTED. Never invent NEGATIVE from absence.

## Frozen resolution states
`INDIVIDUAL_PEPTIDE | MUTATION_LONG_PEPTIDE | POOL | MUTATION_TCR | UNKNOWN`
Non-blank pool id → POOL; long-peptide block, blank pool, no individual proof → MUTATION_LONG_PEPTIDE; explicit
individual minimal-epitope proof → INDIVIDUAL_PEPTIDE; Hudson expander → MUTATION_TCR; `—`/`NA`/blank without
proof → UNKNOWN. `—`/`NA` never imply individual testing. Pool-positive never propagates to member peptides.

## Frozen tables (column order fixed; see spec §2)
- `variant_catalog.csv` — 1/variant (182); index counts cross-checked against page pills.
- `peptide_inventory.csv` — 1/peptide block; `parsed==declared` experiment count, `len==aa` enforced.
- `assay_ledger.csv` — 1/ELISPOT experiment row; raw class **and** text preserved; SHA1 `experiment_key` unique.
- `hudson_tcr_labels.csv` — 1/(timepoint,TRB,mutation); `MUTATION_TCR`; SEPARATE stream, never merged.
- `reachability_funnel.csv` — per target; joins by `(chrom,pos,ref,alt)`/protein, never gene-only.

## Frozen integrity checks (spec §4)
182/44/14; per-block `parsed==declared` & `len==aa`; unique dedup keys; enum membership; coordinate-level
reachability with MAP2/DYNC1H1 dual-coordinate disambiguation; all 3 Hudson positives traced with a non-null
`first_failure_stage`; contradictions reported (never collapsed); network-free rerun from cache is byte-identical.

## Guardrails
No fit/tune/threshold/ranking-gate this milestone. Frozen Epicurus config untouched. `dd3efd1` single-positive
AUROC is descriptive/superseded. Improvement ideas validated later on multimer/Gartner/IMPROVE/CheckMate — never
on this n=1 patient. Emitted tables are derived public data (committable); raw `site_cache/` stays gitignored.

## Status
Phase A (this design) committed first. Phase B (`variant_catalog/peptide_inventory/assay_ledger/
hudson_tcr_labels/reachability_funnel.csv` + `AUDIT.json` + `REPORT.md` + `PROVENANCE.json`) lands as a separate
scoped commit and is reviewed before any downstream use.
