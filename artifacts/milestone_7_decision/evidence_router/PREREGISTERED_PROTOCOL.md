# PREREGISTERED PROTOCOL — evidence router + route-aware top-20 selection

Frozen paired copy of `docs/superpowers/specs/2026-07-12-evidence-router-and-route-aware-selection-preregistration.md`.
Frozen parameters: `configs/frozen/evidence_router_v1.json`. **Committed before any route-aware
selection result or independent-cohort comparison is computed.** The already-read Sid structural
audit motivated the router; Sid is not presented as independent validation of the design.

## Fixed decisions (pinned)

1. **Hard-remove only route-verifiable impossibilities** (first-wins): `MALFORMED_AA` (non-empty
   peptide with a non-standard AA), `BAD_CLASS_I_LENGTH` (class I length ∉ [8,14]), `MUT_NOT_IN_PEPTIDE`
   (missense residue absent from a missense peptide only), `DUP_CANDIDATE`, `HLA_LOH_LOST_ALLELE`.
2. **Never hard-remove for**: `expression_tpm==0`, `expression_call=N`, `rna_mutant_reads==0`, missing
   RNA, frameshift/indel/splice, single-caller. These become orthogonal flags/routes.
3. **Empty peptide ≠ MALFORMED_AA.** Empty peptide (or missing HLA) → `NEEDS_PEPTIDE_GENERATION`,
   `rankable=False`, upstream gap — never removed, never charged to the ranker.
4. **Orthogonal flags** emitted for every candidate; **exactly one `primary_route`** with precedence
   `IMPOSSIBLE → RESCUE → LONGITUDINAL → UNCERTAIN → CORE`.
   - RESCUE = (atypical class OR weak/absent RNA) AND (has presentation OR multi-source support).
   - LONGITUDINAL = multi-source (caller/timepoint/region) recovery, not already RESCUE.
   - UNCERTAIN = missing/conflicting evidence, needs-peptide, or no presentation.
   - CORE = conventional supported (missense, presented, RNA-supported, no conflict).
5. **Candidate generation reachability ≠ reranking.** `rankable` requires peptide+HLA. Missing-peptide
   is an upstream generation gap (`NEEDS_PEPTIDE_GENERATION`), reported in the funnel, never a PRIME miss.
6. **Multi-source union**, patient-scoped when available, by `(genome_build,chrom,pos,ref,alt)` else
   `(gene, exact protein/mutation)` — **never gene-only**;
   populated peptide/HLA fields extend the key so distinct candidate routes from one variant remain
   distinct;
   provenance (callers/timepoints/regions/sources) aggregated; representation conflicts preserved &
   flagged; MAP2's two 4-bp-apart coordinates stay distinct.
7. **Constrained route-aware top-20**: k=20; modest exploration reserve of `reserve_per_route=1` per
   present non-CORE route capped at `max_reserve=3`; remaining slots by frozen incumbent score; diversity
   caps `max_per_mutation=2 / max_per_gene=4 / max_per_hla=None`; graceful backfill; deterministic
   `md5(peptide|hla)` tie-break. No immunogenicity claim without set-level labels.
8. **Funnel** separates `generated → valid → rankable → selected` (machine-readable JSON + Markdown + CSV).
9. **Acceptance metrics**: candidate-generation recall (only where a raw denominator exists), ranker
   hits@20 conditional on reachability, end-to-end hits@20, route composition, no-regression on measured
   cohorts, bootstrap CIs where labels permit. Mechanical reachability reported separately from any
   learned-superiority claim.
10. **Dataset allocation**: Sid = hypothesis-generating reachability/ledger diagnostic only (no policy
    tuning after freeze; not independent validation); multimer /
    Gartner / IMPROVE / CheckMate = conditional reranker + no-regression within their denominator limits;
    RTTP = label-free deployment. No fabricated untouched cohort. No superiority claim over PRIME here.

## Guardrails
Never UNTESTED→negative, pool-positive→peptide-positive, vaccine-inclusion→negative. No fit/tune on Sid.
No edit to `epicurus_v0_1.json`. Legacy gate preserved.

## Deviations
- **D1** — integration delivered as a new module composing product.py by import (product.py/gates.py/
  portfolio_selection.py are untracked in this tree and are neither edited nor committed here). See spec §11.
- **D2** — union identities are explicitly patient-scoped; protein/mutation fallback includes gene.
  This closes cross-patient hotspot and cross-gene notation collisions before implementation/results.
- **D3** — candidate rows extend the base variant key by exact peptide–HLA identity, preventing
  coordinate-level collapse of multiple candidates generated from one mutation.
