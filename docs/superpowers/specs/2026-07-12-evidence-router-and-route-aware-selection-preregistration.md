# Evidence router + route-aware top-20 selection — preregistration

**Design source of truth for the inference-time candidate *evidence router* and the *constrained,
route-aware top-20 selection* it feeds.** This is a **pre-fit policy design**: it freezes the
inference-time routing rules, the route precedence, the multi-source union identity rules, and the
top-20 selection policy **before route-aware selection is evaluated or compared on any cohort**.
The route hypothesis is explicitly informed by the already-read Sid structural audit, so Sid is not
an independent validation cohort for the conception of this router. Paired frozen copy:
`artifacts/milestone_7_decision/evidence_router/PREREGISTERED_PROTOCOL.md`; frozen machine-readable
parameters: `configs/frozen/evidence_router_v1.json`.

> **No route-aware selection metric or independent-cohort comparison is computed until this policy is
> committed.** Sid's already-known structural losses motivate the design but are not used to choose the
> frozen reserve count, precedence, or score. The frozen Epicurus v0.1 config
> (`configs/frozen/epicurus_v0_1.json`) is not read or modified. The legacy deterministic gate
> (`src/epicurus_neo/gates.py::apply_deterministic_gate`) is **preserved unchanged** for backward
> compatibility; the router is an additive v2 path.

North Star (unchanged): maximize the chance a genuinely recognized neoantigen survives
`WES/RNA/HLA → candidate generation → eligibility → top-20`, compared to genuine PRIME where a fair
comparison exists. PRIME is a comparator/ranker component, never a recall mechanism: it cannot recover
a candidate filtered before it can score.

---

## 0. What this changes and why (the Sid reachability finding, commit `6317fc8`)

The reconstructed osteosarc.com/Sid ledger shows that of the 3 IFNγ/TCR-recognized targets, **2 (ASPM,
MAP2) are lost at candidate *generation*** — they are called by DRAGEN/Sarek/oncoanalyser but are
**absent from the single pVACtools 2025.01 candidate universe** — and only DYNC1H1 reached the
shortlist. The dominant, recoverable loss is **candidate recall upstream of any ranker**, not ranking.

Two failure modes in the legacy path conflate *upstream generation gaps* and *cross-sectional RNA
absence* with *biological impossibility*, and would silently drop recoverable targets:

1. `apply_deterministic_gate` **hard-removes `GENE_NOT_EXPRESSED`** (`expression_call ∈ {N,…}`) and the
   downstream `_exclusion_reason` marks `NO_RNA_EXPRESSION` (`expression_tpm ≤ 0`) /
   `NO_MUTANT_RNA_SUPPORT` ineligible. Cross-sectional RNA absence at one timepoint/region is **not**
   biological impossibility (a variant can be expressed in tumor and absent from the sampled RNA).
2. The legacy `MALFORMED_AA` rule treats an **empty peptide** as malformed and removes it, which would
   count an upstream candidate-generation gap (peptide never generated) as a validity impossibility.

The router fixes both: it hard-removes **only route-verifiable impossibilities**, keeps RNA-absent /
atypical-class / single-caller candidates **eligible but flagged**, and reports "no peptide/HLA" as a
distinct **`NEEDS_PEPTIDE_GENERATION`** upstream status — never as a ranker failure.

---

## 1. Scope and non-goals

**In scope (this preregistration + the implementation it licenses):** an additive inference-time
evidence router; a generic multi-source variant/candidate union helper; a constrained route-aware
top-20 selection; a machine-readable per-stage funnel; and their tests. Evaluation is licensed only
after the code commit (§7–§8).

**Non-goals / prohibited:** no model fit, tune, retrain, or threshold search on any cohort; **no tuning
of any policy constant on Sid outcomes** (Sid is locked diagnostic evaluation); no edit to
`epicurus_v0_1.json`; no mutation of legacy v1 (`apply_deterministic_gate`, `score_product_candidates`)
outputs; no claim that route-aware selection improves immunogenicity/response absent set-level labels.

---

## 2. Legacy gate preserved (backward compatibility)

`src/epicurus_neo/gates.py::apply_deterministic_gate` and `src/epicurus_neo/product.py`'s v1 scoring
remain byte-identical and continue to emit `GENE_NOT_EXPRESSED` / `NO_RNA_EXPRESSION`. The router is a
**new module** (`src/epicurus_neo/evidence_router.py`) that *composes* the v1 scorer's public output
(it consumes `epicurus_lower_evidence_score` and the normalized candidate columns) and adds v2 columns
(`primary_route`, orthogonal `flag_*`, `rankable`, `router_removed_reason`). v1 columns are never
overwritten. Existing tests for v1 must stay green.

---

## 3. The evidence router (frozen inference-time rules)

Applied to normalized product candidates (post `normalize_product_candidates`). All rules are
deterministic and first-wins where an order is specified.

### 3.1 IMPOSSIBLE — hard-remove, route-verifiable only (first-wins order)

A candidate is removed (`primary_route = IMPOSSIBLE`, `eligible = False`) iff exactly one of, in order:

1. `MALFORMED_AA` — `mutant_peptide` is **non-empty** and contains a non-standard amino acid
   (outside `ACDEFGHIKLMNPQRSTVWY`). *An empty peptide is NOT malformed* — it is
   `NEEDS_PEPTIDE_GENERATION` (§3.5), never removed.
2. `BAD_CLASS_I_LENGTH` — `mhc_class` is class I **and** `len(mutant_peptide) ∉ [8, 14]`
   (`CLASS_I_MIN_LEN=8`, `CLASS_I_MAX_LEN=14`, identical to the legacy gate; only genuinely impossible
   lengths, keeping documented bulged 12–14mers).
3. `MUT_NOT_IN_PEPTIDE` — `source_variant_type ∈ {SNV, SNP, MISSENSE}` **and** the parsed missense
   mutant residue (`p.<wt><pos><mut>`) is absent from `mutant_peptide`. Runs **only** for missense; a
   frameshift/indel/splice peptide legitimately need not contain a single substituted residue and is
   never removed by this rule.
4. `DUP_CANDIDATE` — duplicate identity on `(patient_id, mutation_id, mutant_peptide, hla_allele)`
   (when ≥3 of those columns are present), `keep=first`.
5. `HLA_LOH_LOST_ALLELE` — `hla_loh_call ∈ {Y, YES, TRUE, 1, LOST}`: the peptide→HLA presentation route
   is through a definitely lost allele. Fires only when an allele is present.

**This IMPOSSIBLE set is the legacy gate's rules minus `GENE_NOT_EXPRESSED`, minus the empty-peptide
arm of `MALFORMED_AA`.** Those two removals are demoted to flags/status below.

### 3.2 NEVER hard-removed (demoted to flags/routes)

None of the following, on its own, ever removes a candidate: `expression_tpm == 0`;
`expression_call ∈ {N,…}`; `rna_mutant_reads == 0`; missing RNA entirely; variant class
frameshift/indel/splice/inframe; single-caller status. Cross-sectional RNA absence and atypical class
are **evidence flags**, not impossibilities.

### 3.3 Orthogonal flags (independent booleans; all emitted)

- `flag_atypical_variant_class` — `source_variant_type` present and **not** in `{SNV, SNP, MISSENSE}`
  (frameshift, indel, ins, del, splice, inframe, …).
- `flag_weak_or_absent_rna` — `(expression_tpm == 0)` **or** `(expression_call ∈ {N,…})` **or**
  `(rna_mutant_reads == 0)`, for the columns that are present.
- `flag_missing_rna` — none of `{expression_tpm, rna_vaf, rna_mutant_reads, expression_call}` is
  present/populated.
- `flag_single_caller` — union provenance shows exactly one caller/pipeline.
- `flag_multi_source_support` — union provenance corroborates across ≥2 distinct
  callers/timepoints/regions (secondary-evidence recovery).
- `flag_has_presentation` — a presentation signal is available (`presentation_score`,
  `binding_percentile_rank`, or `binding_affinity_nm`).
- `flag_conflicting_evidence` — a representation/label conflict recorded by the union (same identity,
  differing annotation) or contradictory recognition markers.
- `flag_needs_peptide_generation` — no `mutant_peptide` **or** no `hla_allele` (see §3.5).

### 3.4 Primary route + precedence (deterministic, first-wins)

Exactly one `primary_route` per candidate. IMPOSSIBLE is decided first (§3.1). Among **eligible**
candidates the precedence is:

1. **RESCUE** — `(flag_atypical_variant_class OR flag_weak_or_absent_rna)` **AND**
   `(flag_has_presentation OR flag_multi_source_support)`. *Atypical variant class or weak/absent RNA
   with other genomic/presentation support* — the candidates the legacy gate would have dropped but
   that carry independent support. Highest non-IMPOSSIBLE precedence so they are surfaced.
2. **LONGITUDINAL** — `flag_multi_source_support` and not already RESCUE. *Recovery/support from a
   secondary timepoint/region/caller.*
3. **UNCERTAIN** — missing or conflicting core evidence: `flag_missing_rna` **or**
   `flag_conflicting_evidence` **or** `flag_needs_peptide_generation` **or** `NOT flag_has_presentation`.
   *We cannot confidently place it as CORE.*
4. **CORE** — default: conventional supported candidate (missense, has presentation, not weak/absent
   RNA, no conflict).

Rationale for precedence: RESCUE > LONGITUDINAL because "would-have-been-dropped-but-supported" is the
most actionable class to protect; UNCERTAIN > CORE because absence of presentation or presence of
conflict disqualifies the CORE label. Flags are orthogonal and all retained, so a RESCUE candidate that
is also multi-caller still exposes `flag_multi_source_support`.

### 3.5 Candidate-generation reachability vs reranking

A candidate is **`rankable`** iff it is eligible **and** has both a non-empty `mutant_peptide` and a
non-empty `hla_allele`. If either is missing, `flag_needs_peptide_generation = True`,
`rankable = False`, `router_status = NEEDS_PEPTIDE_GENERATION`. Such a candidate is **an upstream
candidate-generation gap** and is **never** counted as a PRIME/ranker miss: PRIME may score only
rankable candidates. The funnel (§6) reports these as `generated` and `valid` but not `rankable`.

---

## 4. Multi-source variant/candidate union (generic helper)

`src/epicurus_neo/variant_union.py::union_variants` merges candidate/variant rows from multiple
sources into one deduplicated frame:

- **Identity key** (in priority order, never gene-only and always patient-scoped when a patient column
  exists): `(patient_id, genome_build, chrom, pos, ref, alt)` if the genomic fields are present
  (`genome_build` included when available), else `(patient_id, gene_symbol, exact normalized
  protein_change/mutation_id)`. If only a gene symbol is available, rows are **not**
  merged (kept distinct) — a gene-only union is provably wrong (MAP2 and DYNC1H1 each carry two
  distinct coordinates; MAP2's two frameshift coordinates 4 bp apart are *distinct* keys and stay
  separate rows).
- **Provenance preserved and aggregated**: `callers`, `timepoints`, `regions`, `sources` become sorted
  unique sets; `n_callers`/`n_timepoints`/`n_regions` counts drive `flag_single_caller` /
  `flag_multi_source_support`.
- **Representation conflicts preserved**: if the same identity key carries differing annotations
  (e.g. `protein_change`, `consequence`), the distinct values are recorded in a `representation_conflicts`
  field and never silently collapsed; the row is flagged `flag_conflicting_evidence`.
- **Peptide/HLA availability**: if a unioned candidate lacks a peptide or HLA, its status is
  `NEEDS_PEPTIDE_GENERATION` (§3.5) — the helper never fabricates a peptide or pretends the row can be
  ranked.

The union feeds the router (provenance → multi-source flags) and the reachability funnel (the
multi-caller raw variant union is the candidate-generation-recall denominator, §7).

---

## 5. Constrained, route-aware top-20 selection (frozen)

`src/epicurus_neo/evidence_router.py::select_route_aware_topk`. Selects up to `k = 20` from the
**eligible & rankable** candidates of one patient, deterministically.

**Policy (frozen; `configs/frozen/evidence_router_v1.json`):**
- `k = 20`.
- **Exploration reserve (modest):** for each non-CORE route in `{RESCUE, LONGITUDINAL, UNCERTAIN}`
  that has ≥1 eligible+rankable candidate, reserve `reserve_per_route = 1` slot, capped at
  `max_reserve = 3` total and always `< k` (so CORE keeps ≥ `k − max_reserve = 17` of the slots). This
  **guarantees representation of every available non-CORE evidence route** without letting exploration
  dominate.
- **Fill order:** (a) fill each reserved route with its highest-incumbent-score eligible+rankable
  candidate, respecting diversity caps; (b) fill all remaining slots by descending **frozen incumbent
  score** across all eligible+rankable candidates, respecting diversity caps and never re-selecting.
- **Incumbent score:** a caller-supplied score column (never re-fit here) — the frozen ranker used for
  that cohort: **genuine PRIME where available**, else the frozen Epicurus product score
  (`epicurus_lower_evidence_score`). The selection is score-agnostic; it only orders and reserves.
- **Diversity caps (existing):** `max_per_mutation = 2`, `max_per_gene = 4`, `max_per_hla = None`.
- **Graceful backfill:** if a reserved route is absent (or its only candidates are exhausted by caps),
  the freed slot returns to the score-fill pool — the final set always has `min(k, n_rankable)` rows.
- **Determinism:** ties broken by `md5(mutant_peptide | hla_allele)`, `mergesort`, stable across input
  permutation.

**Guardrail:** this reserves *representation*, spreading selection across evidence routes and diversity
axes. Absent set-level outcome labels it makes **no claim** to improve immunogenicity or response — it
only guarantees non-CORE routes are not starved by a pure score-sorted top-20.

---

## 6. Machine-readable per-stage funnel

`evidence_router` emits a per-patient funnel separating four stages, in order:

`generated` (all input rows) → `valid` (eligible; passed IMPOSSIBLE) → `rankable` (valid **and** has
peptide+HLA) → `selected` (in the route-aware top-k).

Per patient the report records: the four stage counts; `router_removed_reason` counts (IMPOSSIBLE
breakdown); `needs_peptide_generation` count; `route_composition_valid` and `route_composition_selected`
(counts by `primary_route`); and the policy id. Written as JSON + Markdown, alongside a CSV of the
routed candidates. The stage semantics make it impossible to read an upstream generation gap as a
ranker miss.

---

## 7. Acceptance metrics (frozen definitions)

- **Candidate-generation recall of known positives** = fraction of known recognized positives present
  in the multi-caller **raw variant union** (before peptide generation). Reported only where the cohort
  actually carries a raw multi-caller denominator (Sid does; the reranker cohorts largely do not — see
  §8). Never claim generation recall where the raw denominator is absent.
- **Ranker hits@20 conditional on reachability** = hits@20 computed over the subset of positives that
  are **rankable** (peptide+HLA present), so a generation gap is not charged to the ranker.
- **End-to-end hits@20** = hits@20 over all known positives, generation gaps included (the honest
  deployment number; will be ≤ conditional).
- **Route composition** = distribution of `primary_route` over valid and over selected candidates.
- **No-regression** = on each independently measured cohort, router+route-aware-selection top-20 vs the
  static legacy-gate / pure-score top-20 must not lose hits@20 beyond noise (bootstrap CI); the router
  is a recall-preserving change and must not reduce measured hits.
- **Bootstrap CIs** = patient-level (or positive-level) bootstrap CIs on paired hits@20 deltas wherever
  patient-level labels permit; report `CONSISTENT_WITH_NO_EFFECT` when a CI spans 0, never as
  equivalence.

Mechanical reachability (a target present in the union / rankable / selected) is reported **separately**
from any learned-recognition-superiority claim, and neither is used to argue the other.

---

## 8. Dataset allocation (frozen)

- **osteosarc.com / Sid** — **hypothesis-generating reachability + ledger diagnostic ONLY.** Its
  structural audit informed the router, so it is not independent validation. After the code commit it is
  replayed without changing any policy constant/quota. Reports candidate-generation recall of
  the 3 Hudson-recognized targets in the multi-caller raw union; which are peptide-generated/rankable;
  which are selected. It is expected and acceptable to report ASPM/MAP2 as `NEEDS_PEPTIDE_GENERATION`
  until real peptide generation is run.
- **CD8 multimer / Gartner NCI / IMPROVE / CheckMate 153** — conditional reranker + no-regression, each
  respecting its known denominator limits (multimer = fully assayed tested-negatives; Gartner = tested
  vs untested three-state, tested-only; IMPROVE = pre-screened candidate denominator, not full somatic;
  CheckMate = small/external). Use existing local genuine-PRIME artifacts; no new fit.
- **RTTP** — label-free deployment/demonstration only (predicted immunogenicity, no measured label).
- **No new untouched cohort is fabricated.** Superiority over PRIME is claimed only if a frozen model
  wins on an untouched cohort — not licensed here.

---

## 9. Tests (frozen list; TDD, written before implementation)

`tests/test_evidence_router.py` + `tests/test_variant_union.py`:

1. **Sid-like ASPM/MAP2 rescue routing** — an atypical-class / weak-RNA candidate *with* presentation
   or multi-caller support routes `RESCUE`, is retained (eligible), not removed.
2. **DYNC reachability** — a conventional supported SNV with peptide+HLA is `CORE`, `rankable`, and
   selectable.
3. **Expression-N retained in RESCUE** — `expression_call = N` / `expression_tpm = 0` is **not**
   removed; it is flagged and routed (RESCUE with support, else UNCERTAIN).
4. **Lost allele → IMPOSSIBLE** — `hla_loh_call = Y` removes the candidate (`HLA_LOH_LOST_ALLELE`).
5. **Zero mutant RNA reads not excluded** — `rna_mutant_reads = 0` flags but never removes.
6. **Empty peptide is NEEDS_PEPTIDE_GENERATION, not MALFORMED_AA** — empty peptide → not removed,
   `rankable = False`, upstream status; a genuinely bad-AA non-empty peptide → removed.
7. **Coordinate-level union** — rows sharing `(chrom,pos,ref,alt)` merge with provenance aggregated;
   MAP2's two coordinates 4 bp apart stay **distinct**; gene-only rows are not merged.
8. **Route reserves / backfill / determinism** — reserve guarantees ≥1 slot per present non-CORE route;
   absent routes backfill by score; selection is permutation-invariant and deterministic; CORE keeps
   ≥ `k − max_reserve`.
9. **No-regression harness sanity** — with only CORE candidates, route-aware selection reproduces the
   pure-score top-k exactly (reserve is a no-op when non-CORE routes are absent).

Plus: v1 product/gate tests remain green (the router is additive).

---

## 10. Guardrails / out of scope

Never convert UNTESTED→negative, pool-positive→peptide-positive, or vaccine-inclusion→negative. Never
fit on Sid or choose quotas after seeing its labels. Never claim Epicurus beats PRIME unless a frozen
model is externally established on an untouched cohort. Do not touch `epicurus_v0_1.json` or unrelated
dirty files. The legacy gate stays for backward compatibility.

## 11. Registered deviations
_(appended chronologically; the design is fixed above before any route-aware evaluation is run.)_

- **D1 (integration form).** The instruction asked for "small integration into product.py." In this
  working tree `src/epicurus_neo/product.py`, `gates.py`, and `portfolio_selection.py` are **untracked**
  (not in HEAD) alongside ~30 other untracked src files. To honor "do not mutate/falsify legacy v1
  outputs" and "do not touch unrelated dirty files," the integration is delivered as a **new module that
  composes product.py's public API by import**, and product.py/gates.py/portfolio_selection.py are
  **neither edited nor committed** by this milestone. Behavior is identical to a thin product.py hook;
  only the file scope of the commit differs. Recorded here for transparency.
- **D2 (union-key safety clarification, before implementation/results).** Every union identity is
  patient-scoped when `patient_id` exists, and the non-coordinate fallback includes both gene and exact
  protein/mutation identity. This prevents cross-patient hotspot merging and cross-gene `p.V600E`-style
  collisions. It narrows the originally ambiguous phrase "exact protein/mutation key" without changing
  any route or selection constant.
