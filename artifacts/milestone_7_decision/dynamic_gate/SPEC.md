# Dynamic upstream gate — specification (Milestone 7)

_Status: pre-registered design. Fitting/selection happens only on DEV under leave-one-cohort-out;
one candidate is frozen before the LOCKED external test. Companion: `FEASIBILITY.md` (empirical
go-signal), `FALSIFICATION_LEDGER.md` (failed variants)._

> **Post-review scoping correction.** This spec builds a v1 gate whose veto axes are {EL, PRIME}. Those
> are the dominant inputs to the downstream rankers, so v1 cannot change a top-20 they produce — a
> **structural tautology** (`CIRCULARITY_AUDIT.md`). v1's null downstream result is therefore scoped to
> **same-feature presentation gating**; it does NOT falsify the general orthogonal-feature gate. The
> orthogonal challenger — which targets the high-EL/high-PRIME decoy stratum using features NOT in the
> downstream rank — is specified in `V2_CONTRACT.md` and pre-registered in `V2_PREREGISTRATION.md`, and is
> **data-blocked** on WES/RNA reconstruction (Miller IPV `PRJNA980652`, Gartner).

## 0. Problem framing (what this is, and is NOT)

Build a **label-blind upstream gate** that removes a large fraction of tested-negative neoantigen
candidates while retaining essentially all *recognized positives*, then hands the survivors **unchanged**
to the existing rankers (genuine GfellerLab PRIME and frozen Epicurus v0.1). It is a **selective-prediction
/ safe-rejection** problem under extreme asymmetric cost: a false veto (dropping a real positive) is far
worse than retaining a negative. **Missing evidence defaults to KEEP.**

It is **NOT** a classifier, **NOT** a reranker, and it never reweights the downstream ranker (the frozen
expression policy — expression is confidence-only, never a rank penalty — is preserved: the gate removes
candidates but never changes the order of survivors, and expression only ever *rescues*, never vetoes on
its own).

Motivation (`FEASIBILITY.md`, pool-size diagnostic): oracle pruning to 25%-negative pools lifts Epicurus
hits@20 Gartner 0.808→1.652 / IMPROVE 1.230→3.214, but the incumbent pure-EL percentile gate retained only
66.7% / 35.7% of positives there. The feasibility probe shows an **AND-of-independent-vetoes (missing→KEEP)
Pareto-dominates the pure-EL gate on positive retention at matched negative-removal** — so there is real
headroom to close the oracle gap honestly.

## 1. Feature availability (verified; never fabricate a field)

Unified per-candidate features usable label-blind across the eval cohorts. **All within-patient
percentiles** (oriented higher = better), which harmonizes TPM vs decile encodings and removes cohort/study
scale as an input.

| feature | gartner | improve | multimer (in-sample) | checkmate153 (LOCKED) | role in gate |
|---|---|---|---|---|---|
| `el` (NetMHCpan-EL / MHCflurry-pres %) | Y | Y (RankEL) | Y | Y (MHCflurry pres) | veto axis |
| `prime` (genuine GfellerLab PRIME %rank) | Y | Y | Y | Y | veto axis |
| `expr` (RNA; decile/TPM/counts) | Y (decile) | Y | Y (TPM) | Y (counts) | **rescue-only** |
| predictor set for disagreement | 5 (EL,BA,MixMHC,MHCflurry,HLAthena) | 2 (EL,BA) | 2 (EL,BA) | 2 (MHCflurry,MixMHC) | uncertainty rescue |
| tumor VAF | decile only | N | N | Y(raw) | **absent in dev** → extension |
| mutant-RNA VAF / read support / ASE | N | N | N | N | **absent everywhere** → extension |
| agretopicity / foreignness / processing | N | N | Y | N | multimer-only; not cross-cohort |
| `expression_call=N`, `hla_loh_call` (Layer-0 biology) | N | N | N | N | **no cohort populates → silent no-op here** |

**Consequences, stated honestly:**
- The gate this milestone builds is an **availability-aware, peptide-/presentation-only gate**. It does
  **not** learn genomics; there is no raw WES depth, mutant-allele RNA VAF, purity/CCF/clonality, or
  proteasomal processing available across the eval cohorts.
- Layer-0's biological impossibility rules (`GENE_NOT_EXPRESSED`, `HLA_LOH_LOST_ALLELE`) are **inert** on
  every research cohort because none populate `expression_call`/`hla_loh_call`; they remain wired for the
  SHERPA/product path and become active the moment a real WES/RNA patient (Miller/Sid/RTTP) is passed.
- The **WES/RNA extension contract** (§7) specifies exactly which new fields turn IMPROVE's ~1/4
  unrescuable positives into rescuable ones.

## 2. Layered architecture

Every candidate flows through four layers; a KEEP at any protective layer is final (fail-open).

- **Layer 0 — deterministic impossibility gate** (`src/epicurus_neo/gates.py::apply_deterministic_gate`,
  reused unchanged). Removes only rule-verifiable impossibilities (malformed AA, impossible class-I length
  8–14, exact duplicate, mutation-not-in-peptide, lost allele, vendor-called unexpressed) with auditable
  reasons. Fail-open by design; inert where inputs absent.
- **Layer 1 — high-confidence negative veto** (AND-of-independent-vetoes). A candidate is *veto-eligible*
  iff **every present** presentation/PRIME axis is below the operating threshold `t`:
  `veto = (el present ∧ s_el < t) ∧ (prime present ∧ s_prime < t) ∧ (expr present ∧ s_expr < t)`.
  A **missing** axis abstains (does not vote to remove) ⇒ any missing feature yields KEEP. Because it is an
  AND, a candidate strong on *any* axis survives ⇒ high recall by construction.
- **Layer 2 — uncertainty / discordance rescue** (overrides a Layer-1 veto → KEEP):
  - **predictor disagreement**: spread across the available predictor percentiles `> d` ⇒ models disagree ⇒ keep;
  - **near-boundary**: `min_axis |s − t| < m` ⇒ keep (don't veto a coin-flip);
  - **coverage/OOD**: candidate missing any core axis, or patient below coverage/pool floors ⇒ keep.
- **Layer 3 — patient-adaptive controller**. Picks the most aggressive `t` that still satisfies a
  **calibrated positive-retention lower bound** (§3), with hard per-patient safety rails:
  - never veto a patient's **top-`M_floor`-by-EL** (presentation-best always survive);
  - if the veto would remove `> cap` of a patient's pool, or the pool `< pool_floor` candidates, or feature
    coverage `< cov_floor` ⇒ **fall back to permissive / keep-all** for that patient.

Cohort/study identity is **never** an input at any layer.

## 3. Calibration (conformal / Neyman–Pearson, distribution-light)

The single learned quantity is the operating threshold `t` (and optional `d`, `m`). Selection uses **only**
calibration positives, never eval-cohort positives:

1. Sweep `t` over a grid. For each `t`, on the calibration positives compute empirical retention `r(t)` and
   its **Clopper–Pearson lower bound** `LB(t)` at 95% confidence (n = #calibration positives).
2. `t* = max{ t : LB(t) ≥ target }` for target ∈ {0.90, 0.95, 0.975, 0.99} (Pareto operating points).
   This controls the false-veto rate on positives at level `1 − target` (the "risk"), maximizing negative
   removal (the "power") — a Neyman–Pearson constraint realized by conformal thresholding.
3. **Leave-one-cohort-out (LOCO)**: to evaluate cohort X, calibrate `t*` on the *other* dev cohorts only,
   then apply to X. This stresses transfer (positives are not exchangeable across cohorts, so LOCO is
   conservative). Also report leave-one-study-out within IMPROVE (RH-/Neye-/BC- sub-cohorts).
4. **Freeze**: pick ONE deployment config (target = 0.95, LOCO-calibrated on all three dev cohorts) →
   `configs/frozen/dynamic_gate_v1.json`, before touching the LOCKED test.

## 4. Data discipline

- **DEV (fit + LOCO selection + reporting):** gartner, improve, multimer. Multimer is frozen-Epicurus'
  training cohort ⇒ **in-sample, flagged**; used for gate calibration but never for a headline.
- **LOCKED EXTERNAL TEST:** CheckMate 153 (14 pts / 162 pos / ~1035 tested-neg; genuine PRIME + genuine
  tested-negatives; in no training or dev split). Scored **once** with the frozen config.
- Patient-level and leave-one-cohort/study-out evaluation only; cohorts never pooled into one number.
- **Leakage:** exact + near-duplicate (`near_duplicate`, threshold 0.8) peptide removal between calibration
  and eval; PRIME-training leakage mask on any PRIME-derived signal; repeated-antigen audit
  (`quarantined_recurrent_peptides` from the frozen split). Gate calibration uses no peptide that leaks into
  an eval cohort.
- **Labels:** only explicit POSITIVE / TESTED_NEGATIVE. PU cohorts (Müller, Sid, Zhao) are excluded from
  negative-removal claims; Sid is descriptive only, never tuned on.
- **Oracle** positive-retention (100% by construction) is used **only** as a ceiling.

## 5. Metrics

**Primary gate metrics** (per cohort; overall, per-patient, worst-patient):
- positive retention (mean / min / worst-patient) + **Clopper–Pearson lower bound**;
- negative removal fraction; fraction of patients losing ≥1 positive; bootstrap CI on retention.

**Pareto frontier:** achievable negative removal at retention targets 90 / 95 / 97.5 / 99 / 100%.

**Downstream consequence** (the real test — gate → *unchanged* ranker):
- paired recognized **hits@20** and **recall@20** for genuine PRIME and frozen Epicurus, LARGE (ungated) vs
  gated, via `benchmark.stats.paired_bootstrap` (20k) → `pre_registered_verdict`;
- saturation flags (pool ≤ 20 after gating makes top-20 trivially saturate — reported, not hidden).

**Baselines:** (a) incumbent pure-EL percentile gate at matched removal; (b) deterministic-only (Layer 0);
(c) single global EL threshold; (d) random removal at matched count; (e) full-positive-retention oracle
ceiling; (f) keep-all (ungated).

## 6. Predeclared safety bar (deployment candidate must satisfy ALL)

1. Worst-cohort **positive-retention CP lower bound ≥ 0.95** on DEV under LOCO (external cohorts;
   multimer in-sample excluded from the bar).
2. **No downstream regression:** gated recognized hits@20 CI upper bound for both PRIME and frozen Epicurus
   is ≥ the ungated value (gating must not *hurt* the ranker) on every external cohort.
3. Negative removal materially > 0 where the bar allows (report per cohort; a patient that cannot support
   removal keeps everything — that is a pass, not a failure).
4. On the LOCKED test, the frozen config reproduces (1)–(2) with no re-tuning.

## 7. WES/RNA extension contract (what new data unlocks)

The IMPROVE floor (≈1/4 of positives jointly low on EL, PRIME, expr) is **not** closable with peptide-only
features. To make those positives rescuable, the gate needs — as **rescue-only** axes inside Layer 1's AND
(never standalone vetoes) — per candidate:
- **mutant-allele RNA VAF** and **RNA read support** (allele-specific expression of the *mutant* transcript,
  not gene-level TPM) — the single highest-value addition;
- **tumor DNA VAF + read depth + purity/CCF** (clonality) — a clonal, well-covered mutation should never be
  vetoed on weak presentation alone;
- **proteasomal processing / cleavage / stability** (available already in multimer; needed cross-cohort);
- **agretopicity** (mut-vs-WT binding ratio) as an uncertainty-rescue signal.

Contract for Miller IPV / Gartner reconstruction (open WES+RNA): emit these fields per candidate in the
unified schema with an explicit presence mask; the gate consumes them additively (each is a KEEP-only
axis). No field is ever imputed to a value that would *enable* a veto.

## 8. Deliverables & iteration protocol

- `src/event_b/dynamic_gate.py` — gate + calibration + eval (pure, testable).
- `tests/test_dynamic_gate.py` — safety invariants (missing→KEEP, AND-structure, monotonic removal in `t`,
  CP-LB correctness, cohort-not-an-input, top-`M_floor` floor honored).
- `scripts/dynamic_gate.py` — runner: load cohorts, LOCO calibrate, freeze, evaluate, write artifacts.
- `configs/frozen/dynamic_gate_v1.json` — the one frozen config.
- artifacts: `REPORT.md`, `pareto.csv`, `per_patient.csv`, `dynamic_gate.json`, `FALSIFICATION_LEDGER.md`.
- CLI: `epicurus-neo dynamic-gate` subcommand.
- Verdict at the end: **A** clears safety + improves downstream top20 / **B** promising but
  data-insufficient / **C** falsified — with the exact WES/RNA data that would unlock the next level.
