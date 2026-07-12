"""Four-arm generation x scorer patient-level top-k benchmark harness.

Factors the neoantigen pipeline into a **generation** stage and a **scoring/selection** stage
and compares four arms per patient:

  ==================  ======================  ==============  ==================
  arm_id              generation              scorer          selection
  ==================  ======================  ==============  ==================
  pvac_prime          standard pVAC           genuine PRIME   plain top-k
  lossless_prime      lossless-gen union      genuine PRIME   plain top-k
  lossless_epicurus   lossless-gen union      Epicurus        plain top-k
  full_epicurus       lossless-gen union      Epicurus        route-aware
  ==================  ======================  ==============  ==================

Design invariants (all enforced by tests):

* **Strict label isolation.** Measured positives are passed as a *separate* set of mutation ids
  and only ever consulted to score coverage AFTER a ranking is fixed. The scoring frame never
  carries a label column.
* **Explicit NOT_EVALUABLE.** Each arm declares its input requirements. When a requirement is
  missing for a patient/cohort the arm returns ``evaluable=False`` with the missing keys, never a
  silent number. Attribution is only computed when all four arms are evaluable.
* **Additive stage attribution.** ``full_epicurus - pvac_prime`` decomposes exactly into
  generation (arm2-arm1) + scorer (arm3-arm2) + selection (arm4-arm3).
* **Deterministic.** Ties break on ``md5(mutant_peptide|hla_allele)`` (shared with the frozen
  evidence router) so results are permutation-invariant.

Metric granularity is the **mutation**: a positive counts as a top-k hit when at least one of its
candidate rows (peptide x HLA) is selected into the top-k. This matches the measured recognition
labels, which are per-mutation, and the frozen recovery diagnostic.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field

import pandas as pd

from epicurus_neo.evidence_router import (
    DEFAULT_ROUTER_POLICY,
    RouterPolicy,
    route_candidates,
    select_route_aware_topk,
)

# ---------------------------------------------------------------------------
# Requirement keys (a cohort/patient "has" a subset of these)
# ---------------------------------------------------------------------------
REQ_LABELS = "measured_labels"
REQ_PVAC = "pvac_candidates"
REQ_LOSSLESS = "lossless_generation"
REQ_PRIME = "genuine_prime"
REQ_EPICURUS = "epicurus_features"
REQ_ROUTER = "router_features"

DEFAULT_K = 20


@dataclass(frozen=True)
class ArmSpec:
    arm_id: str
    generation: str  # "pvac" | "lossless_union"
    scorer: str  # "genuine_prime" | "epicurus"
    selection: str  # "plain_topk" | "route_aware"
    label: str


FOUR_ARMS: list[ArmSpec] = [
    ArmSpec("pvac_prime", "pvac", "genuine_prime", "plain_topk",
            "Standard pVAC candidates + genuine PRIME (incumbent baseline)"),
    ArmSpec("lossless_prime", "lossless_union", "genuine_prime", "plain_topk",
            "Lossless-generation union + genuine PRIME (isolates the generation/recall stage)"),
    ArmSpec("lossless_epicurus", "lossless_union", "epicurus", "plain_topk",
            "Lossless-generation union + Epicurus (isolates the scorer swap)"),
    ArmSpec("full_epicurus", "lossless_union", "epicurus", "route_aware",
            "Lossless-generation union + Epicurus + route-aware selection (full stack)"),
]


def arm_requirements(arm: ArmSpec) -> tuple[str, ...]:
    """Input requirement keys for an arm. Missing any -> NOT_EVALUABLE."""
    reqs = [REQ_LABELS, REQ_PVAC]
    if arm.generation == "lossless_union":
        reqs.append(REQ_LOSSLESS)
    reqs.append(REQ_PRIME if arm.scorer == "genuine_prime" else REQ_EPICURUS)
    if arm.selection == "route_aware":
        reqs.append(REQ_ROUTER)
    return tuple(reqs)


@dataclass(frozen=True)
class ArmEligibility:
    arm_id: str
    evaluable: bool
    missing: list[str]


def evaluate_eligibility(
    available: set[str], arms: list[ArmSpec] = FOUR_ARMS
) -> dict[str, ArmEligibility]:
    """Which arms are evaluable given the set of requirement keys a cohort/patient satisfies."""
    out: dict[str, ArmEligibility] = {}
    for arm in arms:
        missing = [r for r in arm_requirements(arm) if r not in available]
        out[arm.arm_id] = ArmEligibility(arm.arm_id, not missing, missing)
    return out


# ---------------------------------------------------------------------------
# Coverage container (n of `of`, with the covered ids)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class Coverage:
    n: int
    of: int
    ids: list[str]

    @property
    def frac(self) -> float | None:
        return round(self.n / self.of, 4) if self.of else None


@dataclass
class ArmResult:
    arm_id: str
    evaluable: bool
    missing: list[str] = field(default_factory=list)
    n_positives: int | None = None
    generation_recall: Coverage | None = None
    rankable_recall: Coverage | None = None
    top_k: Coverage | None = None
    n_selected: int | None = None

    # convenience accessors used by attribution / reporting
    @property
    def hits_at_k(self) -> int | None:
        return self.top_k.n if self.top_k is not None else None

    @property
    def recall_at_k(self) -> float | None:
        return self.top_k.frac if self.top_k is not None else None


# ---------------------------------------------------------------------------
# Epicurus scoring (frozen v0.1 formula) attached to a universe frame
# ---------------------------------------------------------------------------
def attach_epicurus_score(
    universe: pd.DataFrame,
    *,
    prime_col: str = "prime",
    el_col: str = "el",
    expr_col: str = "expr",
    out_col: str = "epicurus",
    spec: dict | None = None,
) -> pd.DataFrame:
    """Attach the frozen Epicurus v0.1 residual score (higher = better) to ``universe``.

    Reuses the immutable ``configs/frozen/epicurus_v0_1.json`` formula via
    ``event_b.prime_transfer.score_with_frozen`` (per-patient percentile of prime/el/expr; NaN
    features fall back to the 0.5 percentile per the frozen policy). No retraining occurs.
    """
    from event_b.prime_transfer import score_with_frozen

    frame = universe.rename(
        columns={prime_col: "prime", el_col: "el", expr_col: "expr"}
    )[["patient_id", "prime", "el", "expr"]].copy()
    out = universe.copy()
    out[out_col] = score_with_frozen(frame, spec)
    return out


# ---------------------------------------------------------------------------
# Generation-source classification
# ---------------------------------------------------------------------------
def is_pvac_source(value: object) -> bool:
    """True for standard pVAC candidate sources (``pvac`` / ``pvactools*``)."""
    text = str(value).strip().lower()
    return text == "pvac" or text.startswith("pvac")


def _generation_rows(universe: pd.DataFrame, generation: str) -> pd.DataFrame:
    if generation == "pvac":
        return universe[universe["candidate_source"].map(is_pvac_source)]
    return universe  # lossless_union: pVAC + recovered


def _tie_key(peptide: object, hla: object) -> str:
    return hashlib.md5(f"{peptide}|{hla}".encode()).hexdigest()


def _coverage(mutation_ids, positives: set[str]) -> Coverage:
    hit = sorted(set(mutation_ids) & positives)
    return Coverage(len(hit), len(positives), hit)


def _rankable_pool(rows: pd.DataFrame, policy: RouterPolicy) -> pd.DataFrame:
    """Router-eligible, rankable rows (peptide+HLA present, not IMPOSSIBLE)."""
    routed = route_candidates(rows, policy)
    return routed[routed["router_eligible"].astype(bool) & routed["rankable"].astype(bool)]


def _plain_topk(rows: pd.DataFrame, score_col: str, k: int, policy: RouterPolicy) -> pd.DataFrame:
    pool = _rankable_pool(rows, policy).copy()
    pool = pool[pd.to_numeric(pool[score_col], errors="coerce").notna()]
    if pool.empty:
        return pool
    pool["_tie_key"] = [
        _tie_key(p, h)
        for p, h in zip(pool.get("mutant_peptide", ""), pool.get("hla_allele", ""))
    ]
    ordered = pool.sort_values(
        [score_col, "_tie_key"], ascending=[False, True], kind="mergesort"
    )
    return ordered.head(k)


def _route_aware_topk(rows: pd.DataFrame, score_col: str, policy: RouterPolicy) -> pd.DataFrame:
    routed = route_candidates(rows, policy)
    selected = select_route_aware_topk(routed, score_column=score_col, policy=policy)
    return selected[selected["route_selected"].astype(bool)]


def run_arm(
    universe: pd.DataFrame,
    positives: set[str],
    arm: ArmSpec,
    *,
    k: int = DEFAULT_K,
    policy: RouterPolicy = DEFAULT_ROUTER_POLICY,
) -> ArmResult:
    """Compute one arm's per-patient stage metrics on a single-patient candidate universe.

    ``positives`` (measured recognized mutation ids) is consulted ONLY to score coverage, after
    the ranking is fixed. Returns generation recall, rankable recall, and top-k mutation coverage.
    """
    positives = set(positives)
    gen_rows = _generation_rows(universe, arm.generation)
    generation_recall = _coverage(gen_rows["mutation_id"], positives)

    rankable = _rankable_pool(gen_rows, policy)
    rankable_recall = _coverage(rankable["mutation_id"], positives)

    if arm.selection == "route_aware":
        selected = _route_aware_topk(gen_rows, arm.scorer, policy)
    else:
        selected = _plain_topk(gen_rows, arm.scorer, k, policy)
    top_k = _coverage(selected["mutation_id"], positives)

    return ArmResult(
        arm_id=arm.arm_id,
        evaluable=True,
        n_positives=len(positives),
        generation_recall=generation_recall,
        rankable_recall=rankable_recall,
        top_k=top_k,
        n_selected=int(len(selected)),
    )


# ---------------------------------------------------------------------------
# Auto-detected availability from a universe frame
# ---------------------------------------------------------------------------
def detect_available(universe: pd.DataFrame, positives: set[str]) -> set[str]:
    available: set[str] = set()
    if positives:
        available.add(REQ_LABELS)
    source = universe.get("candidate_source")
    if source is not None and source.map(is_pvac_source).any():
        available.add(REQ_PVAC)
    if source is not None and (~source.map(is_pvac_source)).any():
        available.add(REQ_LOSSLESS)
    if "genuine_prime" in universe and pd.to_numeric(universe["genuine_prime"], errors="coerce").notna().any():
        available.add(REQ_PRIME)
    if "epicurus" in universe and pd.to_numeric(universe["epicurus"], errors="coerce").notna().any():
        available.add(REQ_EPICURUS)
    if {"mutant_peptide", "hla_allele"} <= set(universe.columns):
        available.add(REQ_ROUTER)
    return available


def run_patient(
    universe: pd.DataFrame,
    positives: set[str],
    *,
    available: set[str] | None = None,
    arms: list[ArmSpec] = FOUR_ARMS,
    k: int = DEFAULT_K,
    policy: RouterPolicy = DEFAULT_ROUTER_POLICY,
) -> dict:
    """Run all four arms for one patient, honoring per-arm eligibility.

    ``available`` may be supplied explicitly (e.g. from a cohort audit) or is auto-detected from
    the universe + labels. NOT_EVALUABLE arms return a status-only :class:`ArmResult`.
    """
    positives = set(positives)
    if available is None:
        available = detect_available(universe, positives)
    elig = evaluate_eligibility(available, arms)

    results: dict[str, ArmResult] = {}
    for arm in arms:
        if elig[arm.arm_id].evaluable:
            results[arm.arm_id] = run_arm(universe, positives, arm, k=k, policy=policy)
        else:
            results[arm.arm_id] = ArmResult(
                arm_id=arm.arm_id, evaluable=False, missing=elig[arm.arm_id].missing
            )
    return {"available": sorted(available), "arms": results}


def stage_attribution(results: dict[str, ArmResult]) -> dict:
    """Additive decomposition of the full-stack top-k gain, when all four arms are evaluable."""
    order = ["pvac_prime", "lossless_prime", "lossless_epicurus", "full_epicurus"]
    if not all(a in results and results[a].evaluable for a in order):
        return {"evaluable": False}
    h = {a: results[a].hits_at_k for a in order}
    generation = h["lossless_prime"] - h["pvac_prime"]
    scorer = h["lossless_epicurus"] - h["lossless_prime"]
    selection = h["full_epicurus"] - h["lossless_epicurus"]
    total = h["full_epicurus"] - h["pvac_prime"]
    return {
        "evaluable": True,
        "generation": generation,
        "scorer": scorer,
        "selection": selection,
        "total": total,
        "hits": h,
    }
