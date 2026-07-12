"""Three-level benchmark cohort eligibility audit.

Per the product decision (2026-07-12, `NORTH_STAR_HISTORY.md`), the benchmark is a hierarchy of three
DISTINCT tasks — never pooled into one headline metric:

  1. **REACHABILITY** — raw variants/WES/RNA/HLA *through peptide generation*, with stage-loss
     attribution (which recognized mutations are lost at generation vs rankable vs top-k).
  2. **CONDITIONAL RANKING** — ranking ONLY among candidates actually generated/rankable, each cohort
     interpreted strictly within its own assay/denominator.
  3. **END-TO-END PATIENT UTILITY** (the PRIMARY north star) — recognized mutations in the final top-20
     from *common raw inputs*, compared against a standard pVAC-style pipeline + genuine PRIME.

Each cohort has a FIXED role and is eligible for a subset of levels. The four-arm generation×scorer
harness (`benchmark.four_arm`) is reusable INFRASTRUCTURE; it only produces an end-to-end headline where
a cohort is Level-3 eligible (currently `osteosarc_sid` only, post-hoc).

Availability reflects documented facts in the ledger (`NORTH_STAR_HISTORY.md` / `memory/`), not a live
file probe — much raw input is dbGaP/DUA-gated. ``REQ_LOSSLESS`` (a raw multi-caller GRCh38 callset
carried to peptide generation) is present ONLY for cohorts that actually carry raw variants to
generation; a pre-generated peptide×HLA list does NOT qualify.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from benchmark.four_arm import (
    FOUR_ARMS,
    REQ_EPICURUS,
    REQ_LABELS,
    REQ_LOSSLESS,
    REQ_PRIME,
    REQ_PVAC,
    REQ_ROUTER,
    evaluate_eligibility,
)

ALL_REQS = (REQ_LABELS, REQ_PVAC, REQ_LOSSLESS, REQ_PRIME, REQ_EPICURUS, REQ_ROUTER)
ARM_IDS = [a.arm_id for a in FOUR_ARMS]

# ---------------------------------------------------------------------------
# The three benchmark levels
# ---------------------------------------------------------------------------
LEVEL_REACHABILITY = 1
LEVEL_CONDITIONAL_RANKING = 2
LEVEL_END_TO_END = 3

LEVEL_NAMES = {
    LEVEL_REACHABILITY: "reachability",
    LEVEL_CONDITIONAL_RANKING: "conditional_ranking",
    LEVEL_END_TO_END: "end_to_end_patient_utility",
}

LEVEL_DESCRIPTIONS = {
    LEVEL_REACHABILITY: "raw variants/WES/RNA/HLA through peptide generation, with stage-loss "
                        "attribution (recognized mutations lost at generation vs rankable vs top-k)",
    LEVEL_CONDITIONAL_RANKING: "ranking only among candidates actually generated/rankable, interpreted "
                               "strictly within this cohort's own assay/denominator (never pooled)",
    LEVEL_END_TO_END: "PRIMARY north star: recognized mutations in the final top-20 from common raw "
                      "inputs, vs a standard pVAC-style pipeline + genuine PRIME",
}

# Requirement keys per level. Reachability needs raw->generation AND a recognized set (labels) to
# attribute stage loss against the incumbent (pVAC) generation. Conditional ranking needs rankable
# candidates + the genuine-PRIME incumbent + labels, but NOT raw generation. End-to-end needs all of
# raw->generation, the pVAC incumbent, genuine PRIME, and labels.
_LEVEL_REQUIREMENTS = {
    LEVEL_REACHABILITY: frozenset({REQ_LOSSLESS, REQ_PVAC, REQ_LABELS}),
    LEVEL_CONDITIONAL_RANKING: frozenset({REQ_PVAC, REQ_PRIME, REQ_LABELS}),
    LEVEL_END_TO_END: frozenset({REQ_LOSSLESS, REQ_PVAC, REQ_PRIME, REQ_LABELS}),
}


def level_requirements(level: int) -> frozenset[str]:
    return _LEVEL_REQUIREMENTS[level]


NO_POOLING = (
    "These cohorts serve DIFFERENT tasks over DIFFERENT denominators (broad-somatic Gartner, "
    "prefiltered IMPROVE, presentation/T-cell multimer, external CheckMate, end-to-end Sid) and are "
    "NEVER pooled into one headline metric. A single cross-cohort reranker number would be a category "
    "error. Each level and each cohort is reported and interpreted on its own."
)


@dataclass(frozen=True)
class Cohort:
    cohort_id: str
    available: set[str]
    role: str
    denominator: str
    labels: str
    leakage_blocked: dict[str, str] = field(default_factory=dict)  # arm_id -> reason
    note: str = ""


COHORTS: list[Cohort] = [
    Cohort(
        "osteosarc_sid",
        available={REQ_LABELS, REQ_PVAC, REQ_LOSSLESS, REQ_PRIME, REQ_EPICURUS, REQ_ROUTER},
        role="End-to-end diagnostic (all 3 levels; post-hoc, n=3, single patient)",
        denominator="per-patient somatic pVAC universe (21 curated muts) + input-only lossless recovery",
        labels="Hudson IFNg peptide-expansion; 3 recognized mutations (only 1 in the pVAC universe)",
        note="ONLY cohort with a raw callset carried to generation AND a measured label -> the only "
             "Level-3 (end-to-end) instance. POST-HOC, n=3, 1 patient: diagnostic, not a powered/blinded "
             "gate. Frozen Epicurus is OUT-OF-SAMPLE here (trained on cd8_multimer).",
    ),
    Cohort(
        "gartner_nci",
        available={REQ_LABELS, REQ_PVAC, REQ_PRIME, REQ_EPICURUS, REQ_ROUTER},
        role="Conditional broad-denominator ranking (Level 2)",
        denominator="Gartner TEST minimal peptide-HLA list (SEMI-CONSUMED holdout); no raw callset->gen",
        labels="TIL/ELISpot screening positives + tested negatives",
        note="Level-2 only. No lossless-generation arm: the stored artifact is a pre-generated peptide "
             "list, not a raw multi-caller callset. Holdout semi-consumed -> not a fresh test. Broadest "
             "somatic denominator of the ranking cohorts; interpret ONLY within itself.",
    ),
    Cohort(
        "improve_srhgroup",
        available={REQ_LABELS, REQ_PVAC, REQ_PRIME, REQ_EPICURUS, REQ_ROUTER},
        role="Conditional prefiltered-subset ranking (Level 2)",
        denominator="pre-screened ~200 candidates/patient — NOT the full somatic universe",
        labels="functional response POSITIVE vs TESTED_NEGATIVE",
        note="Level-2 only; denominator is range-restricted (prefiltered) so top-20 recall is optimistic "
             "and NOT comparable to Gartner's broad denominator. No raw callset -> no lossless arm.",
    ),
    Cohort(
        "checkmate153",
        available={REQ_LABELS, REQ_PVAC, REQ_PRIME, REQ_EPICURUS, REQ_ROUTER},
        role="External conditional ranking (Level 2; small, PRIME-untouched)",
        denominator="HLA-resolved 9-mer candidate list (14 patients); raw WES/RNA behind dbGaP",
        labels="combinatorial tetramer+ vs tested-negative",
        note="Level-2 only, underpowered (14 pts). Acquiring dbGaP WES/RNA would upgrade it toward "
             "Level 1/3. Independent + PRIME-untouched but range-restricted.",
    ),
    Cohort(
        "cd8_multimer",
        available={REQ_LABELS, REQ_PVAC, REQ_PRIME, REQ_EPICURUS, REQ_ROUTER},
        role="Presentation / T-cell compatibility asset (and Epicurus v0.1 training set)",
        denominator="pMHC-multimer peptide-HLA candidate list",
        labels="multimer POSITIVE vs TESTED_NEGATIVE",
        leakage_blocked={
            "lossless_epicurus": "epicurus_v0_1_training_cohort",
            "full_epicurus": "epicurus_v0_1_training_cohort",
        },
        note="Presentation/T-cell-compatibility role. It IS the frozen Epicurus v0.1 training cohort -> "
             "the Epicurus arms are leakage-INVALID here even with inputs present; only the genuine-PRIME "
             "scorer is honest. No raw callset -> no lossless arm. Level-2 for genuine PRIME only.",
    ),
    Cohort(
        "cedar_tcell",
        available={REQ_LABELS, REQ_PVAC},
        role="Training / recognition-prior asset",
        denominator="mutation-derived cancer T-cell recognition rows (not a per-patient ranking denom.)",
        labels="curated recognition POSITIVE / explicit NEGATIVE rows",
        note="Recognition-PRIOR / training asset, not a ranking cohort: no per-patient somatic "
             "denominator, no genuine-PRIME incumbent, no raw callset. Feeds priors/representation, never "
             "a benchmark headline. Eligible for none of the three levels.",
    ),
    Cohort(
        "zhao_dc_2026",
        available={REQ_LABELS, REQ_PVAC, REQ_ROUTER},
        role="Training / recognition-prior asset (genuine PRIME blocked)",
        denominator="352 pts / 2317 SNV peptide-HLA list",
        labels="DC-presentation POSITIVE vs explicit NEGATIVE",
        note="Recognition-prior/training asset. genuine PRIME is BLOCKED (only a MixMHCpred proxy; the "
             "incumbent guard forbids labeling a proxy PRIME) -> no Level-2/3 incumbent, no genuine-PRIME "
             "feature -> no Epicurus arm, no raw callset -> no Level 1. Eligible for none of the levels.",
    ),
    Cohort(
        "rttp_sr24_58221",
        available={REQ_PVAC, REQ_LOSSLESS, REQ_PRIME, REQ_EPICURUS, REQ_ROUTER},
        role="Deployment only (no measured recognition label)",
        denominator="complete North-Star INPUT (WES+RNA+HLA+candidate universe), 1 patient",
        labels="NONE measured (Personalis immunogenicity score is PREDICTED, not an assay)",
        note="Deployment/DEMO asset, not a benchmark: no measured recognition label -> no denominator -> "
             "eligible for none of the three levels. Inputs would otherwise support all three.",
    ),
]


def classify_levels(cohort: Cohort) -> dict[int, dict]:
    """Per-level eligibility for a cohort: eligible + missing requirement keys."""
    out: dict[int, dict] = {}
    for level in (LEVEL_REACHABILITY, LEVEL_CONDITIONAL_RANKING, LEVEL_END_TO_END):
        missing = sorted(level_requirements(level) - cohort.available)
        out[level] = {"eligible": not missing, "missing": missing}
    return out


def audit_cohort(cohort: Cohort) -> dict:
    """Per-cohort audit: role, three-level eligibility, and the four-arm infra evaluability."""
    elig = evaluate_eligibility(cohort.available)
    arms = {}
    for arm in FOUR_ARMS:
        missing = list(elig[arm.arm_id].missing)
        if arm.arm_id in cohort.leakage_blocked:
            missing.append(f"LEAKAGE:{cohort.leakage_blocked[arm.arm_id]}")
        arms[arm.arm_id] = {"evaluable": not missing, "missing": missing}
    levels = {
        LEVEL_NAMES[level]: verdict for level, verdict in classify_levels(cohort).items()
    }
    return {
        "cohort_id": cohort.cohort_id,
        "role": cohort.role,
        "available_inputs": sorted(cohort.available),
        "missing_inputs": sorted(set(ALL_REQS) - cohort.available),
        "denominator": cohort.denominator,
        "labels": cohort.labels,
        "levels": levels,
        "arms": arms,
        "note": cohort.note,
    }


def run_cohort_audit(cohorts: list[Cohort] = COHORTS) -> dict:
    audits = [audit_cohort(c) for c in cohorts]
    n_end_to_end = sum(a["levels"][LEVEL_NAMES[LEVEL_END_TO_END]]["eligible"] for a in audits)
    return {
        "policy": "three-level-benchmark-1.0.0",
        "levels": [
            {"level": level, "name": LEVEL_NAMES[level], "description": LEVEL_DESCRIPTIONS[level],
             "requirements": sorted(level_requirements(level))}
            for level in (LEVEL_REACHABILITY, LEVEL_CONDITIONAL_RANKING, LEVEL_END_TO_END)
        ],
        "no_pooling": NO_POOLING,
        "four_arm_harness_role": "reusable infrastructure; produces an end-to-end headline ONLY for "
                                 "Level-3-eligible cohorts (currently osteosarc_sid, post-hoc).",
        "arms": [{"arm_id": a.arm_id, "generation": a.generation, "scorer": a.scorer,
                  "selection": a.selection, "label": a.label} for a in FOUR_ARMS],
        "requirement_keys": list(ALL_REQS),
        "n_cohorts": len(audits),
        "n_end_to_end_eligible": int(n_end_to_end),
        "cohorts": audits,
        "interpretation": (
            "Exactly one cohort (osteosarc_sid) is END-TO-END (Level-3) eligible, and only post-hoc "
            "with n=3. Conditional-ranking (Level-2) cohorts each stand alone within their own "
            "denominator and are never pooled. Reachability (Level-1) needs a raw callset carried to "
            "generation, which only Sid provides among labelled cohorts. The binding constraint is a "
            "DENSE-denominator, PRIME-untouched, end-to-end patient — acquire/identify one; do not "
            "manufacture a pooled reranker headline."
        ),
    }
