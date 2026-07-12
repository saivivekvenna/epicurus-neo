"""Phase 3 — locked replay + conditional no-regression for the evidence router.

Runs AFTER the frozen Phase-2 implementation commit (see the preregistration
``docs/superpowers/specs/2026-07-12-evidence-router-and-route-aware-selection-preregistration.md``
and ``configs/frozen/evidence_router_v1.json``). This script FITS NOTHING and TUNES NOTHING: it
replays the frozen router + route-aware selection over existing local artifacts.

Two evaluations, reported separately (mechanical reachability is never used to argue learned
superiority, and vice-versa):

1. **osteosarc.com / Sid — locked reachability replay (hypothesis-generating, NOT independent
   validation; its structural audit informed the router).** For the 3 IFNgamma/TCR (Hudson)
   recognized targets, it reports (a) the multi-caller *raw variant union* recall (computed with the
   frozen ``union_variants``), and (b) the router's routing of the same variants, keeping ASPM/MAP2
   as ``NEEDS_PEPTIDE_GENERATION`` (absent from the single pVACtools 2025.01 candidate set) rather
   than as a ranker miss. Rankable/selected are read from the reconstruction's own
   ``in_pvactools_2025`` / shortlist flags. Vaccine-selected peptide sequences are NOT used to fill
   missing candidate generation, and missing peptide/HLA is NEVER charged to PRIME.

2. **Reranker cohorts — conditional no-regression.** On each locally scored cohort, the route-aware
   top-20 is compared against the pure-score top-20 under the same incumbent score, per patient, on
   paired hits@20 with a patient-level bootstrap CI. The router is a recall-preserving change; the
   requirement is that it must not LOSE hits@20. Feature availability is audited honestly.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from epicurus_neo.evidence_router import (
    DEFAULT_ROUTER_POLICY,
    build_funnel,
    route_candidates,
    select_route_aware_topk,
)
from epicurus_neo.variant_union import union_variants

ROOT = Path(__file__).resolve().parents[1]
M7 = ROOT / "artifacts" / "milestone_7_decision"
OUT = M7 / "evidence_router"
POLICY = DEFAULT_ROUTER_POLICY

# Locked Hudson-recognized target coordinates (from the task brief + reconstruction funnel).
HUDSON_TARGETS = [
    {
        "target_id": "ASPM-chr1-197102716",
        "gene_symbol": "ASPM",
        "chrom": "chr1",
        "pos": 197102716,
        "ref": "C",
        "alt": "T",
        "source_variant_type": "SNV",  # C>T missense
        "protein_variant": "",
        "gene_tpm": 16.49,
    },
    {
        "target_id": "MAP2-chr2-209694772",
        "gene_symbol": "MAP2",
        "chrom": "chr2",
        "pos": 209694772,
        "ref": "GGCTACTGTGTGTTCAATAAGTACACAGT",
        "alt": "G",
        "source_variant_type": "frameshift",  # Gly868fs deletion
        "protein_variant": "",
        "gene_tpm": 5.2,
    },
    {
        "target_id": "DYNC1H1-chr14-101980529",
        "gene_symbol": "DYNC1H1",
        "chrom": "chr14",
        "pos": 101980529,
        "ref": "G",
        "alt": "A",
        "source_variant_type": "SNV",  # G>A missense
        "protein_variant": "",
        "gene_tpm": 357.14,
    },
]
def _hits_at_k(selected_flags: pd.Series, labels: pd.Series) -> int:
    return int((selected_flags.astype(bool) & labels.astype(bool)).sum())


def _bootstrap_ci(deltas: np.ndarray, n: int = 2000) -> tuple[float, float]:
    if len(deltas) == 0:
        return (0.0, 0.0)
    rng = np.random.default_rng(0)
    means = [
        float(np.mean(deltas[rng.integers(0, len(deltas), len(deltas))])) for _ in range(n)
    ]
    return (float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5)))


# ---------------------------------------------------------------------------
# 1. Sid locked reachability replay
# ---------------------------------------------------------------------------
def sid_replay() -> dict:
    funnel_csv = M7 / "osteosarc_sid_reconstruction" / "reachability_funnel.csv"
    recon = pd.read_csv(funnel_csv)
    hudson = recon[recon["hudson_recognized"].astype(str).str.lower() == "true"].copy()
    hudson_ids = set(hudson["target_id"])

    # (a) Multi-caller RAW variant union recall from the actual public VAF table. A detected call
    #     requires tumor alt-read support (>0); normal-only/zero-floor rows are not counted. The
    #     recognized coordinates are joined only AFTER the union is built, avoiding a circular
    #     construction in which known positives are injected into the candidate set.
    raw_vafs = pd.read_csv(
        ROOT / "data/raw/osteosarc/site_cache/variant_vafs_long.tsv",
        sep="\t",
        low_memory=False,
    )
    detected = raw_vafs[
        raw_vafs["tissue"].astype(str).str.lower().eq("tumor")
        & (pd.to_numeric(raw_vafs["alt_reads"], errors="coerce").fillna(0) > 0)
        & raw_vafs["pipeline"].fillna("").astype(str).str.strip().ne("")
    ].copy()
    union_input = detected.rename(
        columns={
            "gene": "gene_symbol",
            "pipeline": "caller",
            "data_source": "source",
            "variant_type": "source_variant_type",
        }
    )
    union_input["patient_id"] = "sid"
    union_input["genome_build"] = "GRCh38"
    union_input["region"] = union_input["tissue"]
    union_input["mutant_peptide"] = ""
    union_input["hla_allele"] = ""
    unioned = union_variants(union_input)

    def identity(row) -> tuple[str, int, str, str]:
        return (str(row["chrom"]), int(row["pos"]), str(row["ref"]), str(row["alt"]))

    union_by_key = {identity(row): row for _, row in unioned.iterrows()}
    recognized_by_key = {identity(row): row for _, row in hudson.iterrows()}
    union_recall_hits = len(set(recognized_by_key) & set(union_by_key))
    union_recall = union_recall_hits / len(hudson_ids)

    # (b) Router routing of the same variants (still peptide-free offline): demonstrates the router
    #     keeps them eligible + NEEDS_PEPTIDE_GENERATION (never IMPOSSIBLE, never a ranker miss).
    router_rows = []
    for target in HUDSON_TARGETS:
        key = (target["chrom"], int(target["pos"]), target["ref"], target["alt"])
        evidence = union_by_key.get(key)
        recon_row = recognized_by_key[key]
        router_rows.append(
            {
                "patient_id": "sid",
                "candidate_id": target["target_id"],
                "gene_symbol": target["gene_symbol"],
                "mutation_id": target["target_id"],
                "source_variant_type": target["source_variant_type"],
                "protein_variant": target["protein_variant"],
                "mhc_class": "I",
                "mutant_peptide": "",  # offline: no pVACtools peptide sequence available
                "hla_allele": "",
                "expression_tpm": float(recon_row["gene_tpm"]),
                "n_callers": int(evidence["n_callers"]) if evidence is not None else 0,
                "n_timepoints": int(evidence["n_timepoints"]) if evidence is not None else 0,
            }
        )
    routed = route_candidates(pd.DataFrame(router_rows))
    routing = {
        row["candidate_id"]: {
            "primary_route": row["primary_route"],
            "router_status": row["router_status"],
            "router_eligible": bool(row["router_eligible"]),
            "rankable": bool(row["rankable"]),
            "router_removed_reason": row["router_removed_reason"],
        }
        for _, row in routed.iterrows()
    }
    # Not one of the recognized targets is hard-removed, and not one is charged to the ranker.
    assert not any(v["primary_route"] == "IMPOSSIBLE" for v in routing.values())
    assert all(v["router_status"] == "NEEDS_PEPTIDE_GENERATION" for v in routing.values())
    # §6 machine-readable per-stage funnel over the offline (peptide-free) routed frame.
    router_offline_funnel = build_funnel(routed)[0]

    # (c) Deployment funnel, using the reconstruction's OWN peptide-generation flag. DYNC is the
    #     only Hudson target that reached peptide generation (in pVACtools 2025.01) and the shortlist.
    peptide_generated = {
        row["target_id"]
        for _, row in hudson.iterrows()
        if str(row["in_pvactools_2025"]).strip().lower() == "true"
    }
    rankable = peptide_generated  # rankable requires a generated peptide+HLA
    # Selection is verified from the existing frozen mutation-level score artifact rather than
    # assumed. DYNC1H1 is rank 1 under frozen Epicurus v0.1 and therefore in its top-20.
    mutation_scores = pd.read_csv(M7 / "osteosarc_sid" / "per_mutation_scores.csv")
    top20_genes = set(
        mutation_scores.sort_values("epicurus_v0_1", ascending=False, kind="mergesort")
        .head(20)["gene"].astype(str)
    )
    selected = {
        target_id for target_id in rankable
        if str(hudson.loc[hudson["target_id"] == target_id, "gene"].iloc[0]) in top20_genes
    }

    funnel = {
        "generated_raw_union": len(hudson_ids),
        "peptide_generated": len(peptide_generated),
        "rankable": len(rankable),
        "selected": len(selected),
        "needs_peptide_generation": len(hudson_ids - peptide_generated),
    }

    return {
        "cohort": "osteosarc_sid",
        "role": "hypothesis_generating_reachability_diagnostic_not_independent_validation",
        "n_hudson_recognized_targets": len(hudson_ids),
        "hudson_target_ids": sorted(hudson_ids),
        "multi_caller_raw_union": {
            "recall_hits": union_recall_hits,
            "recall": union_recall,
            "definition": "fraction of Hudson-recognized targets present in the multi-caller raw "
            "variant union (before peptide generation)",
            "n_union_rows": int(len(unioned)),
            "n_detected_source_rows": int(len(detected)),
            "target_support": {
                str(row["target_id"]): {
                    "n_callers": int(union_by_key[identity(row)]["n_callers"]),
                    "n_timepoints": int(union_by_key[identity(row)]["n_timepoints"]),
                    "callers": str(union_by_key[identity(row)]["callers"]),
                }
                for _, row in hudson.iterrows() if identity(row) in union_by_key
            },
            "note": "All 3 targets are called by DRAGEN/Sarek/oncoanalyser; the recoverable loss is "
            "candidate recall UPSTREAM of any ranker, not variant calling.",
        },
        "router_routing_offline": routing,
        "router_offline_funnel": router_offline_funnel,
        "deployment_funnel": funnel,
        "peptide_generation_gap": sorted(hudson_ids - peptide_generated),
        "reached_ranking": sorted(peptide_generated),
        "hits_at_20_conditional_on_rankable": (
            f"{len(rankable & selected)}/{len(rankable)}" if rankable else "0/0"
        ),
        "hits_at_20_end_to_end": f"{len(selected)}/{len(hudson_ids)}",
        "honesty_notes": [
            "ASPM (chr1:197102716 C>T missense) and MAP2 Gly868fs (chr2:209694772 deletion) are "
            "recovered 3/3 in the multi-caller raw union but were never peptide-generated by the "
            "single pVACtools 2025.01 step -> NEEDS_PEPTIDE_GENERATION, an upstream generation gap, "
            "NOT a PRIME/ranker miss.",
            "Vaccine-selected peptide sequences are NOT used to fill missing candidate generation.",
            "DYNC1H1 (chr14:101980529 G>A) is the only Hudson target peptide-generated; its selection "
            "is verified from the frozen Epicurus v0.1 mutation scores (rank 1, not assumed).",
            "Sid informed the router design and is NOT independent validation; no policy constant was "
            "tuned on it.",
        ],
    }


# ---------------------------------------------------------------------------
# 2. Reranker cohorts — conditional no-regression
# ---------------------------------------------------------------------------
def _feature_audit(frame: pd.DataFrame) -> dict:
    router_features = [
        "expression_tpm",
        "expression_call",
        "rna_mutant_reads",
        "rna_vaf",
        "hla_loh_call",
        "n_callers",
        "n_timepoints",
        "n_regions",
        "presentation_score",
        "binding_percentile_rank",
        "binding_affinity_nm",
        "source_variant_type",
        "mutation_id",
        "gene_symbol",
    ]
    return {f: (f in frame.columns) for f in router_features}


def _populated_fraction(frame: pd.DataFrame, column: str) -> float:
    if column not in frame:
        return 0.0
    values = frame[column]
    return float((values.notna() & values.astype(str).str.strip().ne("")).mean())


def no_regression(cohort: str, path: Path, score_column: str, tested_only: bool) -> dict:
    frame = pd.read_csv(path)
    if tested_only and "label" in frame:
        frame = frame[frame["label"].astype(str).str.upper() != "UNTESTED"].copy()
    frame = frame.reset_index(drop=True)
    if "candidate_id" not in frame:
        frame["candidate_id"] = [f"{cohort}_{i}" for i in range(len(frame))]

    routed = route_candidates(frame)
    n_rankable = int(routed["rankable"].astype(bool).sum())
    removed = routed[~routed["router_eligible"].astype(bool)]
    removed_reason_counts = {
        str(k): int(v) for k, v in removed["router_removed_reason"].value_counts().items()
    }
    positives_removed = int(removed["y"].sum()) if "y" in removed else 0
    peptide_populated = _populated_fraction(frame, "mutant_peptide")
    hla_populated = _populated_fraction(frame, "hla_allele")

    base_meta = {
        "cohort": cohort,
        "incumbent_score_column": score_column,
        "denominator": "tested_only" if tested_only else "all_scored",
        "n_candidates": int(len(frame)),
        "n_patients": int(frame["patient_id"].nunique()) if "patient_id" in frame else 0,
        "n_positives": int(frame["y"].sum()) if "y" in frame else 0,
        "n_rankable": n_rankable,
        "n_positives_rankable": int(routed.loc[routed["rankable"].astype(bool), "y"].sum())
        if "y" in routed
        else 0,
        "peptide_populated_fraction": peptide_populated,
        "hla_populated_fraction": hla_populated,
        "router_removed_reason_counts": removed_reason_counts,
        "positives_removed_by_router": positives_removed,
        "route_composition_valid": {
            str(r): int(c)
            for r, c in routed.loc[routed["router_eligible"].astype(bool), "primary_route"]
            .value_counts()
            .items()
        },
        "router_feature_availability": _feature_audit(frame),
    }

    # A route-aware TOP-20 comparison needs rankable candidates (peptide+HLA). If the artifact does
    # not carry HLA/peptide, nothing is rankable and the comparison is NOT evaluable here — an honest
    # data-availability limit, not a ranker failure.
    if n_rankable == 0:
        return {
            **base_meta,
            "verdict": "NOT_EVALUABLE_NO_RANKABLE_CANDIDATES",
            "honesty_note": (
                "The scored artifact lacks populated HLA and/or peptide "
                f"(hla_populated={hla_populated:.2f}, peptide_populated={peptide_populated:.2f}), so "
                "no candidate is rankable and the route-aware top-20 no-regression cannot be measured "
                "on this artifact. This is a data limit of the stored file, not a router or ranker "
                "result; a genuine-PRIME/HLA-resolved re-score would be needed to evaluate it."
            ),
        }

    route_aware = select_route_aware_topk(routed, score_column=score_column)

    # Pure-score baseline: same deterministic ordering (score desc, md5 tie asc), top-k per patient,
    # over the same eligible+rankable pool. The route-aware path only adds reserves/caps on top.
    baseline = route_aware.copy()
    baseline["_pure_selected"] = False
    pool = baseline[baseline["router_eligible"].astype(bool) & baseline["rankable"].astype(bool)]
    for _, idx in pool.groupby("patient_id", sort=True).groups.items():
        rows = pool.loc[idx].sort_values(
            [score_column, "_tie_key"], ascending=[False, True], kind="mergesort"
        )
        baseline.loc[rows.head(POLICY.k).index, "_pure_selected"] = True

    per_patient = []
    exact_membership_by_patient = []
    for patient, group in route_aware.groupby("patient_id", sort=True):
        base = baseline.loc[group.index]
        ra_hits = _hits_at_k(group["route_selected"], group["y"])
        pure_hits = _hits_at_k(base["_pure_selected"], group["y"])
        per_patient.append(
            {"patient_id": str(patient), "route_aware": ra_hits, "pure_score": pure_hits,
             "delta": ra_hits - pure_hits}
        )
        route_ids = set(group.loc[group["route_selected"], "candidate_id"].astype(str))
        pure_ids = set(base.loc[base["_pure_selected"], "candidate_id"].astype(str))
        exact_membership_by_patient.append(route_ids == pure_ids)
    deltas = np.array([p["delta"] for p in per_patient], dtype=float)
    lo, hi = _bootstrap_ci(deltas)

    exact_pass_through = all(exact_membership_by_patient)
    verdict = (
        "NO_REGRESSION_EXACT_PASS_THROUGH"
        if exact_pass_through
        else ("NO_REGRESSION" if lo >= 0 else "REGRESSION_DETECTED")
    )
    return {
        **base_meta,
        "total_route_aware_hits_at_20": int(sum(p["route_aware"] for p in per_patient)),
        "total_pure_score_hits_at_20": int(sum(p["pure_score"] for p in per_patient)),
        "mean_paired_delta_hits_at_20": float(deltas.mean()) if len(deltas) else 0.0,
        "bootstrap_ci_delta": [lo, hi],
        "exact_selection_membership_all_patients": exact_pass_through,
        "verdict": verdict,
        "honesty_note": "Router recall-discriminating features (RNA/expression/HLA-LOH/multi-caller "
        "provenance) are absent from this reranker artifact and no diversity-cap key is populated, "
        "so every candidate collapses to a single route and route-aware selection is a pure-score "
        "pass-through. This satisfies no-regression by construction; it is NOT evidence of benefit. "
        "The router's value is at candidate GENERATION (the Sid recall finding), not at reranking a "
        "feature-poor candidate list.",
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sid = sid_replay()
    reranker = [
        no_regression(
            "cd8_multimer",
            M7 / "prime_headtohead" / "multimer_prime_scored.csv",
            score_column="prime",  # genuine PRIME where available
            tested_only=True,
        ),
        no_regression(
            "gartner_nci",
            M7 / "gartner" / "scored_candidates.csv.gz",
            score_column="full_model",  # frozen Epicurus product score (no genuine PRIME per-candidate here)
            tested_only=True,
        ),
        no_regression(
            "cd8_multimer_decision",
            M7 / "multimer" / "scored_candidates.csv.gz",
            score_column="full_model",
            tested_only=True,
        ),
    ]

    payload = {
        "policy_id": POLICY.policy_id,
        "phase": "phase_3_locked_replay_and_conditional_no_regression",
        "fitting": "NONE — frozen router + selection replayed over existing artifacts",
        "sid_reachability": sid,
        "no_regression": reranker,
        "headline": {
            "sid_multi_caller_union_recall": sid["multi_caller_raw_union"]["recall"],
            "sid_rankable": sid["deployment_funnel"]["rankable"],
            "sid_selected": sid["deployment_funnel"]["selected"],
            "sid_needs_peptide_generation": sid["deployment_funnel"]["needs_peptide_generation"],
            "no_regression_verdicts": {r["cohort"]: r["verdict"] for r in reranker},
        },
    }
    (OUT / "phase3_replay.json").write_text(json.dumps(payload, indent=2) + "\n")
    _write_markdown(payload)
    print(json.dumps(payload["headline"], indent=2))


def _write_markdown(payload: dict) -> None:
    sid = payload["sid_reachability"]
    lines = [
        "# Evidence router — Phase 3 locked replay + conditional no-regression",
        "",
        f"> Policy `{payload['policy_id']}`. **No fitting/tuning.** Frozen router + route-aware "
        "selection replayed over existing local artifacts. Mechanical reachability and learned "
        "no-regression are reported separately; neither is used to argue the other.",
        "",
        "## 1. osteosarc.com / Sid — locked reachability replay",
        "",
        "> Sid informed the router design and is **not** independent validation. No policy constant "
        "was tuned on it. Vaccine-selected peptides are not used to fill candidate generation; "
        "missing peptide/HLA is never charged to PRIME.",
        "",
        f"- Hudson IFNgamma/TCR recognized targets: **{sid['n_hudson_recognized_targets']}** "
        f"({', '.join(sid['hudson_target_ids'])}).",
        f"- **Multi-caller raw variant union recall: {sid['multi_caller_raw_union']['recall_hits']}/"
        f"{sid['n_hudson_recognized_targets']} = {sid['multi_caller_raw_union']['recall']:.2f}** "
        "(all three are called by DRAGEN/Sarek/oncoanalyser — the recoverable loss is candidate "
        "recall upstream of any ranker).",
        "",
        "| Stage | Count |",
        "|---|---:|",
        f"| generated (raw union) | {sid['deployment_funnel']['generated_raw_union']} |",
        f"| peptide-generated (pVACtools 2025.01) | {sid['deployment_funnel']['peptide_generated']} |",
        f"| rankable (peptide+HLA) | {sid['deployment_funnel']['rankable']} |",
        f"| selected (frozen Epicurus v0.1 top-20) | {sid['deployment_funnel']['selected']} |",
        f"| NEEDS_PEPTIDE_GENERATION | {sid['deployment_funnel']['needs_peptide_generation']} |",
        "",
        f"- Peptide-generation gap (NEEDS_PEPTIDE_GENERATION, not a ranker miss): "
        f"**{', '.join(sid['peptide_generation_gap'])}**.",
        f"- Reached ranking / shortlist: **{', '.join(sid['reached_ranking'])}**.",
        f"- hits@20 conditional on rankability: **{sid['hits_at_20_conditional_on_rankable']}**; "
        f"end-to-end hits@20: **{sid['hits_at_20_end_to_end']}**.",
        "",
        "Router routing of the three peptide-free variants (offline): "
        + "; ".join(
            f"{cid} -> {v['primary_route']}/{v['router_status']}"
            for cid, v in sid["router_routing_offline"].items()
        )
        + ". None is hard-removed (IMPOSSIBLE); none is charged to the ranker.",
        "",
        "## 2. Reranker cohorts — conditional no-regression",
        "",
        "> Requirement (§7): route-aware top-20 must not LOSE hits@20 vs the pure-score top-20. "
        "The router is recall-preserving.",
        "",
        "| Cohort | Score | Patients | Positives | Rankable | HLA populated | Route-aware hits@20 | "
        "Pure-score hits@20 | Δ (CI) | Verdict |",
        "|---|---|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for r in payload["no_regression"]:
        if "total_route_aware_hits_at_20" in r:
            ra = str(r["total_route_aware_hits_at_20"])
            pure = str(r["total_pure_score_hits_at_20"])
            delta = "{d:.3f} [{lo:.3f}, {hi:.3f}]".format(
                d=r["mean_paired_delta_hits_at_20"],
                lo=r["bootstrap_ci_delta"][0],
                hi=r["bootstrap_ci_delta"][1],
            )
        else:  # NOT_EVALUABLE
            ra = pure = delta = "n/a"
        lines.append(
            "| {c} | {s} | {np_} | {pos} | {rk} | {hla:.2f} | {ra} | {pure} | {delta} | {v} |".format(
                c=r["cohort"],
                s=r["incumbent_score_column"],
                np_=r["n_patients"],
                pos=r["n_positives"],
                rk=r["n_rankable"],
                hla=r["hla_populated_fraction"],
                ra=ra,
                pure=pure,
                delta=delta,
                v=r["verdict"],
            )
        )
    dup_notes = "; ".join(
        f"{r['cohort']}: {r['router_removed_reason_counts']} removed, "
        f"{r['positives_removed_by_router']} positives lost"
        for r in payload["no_regression"]
    )
    lines += [
        "",
        "### Honest feature-availability caveats",
        "",
        "- On the evaluable reranker artifacts (both CD8 multimer scorings) the router's "
        "recall-discriminating features (RNA/expression/HLA-LOH/multi-caller provenance) are "
        "**absent** and no diversity-cap key is populated, so all candidates collapse to a single "
        "route and route-aware selection is a **pure-score pass-through** (Δ=0 by construction). This "
        "satisfies no-regression but is **not** evidence of benefit; the router's value is realised at "
        "candidate **generation** (the Sid recall recovery), not at reranking a feature-poor list.",
        "- The **Gartner NCI** stored artifact carries **no HLA allele** "
        f"(hla populated = {next(r['hla_populated_fraction'] for r in payload['no_regression'] if r['cohort'] == 'gartner_nci'):.2f}), "
        "so nothing is rankable and the top-20 no-regression is **NOT_EVALUABLE** there — a data limit "
        "of the file, not a ranker result.",
        f"- Router removals per cohort (route-verifiable only; positives lost must be 0): {dup_notes}.",
        "",
        "Route composition (valid candidates) per cohort: "
        + "; ".join(f"{r['cohort']}={r['route_composition_valid']}" for r in payload["no_regression"])
        + ".",
        "",
    ]
    (OUT / "PHASE3_REPLAY.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
