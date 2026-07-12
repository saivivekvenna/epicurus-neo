"""Sid full Epicurus filter-stack benchmark, with labels joined only after selections freeze.

Consumes the label-free, label-blind generated candidate table produced by
``python -m scripts.sid_benchmark_generate --offline`` and evaluates every currently runnable product
layer from the same candidate universe:

* deterministic validity gate;
* frozen dynamic safe-rejection gate;
* genuine PRIME and frozen Epicurus v0.1 rankers;
* frozen evidence router + route-aware selection;
* current product-v1 evidence policy and complete composed stack.

The three Hudson-recognized mutations are not imported until ``frozen_selections.json`` has been written.
This is post-hoc n=1/3 and the upstream generator covers 130/147 eligible variants, so it is a diagnostic,
not a superiority claim.
"""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

import pandas as pd

from epicurus_neo.evidence_router import route_candidates, select_route_aware_topk
from epicurus_neo.gates import apply_deterministic_gate, summarize_gate
from epicurus_neo.product import InferenceConfig, score_product_candidates
from event_b.dynamic_gate import GateConfig, apply_gate
from event_b.sid_full_pipeline import (
    assert_label_blind,
    evaluate_frozen,
    freeze_mutation_topk,
    freeze_portfolio,
    prepare_sid_gate_frame,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "artifacts/milestone_7_decision/sid_benchmark/scored_candidates.csv.gz"
VAF = ROOT / "data/raw/osteosarc/site_cache/variant_vafs_long.tsv"
DYNAMIC_CONFIG = ROOT / "configs/frozen/dynamic_gate_v1.json"
OUT = ROOT / "artifacts/milestone_7_decision/sid_full_pipeline"
K = 20


def _dynamic_config() -> GateConfig:
    return GateConfig.from_json(json.loads(DYNAMIC_CONFIG.read_text())["config"])


def _freeze_route_arm(frame: pd.DataFrame, score: str) -> tuple[dict, pd.DataFrame]:
    routed = route_candidates(frame)
    selected = select_route_aware_topk(routed, score_column=score)
    return freeze_portfolio(selected, "route_selected", "route_rank", k=K), selected


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    base = prepare_sid_gate_frame(SOURCE, VAF)
    assert_label_blind(base)
    cfg = _dynamic_config()

    # Layer 0 — deterministic validity. Sid HLA-LOH is unavailable, so that rule abstains.
    deterministic = apply_deterministic_gate(base)
    det_pass = deterministic[deterministic["deterministic_gate_pass"]].copy()

    # Layers 1–3 — frozen dynamic safe rejection, both alone and after deterministic validity.
    dynamic = apply_gate(base, cfg)
    dyn_pass = dynamic[dynamic["dyn_gate_keep"]].copy()
    det_dynamic = apply_gate(det_pass, cfg)
    det_dyn_pass = det_dynamic[det_dynamic["dyn_gate_keep"]].copy()

    frozen: dict[str, dict] = {}
    frozen["baseline_genuine_prime"] = freeze_mutation_topk(base, "prime", k=K, ascending=True)
    frozen["baseline_frozen_epicurus_v0_1"] = freeze_mutation_topk(
        base, "arm_frozen_epicurus_v0_1", k=K)
    frozen["deterministic_then_prime"] = freeze_mutation_topk(det_pass, "prime", k=K, ascending=True)
    frozen["dynamic_then_prime"] = freeze_mutation_topk(dyn_pass, "prime", k=K, ascending=True)
    frozen["deterministic_dynamic_then_prime"] = freeze_mutation_topk(
        det_dyn_pass, "prime", k=K, ascending=True)
    frozen["deterministic_dynamic_then_frozen_epicurus"] = freeze_mutation_topk(
        det_dyn_pass, "arm_frozen_epicurus_v0_1", k=K)

    # Unit-corrected sensitivity: the frozen dynamic policy was calibrated on candidate rows, but the
    # clinical decision is mutation-level. Aggregate each veto axis to its best available route and apply
    # the unchanged frozen threshold once per mutation. This is diagnostic, not a new validated policy.
    mutation_gate_input = (
        base.sort_values("prime", kind="mergesort").drop_duplicates("mutation_id", keep="first").copy()
    )
    best_axes = base.groupby("mutation_id").agg(el=("el", "min"), prime=("prime", "min"), expr=("expr", "first"))
    mutation_gate_input = mutation_gate_input.drop(columns=["el", "prime", "expr"]).join(
        best_axes, on="mutation_id")
    mutation_dynamic = apply_gate(mutation_gate_input, cfg)
    mutation_dyn_pass = mutation_dynamic[mutation_dynamic["dyn_gate_keep"]].copy()
    frozen["mutation_level_dynamic_sensitivity_then_prime"] = freeze_mutation_topk(
        mutation_dyn_pass, "prime", k=K, ascending=True)

    frozen["deterministic_dynamic_router_prime"], route_prime = _freeze_route_arm(
        det_dyn_pass, "genuine_prime_score")
    frozen["deterministic_dynamic_router_frozen_epicurus"], route_epicurus = _freeze_route_arm(
        det_dyn_pass, "arm_frozen_epicurus_v0_1")

    # Current product-v1 score maps percentile ranks monotonically to [0,1], then uses measured TPM and
    # longitudinal DNA VAF. No missing RNA-VAF/read value is invented.
    product = score_product_candidates(base, InferenceConfig())
    frozen["product_v1_as_shipped"] = freeze_portfolio(product, "selected", "rank", k=K)

    # Complete composed stack: product validity/evidence eligibility -> frozen dynamic gate -> frozen
    # evidence router/route-aware selection, using the product's lower-confidence evidence score.
    product_eligible = product[product["eligible"]].copy()
    product_dynamic = apply_gate(product_eligible, cfg)
    product_dyn_pass = product_dynamic[product_dynamic["dyn_gate_keep"]].copy()
    frozen["full_stack_mutation_level_fair"] = freeze_mutation_topk(
        product_dyn_pass, "epicurus_lower_evidence_score", k=K)
    frozen["full_current_epicurus_stack"], full_selected = _freeze_route_arm(
        product_dyn_pass, "epicurus_lower_evidence_score")

    # HARD LABEL BARRIER: serialize every selection before importing evaluation labels.
    frozen_payload = {
        "status": "FROZEN_BEFORE_LABEL_JOIN",
        "source": str(SOURCE.relative_to(ROOT)),
        "source_has_evaluation_columns": False,
        "k": K,
        "arms": frozen,
    }
    frozen_path = OUT / "frozen_selections.json"
    frozen_path.write_text(json.dumps(frozen_payload, indent=2, sort_keys=True) + "\n")

    # EVALUATION ONLY — import exact Hudson mutation labels after the on-disk freeze above.
    from event_b.sid_benchmark import hudson_positive_variant_ids

    positives = hudson_positive_variant_ids()
    evaluation = {name: evaluate_frozen(selection, positives) for name, selection in frozen.items()}

    det_summary = summarize_gate(deterministic)
    dynamic_removed = ~dynamic["dyn_gate_keep"].astype(bool)
    det_dynamic_removed = ~det_dynamic["dyn_gate_keep"].astype(bool)
    product_excluded = ~product["eligible"].astype(bool)
    report = {
        "experiment": "sid_full_filter_stack_label_blind",
        "status": "POST_HOC_N1_3_DIAGNOSTIC_GENERATOR_COVERAGE_88.4_PERCENT",
        "label_barrier": {
            "selections_frozen_before_label_import": True,
            "frozen_path": str(frozen_path.relative_to(ROOT)),
            "positive_ids_joined_after_freeze": sorted(positives),
        },
        "input": {
            "candidate_rows": int(len(base)),
            "generated_mutations": int(base["mutation_id"].nunique()),
            "eligible_upstream_mutations": 147,
            "generator_coverage": round(base["mutation_id"].nunique() / 147, 4),
        },
        "feature_availability": {
            "peptide_sequence": "COMPLETE_ON_GENERATED",
            "patient_hla": "COMPLETE_MATCHED_SID_A01_B08_B27_C01_C07",
            "tumor_gene_tpm": "COMPLETE_ON_GENERATED_MATCHED_T2_RSEM",
            "tumor_dna_vaf": "COMPLETE_ON_GENERATED_LONGITUDINAL_WES_WGS",
            "multi_caller_timepoint_support": "COMPLETE_ON_GENERATED",
            "mixmhcpred_and_genuine_prime": "COMPLETE_ON_GENERATED",
            "hla_loh": "NOT_ASSESSED_FOR_SID",
            "rna_vaf_and_mutant_reads": "COMPLETE_ON_GENERATED_MATCHED_T2_BULK_RNA",
            "longitudinal_rna_support": "COMPLETE_ON_GENERATED_PREDECISION_T0_T1_T2_BULK_AND_SINGLE_CELL",
        },
        "filter_effects_before_label_join": {
            "deterministic": {
                **asdict(det_summary),
                "evidence_status": "PARTIAL: sequence/duplicate rules runnable; Sid HLA-LOH unavailable",
                "mutations_removed": int(base["mutation_id"].nunique() - det_pass["mutation_id"].nunique()),
            },
            "dynamic_v1": {
                "candidate_rows_removed": int(dynamic_removed.sum()),
                "candidate_row_removal_fraction": round(float(dynamic_removed.mean()), 4),
                "mutations_removed": int(base["mutation_id"].nunique() - dyn_pass["mutation_id"].nunique()),
                "reason_counts": dynamic.loc[dynamic_removed, "dyn_gate_reason"].value_counts().to_dict(),
                "config": cfg.to_json(),
            },
            "deterministic_plus_dynamic": {
                "candidate_rows_removed_after_deterministic": int(det_dynamic_removed.sum()),
                "mutations_removed": int(det_pass["mutation_id"].nunique() - det_dyn_pass["mutation_id"].nunique()),
            },
            "mutation_level_dynamic_sensitivity": {
                "status": "POST_HOC_UNIT_CORRECTION_WITH_FROZEN_THRESHOLD_NOT_A_VALIDATED_POLICY",
                "mutations_input": int(len(mutation_dynamic)),
                "mutations_removed": int((~mutation_dynamic["dyn_gate_keep"].astype(bool)).sum()),
            },
            "product_v1_expression_evidence_policy": {
                "candidate_rows_excluded": int(product_excluded.sum()),
                "mutations_remaining": int(product.loc[~product_excluded, "mutation_id"].nunique()),
                "exclusion_reasons": product.loc[product_excluded, "exclusion_reason"].value_counts().to_dict(),
            },
            "route_distributions": {
                "prime": route_prime["primary_route"].value_counts().to_dict(),
                "frozen_epicurus": route_epicurus["primary_route"].value_counts().to_dict(),
                "full_stack": full_selected["primary_route"].value_counts().to_dict(),
            },
        },
        "evaluation_only": evaluation,
        "verdict": {
            "full_stack_hits_at_20": evaluation["full_current_epicurus_stack"]["hits_at_20"],
            "genuine_prime_hits_at_20": evaluation["baseline_genuine_prime"]["hits_at_20"],
            "filters_improve_prime": (
                evaluation["deterministic_dynamic_router_prime"]["hits_at_20"]
                > evaluation["baseline_genuine_prime"]["hits_at_20"]
            ),
            "full_stack_beats_prime": (
                evaluation["full_current_epicurus_stack"]["hits_at_20"]
                > evaluation["baseline_genuine_prime"]["hits_at_20"]
            ),
        },
    }
    (OUT / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    (OUT / "REPORT.md").write_text(_markdown(report))
    _write_selected_csv(base, frozen, positives)
    print(json.dumps(report["verdict"], indent=2))
    return 0


def _write_selected_csv(base: pd.DataFrame, frozen: dict[str, dict], positives: set[str]) -> None:
    identity = base.set_index("candidate_id")
    rows = []
    for arm, selection in frozen.items():
        for rank, candidate_id in enumerate(selection["selected_candidate_ids"], start=1):
            row = identity.loc[candidate_id]
            rows.append({
                "arm": arm,
                "selected_rank": rank,
                "candidate_id": candidate_id,
                "mutation_id": row["mutation_id"],
                "gene_symbol": row["gene_symbol"],
                "mutant_peptide": row["mutant_peptide"],
                "hla_allele": row["hla_allele"],
                "is_recognized_evaluation_only": row["mutation_id"] in positives,
            })
    pd.DataFrame(rows).to_csv(OUT / "selected_candidates.csv", index=False)


def _markdown(r: dict) -> str:
    effects = r["filter_effects_before_label_join"]
    ev = r["evaluation_only"]
    lines = [
        "# Sid full Epicurus filter-stack benchmark",
        "",
        "> Post-hoc n=1 patient / 3 recognized mutations. All selections were serialized before the exact "
        "Hudson labels were imported. Upstream generation covers 130/147 eligible mutations (88.4%), so "
        "this is a diagnostic, not a general superiority claim.",
        "",
        "## Filter effects (label-blind)",
        "",
        f"- Deterministic validity: removed **{effects['deterministic']['removed_count']}** candidate rows "
        f"and **{effects['deterministic']['mutations_removed']}** mutations. Sid HLA-LOH is unavailable.",
        f"- Dynamic gate v1: removed **{effects['dynamic_v1']['candidate_rows_removed']:,}** rows "
        f"({100 * effects['dynamic_v1']['candidate_row_removal_fraction']:.1f}%) but "
        f"**{effects['dynamic_v1']['mutations_removed']} mutations**.",
        f"- Product RNA policy: excluded **{effects['product_v1_expression_evidence_policy']['candidate_rows_excluded']:,}** "
        f"rows; {effects['product_v1_expression_evidence_policy']['mutations_remaining']} mutations remained.",
        "",
        "## Recognized mutations in the selected top 20",
        "",
        "| Arm | Hits / 3 | Selected routes | Unique mutations | Hit mutations |",
        "|---|---:|---:|---:|---|",
    ]
    for name, result in ev.items():
        lines.append(
            f"| `{name}` | **{result['hits_at_20']}/3** | {result['selected_routes']} | "
            f"{result['unique_selected_mutations']} | {', '.join(result['hit_variant_ids']) or 'none'} |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        f"- Genuine PRIME baseline: **{r['verdict']['genuine_prime_hits_at_20']}/3**.",
        f"- Complete currently runnable Epicurus stack: **{r['verdict']['full_stack_hits_at_20']}/3**.",
        f"- Filters improve PRIME on this patient: **{r['verdict']['filters_improve_prime']}**.",
        f"- Full stack beats PRIME on this patient: **{r['verdict']['full_stack_beats_prime']}**.",
        "",
        "The deterministic gate has nothing impossible to remove in the generated pool. The dynamic gate "
        "removes low-scoring peptide×HLA routes but no whole mutations and does not change PRIME's "
        "mutation-level top-20. With the newly recovered matched T2 RNA counts, the legacy product policy "
        "hard-excludes zero-mutant-read mutations and falls to 1/3. That policy is over-aggressive here: "
        "MAP2 is experimentally recognized despite zero observed mutant-allele RNA reads.",
        "",
        "Matched T2 RNA depth, mutant reads, and RNA VAF are now evaluated for every generated mutation. "
        "Sid-specific HLA-LOH remains unavailable and is not silently imputed. A separate absence-safe RNA "
        "gate is reported in sid_benchmark/rna_support_arm.json.",
    ])
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())
