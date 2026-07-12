"""Tests for the additive inference-time evidence router + route-aware top-k selection.

Frozen preregistration: docs/superpowers/specs/
2026-07-12-evidence-router-and-route-aware-selection-preregistration.md (§3, §5, §6, §9).
Frozen policy: configs/frozen/evidence_router_v1.json.

The router hard-removes ONLY route-verifiable impossibilities, keeps cross-sectional-RNA-absent /
atypical-class / single-caller candidates eligible-but-flagged, reports "no peptide/HLA" as the
distinct upstream status NEEDS_PEPTIDE_GENERATION (never a ranker miss), and feeds a constrained
route-aware top-k with modest exploration reserves.
"""

import pandas as pd

from epicurus_neo.evidence_router import (
    DEFAULT_ROUTER_POLICY,
    build_funnel,
    route_candidates,
    select_route_aware_topk,
)


def _base_row(**overrides) -> dict:
    row = {
        "patient_id": "sid",
        "candidate_id": "c0",
        "gene_symbol": "GENE",
        "mutation_id": "m0",
        "protein_variant": "",
        "source_variant_type": "SNV",
        "mhc_class": "I",
        "mutant_peptide": "GVSVEIALK",
        "hla_allele": "HLA-A*02:01",
        "hla_loh_call": "N",
        "expression_call": "Y",
        "expression_tpm": 40.0,
        "rna_mutant_reads": 20.0,
        "presentation_score": 0.6,
    }
    row.update(overrides)
    return row


# ---------------------------------------------------------------------------
# §9.1 Sid-like ASPM/MAP2 rescue routing
# ---------------------------------------------------------------------------
def test_atypical_weak_rna_with_support_routes_rescue_and_is_retained():
    # Atypical variant class (frameshift) + weak/absent RNA (expression_call=N), but with
    # presentation AND multi-caller support -> RESCUE, eligible, never removed.
    frame = pd.DataFrame(
        [
            _base_row(
                gene_symbol="MAP2",
                source_variant_type="frameshift",
                expression_call="N",
                expression_tpm=0.0,
                n_callers=3,
            )
        ]
    )
    routed = route_candidates(frame)
    row = routed.iloc[0]
    assert row["primary_route"] == "RESCUE"
    assert bool(row["router_eligible"]) is True
    assert row["router_removed_reason"] == ""
    assert bool(row["flag_atypical_variant_class"]) is True
    assert bool(row["flag_weak_or_absent_rna"]) is True


# ---------------------------------------------------------------------------
# §9.2 DYNC reachability — conventional supported SNV is CORE, rankable, selectable
# ---------------------------------------------------------------------------
def test_conventional_snv_with_peptide_hla_is_core_rankable_and_selectable():
    frame = pd.DataFrame([_base_row(gene_symbol="DYNC1H1", expression_tpm=357.14)])
    routed = route_candidates(frame)
    row = routed.iloc[0]
    assert row["primary_route"] == "CORE"
    assert bool(row["rankable"]) is True
    assert row["router_status"] == "RANKABLE"

    selected = select_route_aware_topk(
        routed.assign(epicurus_lower_evidence_score=0.9)
    )
    assert bool(selected.iloc[0]["route_selected"]) is True


# ---------------------------------------------------------------------------
# §9.3 Expression-N retained in a route, not removed
# ---------------------------------------------------------------------------
def test_expression_n_retained_and_routed_not_removed():
    # expression_call=N / expression_tpm=0 with no support: flagged, routed UNCERTAIN, NOT removed.
    frame = pd.DataFrame(
        [
            _base_row(
                expression_call="N",
                expression_tpm=0.0,
                presentation_score=float("nan"),
                binding_percentile_rank=float("nan"),
                binding_affinity_nm=float("nan"),
            )
        ]
    )
    routed = route_candidates(frame)
    row = routed.iloc[0]
    assert bool(row["router_eligible"]) is True
    assert row["router_removed_reason"] == ""
    assert bool(row["flag_weak_or_absent_rna"]) is True
    assert row["primary_route"] == "UNCERTAIN"


# ---------------------------------------------------------------------------
# §9.4 Lost allele -> IMPOSSIBLE
# ---------------------------------------------------------------------------
def test_lost_allele_is_impossible_and_removed():
    frame = pd.DataFrame([_base_row(hla_loh_call="Y")])
    routed = route_candidates(frame)
    row = routed.iloc[0]
    assert row["primary_route"] == "IMPOSSIBLE"
    assert bool(row["router_eligible"]) is False
    assert row["router_removed_reason"] == "HLA_LOH_LOST_ALLELE"


def test_frozen_impossible_rule_precedence_puts_duplicate_before_lost_allele():
    first = _base_row(candidate_id="first", hla_loh_call="Y")
    duplicate = _base_row(candidate_id="duplicate", hla_loh_call="Y")
    routed = route_candidates(pd.DataFrame([first, duplicate])).set_index("candidate_id")
    assert routed.loc["first", "router_removed_reason"] == "HLA_LOH_LOST_ALLELE"
    assert routed.loc["duplicate", "router_removed_reason"] == "DUP_CANDIDATE"


# ---------------------------------------------------------------------------
# §9.5 Zero mutant RNA reads not excluded
# ---------------------------------------------------------------------------
def test_zero_mutant_rna_reads_flags_but_never_removes():
    frame = pd.DataFrame([_base_row(rna_mutant_reads=0.0)])
    routed = route_candidates(frame)
    row = routed.iloc[0]
    assert bool(row["router_eligible"]) is True
    assert row["router_removed_reason"] == ""
    assert bool(row["flag_weak_or_absent_rna"]) is True


# ---------------------------------------------------------------------------
# §9.6 Empty peptide is NEEDS_PEPTIDE_GENERATION, not MALFORMED_AA
# ---------------------------------------------------------------------------
def test_empty_peptide_is_needs_peptide_generation_and_bad_aa_is_removed():
    frame = pd.DataFrame(
        [
            _base_row(candidate_id="empty", mutant_peptide=""),  # upstream generation gap
            _base_row(candidate_id="badaa", mutant_peptide="XBJOUZ"),  # genuine bad AA
        ]
    )
    routed = route_candidates(frame).set_index("candidate_id")

    empty = routed.loc["empty"]
    assert empty["primary_route"] != "IMPOSSIBLE"
    assert bool(empty["router_eligible"]) is True
    assert bool(empty["rankable"]) is False
    assert bool(empty["flag_needs_peptide_generation"]) is True
    assert empty["router_status"] == "NEEDS_PEPTIDE_GENERATION"

    bad = routed.loc["badaa"]
    assert bad["primary_route"] == "IMPOSSIBLE"
    assert bad["router_removed_reason"] == "MALFORMED_AA"


# ---------------------------------------------------------------------------
# §9.8 Route reserves / backfill / determinism
# ---------------------------------------------------------------------------
def _routed_fixture() -> pd.DataFrame:
    rows = []
    # 20 high-scoring CORE candidates (distinct mutations/genes/hlas), scores 0.90..0.71.
    for i in range(20):
        rows.append(
            {
                "patient_id": "p",
                "candidate_id": f"core{i}",
                "gene_symbol": f"G{i}",
                "mutation_id": f"m{i}",
                "hla_allele": f"HLA-A*{i:02d}:01",
                "mutant_peptide": f"PEPTIDE{i:03d}A",
                "primary_route": "CORE",
                "router_eligible": True,
                "rankable": True,
                "epicurus_lower_evidence_score": 0.90 - 0.01 * i,
            }
        )
    # One low-scoring RESCUE and one low-scoring UNCERTAIN; no LONGITUDINAL present.
    rows.append(
        {
            "patient_id": "p",
            "candidate_id": "rescue0",
            "gene_symbol": "RG",
            "mutation_id": "rm",
            "hla_allele": "HLA-B*07:02",
            "mutant_peptide": "RESCUEPEPA",
            "primary_route": "RESCUE",
            "router_eligible": True,
            "rankable": True,
            "epicurus_lower_evidence_score": 0.10,
        }
    )
    rows.append(
        {
            "patient_id": "p",
            "candidate_id": "uncertain0",
            "gene_symbol": "UG",
            "mutation_id": "um",
            "hla_allele": "HLA-B*08:01",
            "mutant_peptide": "UNCERTPEPA",
            "primary_route": "UNCERTAIN",
            "router_eligible": True,
            "rankable": True,
            "epicurus_lower_evidence_score": 0.05,
        }
    )
    return pd.DataFrame(rows)


def test_reserves_guarantee_non_core_routes_and_backfill_and_are_deterministic():
    routed = _routed_fixture()
    selected = select_route_aware_topk(routed)
    chosen = set(selected.loc[selected["route_selected"], "candidate_id"])

    # k=20 slots filled.
    assert len(chosen) == 20
    # Each present non-CORE route reserved a slot despite low scores.
    assert "rescue0" in chosen
    assert "uncertain0" in chosen
    # CORE keeps >= k - max_reserve (=17); exactly 18 here (2 non-CORE reserved).
    core_selected = sum(1 for cid in chosen if cid.startswith("core"))
    assert core_selected == 18
    assert core_selected >= 20 - DEFAULT_ROUTER_POLICY.max_reserve
    # The two lowest-scored CORE candidates are displaced by the reserves.
    assert "core19" not in chosen
    assert "core18" not in chosen

    # Reserve changes membership, not rank semantics: selected output remains incumbent-score ordered.
    ranked_scores = (
        selected.loc[selected["route_selected"]]
        .sort_values("route_rank")["epicurus_lower_evidence_score"]
        .tolist()
    )
    assert ranked_scores == sorted(ranked_scores, reverse=True)

    # Permutation invariance + determinism: shuffle rows, identical selection and rank order.
    shuffled = routed.sample(frac=1.0, random_state=7).reset_index(drop=True)
    selected_shuf = select_route_aware_topk(shuffled)
    order = (
        selected.loc[selected["route_selected"]]
        .sort_values("route_rank")["candidate_id"]
        .tolist()
    )
    order_shuf = (
        selected_shuf.loc[selected_shuf["route_selected"]]
        .sort_values("route_rank")["candidate_id"]
        .tolist()
    )
    assert order == order_shuf


# ---------------------------------------------------------------------------
# §9.9 No-regression sanity — only CORE reproduces the pure-score top-k exactly
# ---------------------------------------------------------------------------
def test_only_core_reproduces_pure_score_topk():
    rows = []
    for i in range(25):
        rows.append(
            {
                "patient_id": "p",
                "candidate_id": f"core{i:02d}",
                "gene_symbol": f"G{i}",
                "mutation_id": f"m{i}",
                "hla_allele": f"HLA-A*{i:02d}:01",
                "mutant_peptide": f"PEPTIDE{i:03d}A",
                "primary_route": "CORE",
                "router_eligible": True,
                "rankable": True,
                "epicurus_lower_evidence_score": 1.0 - 0.01 * i,
            }
        )
    routed = pd.DataFrame(rows)
    selected = select_route_aware_topk(routed)
    chosen = set(selected.loc[selected["route_selected"], "candidate_id"])

    pure_score_top20 = set(
        routed.sort_values("epicurus_lower_evidence_score", ascending=False)
        .head(20)["candidate_id"]
    )
    assert chosen == pure_score_top20
    # Rank order matches pure descending score (reserve is a no-op with no non-CORE routes).
    order = (
        selected.loc[selected["route_selected"]]
        .sort_values("route_rank")["candidate_id"]
        .tolist()
    )
    expected = (
        routed.sort_values("epicurus_lower_evidence_score", ascending=False)
        .head(20)["candidate_id"]
        .tolist()
    )
    assert order == expected


# ---------------------------------------------------------------------------
# §6 Funnel semantics — generated -> valid -> rankable -> selected
# ---------------------------------------------------------------------------
def test_funnel_separates_generation_gap_from_ranker_miss():
    frame = pd.DataFrame(
        [
            _base_row(candidate_id="core", gene_symbol="DYNC1H1"),
            _base_row(candidate_id="gap", mutant_peptide=""),  # NEEDS_PEPTIDE_GENERATION
            _base_row(candidate_id="lost", hla_loh_call="Y"),  # IMPOSSIBLE
        ]
    )
    routed = select_route_aware_topk(
        route_candidates(frame).assign(epicurus_lower_evidence_score=0.9)
    )
    funnel = build_funnel(routed)
    entry = next(f for f in funnel if f["patient_id"] == "sid")
    assert entry["generated"] == 3
    assert entry["valid"] == 2  # lost removed
    assert entry["rankable"] == 1  # gap is valid but not rankable
    assert entry["needs_peptide_generation"] == 1
    assert entry["selected"] == 1
