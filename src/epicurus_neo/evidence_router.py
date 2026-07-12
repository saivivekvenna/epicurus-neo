"""Additive inference-time candidate evidence router + constrained route-aware top-k selection.

Frozen preregistration: ``docs/superpowers/specs/
2026-07-12-evidence-router-and-route-aware-selection-preregistration.md``.
Frozen policy: ``configs/frozen/evidence_router_v1.json``.

This module is **additive and non-destructive**. It composes the v1 product scorer's public output
by import (it reads the normalized candidate columns and, for selection, an incumbent score column
such as ``epicurus_lower_evidence_score``) and adds v2 columns
(``primary_route``, orthogonal ``flag_*``, ``rankable``, ``router_removed_reason`` …). It never
overwrites v1 columns and never mutates ``gates.py`` / ``product.py``. The legacy deterministic gate
(``apply_deterministic_gate``) is preserved unchanged for backward compatibility.

Key corrections vs. the legacy gate (motivated by the already-read Sid structural audit, which is
therefore NOT independent validation):

* ``GENE_NOT_EXPRESSED`` is **not** a router impossibility — cross-sectional RNA absence at one
  timepoint/region is not biological impossibility. It is demoted to ``flag_weak_or_absent_rna``.
* An **empty peptide** is **not** ``MALFORMED_AA`` — it is an upstream candidate-generation gap,
  reported as ``NEEDS_PEPTIDE_GENERATION`` and never removed and never charged to the ranker.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd

_POLICY_PATH = (
    Path(__file__).resolve().parents[2] / "configs" / "frozen" / "evidence_router_v1.json"
)


@dataclass(frozen=True)
class RouterPolicy:
    policy_id: str = "epicurus-evidence-router-1.0.0"
    amino_acids: frozenset[str] = frozenset("ACDEFGHIKLMNPQRSTVWY")
    class_i_min_len: int = 8
    class_i_max_len: int = 14
    missense_variant_types: frozenset[str] = frozenset({"SNV", "SNP", "MISSENSE"})
    not_expressed_tokens: frozenset[str] = frozenset({"N", "NO", "FALSE", "0", "NOT_EXPRESSED"})
    lost_allele_tokens: frozenset[str] = frozenset({"Y", "YES", "TRUE", "1", "LOST"})
    k: int = 20
    reserve_per_route: int = 1
    max_reserve: int = 3
    reserve_routes: tuple[str, ...] = ("RESCUE", "LONGITUDINAL", "UNCERTAIN")
    max_per_mutation: int | None = 2
    max_per_gene: int | None = 4
    max_per_hla: int | None = None
    incumbent_score_default_column: str = "epicurus_lower_evidence_score"


def load_router_policy(path: str | Path | None = None) -> RouterPolicy:
    """Load the frozen router policy from JSON (defaults to the committed v1 policy)."""
    source = Path(path) if path is not None else _POLICY_PATH
    data = json.loads(Path(source).read_text())
    selection = data.get("selection", {})
    return RouterPolicy(
        policy_id=data["policy_id"],
        amino_acids=frozenset(data["amino_acids"]),
        class_i_min_len=int(data["class_i_len"]["min"]),
        class_i_max_len=int(data["class_i_len"]["max"]),
        missense_variant_types=frozenset(data["missense_variant_types"]),
        not_expressed_tokens=frozenset(data["not_expressed_tokens"]),
        lost_allele_tokens=frozenset(data["lost_allele_tokens"]),
        k=int(selection["k"]),
        reserve_per_route=int(selection["reserve_per_route"]),
        max_reserve=int(selection["max_reserve"]),
        reserve_routes=tuple(selection["reserve_routes"]),
        max_per_mutation=selection["max_per_mutation"],
        max_per_gene=selection["max_per_gene"],
        max_per_hla=selection["max_per_hla"],
        incumbent_score_default_column=selection["incumbent_score_default_column"],
    )


try:  # Prefer the frozen JSON so code and config never drift; fall back to dataclass defaults.
    DEFAULT_ROUTER_POLICY = load_router_policy()
except (FileNotFoundError, KeyError, ValueError):  # pragma: no cover - defensive
    DEFAULT_ROUTER_POLICY = RouterPolicy()

POLICY_ID = DEFAULT_ROUTER_POLICY.policy_id


# ---------------------------------------------------------------------------
# Cell-level helpers (present/populated semantics, never inventing evidence)
# ---------------------------------------------------------------------------
def _text(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


def _upper(value: object) -> str:
    return _text(value).upper()


def _has_number(row: pd.Series, column: str) -> bool:
    return column in row and pd.notna(row.get(column))


def _number(row: pd.Series, column: str) -> float | None:
    if not _has_number(row, column):
        return None
    try:
        return float(row.get(column))
    except (TypeError, ValueError):
        return None


def _md5(value: str) -> str:
    return hashlib.md5(value.encode()).hexdigest()


# ---------------------------------------------------------------------------
# §3 Router: flags, IMPOSSIBLE rules, primary route, reachability
# ---------------------------------------------------------------------------
def _flags(row: pd.Series, policy: RouterPolicy) -> dict[str, bool]:
    peptide = _upper(row.get("mutant_peptide"))
    hla = _text(row.get("hla_allele"))
    variant_type = _upper(row.get("source_variant_type"))

    atypical = bool(variant_type) and variant_type not in policy.missense_variant_types

    weak_rna = False
    tpm = _number(row, "expression_tpm")
    if tpm is not None and tpm == 0:
        weak_rna = True
    if _text(row.get("expression_call")) and _upper(row.get("expression_call")) in policy.not_expressed_tokens:
        weak_rna = True
    reads = _number(row, "rna_mutant_reads")
    if reads is not None and reads == 0:
        weak_rna = True

    rna_columns = ("expression_tpm", "rna_vaf", "rna_mutant_reads", "expression_call")
    missing_rna = not any(
        (_has_number(row, c) if c != "expression_call" else bool(_text(row.get("expression_call"))))
        for c in rna_columns
    )

    # Provenance counts drive single/multi-source flags; unknown provenance asserts neither.
    counts = [int(row[c]) for c in ("n_callers", "n_timepoints", "n_regions") if _has_number(row, c)]
    single_caller = _has_number(row, "n_callers") and int(row["n_callers"]) == 1
    multi_source = any(c >= 2 for c in counts)

    has_presentation = any(
        _has_number(row, c)
        for c in ("presentation_score", "binding_percentile_rank", "binding_affinity_nm")
    )

    conflicting = bool(_text(row.get("representation_conflicts"))) or bool(
        _text(row.get("recognition_conflict"))
    )

    needs_peptide_generation = (peptide == "") or (hla == "")

    return {
        "flag_atypical_variant_class": atypical,
        "flag_weak_or_absent_rna": weak_rna,
        "flag_missing_rna": missing_rna,
        "flag_single_caller": bool(single_caller),
        "flag_multi_source_support": multi_source,
        "flag_has_presentation": has_presentation,
        "flag_conflicting_evidence": conflicting,
        "flag_needs_peptide_generation": needs_peptide_generation,
    }


def _impossible_reason(row: pd.Series, is_duplicate: bool, policy: RouterPolicy) -> str:
    """First-wins IMPOSSIBLE reason, in the frozen ``impossible_rule_order``.

    Order (route-verifiable only): MALFORMED_AA, BAD_CLASS_I_LENGTH, MUT_NOT_IN_PEPTIDE,
    DUP_CANDIDATE, HLA_LOH_LOST_ALLELE. ``GENE_NOT_EXPRESSED`` and the empty-peptide arm of
    ``MALFORMED_AA`` are intentionally absent (demoted to flags / NEEDS_PEPTIDE_GENERATION).
    """
    peptide = _upper(row.get("mutant_peptide"))

    # 1. MALFORMED_AA — non-empty peptide with a non-standard residue. Empty peptide is NEVER
    #    malformed (it is NEEDS_PEPTIDE_GENERATION, handled in reachability).
    if peptide and set(peptide) - policy.amino_acids:
        return "MALFORMED_AA"

    # 2. BAD_CLASS_I_LENGTH — class I with a genuinely impossible length. Only for a non-empty
    #    peptide (an empty peptide is a generation gap, not a bad length).
    mhc_class = _upper(row.get("mhc_class"))
    if peptide and mhc_class in {"I", "CLASS_I", "CLASS I"}:
        if not (policy.class_i_min_len <= len(peptide) <= policy.class_i_max_len):
            return "BAD_CLASS_I_LENGTH"

    # 3. MUT_NOT_IN_PEPTIDE — missense only, and only when a peptide exists to check.
    variant_type = _upper(row.get("source_variant_type"))
    if peptide and variant_type in policy.missense_variant_types:
        residue = _missense_mutant_residue(row.get("protein_variant"))
        if residue and residue not in peptide:
            return "MUT_NOT_IN_PEPTIDE"

    # 4. DUP_CANDIDATE — decided frame-wide (keep=first); precedes HLA_LOH per the frozen order.
    if is_duplicate:
        return "DUP_CANDIDATE"

    # 5. HLA_LOH_LOST_ALLELE — fires only for a specified peptide-HLA route.
    if _text(row.get("hla_allele")) and _upper(row.get("hla_loh_call")) in policy.lost_allele_tokens:
        return "HLA_LOH_LOST_ALLELE"

    return ""


def _missense_mutant_residue(value: object) -> str | None:
    match = re.fullmatch(r"P\.([A-Z])\d+([A-Z])", _upper(value))
    return match.group(2) if match else None


def _duplicate_mask(frame: pd.DataFrame) -> pd.Series:
    identity = [
        c for c in ("patient_id", "mutation_id", "mutant_peptide", "hla_allele") if c in frame
    ]
    if len(identity) >= 3:
        return frame.duplicated(identity, keep="first")
    return pd.Series(False, index=frame.index)


def _primary_route(flags: dict[str, bool], removed_reason: str) -> str:
    if removed_reason:
        return "IMPOSSIBLE"
    # RESCUE: would-have-been-dropped-but-supported (highest non-IMPOSSIBLE precedence).
    if (flags["flag_atypical_variant_class"] or flags["flag_weak_or_absent_rna"]) and (
        flags["flag_has_presentation"] or flags["flag_multi_source_support"]
    ):
        return "RESCUE"
    # LONGITUDINAL: secondary timepoint/region/caller support, not already RESCUE.
    if flags["flag_multi_source_support"]:
        return "LONGITUDINAL"
    # UNCERTAIN: missing/conflicting core evidence or no presentation.
    if (
        flags["flag_missing_rna"]
        or flags["flag_conflicting_evidence"]
        or flags["flag_needs_peptide_generation"]
        or not flags["flag_has_presentation"]
    ):
        return "UNCERTAIN"
    return "CORE"


def route_candidates(
    frame: pd.DataFrame, policy: RouterPolicy = DEFAULT_ROUTER_POLICY
) -> pd.DataFrame:
    """Route normalized candidates, adding additive v2 columns without overwriting v1 columns.

    Adds: ``primary_route``, ``router_eligible``, ``router_removed_reason``, ``router_status``,
    ``rankable``, the eight orthogonal ``flag_*`` booleans, and ``router_policy_id``.
    """
    out = frame.reset_index(drop=True).copy()
    duplicate = _duplicate_mask(out)

    routes, eligible, removed, status, rankable = [], [], [], [], []
    flag_records: list[dict[str, bool]] = []
    for position, (_, row) in enumerate(out.iterrows()):
        flags = _flags(row, policy)
        reason = _impossible_reason(row, bool(duplicate.iloc[position]), policy)

        route = _primary_route(flags, reason)
        is_eligible = reason == ""
        is_rankable = is_eligible and not flags["flag_needs_peptide_generation"]
        if not is_eligible:
            row_status = "IMPOSSIBLE"
        elif flags["flag_needs_peptide_generation"]:
            row_status = "NEEDS_PEPTIDE_GENERATION"
        else:
            row_status = "RANKABLE"

        routes.append(route)
        eligible.append(is_eligible)
        removed.append(reason)
        status.append(row_status)
        rankable.append(is_rankable)
        flag_records.append(flags)

    for name in flag_records[0] if flag_records else []:
        out[name] = [record[name] for record in flag_records]
    out["primary_route"] = routes
    out["router_eligible"] = eligible
    out["router_removed_reason"] = removed
    out["router_status"] = status
    out["rankable"] = rankable
    out["router_policy_id"] = policy.policy_id
    return out


# ---------------------------------------------------------------------------
# §5 Constrained, route-aware top-k selection
# ---------------------------------------------------------------------------
def select_route_aware_topk(
    routed: pd.DataFrame,
    *,
    score_column: str | None = None,
    policy: RouterPolicy = DEFAULT_ROUTER_POLICY,
    patient_column: str = "patient_id",
) -> pd.DataFrame:
    """Select up to ``k`` eligible+rankable candidates per patient, deterministically.

    Modest exploration reserves guarantee representation of every present non-CORE route (capped at
    ``max_reserve`` total), diversity caps bound per-mutation/gene/HLA, freed reserves backfill by
    score, and ties break on ``md5(mutant_peptide|hla_allele)`` for permutation-invariant output.
    The selection is score-agnostic: it only orders and reserves an incumbent score never fit here.
    """
    column = score_column or policy.incumbent_score_default_column
    if column not in routed.columns:
        raise ValueError(
            f"route-aware selection needs an incumbent score column '{column}' "
            "(genuine PRIME where available, else epicurus_lower_evidence_score); none found"
        )

    out = routed.reset_index(drop=True).copy()
    empty = pd.Series("", index=out.index)
    peptide = out["mutant_peptide"] if "mutant_peptide" in out else empty
    hla = out["hla_allele"] if "hla_allele" in out else empty
    # NaN-safe, dtype-agnostic (real CSVs may carry arrow-backed strings / NaN peptides).
    out["_tie_key"] = [_md5(f"{_text(p)}|{_text(h)}") for p, h in zip(peptide, hla)]
    out["route_selected"] = False
    out["route_rank"] = pd.Series(pd.NA, index=out.index, dtype="Int64")
    out["route_selection_kind"] = ""

    patients = out[patient_column] if patient_column in out else pd.Series("_all", index=out.index)
    for _, index in out.groupby(patients, sort=True).groups.items():
        selected = _select_one_patient(out.loc[index], column, policy)
        for rank, (idx, kind) in enumerate(selected, start=1):
            out.loc[idx, "route_selected"] = True
            out.loc[idx, "route_rank"] = rank
            out.loc[idx, "route_selection_kind"] = kind

    return out


def _select_one_patient(
    rows: pd.DataFrame, column: str, policy: RouterPolicy
) -> list[tuple[int, str]]:
    pool = rows[rows["router_eligible"].astype(bool) & rows["rankable"].astype(bool)]
    if pool.empty:
        return []
    ordered = pool.sort_values(
        [column, "_tie_key"], ascending=[False, True], kind="mergesort"
    )

    selected: list[tuple[int, str]] = []
    chosen: set[int] = set()
    counts: dict[tuple[str, str], int] = {}

    def _caps_ok(row: pd.Series) -> bool:
        limits = (
            ("mutation", _text(row.get("mutation_id")), policy.max_per_mutation),
            ("gene", _text(row.get("gene_symbol")), policy.max_per_gene),
            ("hla", _text(row.get("hla_allele")), policy.max_per_hla),
        )
        return not any(
            value and limit is not None and counts.get((kind, value), 0) >= limit
            for kind, value, limit in limits
        )

    def _take(idx: int, row: pd.Series, kind: str) -> None:
        selected.append((idx, kind))
        chosen.add(idx)
        for axis, value in (
            ("mutation", _text(row.get("mutation_id"))),
            ("gene", _text(row.get("gene_symbol"))),
            ("hla", _text(row.get("hla_allele"))),
        ):
            if value:
                counts[(axis, value)] = counts.get((axis, value), 0) + 1

    # (a) Exploration reserves: >=1 slot per present non-CORE route, capped at max_reserve total.
    reserves_used = 0
    for route in policy.reserve_routes:
        if reserves_used >= policy.max_reserve or len(selected) >= policy.k:
            break
        route_rows = ordered[(ordered["primary_route"] == route) & ~ordered.index.isin(chosen)]
        for idx, row in route_rows.iterrows():
            if _caps_ok(row):
                _take(idx, row, f"RESERVE:{route}")
                reserves_used += 1
                break
        # If no candidate of this route can be placed (absent or caps-exhausted), the slot is freed
        # back to the score-fill pool below.

    # (b) Score-fill: remaining slots by descending incumbent score across all eligible+rankable.
    for idx, row in ordered.iterrows():
        if len(selected) >= policy.k:
            break
        if idx in chosen or not _caps_ok(row):
            continue
        _take(idx, row, "SCORE")

    # Reserves change set membership, not rank semantics. Return selected rows in incumbent-score
    # order so a low-score exploratory reserve does not masquerade as rank #1.
    kind_by_index = dict(selected)
    return [(idx, kind_by_index[idx]) for idx in ordered.index if idx in chosen]


# ---------------------------------------------------------------------------
# §6 Machine-readable per-stage funnel
# ---------------------------------------------------------------------------
def build_funnel(
    routed: pd.DataFrame,
    *,
    patient_column: str = "patient_id",
    policy: RouterPolicy = DEFAULT_ROUTER_POLICY,
) -> list[dict[str, Any]]:
    """Per-patient ``generated -> valid -> rankable -> selected`` funnel.

    Stage semantics make it impossible to read an upstream generation gap as a ranker miss:
    ``needs_peptide_generation`` rows are ``valid`` but not ``rankable``.
    """
    has_selection = "route_selected" in routed.columns
    patients = routed[patient_column] if patient_column in routed else pd.Series(
        "_all", index=routed.index
    )
    entries: list[dict[str, Any]] = []
    for patient, group in routed.groupby(patients, sort=True):
        valid = group[group["router_eligible"].astype(bool)]
        rankable = valid[valid["rankable"].astype(bool)]
        removed = group[~group["router_eligible"].astype(bool)]
        entries.append(
            {
                "patient_id": str(patient),
                "policy_id": policy.policy_id,
                "generated": int(len(group)),
                "valid": int(len(valid)),
                "rankable": int(len(rankable)),
                "selected": int(group["route_selected"].astype(bool).sum()) if has_selection else 0,
                "needs_peptide_generation": int(
                    (group["router_status"] == "NEEDS_PEPTIDE_GENERATION").sum()
                ),
                "router_removed_reason_counts": {
                    str(reason): int(count)
                    for reason, count in removed["router_removed_reason"].value_counts().items()
                },
                "route_composition_valid": {
                    str(route): int(count)
                    for route, count in valid["primary_route"].value_counts().items()
                },
                "route_composition_selected": (
                    {
                        str(route): int(count)
                        for route, count in group.loc[
                            group["route_selected"].astype(bool), "primary_route"
                        ].value_counts().items()
                    }
                    if has_selection
                    else {}
                ),
            }
        )
    return entries


def write_router_report(
    routed: pd.DataFrame,
    output_dir: str | Path,
    *,
    policy: RouterPolicy = DEFAULT_ROUTER_POLICY,
    patient_column: str = "patient_id",
) -> dict[str, str]:
    """Write the routed candidates (CSV) plus a JSON + Markdown per-stage funnel."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "routed_candidates.csv"
    json_path = output / "router_funnel.json"
    markdown_path = output / "router_funnel.md"

    routed.to_csv(csv_path, index=False)
    funnel = build_funnel(routed, patient_column=patient_column, policy=policy)
    payload = {
        "policy_id": policy.policy_id,
        "funnel_stages": ["generated", "valid", "rankable", "selected"],
        "patients": funnel,
    }
    json_path.write_text(json.dumps(payload, indent=2) + "\n")

    lines = [
        "# Evidence router funnel",
        "",
        f"> Policy `{policy.policy_id}`. Stages: generated → valid → rankable → selected.",
        "> `NEEDS_PEPTIDE_GENERATION` rows are valid but not rankable — an upstream generation gap,",
        "> never a ranker miss.",
        "",
        "| Patient | Generated | Valid | Rankable | Selected | NeedsPeptideGen | Removed |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for entry in funnel:
        lines.append(
            "| {patient} | {generated} | {valid} | {rankable} | {selected} | {npg} | {removed} |".format(
                patient=entry["patient_id"],
                generated=entry["generated"],
                valid=entry["valid"],
                rankable=entry["rankable"],
                selected=entry["selected"],
                npg=entry["needs_peptide_generation"],
                removed=entry["router_removed_reason_counts"] or {},
            )
        )
    markdown_path.write_text("\n".join(lines) + "\n")
    return {"csv": str(csv_path), "json": str(json_path), "markdown": str(markdown_path)}
