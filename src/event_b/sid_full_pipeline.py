"""Leakage-safe helpers for the Sid full-filter-stack benchmark.

This module deliberately knows nothing about the Hudson recognition labels.  It prepares the matched
Sid candidate/evidence frame and freezes mutation-level or portfolio selections.  Evaluation labels are
joined by the runner only after the selections have been serialized.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


FORBIDDEN_EVALUATION_COLUMNS = frozenset({"label", "is_recognized", "recognized_label", "hudson_label"})


def assert_label_blind(frame: pd.DataFrame) -> None:
    """Fail closed if an evaluation-label column reaches a generation/gating/scoring frame."""
    present = FORBIDDEN_EVALUATION_COLUMNS & set(frame.columns)
    if present:
        raise ValueError(f"evaluation label column(s) reached the pipeline: {sorted(present)}")


def prepare_sid_gate_frame(
    scored_candidates_path: str | Path,
    vaf_path: str | Path,
) -> pd.DataFrame:
    """Attach label-blind longitudinal WES provenance and canonical gate/product feature aliases."""
    candidates = pd.read_csv(scored_candidates_path)
    assert_label_blind(candidates)

    vaf = pd.read_csv(vaf_path, sep="\t", low_memory=False)
    tumor = vaf[
        vaf["tissue"].astype(str).str.lower().eq("tumor")
        & (pd.to_numeric(vaf["alt_reads"], errors="coerce").fillna(0) > 0)
    ].copy()
    # DNA provenance must be DNA-only.  Earlier benchmark code accidentally allowed RNA/scRNA VAFs
    # into ``dna_vaf`` by aggregating every tumour assay together.
    tumor_dna = tumor[tumor["assay_type"].isin({"WES", "WGS"})].copy()
    provenance = tumor_dna.groupby("variant_id").agg(
        n_callers=("pipeline", lambda x: x.dropna().astype(str).nunique()),
        n_timepoints=("timepoint", lambda x: x.dropna().astype(str).nunique()),
        dna_vaf=("vaf", "max"),
    ).reset_index()
    out = candidates.merge(provenance, left_on="mutation_id", right_on="variant_id", how="left")
    out = out.drop(columns=["variant_id"])

    # Matched decision-time RNA: T2 bulk RNA from UCLA, the same tumour/timepoint as the RSEM TPM and
    # pVAC package.  Multiple rows are alternative count derivations of the same public data; max depth,
    # alt reads, and VAF preserve observed support without summing duplicates.  Longitudinal features are
    # kept separately as rescue/provenance and never substituted for the matched T2 values.
    t2_rna = vaf[
        vaf["tissue"].astype(str).str.lower().eq("tumor")
        & vaf["assay_type"].eq("RNA")
        & vaf["timepoint"].eq("T2")
        & vaf["data_source"].eq("UCLA")
    ].copy()
    t2_evidence = t2_rna.groupby("variant_id").agg(
        rna_depth=("total_reads", "max"),
        rna_mutant_reads=("alt_reads", "max"),
        rna_vaf=("vaf", "max"),
    ).reset_index()
    out = out.merge(t2_evidence, left_on="mutation_id", right_on="variant_id", how="left").drop(
        columns=["variant_id"])

    all_rna = vaf[
        vaf["tissue"].astype(str).str.lower().eq("tumor")
        & vaf["assay_type"].isin({"RNA", "scRNA", "scRNA_ONT"})
        & vaf["timepoint"].isin({"T0", "T1", "T2"})
    ].copy()
    all_rna["_supported_timepoint"] = all_rna["timepoint"].where(
        pd.to_numeric(all_rna["alt_reads"], errors="coerce").fillna(0) > 0)
    longitudinal = all_rna.groupby("variant_id").agg(
        longitudinal_rna_max_alt_reads=("alt_reads", "max"),
        longitudinal_rna_max_vaf=("vaf", "max"),
        longitudinal_rna_positive_timepoints=("_supported_timepoint", lambda x: x.dropna().astype(str).nunique()),
    ).reset_index()
    out = out.merge(longitudinal, left_on="mutation_id", right_on="variant_id", how="left").drop(
        columns=["variant_id"])

    # Stable candidate identity is independent of row order and contains no assay label.
    identity = (
        out["mutation_id"].astype(str) + "|" + out["mutant_peptide"].astype(str) + "|"
        + out["hla_allele"].astype(str)
    )
    out["candidate_id"] = identity.map(lambda value: f"sid:{hashlib.sha256(value.encode()).hexdigest()[:20]}")
    if out["candidate_id"].duplicated().any():
        raise ValueError("candidate identity is not unique after generation deduplication")

    out["el"] = pd.to_numeric(out["mixmhcpred_rank"], errors="coerce")
    out["prime"] = pd.to_numeric(out["prime_rank"], errors="coerce")
    out["expr"] = pd.to_numeric(out["expression_tpm"], errors="coerce")
    out["genuine_prime_score"] = -out["prime"]
    out["binding_percentile_rank"] = out["el"]
    out["presentation_score"] = (1.0 - out["el"] / 100.0).clip(0.0, 1.0)
    out["recognition_score"] = (1.0 - out["prime"] / 100.0).clip(0.0, 1.0)
    out["hla_loh_call"] = ""  # unavailable for Sid; explicit missing evidence, never imputed
    out["expression_call"] = ""  # use measured TPM; no vendor binary call at this boundary
    assert_label_blind(out)
    return out


def freeze_mutation_topk(
    frame: pd.DataFrame,
    score_column: str,
    *,
    k: int = 20,
    ascending: bool = False,
) -> dict:
    """Freeze a mutation-level top-k after choosing the best peptide×HLA route per mutation."""
    assert_label_blind(frame)
    ranked = (
        frame.dropna(subset=[score_column])
        .sort_values(score_column, ascending=ascending, kind="mergesort")
        .drop_duplicates("mutation_id", keep="first")
        .reset_index(drop=True)
    )
    ranked["mutation_rank"] = np.arange(1, len(ranked) + 1)
    top = ranked.head(k)
    return {
        "selection_unit": "mutation",
        "k": k,
        "n_candidate_rows": int(len(frame)),
        "n_ranked_mutations": int(len(ranked)),
        "selected_mutation_ids": top["mutation_id"].astype(str).tolist(),
        "selected_candidate_ids": top["candidate_id"].astype(str).tolist(),
        "all_mutation_ranks": dict(zip(ranked["mutation_id"].astype(str), ranked["mutation_rank"].astype(int), strict=False)),
    }


def freeze_portfolio(selection: pd.DataFrame, selected_column: str, rank_column: str, *, k: int = 20) -> dict:
    """Freeze an actual peptide×HLA portfolio, preserving duplicate mutation routes when policy allows."""
    assert_label_blind(selection)
    selected = selection[selection[selected_column].astype(bool)].copy()
    selected = selected.sort_values(rank_column, kind="mergesort").head(k)
    return {
        "selection_unit": "peptide_hla_route",
        "k": k,
        "n_candidate_rows": int(len(selection)),
        "n_selected_routes": int(len(selected)),
        "n_unique_selected_mutations": int(selected["mutation_id"].nunique()),
        "selected_mutation_ids": selected["mutation_id"].astype(str).tolist(),
        "selected_candidate_ids": selected["candidate_id"].astype(str).tolist(),
        "selected_route_ranks": selected[rank_column].astype(int).tolist(),
    }


def evaluate_frozen(selection: dict, positive_ids: set[str]) -> dict:
    """Evaluation-only join: called after the selection has been frozen to disk."""
    selected = selection["selected_mutation_ids"]
    unique = list(dict.fromkeys(selected))
    hits = sorted(set(unique) & set(positive_ids))
    all_ranks = selection.get("all_mutation_ranks", {})
    return {
        "hits_at_20": len(hits),
        "recall_at_20": round(len(hits) / len(positive_ids), 4) if positive_ids else None,
        "hit_variant_ids": hits,
        "positive_ranks": {pid: int(all_ranks[pid]) for pid in sorted(positive_ids) if pid in all_ranks},
        "selected_routes": len(selected),
        "unique_selected_mutations": len(unique),
    }
