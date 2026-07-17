from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

from epicurus_neo.contracts import (
    CANDIDATE_SCHEMA_VERSION,
    RANKED_SCHEMA_VERSION,
    validate_candidate_contract,
    validate_ranked_contract,
)
from epicurus_neo.gates import apply_deterministic_gate, summarize_gate
from epicurus_neo.normalize import read_table
from epicurus_neo.schema import normalize_hla, normalize_peptide


_ALIASES: dict[str, tuple[str, ...]] = {
    "candidate_id": ("candidate_id", "Index", "ID"),
    "patient_id": ("patient_id", "Patient", "Sample", "sample_name"),
    "gene_symbol": ("gene_symbol", "Gene Name", "Gene", "SYMBOL", "Gene Symbol"),
    "transcript_id": ("transcript_id", "Transcript", "Ensembl Transcript ID"),
    "mutation_id": ("mutation_id", "Mutation", "variant_id", "Variant"),
    "protein_variant": ("protein_variant", "Protein Variant"),
    "source_variant_type": ("source_variant_type", "Source Variant Type", "Variant Type"),
    "mhc_class": ("mhc_class", "MHC Class"),
    "hla_allele": ("hla_allele", "HLA Allele", "Allele", "Best Allele", "HLA"),
    "hla_loh_call": ("hla_loh_call", "HLA LOH"),
    "expression_call": ("expression_call", "Expressed"),
    "mutant_peptide": (
        "mutant_peptide",
        "MT Epitope Seq",
        "Mut Epitope",
        "peptide",
        "Peptide",
    ),
    "wildtype_peptide": ("wildtype_peptide", "WT Epitope Seq", "Wt Epitope"),
    "dna_vaf": ("dna_vaf", "Tumor DNA VAF", "DNA VAF", "DNA Allelic Fraction"),
    "rna_depth": ("rna_depth", "Tumor RNA Depth", "RNA Depth"),
    "rna_vaf": ("rna_vaf", "Tumor RNA VAF", "RNA VAF", "RNA Allelic Fraction"),
    "rna_mutant_reads": ("rna_mutant_reads", "Tumor RNA Alt Read Support", "RNA Mutant Reads"),
    "expression_tpm": (
        "expression_tpm",
        "Gene Expression",
        "Transcript Expression",
        "RNA expression (TPM)",
        "Gene Level Expression TPM",
    ),
    "clonality_ccf": ("clonality_ccf", "CCF", "Cancer Cell Fraction"),
    "binding_affinity_nm": (
        "binding_affinity_nm",
        "Best MT Score",
        "Median MT IC50 Score",
        "IC50 MT",
    ),
    "binding_percentile_rank": (
        "binding_percentile_rank",
        "Best MT Percentile",
        "Median MT Percentile",
        "MT Percentile",
        "SHERPA Presentation Rank",
    ),
    "wildtype_binding_affinity_nm": (
        "wildtype_binding_affinity_nm",
        "Corresponding WT Score",
        "Median WT IC50 Score",
    ),
    "presentation_score": (
        "presentation_score",
        "mhcflurry_presentation_score",
        "bigmhc_el_score",
    ),
    "recognition_score": (
        "recognition_score",
        "prime_source_score",
        "foreignness_score",
    ),
}

_NUMERIC_COLUMNS = {
    "dna_vaf",
    "rna_depth",
    "rna_vaf",
    "rna_mutant_reads",
    "expression_tpm",
    "clonality_ccf",
    "binding_affinity_nm",
    "binding_percentile_rank",
    "wildtype_binding_affinity_nm",
    "presentation_score",
    "recognition_score",
}


@dataclass(frozen=True)
class InferenceConfig:
    k: int = 20
    max_per_mutation: int = 1
    max_per_gene: int = 4
    max_per_hla: int | None = None
    core_threshold: float = 0.55
    supporting_threshold: float = 0.35
    max_core_uncertainty: float = 0.35
    apply_validity_gate: bool = True


@dataclass(frozen=True)
class PatientSummary:
    patient_id: str
    candidate_count: int
    eligible_count: int
    selected_count: int
    core_count: int
    supporting_count: int
    filler_count: int
    abstained: bool
    abstention_reason: str | None


def _first(frame: pd.DataFrame, aliases: Iterable[str]) -> str | None:
    return next((column for column in aliases if column in frame.columns), None)


def _stable_candidate_id(row: pd.Series, patient_id: str) -> str:
    identity = "|".join(
        str(row.get(column, ""))
        for column in ("mutation_id", "transcript_id", "mutant_peptide", "hla_allele")
    )
    digest = hashlib.sha256(f"{patient_id}|{identity}".encode()).hexdigest()[:16]
    return f"{patient_id}:{digest}"


def normalize_product_candidates(
    frame: pd.DataFrame,
    *,
    patient_id: str | None = None,
    source_name: str = "pvacseq",
) -> pd.DataFrame:
    """Normalize pVACseq or canonical product rows without inventing assay labels."""
    out = pd.DataFrame(index=frame.index)
    for canonical, aliases in _ALIASES.items():
        source = _first(frame, aliases)
        if source is not None:
            out[canonical] = frame[source]

    # Common lossless-generation / PRIME adapter columns. Percentile ranks are
    # lower-is-better, while the product contract requires higher-is-better
    # evidence scores in [0, 1]. Derive these only when an explicit canonical
    # value was not supplied; never overwrite caller-provided evidence.
    if "dna_vaf" not in out and "tumor_vaf" in frame:
        out["dna_vaf"] = frame["tumor_vaf"]
    if "presentation_score" not in out and "mixmhcpred_rank" in frame:
        rank = pd.to_numeric(frame["mixmhcpred_rank"], errors="coerce")
        out["presentation_score"] = (1.0 - rank / 100.0).clip(0.0, 1.0)
    if "recognition_score" not in out and "prime_rank" in frame:
        rank = pd.to_numeric(frame["prime_rank"], errors="coerce")
        out["recognition_score"] = (1.0 - rank / 100.0).clip(0.0, 1.0)

    if patient_id is not None:
        out["patient_id"] = patient_id
    elif "patient_id" not in out:
        raise ValueError("patient_id is required when the input table has no patient column")

    for column in (
        "gene_symbol",
        "transcript_id",
        "mutation_id",
        "protein_variant",
        "source_variant_type",
        "mhc_class",
        "hla_allele",
        "hla_loh_call",
        "expression_call",
        "wildtype_peptide",
    ):
        if column not in out:
            out[column] = ""

    if "mutant_peptide" not in out:
        raise ValueError("input has no recognized mutant peptide column")

    out["patient_id"] = out["patient_id"].fillna("").astype(str).str.strip()
    out["mutant_peptide"] = out["mutant_peptide"].map(normalize_peptide)
    out["wildtype_peptide"] = out["wildtype_peptide"].map(normalize_peptide)
    out["hla_allele"] = out["hla_allele"].map(normalize_hla)
    for column in _NUMERIC_COLUMNS & set(out.columns):
        out[column] = pd.to_numeric(out[column], errors="coerce")

    if "candidate_id" not in out:
        base_ids = [
            _stable_candidate_id(row, str(row["patient_id"])) for _, row in out.iterrows()
        ]
        occurrences: dict[str, int] = {}
        candidate_ids = []
        for base_id in base_ids:
            occurrence = occurrences.get(base_id, 0)
            occurrences[base_id] = occurrence + 1
            candidate_ids.append(base_id if occurrence == 0 else f"{base_id}:dup{occurrence}")
        out["candidate_id"] = candidate_ids
    else:
        missing_id = out["candidate_id"].isna() | out["candidate_id"].astype(str).str.strip().eq("")
        out["candidate_id"] = out["candidate_id"].astype(str)
        out.loc[missing_id, "candidate_id"] = [
            _stable_candidate_id(row, str(row["patient_id"]))
            for _, row in out.loc[missing_id].iterrows()
        ]

    out["source_name"] = source_name
    out["candidate_schema_version"] = CANDIDATE_SCHEMA_VERSION
    report = validate_candidate_contract(out)
    if not report.ok:
        raise ValueError(f"candidate input violates product contract: {report}")
    return out.reset_index(drop=True)


def merge_rna_evidence(candidates: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Merge RNA evidence at the most specific shared identifier level."""
    normalized = evidence.copy()
    rename: dict[str, str] = {}
    for canonical, aliases in _ALIASES.items():
        source = _first(normalized, aliases)
        if source is not None and canonical not in normalized:
            rename[source] = canonical
    normalized = normalized.rename(columns=rename)

    join_priority = (
        ("patient_id", "candidate_id"),
        ("patient_id", "mutation_id", "transcript_id"),
        ("patient_id", "mutation_id"),
        ("patient_id", "transcript_id"),
        ("patient_id", "gene_symbol"),
    )
    join_columns = next(
        (keys for keys in join_priority if set(keys) <= set(candidates) and set(keys) <= set(normalized)),
        None,
    )
    if join_columns is None:
        raise ValueError("RNA evidence has no shared candidate, mutation, transcript, or gene key")

    evidence_columns = [column for column in _NUMERIC_COLUMNS if column in normalized]
    if not evidence_columns:
        raise ValueError("RNA evidence has no recognized numeric evidence columns")
    subset = normalized[list(join_columns) + evidence_columns].copy()
    if subset.duplicated(list(join_columns)).any():
        raise ValueError(f"RNA evidence is not unique at join key {join_columns}")
    for column in evidence_columns:
        subset[column] = pd.to_numeric(subset[column], errors="coerce")

    merged = candidates.merge(subset, on=list(join_columns), how="left", suffixes=("", "_rna"))
    for column in evidence_columns:
        rna_column = f"{column}_rna"
        if rna_column in merged:
            merged[column] = merged[rna_column].combine_first(merged.get(column))
            merged = merged.drop(columns=[rna_column])
    return merged


def load_product_candidates(
    path: str | Path,
    *,
    patient_id: str | None = None,
    rna_evidence_path: str | Path | None = None,
) -> pd.DataFrame:
    candidates = normalize_product_candidates(read_table(path), patient_id=patient_id)
    if rna_evidence_path is not None:
        candidates = merge_rna_evidence(candidates, read_table(rna_evidence_path))
    return candidates


def _bounded(value: Any) -> float | None:
    if value is None or pd.isna(value):
        return None
    return float(np.clip(float(value), 0.0, 1.0))


def _translated_score(row: pd.Series) -> tuple[float, bool]:
    values: list[float] = []
    expression = row.get("expression_tpm")
    if pd.notna(expression):
        values.append(float(max(expression, 0.0) / (max(expression, 0.0) + 10.0)))
    rna_vaf = _bounded(row.get("rna_vaf"))
    if rna_vaf is not None:
        values.append(math.sqrt(rna_vaf))
    reads = row.get("rna_mutant_reads")
    if pd.notna(reads):
        values.append(float(max(reads, 0.0) / (max(reads, 0.0) + 5.0)))
    return (float(np.mean(values)), True) if values else (0.5, False)


def _presentation_score(row: pd.Series) -> tuple[float, bool]:
    direct = _bounded(row.get("presentation_score"))
    if direct is not None:
        return direct, True
    percentile = row.get("binding_percentile_rank")
    if pd.notna(percentile):
        return float(np.clip(1.0 - float(percentile) / 100.0, 0.0, 1.0)), True
    affinity = row.get("binding_affinity_nm")
    if pd.notna(affinity):
        affinity = max(float(affinity), 0.0)
        return float(1.0 / (1.0 + affinity / 500.0)), True
    return 0.5, False


def _recognition_score(row: pd.Series) -> tuple[float, bool]:
    direct = _bounded(row.get("recognition_score"))
    if direct is not None:
        return direct, True
    mutant = row.get("binding_affinity_nm")
    wildtype = row.get("wildtype_binding_affinity_nm")
    if pd.notna(mutant) and pd.notna(wildtype) and float(mutant) > 0 and float(wildtype) > 0:
        delta = math.log(float(wildtype) / float(mutant))
        return float(1.0 / (1.0 + math.exp(-delta))), True
    return 0.5, False


def _coverage_score(row: pd.Series) -> tuple[float, bool]:
    clonality = _bounded(row.get("clonality_ccf"))
    if clonality is not None:
        return clonality, True
    dna_vaf = _bounded(row.get("dna_vaf"))
    if dna_vaf is not None:
        return float(np.clip(dna_vaf / 0.5, 0.0, 1.0)), True
    return 0.5, False


def _exclusion_reason(row: pd.Series) -> str:
    gate_reason = str(row.get("deterministic_gate_reason", "") or "")
    if gate_reason:
        return gate_reason
    expression = row.get("expression_tpm")
    rna_depth = row.get("rna_depth")
    reads = row.get("rna_mutant_reads")
    if pd.notna(expression) and float(expression) <= 0:
        return "NO_RNA_EXPRESSION"
    if pd.notna(rna_depth) and float(rna_depth) >= 10 and pd.notna(reads) and float(reads) <= 0:
        return "NO_MUTANT_RNA_SUPPORT"
    return ""


def score_product_candidates(frame: pd.DataFrame, config: InferenceConfig = InferenceConfig()) -> pd.DataFrame:
    report = validate_candidate_contract(frame)
    if not report.ok:
        raise ValueError(f"candidate input violates product contract: {report}")

    out = apply_deterministic_gate(frame) if config.apply_validity_gate else frame.copy()
    if not config.apply_validity_gate:
        out["deterministic_gate_reason"] = ""
        out["deterministic_gate_pass"] = True
        out["deterministic_gate_policy"] = "DISABLED"
    components = {
        "translated": _translated_score,
        "presented": _presentation_score,
        "recognized": _recognition_score,
        "coverage": _coverage_score,
    }
    weights = {"translated": 0.25, "presented": 0.35, "recognized": 0.25, "coverage": 0.15}
    available_columns: list[str] = []
    for name, scorer in components.items():
        scored = out.apply(scorer, axis=1)
        out[f"{name}_evidence_score"] = scored.map(lambda item: item[0])
        out[f"{name}_evidence_available"] = scored.map(lambda item: item[1])
        available_columns.append(f"{name}_evidence_available")

    log_score = sum(
        weight * np.log(out[f"{name}_evidence_score"].clip(lower=1e-6))
        for name, weight in weights.items()
    )
    out["epicurus_evidence_score"] = np.exp(log_score)
    out["evidence_completeness"] = out[available_columns].mean(axis=1)
    out["evidence_uncertainty"] = 1.0 - out["evidence_completeness"]
    out["epicurus_lower_evidence_score"] = out["epicurus_evidence_score"] * (
        1.0 - 0.35 * out["evidence_uncertainty"]
    )
    out["exclusion_reason"] = out.apply(_exclusion_reason, axis=1)
    out["eligible"] = out["exclusion_reason"].eq("")

    core = (
        (out["epicurus_lower_evidence_score"] >= config.core_threshold)
        & (out["evidence_uncertainty"] <= config.max_core_uncertainty)
    )
    supporting = out["epicurus_lower_evidence_score"] >= config.supporting_threshold
    out["evidence_tier"] = np.select(
        [~out["eligible"], core, supporting],
        ["EXCLUDED", "CORE", "SUPPORTING"],
        default="FILLER",
    )
    out["selection_reason"] = out.apply(_selection_reason, axis=1)
    out["ranked_schema_version"] = RANKED_SCHEMA_VERSION
    return _select_patient_portfolios(out, config)


def _selection_reason(row: pd.Series) -> str:
    if row.get("exclusion_reason"):
        return str(row["exclusion_reason"])
    available = [
        name.upper()
        for name in ("translated", "presented", "recognized", "coverage")
        if bool(row.get(f"{name}_evidence_available", False))
    ]
    missing = [
        name.upper()
        for name in ("translated", "presented", "recognized", "coverage")
        if not bool(row.get(f"{name}_evidence_available", False))
    ]
    text = f"AVAILABLE={','.join(available) or 'NONE'}"
    if missing:
        text += f"; MISSING={','.join(missing)}"
    return text


def _select_patient_portfolios(frame: pd.DataFrame, config: InferenceConfig) -> pd.DataFrame:
    out = frame.copy()
    tie_key = out["mutant_peptide"].astype(str) + "|" + out.get("hla_allele", "").astype(str)
    out["deterministic_tie_key"] = tie_key.map(lambda value: hashlib.md5(value.encode()).hexdigest())
    out["selected"] = False
    out["rank"] = pd.Series(pd.NA, index=out.index, dtype="Int64")

    tier_order = {"CORE": 0, "SUPPORTING": 1, "FILLER": 2, "EXCLUDED": 3}
    for _, patient_rows in out.groupby("patient_id", sort=True):
        ordered = patient_rows.assign(
            _tier=patient_rows["evidence_tier"].map(tier_order),
        ).sort_values(
            ["_tier", "epicurus_lower_evidence_score", "deterministic_tie_key"],
            ascending=[True, False, True],
            kind="mergesort",
        )
        selected: list[int] = []
        counts: dict[tuple[str, str], int] = {}
        for index, row in ordered.iterrows():
            if not row["eligible"] or len(selected) >= config.k:
                continue
            limits = (
                ("mutation", str(row.get("mutation_id", "")), config.max_per_mutation),
                ("gene", str(row.get("gene_symbol", "")), config.max_per_gene),
                ("hla", str(row.get("hla_allele", "")), config.max_per_hla),
            )
            if any(value and limit is not None and counts.get((kind, value), 0) >= limit for kind, value, limit in limits):
                continue
            selected.append(index)
            for kind, value, _ in limits:
                if value:
                    counts[(kind, value)] = counts.get((kind, value), 0) + 1
        out.loc[selected, "selected"] = True
        for rank, index in enumerate(selected, start=1):
            out.loc[index, "rank"] = rank

    final_report = validate_ranked_contract(out)
    if not final_report.ok:
        raise AssertionError(f"ranked output violates product contract: {final_report}")
    return out.sort_values(
        ["patient_id", "selected", "rank", "epicurus_lower_evidence_score"],
        ascending=[True, False, True, False],
        kind="mergesort",
    ).reset_index(drop=True)


def patient_summaries(scored: pd.DataFrame) -> list[PatientSummary]:
    summaries: list[PatientSummary] = []
    for patient_id, rows in scored.groupby("patient_id", sort=True):
        selected = rows[rows["selected"]]
        core_count = int((selected["evidence_tier"] == "CORE").sum())
        abstained = core_count == 0
        summaries.append(
            PatientSummary(
                patient_id=str(patient_id),
                candidate_count=len(rows),
                eligible_count=int(rows["eligible"].sum()),
                selected_count=len(selected),
                core_count=core_count,
                supporting_count=int((selected["evidence_tier"] == "SUPPORTING").sum()),
                filler_count=int((selected["evidence_tier"] == "FILLER").sum()),
                abstained=abstained,
                abstention_reason=(
                    "NO_ELIGIBLE_CANDIDATES"
                    if not rows["eligible"].any()
                    else "NO_CANDIDATE_CLEARS_CORE_EVIDENCE_POLICY"
                    if abstained
                    else None
                ),
            )
        )
    return summaries


def write_product_report(scored: pd.DataFrame, output_dir: str | Path) -> dict[str, str]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "ranked_candidates.csv"
    json_path = output / "report.json"
    markdown_path = output / "report.md"
    scored.to_csv(csv_path, index=False)
    summaries = patient_summaries(scored)
    gate_summary = summarize_gate(scored)
    payload = {
        "schema_version": RANKED_SCHEMA_VERSION,
        "policy": "deterministic_evidence_policy_v1_not_a_validated_response_probability",
        "deterministic_gate": asdict(gate_summary),
        "patients": [asdict(item) for item in summaries],
    }
    json_path.write_text(json.dumps(payload, indent=2) + "\n")

    lines = [
        "# Epicurus Neo patient report",
        "",
        "> Scores are evidence-prioritization scores, not validated probabilities of vaccine response.",
        "",
        "## Deterministic validity gate",
        "",
        f"- Input candidates: {gate_summary.input_count}",
        f"- Survivors: {gate_summary.survivor_count}",
        f"- Removed: {gate_summary.removed_count} ({100 * gate_summary.removed_fraction:.1f}%)",
        f"- Removal reasons: {gate_summary.reason_counts or {}}",
        "",
    ]
    for summary in summaries:
        lines.extend(
            [
                f"## Patient {summary.patient_id}",
                "",
                f"- Candidates: {summary.candidate_count} ({summary.eligible_count} eligible)",
                f"- Selected: {summary.selected_count}",
                f"- Evidence tiers: {summary.core_count} core, {summary.supporting_count} supporting, {summary.filler_count} filler",
                f"- Patient abstention: {'YES — ' + str(summary.abstention_reason) if summary.abstained else 'NO'}",
                "",
                "| Rank | Tier | Gene | Mutation | Peptide | HLA | Score | Lower score | Evidence |",
                "|---:|---|---|---|---|---|---:|---:|---|",
            ]
        )
        selected = scored[(scored["patient_id"] == summary.patient_id) & scored["selected"]]
        selected = selected.sort_values("rank")
        for _, row in selected.iterrows():
            lines.append(
                "| {rank} | {tier} | {gene} | {mutation} | {peptide} | {hla} | {score:.3f} | {lower:.3f} | {reason} |".format(
                    rank=int(row["rank"]),
                    tier=row["evidence_tier"],
                    gene=row.get("gene_symbol", ""),
                    mutation=row.get("mutation_id", ""),
                    peptide=row["mutant_peptide"],
                    hla=row.get("hla_allele", ""),
                    score=row["epicurus_evidence_score"],
                    lower=row["epicurus_lower_evidence_score"],
                    reason=row["selection_reason"],
                )
            )
        lines.append("")
    markdown_path.write_text("\n".join(lines) + "\n")
    return {"csv": str(csv_path), "json": str(json_path), "markdown": str(markdown_path)}


def run_product_inference(
    input_path: str | Path,
    output_dir: str | Path,
    *,
    patient_id: str | None = None,
    rna_evidence_path: str | Path | None = None,
    config: InferenceConfig = InferenceConfig(),
) -> dict[str, str]:
    candidates = load_product_candidates(
        input_path,
        patient_id=patient_id,
        rna_evidence_path=rna_evidence_path,
    )
    return write_product_report(score_product_candidates(candidates, config), output_dir)
