"""Generic multi-source variant/candidate union helper.

Frozen preregistration: ``docs/superpowers/specs/
2026-07-12-evidence-router-and-route-aware-selection-preregistration.md`` (§4).

This helper merges candidate/variant rows drawn from multiple callers, timepoints, regions, or
sources into one deduplicated frame **without ever fabricating a peptide or collapsing distinct
biology**. It is deliberately upstream of ``normalize_product_candidates`` (rows may legitimately
carry an empty peptide/HLA and would not satisfy the candidate contract), and it feeds:

* the evidence router (aggregated provenance -> the multi-source / single-caller flags), and
* the candidate-generation-reachability funnel (the multi-caller raw union is the recall denominator).

Identity rules (frozen, priority order, never gene-only), all scoped to ``patient_id`` (and to
``genome_build`` for coordinates) when those columns are present — the same hotspot in two patients,
or the same chrom/pos under two genome builds, is never merged:

1. ``(chrom, pos, ref, alt)`` when all four are present and non-empty for a row.
2. else exact normalized ``protein_change`` / ``mutation_id``, scoped by gene symbol.
3. else the row is kept **distinct** (a gene-only merge is provably wrong: MAP2 and DYNC1H1 each
   carry two distinct coordinates, and MAP2's two frameshift coordinates 4 bp apart are different
   keys that must stay separate rows).

The ``mutant_peptide`` and ``hla_allele`` extend that variant identity: one variant can generate
several distinct peptide x HLA candidates, and those stay separate rows (never collapsed).
"""

from __future__ import annotations

from typing import Iterable

import pandas as pd


# Columns whose per-row values are aggregated into sorted-unique provenance sets.
_PROVENANCE = {
    "caller": ("callers", "n_callers"),
    "timepoint": ("timepoints", "n_timepoints"),
    "region": ("regions", "n_regions"),
    "source": ("sources", "n_sources"),
}

# Annotation columns whose differing values within one identity key are recorded (never collapsed).
_ANNOTATION_COLUMNS = ("gene_symbol", "protein_change", "mutation_id", "consequence")

# Columns copied through from a representative (first non-empty) row of each identity group. Their
# original values are preserved (dtype intact) rather than the cleaned string form.
_REPRESENTATIVE = (
    "patient_id",
    "genome_build",
    "chrom",
    "pos",
    "ref",
    "alt",
    "gene_symbol",
    "protein_change",
    "mutation_id",
    "consequence",
    "source_variant_type",
    "mhc_class",
)


def _clean(value: object) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    text = str(value).strip()
    return "" if text.lower() in {"", "nan", "none"} else text


def _identity_key(row: pd.Series, index: object) -> tuple:
    # Peptide+HLA extend the variant identity: one variant can generate several distinct
    # peptide x HLA candidates, and those are different rows that must never collapse together.
    patient = _clean(row.get("patient_id"))
    build = _clean(row.get("genome_build"))
    peptide = _clean(row.get("mutant_peptide")).upper()
    hla = _clean(row.get("hla_allele")).upper()
    chrom, pos, ref, alt = (_clean(row.get(c)) for c in ("chrom", "pos", "ref", "alt"))
    if chrom and pos and ref and alt:
        return ("coord", patient, build, chrom, pos, ref, alt, peptide, hla)
    for column in ("protein_change", "mutation_id"):
        value = _clean(row.get(column))
        if value:
            gene = _clean(row.get("gene_symbol"))
            return ("annot", patient, gene.upper(), value.upper(), peptide, hla)
    # Gene-only (or fully unkeyed) rows are never merged -> unique per row.
    return ("row", index)


def _first_nonempty(series: pd.Series) -> str:
    for value in series:
        cleaned = _clean(value)
        if cleaned:
            return cleaned
    return ""


def _first_original(series: pd.Series) -> object:
    """Return the original (dtype-preserved) value of the first non-empty row, else ``""``."""
    for value in series:
        if _clean(value):
            return value
    return ""


def union_variants(frames: pd.DataFrame | Iterable[pd.DataFrame]) -> pd.DataFrame:
    """Merge multi-source variant/candidate rows on the frozen identity rules.

    Accepts one frame or an iterable of frames (concatenated first). Returns one row per identity
    key with aggregated provenance sets/counts, preserved representation conflicts, and a
    ``union_status`` of ``RANKABLE`` or ``NEEDS_PEPTIDE_GENERATION`` (never a fabricated peptide).
    """
    if isinstance(frames, pd.DataFrame):
        combined = frames.copy()
    else:
        parts = [f for f in frames if f is not None and len(f)]
        combined = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()
    if combined.empty:
        return combined.reset_index(drop=True)

    combined = combined.reset_index(drop=True)
    keys = [_identity_key(row, index) for index, row in combined.iterrows()]
    combined = combined.assign(_union_key=keys)

    records: list[dict] = []
    for _, group in combined.groupby("_union_key", sort=False):
        record: dict = {}
        for column in _REPRESENTATIVE:
            if column in group:
                record[column] = _first_original(group[column])

        # Peptide / HLA are copied through only if a source actually carried them; never fabricated.
        peptide = _first_nonempty(group["mutant_peptide"]) if "mutant_peptide" in group else ""
        hla = _first_nonempty(group["hla_allele"]) if "hla_allele" in group else ""
        record["mutant_peptide"] = peptide
        record["hla_allele"] = hla

        for raw, (set_name, count_name) in _PROVENANCE.items():
            if raw in group:
                values = sorted({v for v in (_clean(x) for x in group[raw]) if v})
            elif set_name in group:  # already-aggregated "; "-joined provenance
                values = sorted(
                    {v.strip() for cell in group[set_name] for v in str(cell).split(";") if v.strip()}
                )
            else:
                values = []
            record[set_name] = "; ".join(values)
            record[count_name] = len(values)

        conflicts: dict[str, list[str]] = {}
        for column in _ANNOTATION_COLUMNS:
            if column not in group:
                continue
            distinct = sorted({v for v in (_clean(x) for x in group[column]) if v})
            if len(distinct) > 1:
                conflicts[column] = distinct
        record["representation_conflicts"] = (
            "; ".join(f"{col}={'|'.join(vals)}" for col, vals in conflicts.items()) if conflicts else ""
        )

        record["n_sources_rows"] = len(group)
        record["union_status"] = "RANKABLE" if peptide and hla else "NEEDS_PEPTIDE_GENERATION"
        records.append(record)

    return pd.DataFrame(records).reset_index(drop=True)
