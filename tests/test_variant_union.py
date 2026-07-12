"""Tests for the generic multi-source variant/candidate union helper.

Frozen preregistration: docs/superpowers/specs/
2026-07-12-evidence-router-and-route-aware-selection-preregistration.md (§4, §9.7).
Identity is patient-scoped coordinate+allele first, patient+gene+exact protein/mutation fallback,
never gene-only.
Provenance is aggregated; representation conflicts are preserved; a unioned row that lacks a
peptide or HLA stays NEEDS_PEPTIDE_GENERATION (the helper never fabricates a peptide).
"""

import pandas as pd

from epicurus_neo.variant_union import union_variants


def test_coordinate_level_union_aggregates_provenance():
    # Two source rows for the SAME (chrom,pos,ref,alt) from two callers/timepoints must merge into
    # one row whose provenance sets are aggregated and counted.
    frame = pd.DataFrame(
        {
            "chrom": ["chr14", "chr14"],
            "pos": [101980529, 101980529],
            "ref": ["G", "G"],
            "alt": ["A", "A"],
            "gene_symbol": ["DYNC1H1", "DYNC1H1"],
            "protein_change": ["p.V314I", "p.V314I"],
            "mutant_peptide": ["GVSVEIALK", "GVSVEIALK"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:01"],
            "caller": ["DRAGEN", "Sarek"],
            "timepoint": ["May", "Aug"],
            "region": ["primary", "primary"],
            "source": ["site", "site"],
        }
    )
    unioned = union_variants(frame)
    assert len(unioned) == 1
    row = unioned.iloc[0]
    assert int(row["n_callers"]) == 2
    assert row["callers"] == "DRAGEN; Sarek"
    assert int(row["n_timepoints"]) == 2
    assert row["union_status"] == "RANKABLE"


def test_one_variant_with_multiple_peptide_hla_routes_never_collapses():
    frame = pd.DataFrame(
        {
            "patient_id": ["P1", "P1", "P1"],
            "chrom": ["chr14"] * 3,
            "pos": [101980529] * 3,
            "ref": ["G"] * 3,
            "alt": ["A"] * 3,
            "gene_symbol": ["DYNC1H1"] * 3,
            "protein_change": ["p.V314I"] * 3,
            "mutant_peptide": ["KRFHATISF", "KRFHATISF", "RFHATISF"],
            "hla_allele": ["HLA-B*27:05", "HLA-A*01:01", "HLA-B*27:05"],
            "caller": ["DRAGEN", "DRAGEN", "Sarek"],
        }
    )
    unioned = union_variants(frame)
    assert len(unioned) == 3
    assert set(zip(unioned["mutant_peptide"], unioned["hla_allele"])) == {
        ("KRFHATISF", "HLA-B*27:05"),
        ("KRFHATISF", "HLA-A*01:01"),
        ("RFHATISF", "HLA-B*27:05"),
    }


def test_distinct_coordinates_4bp_apart_stay_separate():
    # MAP2 carries two overlapping frameshift annotations 4 bp apart with different alleles and
    # different neo-frames. They are DISTINCT identity keys and must never collapse into one row.
    frame = pd.DataFrame(
        {
            "chrom": ["chr2", "chr2"],
            "pos": [209694768, 209694772],
            "ref": ["CCTGGGCTACTGTGTGTTCAATA", "GGCTACTGTGTGTTCAATAAGTACACAGT"],
            "alt": ["C", "G"],
            "gene_symbol": ["MAP2", "MAP2"],
            "protein_change": ["p.Leu867fs", "p.Gly868fs"],
            "caller": ["DRAGEN", "DRAGEN"],
        }
    )
    unioned = union_variants(frame)
    assert len(unioned) == 2
    assert set(unioned["pos"]) == {209694768, 209694772}


def test_gene_only_rows_are_not_merged():
    # No coordinates and no protein/mutation identity -> only a gene symbol is shared. A gene-only
    # union is provably wrong, so the rows must remain distinct.
    frame = pd.DataFrame(
        {
            "gene_symbol": ["TP53", "TP53"],
            "caller": ["DRAGEN", "Sarek"],
        }
    )
    unioned = union_variants(frame)
    assert len(unioned) == 2


def test_protein_change_fallback_merges_when_no_coordinates():
    # Without coordinates, exact normalized protein_change is the fallback identity key.
    frame = pd.DataFrame(
        {
            "gene_symbol": ["KRAS", "KRAS"],
            "protein_change": ["p.G12D", "p.G12D"],
            "mutant_peptide": ["KLVVVGADGV", "KLVVVGADGV"],
            "hla_allele": ["HLA-A*02:01", "HLA-A*02:01"],
            "caller": ["DRAGEN", "oncoanalyser"],
        }
    )
    unioned = union_variants(frame)
    assert len(unioned) == 1
    assert int(unioned.iloc[0]["n_callers"]) == 2


def test_same_hotspot_in_two_patients_never_merges():
    frame = pd.DataFrame(
        {
            "patient_id": ["P1", "P2"],
            "genome_build": ["GRCh38", "GRCh38"],
            "chrom": ["chr7", "chr7"],
            "pos": [140753336, 140753336],
            "ref": ["A", "A"],
            "alt": ["T", "T"],
            "gene_symbol": ["BRAF", "BRAF"],
            "protein_change": ["p.V600E", "p.V600E"],
            "caller": ["DRAGEN", "DRAGEN"],
        }
    )
    unioned = union_variants(frame)
    assert len(unioned) == 2
    assert set(unioned["patient_id"]) == {"P1", "P2"}


def test_protein_fallback_is_gene_scoped_and_builds_do_not_collapse():
    protein_rows = pd.DataFrame(
        {
            "patient_id": ["P1", "P1"],
            "gene_symbol": ["BRAF", "GENE2"],
            "protein_change": ["p.V600E", "p.V600E"],
            "caller": ["DRAGEN", "Sarek"],
        }
    )
    assert len(union_variants(protein_rows)) == 2

    coordinate_rows = pd.DataFrame(
        {
            "patient_id": ["P1", "P1"],
            "genome_build": ["GRCh37", "GRCh38"],
            "chrom": ["chr7", "chr7"],
            "pos": [140453136, 140453136],
            "ref": ["A", "A"],
            "alt": ["T", "T"],
            "gene_symbol": ["BRAF", "BRAF"],
            "caller": ["DRAGEN", "Sarek"],
        }
    )
    assert len(union_variants(coordinate_rows)) == 2


def test_representation_conflict_preserved_and_flagged():
    # Same identity key, differing annotation (two protein_change strings) -> the distinct values
    # are recorded and never silently collapsed.
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "pos": [197102716, 197102716],
            "ref": ["C", "C"],
            "alt": ["T", "T"],
            "gene_symbol": ["ASPM", "ASPM"],
            "protein_change": ["p.G2179R", "p.Gly2179Arg"],
            "caller": ["DRAGEN", "Sarek"],
        }
    )
    unioned = union_variants(frame)
    assert len(unioned) == 1
    conflicts = str(unioned.iloc[0]["representation_conflicts"])
    assert "p.G2179R" in conflicts and "p.Gly2179Arg" in conflicts


def test_missing_peptide_or_hla_stays_needs_peptide_generation():
    # A unioned candidate called by multiple pipelines but never peptide-generated must be reported
    # as an upstream generation gap, not a rankable row.
    frame = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "pos": [197102716, 197102716],
            "ref": ["C", "C"],
            "alt": ["T", "T"],
            "gene_symbol": ["ASPM", "ASPM"],
            "protein_change": ["p.G2179R", "p.G2179R"],
            "mutant_peptide": ["", ""],
            "hla_allele": ["", ""],
            "caller": ["DRAGEN", "Sarek"],
        }
    )
    unioned = union_variants(frame)
    assert len(unioned) == 1
    assert unioned.iloc[0]["union_status"] == "NEEDS_PEPTIDE_GENERATION"
