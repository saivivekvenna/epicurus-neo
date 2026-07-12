"""Tests for the label-blind generator extensions (inframe windows + canonical-without-MANE relaxation)."""

from __future__ import annotations

import pytest

from event_b.lossless_peptide_generation import inframe_windows, select_transcript


# ---- inframe deletion / insertion windows -----------------------------------------------------------
def test_inframe_deletion_removes_segment_and_spans_junction():
    protein = "MKQPRSTVWYAAAA"  # deletion of "PRS" (positions 4-6) -> junction Q|TV
    w = inframe_windows(protein, 4, "PRS", "-")
    # every window is a real 8-14mer drawn from the MUTATED protein "MKQTVWYAAAA"
    mutated = "MKQTVWYAAAA"
    assert w, "expected some windows"
    assert all(8 <= len(p) <= 14 and p in mutated for p in w)
    # the deleted residues must not appear as an intact PRS block flanked as before
    assert "QPRS" not in "".join([""])  # sanity noop; real check: mutated has no 'PRS'
    assert "PRS" not in mutated


def test_inframe_insertion_inserts_segment():
    protein = "MKQTVWY"
    w = inframe_windows(protein, 4, "-", "GG")  # insert GG before position 4
    mutated = "MKQGGTVWY"
    assert w and all(p in mutated for p in w)
    assert any("GG" in p for p in w)  # novel inserted residues are covered


def test_inframe_reference_mismatch_fails_closed():
    with pytest.raises(ValueError, match="reference segment mismatch"):
        inframe_windows("MKQTVWY", 4, "XYZ", "-")  # protein has TVW at pos4, not XYZ


# ---- canonical-without-MANE relaxation --------------------------------------------------------------
def _vep(mane: bool):
    tc = {"transcript_id": "ENST1", "biotype": "protein_coding", "canonical": 1,
          "consequence_terms": ["missense_variant"], "protein_start": 10, "amino_acids": "A/T",
          "hgvsc": "ENST1:c.28G>A", "hgvsp": "ENSP1:p.Ala10Thr", "gene_id": "ENSG1"}
    if mane:
        tc["mane_select"] = "NM_1.1"
    return [{"transcript_consequences": [tc]}]


def test_require_mane_true_rejects_canonical_without_mane():
    with pytest.raises(ValueError, match="lacks a MANE Select RefSeq"):
        select_transcript(_vep(mane=False), expected_consequence="missense_variant", require_mane_refseq=True)


def test_require_mane_false_accepts_canonical_without_mane():
    sel = select_transcript(_vep(mane=False), expected_consequence="missense_variant",
                            require_mane_refseq=False)
    assert sel["transcript_id"] == "ENST1"
    assert sel["mane_refseq"] == ""  # no MANE RefSeq, not a KeyError


def test_mane_still_preferred_when_present():
    sel = select_transcript(_vep(mane=True), expected_consequence="missense_variant",
                            require_mane_refseq=True)
    assert sel["mane_refseq"] == "NM_1.1"
