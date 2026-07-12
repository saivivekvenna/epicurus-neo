"""Tests for the input-only lossless variant-to-peptide generator.

Frozen exploratory protocol:
``docs/superpowers/specs/2026-07-12-osteosarc-lossless-peptide-recovery-exploratory-protocol.md``
(paired copy ``artifacts/milestone_7_decision/peptide_recovery/EXPLORATORY_PROTOCOL.md``).

These tests are pure logic, fully offline. Primary-source reference-context slices come from the
committed fixture ``tests/fixtures/peptide_recovery_ref_slices.json`` (short junction context only,
never full third-party transcripts). No assay / vaccine / recognition-label data is read here; the
generator itself is asserted to reference no such input (import/input hygiene test).
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from event_b.lossless_peptide_generation import (
    CacheMiss,
    EnsemblClient,
    MAX_LEN,
    MIN_LEN,
    POLICY_ID,
    STD_AA,
    enumerate_windows_covering,
    frameshift_novel_protein,
    frameshift_windows,
    genomic_hgvs,
    missense_windows,
    read_hla_panel,
    select_transcript,
    translate_to_stop,
    union_candidates,
    verify_transcript,
)
from epicurus_neo.evidence_router import route_candidates, select_route_aware_topk

FIXTURE = Path(__file__).parent / "fixtures" / "peptide_recovery_ref_slices.json"
REF = json.loads(FIXTURE.read_text())


# ---------------------------------------------------------------------------
# 1. Raw GRCh38 allele -> genomic HGVS (SNV, simple deletion, fail-closed)
# ---------------------------------------------------------------------------
def test_genomic_hgvs_snv_strips_chr():
    assert genomic_hgvs("chr1", 197102716, "C", "T") == "1:g.197102716C>T"
    assert genomic_hgvs("chr14", 101980529, "G", "A") == "14:g.101980529G>A"


def test_genomic_hgvs_map2_left_anchored_deletion():
    # ref GGCTA...(29) alt G: delete the 28 nt after the shared anchor, 1-based inclusive.
    hgvs = genomic_hgvs("chr2", 209694772, "GGCTACTGTGTGTTCAATAAGTACACAGT", "G")
    assert hgvs == "2:g.209694773_209694800del"


def test_genomic_hgvs_unsupported_class_fails_closed():
    with pytest.raises(ValueError):
        genomic_hgvs("chr1", 100, "A", "AT")  # insertion — not needed here, must abort


# ---------------------------------------------------------------------------
# 2. MANE/VEP transcript selection + fail-closed field verification
# ---------------------------------------------------------------------------
def _vep_payload(**tc_overrides) -> list:
    mane = {
        "transcript_id": "ENST00000367409",
        "mane_select": "NM_018136.5",
        "biotype": "protein_coding",
        "consequence_terms": ["missense_variant"],
        "protein_start": 2179,
        "amino_acids": "G/R",
        "hgvsc": "ENST00000367409.9:c.6535G>A",
        "hgvsp": "ENSP00000356379.4:p.Gly2179Arg",
        "gene_id": "ENSG00000066279",
        "gene_symbol": "ASPM",
    }
    mane.update(tc_overrides)
    decoy = {
        "transcript_id": "ENST00000000001",
        "biotype": "protein_coding",
        "consequence_terms": ["missense_variant"],
        "protein_start": 999,
        "amino_acids": "A/T",
    }
    return [{"transcript_consequences": [decoy, mane]}]


def test_select_transcript_picks_mane_with_matching_consequence():
    selected = select_transcript(_vep_payload(), expected_consequence="missense_variant")
    assert selected["transcript_id"] == "ENST00000367409"
    assert selected["mane_refseq"] == "NM_018136.5"
    assert selected["protein_start"] == 2179
    assert selected["amino_acids"] == "G/R"
    assert selected["hgvsc"] == "c.6535G>A"  # version prefix stripped


def test_select_transcript_fails_closed_when_no_mane_or_canonical():
    payload = [{"transcript_consequences": [
        {"transcript_id": "ENST00000000001", "biotype": "protein_coding",
         "consequence_terms": ["missense_variant"], "protein_start": 1, "amino_acids": "A/T"},
    ]}]
    with pytest.raises((ValueError, KeyError)):
        select_transcript(payload, expected_consequence="missense_variant")


def test_verify_transcript_aborts_on_drift():
    selected = select_transcript(_vep_payload(), expected_consequence="missense_variant")
    expected = {
        "transcript_id": "ENST00000367409",
        "mane_refseq": "NM_018136.5",
        "protein_start": 2179,
        "amino_acids": "G/R",
        "hgvsc": "c.6535G>A",
    }
    verify_transcript(selected, expected)  # matches: no raise
    drifted = {**expected, "transcript_id": "ENST99999999999"}
    with pytest.raises((ValueError, AssertionError)):
        verify_transcript(selected, drifted)


# ---------------------------------------------------------------------------
# 3. ASPM missense windows — mutant residue present, exact count 77
# ---------------------------------------------------------------------------
def test_aspm_missense_windows_count_and_contain_mutant():
    aspm = REF["aspm_missense"]
    pos_in_slice = aspm["mutant_protein_pos"] - aspm["protein_slice_start_1based"] + 1
    windows = missense_windows(
        aspm["protein_slice"], pos_in_slice, aspm["ref_aa"], aspm["mut_aa"]
    )
    # Interior residue with >=13 residues of flank each side -> sum_{L=8..14} L = 77.
    assert len(windows) == 77
    assert all(MIN_LEN <= len(w) <= MAX_LEN for w in windows)
    assert all(set(w).issubset(STD_AA) for w in windows)
    # Every window spans the substituted residue.
    mutant = aspm["mut_aa"]
    mutated_slice = (
        aspm["protein_slice"][: pos_in_slice - 1] + mutant + aspm["protein_slice"][pos_in_slice:]
    )
    for w in windows:
        start = mutated_slice.index(w)
        assert start + 1 <= pos_in_slice <= start + len(w)
    assert "RRVRVRRTLR" in windows


def test_missense_windows_fail_closed_on_wrong_reference_residue():
    aspm = REF["aspm_missense"]
    pos_in_slice = aspm["mutant_protein_pos"] - aspm["protein_slice_start_1based"] + 1
    with pytest.raises((ValueError, AssertionError)):
        missense_windows(aspm["protein_slice"], pos_in_slice, "A", aspm["mut_aa"])  # slice has G


# ---------------------------------------------------------------------------
# 4. MAP2 frameshift — junction ORF, stop handling, no WT-only windows
# ---------------------------------------------------------------------------
def test_map2_frameshift_junction_and_minimal_epitope():
    m = REF["map2_frameshift"]
    novel = frameshift_novel_protein(
        m["cds_slice"], m["cds_slice_start_1based"], m["del_c_start"], m["del_c_end"]
    )
    assert novel == m["expected_novel_junction_from_861"]
    assert novel.startswith("DSQLEDL")  # shared reference frame up to the junction
    assert m["expected_novel_frame_minimal"] in novel  # RVVPFTKAL
    # Translation stops at the first stop codon (stop symbol never enters the protein).
    assert "*" not in novel


def test_map2_frameshift_windows_all_cover_novel_frame():
    m = REF["map2_frameshift"]
    novel, windows = frameshift_windows(
        m["cds_slice"], m["cds_slice_start_1based"], m["del_c_start"], m["del_c_end"],
        protein_start=m["protein_start"],
    )
    offset_protein = (m["cds_slice_start_1based"] - 1) // 3 + 1  # 861
    novel_local_start = m["protein_start"] - offset_protein + 1  # first novel-frame local pos
    assert all(MIN_LEN <= len(w) <= MAX_LEN for w in windows)
    # No window lies entirely in the shared reference prefix: each spans a novel-frame residue.
    for w in windows:
        start = novel.index(w)  # 0-based
        window_end_1based = start + len(w)
        assert window_end_1based >= novel_local_start
    assert m["expected_novel_frame_minimal"] in windows or any(
        m["expected_novel_frame_minimal"] in w for w in windows
    )


def test_frameshift_window_count_is_novel_positions_times_lengths():
    # 37-residue C-terminal novel frame with ample standard-AA flank -> 37 * 7 = 259.
    novel_frame = "AHCHHLFKTVRIYQGRVVPFTKALMIKFEEIWPQTFH"
    assert len(novel_frame) == 37
    flank = "D" * 867  # WT flank; 'D' does not occur in the novel frame
    protein = flank + novel_frame
    novel_positions = set(range(len(flank) + 1, len(protein) + 1))  # 868..904
    windows = enumerate_windows_covering(protein, novel_positions)
    assert len(windows) == 259
    assert "DDDDDDDD" not in windows  # no purely-WT-flank window is emitted
    for w in windows:
        start = protein.index(w)
        assert start + len(w) >= len(flank) + 1  # spans >=1 novel-frame residue


# ---------------------------------------------------------------------------
# 5. DYNC1H1 input-only positive control reproduces a known pVAC MT epitope
# ---------------------------------------------------------------------------
def test_dync_positive_control_reproduces_known_pvac_epitope():
    d = REF["dync_missense_positive_control"]
    pos_in_slice = d["mutant_protein_pos"] - d["protein_slice_start_1based"] + 1
    windows = missense_windows(d["protein_slice"], pos_in_slice, d["ref_aa"], d["mut_aa"])
    assert d["known_pvac_mt_epitope"] in windows  # ATISFDTDT
    assert "KRFHATISF" in windows


# ---------------------------------------------------------------------------
# 6. translate_to_stop — standard codon table, stop excluded
# ---------------------------------------------------------------------------
def test_translate_to_stop_excludes_stop_and_trailing_partial_codon():
    # ATG AAA TAA (M K stop) then a dangling 'GG' partial codon.
    assert translate_to_stop("ATGAAATAAGG") == "MK"
    assert translate_to_stop("GACAGTCAGTGA") == "DSQ"  # DSQ then TGA stop


# ---------------------------------------------------------------------------
# 7. Ensembl client — offline cache with URL/SHA, fail-closed when uncached
# ---------------------------------------------------------------------------
def test_ensembl_client_caches_online_then_serves_offline(tmp_path):
    calls: list[str] = []

    def fake_fetch(url: str) -> bytes:
        calls.append(url)
        payload = json.dumps(_vep_payload()).encode()
        return payload

    online = EnsemblClient(cache_dir=tmp_path, offline=False, fetcher=fake_fetch)
    first = online.vep_hgvs("1:g.197102716C>T")
    assert first["json"] == _vep_payload()
    assert len(first["sha256"]) == 64
    assert first["url"].endswith("1:g.197102716C>T?mane=1;canonical=1;hgvs=1")
    assert len(calls) == 1

    # A second offline client over the same cache dir serves the cached bytes without fetching.
    offline = EnsemblClient(cache_dir=tmp_path, offline=True, fetcher=None)
    again = offline.vep_hgvs("1:g.197102716C>T")
    assert again["json"] == first["json"]
    assert again["sha256"] == first["sha256"]
    assert len(calls) == 1  # no new fetch

    # Uncached request offline fails closed.
    with pytest.raises(CacheMiss):
        offline.vep_hgvs("2:g.209694773_209694800del")


def test_ensembl_client_sequence_roundtrip_offline(tmp_path):
    def fake_fetch(url: str) -> bytes:
        return json.dumps({"id": "ENST00000367409", "seq": "MKAAA"}).encode()

    online = EnsemblClient(cache_dir=tmp_path, offline=False, fetcher=fake_fetch)
    rec = online.sequence("ENST00000367409", "protein")
    assert rec["seq"] == "MKAAA"
    offline = EnsemblClient(cache_dir=tmp_path, offline=True, fetcher=None)
    assert offline.sequence("ENST00000367409", "protein")["seq"] == "MKAAA"


# ---------------------------------------------------------------------------
# 7b. Patient HLA panel read programmatically from the pVAC table (not hardcoded)
# ---------------------------------------------------------------------------
def test_read_hla_panel_returns_sorted_unique_alleles(tmp_path):
    pvac = tmp_path / "pvac.tsv"
    pvac.write_text(
        "MT Epitope Seq\tHLA Allele\tGene Name\n"
        "AAAAAAAAA\tHLA-C*07:01\tX\n"
        "AAAAAAAAA\tHLA-A*01:01\tX\n"
        "BBBBBBBBB\tHLA-A*01:01\tY\n"
        "BBBBBBBBB\tHLA-B*27:05\tY\n"
    )
    assert read_hla_panel(pvac) == ["HLA-A*01:01", "HLA-B*27:05", "HLA-C*07:01"]


# ---------------------------------------------------------------------------
# 8. Union + frozen route-aware selection: deterministic, caps respected
# ---------------------------------------------------------------------------
def _cand(mutation_id, gene, peptide, hla, score) -> dict:
    return {
        "patient_id": "sid",
        "mutation_id": mutation_id,
        "gene_symbol": gene,
        "source_variant_type": "SNV",
        "mhc_class": "I",
        "mutant_peptide": peptide,
        "hla_allele": hla,
        "genuine_prime": score,
    }


def test_union_dedups_on_genomic_candidate_identity():
    pvac = pd.DataFrame([
        _cand("DYNC1H1-chr14-101980529", "DYNC1H1", "ATISFDTDT", "HLA-A*01:01", -0.5),
    ])
    recovered = pd.DataFrame([
        _cand("DYNC1H1-chr14-101980529", "DYNC1H1", "ATISFDTDT", "HLA-A*01:01", -0.4),  # dup
        _cand("ASPM-chr1-197102716", "ASPM", "RRVRVRRTLR", "HLA-B*27:05", -0.088),  # new
    ])
    union = union_candidates([pvac, recovered])
    assert len(union) == 2  # duplicate (peptide,HLA) collapses; new candidate kept
    dync = union[union["mutation_id"] == "DYNC1H1-chr14-101980529"]
    assert len(dync) == 1
    assert float(dync["genuine_prime"].iloc[0]) == -0.5  # pVAC representative wins (keep=first)


def _valid_peptide(seed: int) -> str:
    aas = "ACDEFGHIKLMNPQRSTVWY"
    return "".join(aas[(seed * 7 + j) % len(aas)] for j in range(9))


def test_route_aware_selection_permutation_invariant_and_capped():
    rows = [
        _cand("ASPM-chr1-197102716", "ASPM", _valid_peptide(i), "HLA-A*01:01", -0.01 * (i + 1))
        for i in range(5)
    ] + [
        _cand("MAP2-chr2-209694772", "MAP2", _valid_peptide(i + 100), "HLA-B*08:01", -0.02 * (i + 1))
        for i in range(5)
    ]
    frame = pd.DataFrame(rows)

    def selected_set(df: pd.DataFrame) -> set:
        routed = route_candidates(df.reset_index(drop=True))
        sel = select_route_aware_topk(routed, score_column="genuine_prime")
        chosen = sel[sel["route_selected"].astype(bool)]
        return set(zip(chosen["mutant_peptide"], chosen["hla_allele"]))

    forward = selected_set(frame)
    shuffled = selected_set(frame.iloc[::-1])
    assert forward == shuffled  # permutation-invariant
    # max_per_mutation cap (2) respected while both mutations are represented.
    routed = route_candidates(frame)
    sel = select_route_aware_topk(routed, score_column="genuine_prime")
    chosen = sel[sel["route_selected"].astype(bool)]
    per_mutation = chosen["mutation_id"].value_counts()
    assert per_mutation.max() <= 2
    assert set(chosen["mutation_id"]) == {"ASPM-chr1-197102716", "MAP2-chr2-209694772"}


# ---------------------------------------------------------------------------
# 9. Generator import/input hygiene — no assay/vaccine/recognition-label paths
# ---------------------------------------------------------------------------
def test_generation_module_reads_no_assay_or_label_inputs():
    import event_b.lossless_peptide_generation as mod

    source = Path(mod.__file__).read_text().lower()
    forbidden = [
        "expander", "_all_expanders", "vaccine", "elispot", "ifng", "ifngamma",
        "hudson", "minimal_epitope", "adjudicat", "osteosarc_sid", "recognized_gene",
        "may_all", "aug_all", "expander",
    ]
    hits = [token for token in forbidden if token in source]
    assert hits == [], f"generation module references forbidden assay/label inputs: {hits}"
    assert POLICY_ID == "lossless-peptide-generation-1.0.0"
