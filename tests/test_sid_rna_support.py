"""Tests for the exploratory Sid RNA-support gate (src/event_b/sid_rna_support.py).

Locks the leakage-safe invariants: absence is never a veto; only positive non-expression evidence removes;
recognized positives are never removed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from event_b.sid_rna_support import (
    confidently_unexpressed,
    load_tumor_rna_support,
    rna_support_gate,
)


def _rna():
    return pd.DataFrame({
        "variant_id": ["EXPRESSED-1", "UNEXPRESSED-1"],
        "rna_max_alt_reads": [10.0, 0.0], "rna_max_vaf": [0.4, 0.0], "rna_assay_present": [True, True],
    })


def test_absence_is_never_a_veto():
    # a mutation with NO RNA row must be KEPT even if TPM is 0
    muts = pd.Series(["NO-RNA-ROW"])
    tpm = pd.Series([0.0])
    flag = confidently_unexpressed(muts, tpm, _rna())
    assert flag[0] == False  # noqa: E712 — absent evidence -> not unexpressed -> keep


def test_positive_rna_evidence_keeps():
    # mutant reads present -> KEEP even if TPM==0
    flag = confidently_unexpressed(pd.Series(["EXPRESSED-1"]), pd.Series([0.0]), _rna())
    assert flag[0] == False  # noqa: E712


def test_confidently_unexpressed_removes():
    # RNA row exists, zero mutant reads, TPM 0 -> unexpressed -> remove
    flag = confidently_unexpressed(pd.Series(["UNEXPRESSED-1"]), pd.Series([0.0]), _rna())
    assert flag[0] == True  # noqa: E712


def test_expression_present_keeps_even_if_no_mutant_reads():
    # zero mutant reads but TPM>0 -> ambiguous -> KEEP (do not veto on TPM alone)
    flag = confidently_unexpressed(pd.Series(["UNEXPRESSED-1"]), pd.Series([3.5]), _rna())
    assert flag[0] == False  # noqa: E712


def test_gate_adds_keep_column():
    cand = pd.DataFrame({"mutation_id": ["EXPRESSED-1", "UNEXPRESSED-1"], "expression_tpm": [12.0, 0.0]})
    out = rna_support_gate(cand, rna=_rna())
    assert list(out["rna_gate_keep"]) == [True, False]


def test_real_sid_positives_are_never_unexpressed():
    """On the real Sid RNA table, none of the 3 recognized positives are flagged unexpressed."""
    rna = load_tumor_rna_support()
    pos = pd.Series(["ASPM-chr1-197102716", "DYNC1H1-chr14-101980529", "MAP2-chr2-209694772"])
    # recognized positives are transcribed (they carry mutant RNA reads) -> never removed, regardless of TPM
    tpm = pd.Series([16.49, 100.0, 5.2])
    flag = confidently_unexpressed(pos, tpm, rna)
    assert not flag.any()
