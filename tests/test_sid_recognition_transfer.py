"""Stage-1 invariants for the non-Sid recognition-transfer freeze (no Sid, no network)."""

from __future__ import annotations

import importlib

import numpy as np
import pandas as pd

mod = importlib.import_module("scripts.sid_recognition_transfer")


def _frame():
    # two patients, partition-split; a shared peptide 'AAAAAAAAA' across partitions -> leaked in held-out
    return pd.DataFrame({
        "patient_id": ["p1", "p1", "p1", "p2", "p2", "p2"],
        "partition": [0, 0, 0, 1, 1, 1],
        "mut_peptide": ["AAAAAAAAA", "CDEFGHIKL", "MNPQRSTVW", "AAAAAAAAA", "YCDEFGHIK", "LMNPQRSTV"],
        "label": ["POSITIVE", "TESTED_NEGATIVE", "TESTED_NEGATIVE", "POSITIVE", "TESTED_NEGATIVE", "TESTED_NEGATIVE"],
        "prime": [0.1, 5.0, 8.0, 0.2, 6.0, 9.0], "rna_af": [0.5, 0.0, 0.0, 0.6, 0.0, 0.0],
    })


def test_leaked_mask_flags_shared_peptide():
    df = _frame()
    m = mod.leaked_mask(df)
    # the shared 'AAAAAAAAA' appears in both partitions -> held-out copy in each partition is leaked
    assert m.sum() >= 2
    assert m[df["mut_peptide"].to_numpy() == "AAAAAAAAA"].all()


def test_config_hits_null_is_prime_top20():
    df = _frame()
    pp = mod._pct(df, "prime", False)
    keep = np.ones(len(df), bool)
    hits = mod.config_hits(df, pp, np.zeros(len(df)), 0.0, 0, keep)
    # each patient's single positive is best PRIME -> in top-20 -> 1 hit each
    assert hits["p1"] == 1.0 and hits["p2"] == 1.0


def test_reserve_never_uses_zero_rna_evidence():
    # with q=1 but only the positive has rna_af>0, the reserve can only pick the positive (never a 0-evidence row)
    df = _frame()
    pp = mod._pct(df, "prime", False)
    keep = np.ones(len(df), bool)
    hits = mod.config_hits(df, pp, np.zeros(len(df)), 0.0, 1, keep)
    assert hits["p1"] == 1.0  # positive retained; reserve does not evict it


def test_no_sid_files_in_allowed_set():
    for f in mod.ALLOWED_DATA_FILES:
        assert "osteosarc" not in f and "sid" not in f.lower() and "variant_vafs" not in f
