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


def test_full_pool_mask_leaked_positive_does_not_count_but_still_competes():
    # one patient: rank1 positive is LEAKED (masked out of hits); a clean positive sits at rank 2.
    df = pd.DataFrame({
        "patient_id": ["p"] * 3,
        "prime": [0.1, 0.2, 9.0], "rna_af": [0.0, 0.0, 0.0],
        "label": ["POSITIVE", "POSITIVE", "TESTED_NEGATIVE"]})
    pp = mod._pct(df, "prime", False)
    clean = np.array([False, True, True])  # rank-1 positive is leaked
    hits = mod.config_hits(df, pp, np.zeros(3), 0.0, 0, clean)
    # the leaked rank-1 positive still occupies a slot (pool not shrunk) but contributes 0; clean rank-2
    # positive counts -> exactly 1 clean hit
    assert hits["p"] == 1.0


def test_null_included_and_tie_break_prefers_null():
    # equal hits -> tiebreak_key must rank the null (alpha=0,q=0) strictly before any acting config
    kn = mod.tiebreak_key(("core_deployable", 1.0, 0.0, 0), 5)
    ka = mod.tiebreak_key(("core_deployable", 1.0, 0.10, 1), 5)
    kq = mod.tiebreak_key(("core_deployable", 1.0, 0.0, 1), 5)
    assert kn < ka and kn < kq                 # null preferred over acting configs at equal hits
    # lower q preferred, then lower alpha, then simpler (core) arm
    assert mod.tiebreak_key(("core_deployable", 1.0, 0.10, 1), 5) < mod.tiebreak_key(("core_deployable", 1.0, 0.10, 2), 5)
    assert mod.tiebreak_key(("core_deployable", 1.0, 0.10, 1), 5) < mod.tiebreak_key(("improve_rich_partial_bridge", 1.0, 0.10, 1), 5)


def test_inner_leaked_mask_independent_of_outer_test():
    # a peptide shared ONLY between an inner-validation partition and the OUTER-TEST partition must NOT be
    # flagged leaked when the inner mask is recomputed on the outer-train subset (outer-test excluded).
    df = pd.DataFrame({
        "patient_id": ["a", "b", "c"], "partition": [0, 1, 2],
        "mut_peptide": ["SHAREDPEP", "OTHERPEPT", "SHAREDPEP"],  # shared between p0 (inner) and p2 (outer-test)
        "label": ["POSITIVE", "TESTED_NEGATIVE", "POSITIVE"],
        "prime": [0.1, 0.2, 0.3], "rna_af": [0.0, 0.0, 0.0]})
    outer_test = 2
    inner = df[df["partition"] != outer_test].reset_index(drop=True)  # partitions 0,1
    m = mod.leaked_mask(inner)  # recomputed within outer-train only
    # 'SHAREDPEP' in inner appears once (p0); its only cross-partition match is the OUTER-TEST copy, which
    # is absent from `inner` -> so within inner it is NOT leaked.
    assert not m.any()
