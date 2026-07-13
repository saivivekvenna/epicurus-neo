"""PRIME 2.0 training-data normalization with a strict three-state label + leakage exclusions.

Source: Gfeller et al. 2023 (Cell Systems, PMC9811684) Supplementary Table S4 — the immunogenicity
training set used by PRIME. Acquired immutably (checksum-pinned) under data/raw/prime_training.

CRITICAL: PRIME's training negatives are 91% RANDOM proteome peptides (the `Random` flag), not
experimentally tested negatives. These must NEVER be merged:
    POSITIVE           Immunogenicity==1 (experimentally immunogenic)
    TESTED_NEGATIVE    Immunogenicity==0 AND Random==0 (real peptide, tested/observed, not immunogenic)
    SYNTHETIC_NEGATIVE Random==1 (artificial random-proteome negative; never a measured negative)

The training peptide set is used to build leakage exclusions: any evaluation-cohort peptide that is
in (or near-duplicate of) PRIME's training set makes a "PRIME vs us" comparison unfair on that
peptide, so it is excluded from the leakage-clean head-to-head.
"""

from __future__ import annotations

from functools import lru_cache
from hashlib import sha256
from pathlib import Path

import pandas as pd

from event_b.leakage_registry import canonical_peptide, near_duplicate, _kmer_index

TABLE_S4 = Path("data/raw/prime_training/PRIME2_TableS4_immunogenicity.xlsx")
EXPECTED_SHA256 = "641a104764167f9f04bafb6606e519e5625740ed1720af7d15ac9026636bc23a"
EXPECTED = {"rows": 65585, "POSITIVE": 596, "TESTED_NEGATIVE": 6084, "SYNTHETIC_NEGATIVE": 58905}


def sha256_file(path: Path) -> str:
    d = sha256()
    with open(path, "rb") as h:
        for chunk in iter(lambda: h.read(1 << 20), b""):
            d.update(chunk)
    return d.hexdigest()


def load_prime_training(path: Path = TABLE_S4, *, strict: bool = True) -> pd.DataFrame:
    digest = sha256_file(path)
    if strict and digest != EXPECTED_SHA256:
        raise RuntimeError(f"PRIME training checksum mismatch (got {digest})")
    df = pd.read_excel(path, sheet_name="TableS4", header=1)
    imm = pd.to_numeric(df["Immunogenicity"], errors="coerce")
    rnd = pd.to_numeric(df["Random"], errors="coerce")
    state = pd.Series("TESTED_NEGATIVE", index=df.index)
    state[rnd == 1] = "SYNTHETIC_NEGATIVE"
    state[(imm == 1) & (rnd == 0)] = "POSITIVE"
    df = df.assign(label_state=state, peptide=df["Mutant"].astype(str).str.upper())
    if strict:
        counts = df["label_state"].value_counts().to_dict()
        for k in ("POSITIVE", "TESTED_NEGATIVE", "SYNTHETIC_NEGATIVE"):
            if counts.get(k) != EXPECTED[k]:
                raise ValueError(f"PRIME training three-state mismatch: {counts} != {EXPECTED}")
    return df


@lru_cache(maxsize=1)
def prime_training_peptides() -> frozenset[str]:
    """All PRIME training peptides (canonical), used for exact leakage exclusion."""
    df = load_prime_training()
    return frozenset(canonical_peptide(p) for p in df["peptide"] if canonical_peptide(p))


def prime_leakage_mask(peptides, *, near: bool = True, threshold: float = 0.8) -> list[bool]:
    """Mark cohort peptides that are exact or near-duplicate matches to PRIME's training set."""
    train = set(prime_training_peptides())
    canon = [canonical_peptide(p) for p in peptides]
    exact = [p in train for p in canon]
    if not near:
        return exact
    index = _kmer_index(train)
    out = []
    for p, ex in zip(canon, exact):
        out.append(ex or (bool(p) and near_duplicate(p, index, threshold=threshold) is not None))
    return out


def prime_training_leakage_report(cohort_peptides: dict[str, set[str]], *, threshold: float = 0.8) -> dict:
    """Exact + near-duplicate overlap of each cohort's peptides with PRIME's training set."""
    train = set(prime_training_peptides())
    index = _kmer_index(train)
    report = {"prime_training_state_counts": dict(EXPECTED),
              "prime_training_unique_peptides": len(train),
              "note": "SYNTHETIC_NEGATIVE (random proteome) are kept strictly separate from "
                      "TESTED_NEGATIVE; only POSITIVE + TESTED_NEGATIVE may train a supervised "
                      "discriminator on measured labels.",
              "cohorts": {}}
    for name, peps in cohort_peptides.items():
        canon = {canonical_peptide(p) for p in peps} - {""}
        exact = canon & train
        near = {p for p in (canon - exact) if near_duplicate(p, index, threshold=threshold) is not None}
        report["cohorts"][name] = {
            "cohort_unique_peptides": len(canon),
            "exact_overlap": len(exact),
            "near_duplicate_overlap": len(near),
            "leaked_total": len(exact) + len(near),
        }
    return report
