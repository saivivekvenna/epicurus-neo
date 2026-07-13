"""CD8 pMHC-multimer neoantigen cohort (2025) — an INDEPENDENT decision-problem corpus.

A second within-patient candidate universe, distinct from Müller NCI in both patients ("PBMC#NN")
and assay (pMHC multimer, not IFN-gamma ELISpot; peptide overlap with Müller is 16/7035 ≈ 0.2%).
Crucially it ships orthogonal recognition features Müller lacks — proteasomal processing and an
explicit foreignness score — so it is the sharpest available external test of whether ANY orthogonal
recognition signal beats presentation at within-patient top-k.

Source: PMC12318345 supplementary mmc2.xlsx (checksum-pinned).
    8,103 candidates / 26 patients / 34 multimer-positive (across 19 patients).

LABEL STATE (source-verified): every candidate in mmc2 WAS multimer-tested. mmc1's per-patient
"# pHLA tested" sums to the full mmc2 row count (8,103) and "# responses detected" to the YES count
(34). So Response=YES -> POSITIVE and Response=NO -> GENUINE TESTED_NEGATIVE (this cohort therefore
DOES support tested-pos-vs-tested-neg AUROC/AP and supervised negatives). This is verified at load
time by ``_TESTED_TOTAL`` reconciliation.

HONESTY: the sheet's "RF classifier score" is the paper's OWN trained recognition model (fit with
knowledge of these labels) — it is NOT a clean feature and is excluded from the benchmark.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_PATH = Path("data/raw/cd8_multimer_2025/files/mmc2.xlsx")
EXPECTED_SHA256 = "dde900d025f1fdb2b2c42efe891f67e9fe4141153fb8e386e3797242ce57b439"
EXPECTED = {"rows": 8103, "patients": 26, "positives": 34}
# Source-verified from mmc1: total pHLA multimer-tested == mmc2 row count; responses == positives.
# This is why Response=NO is a genuine TESTED_NEGATIVE rather than an unknown/untested row.
_TESTED_TOTAL = 8103
_RESPONSES_TOTAL = 34

COLUMN_MAP = {
    "Patient ID": "patient_id",
    "MUT epitope": "mutant_peptide",
    "HLA": "hla_allele",
}

# Oriented (+1 higher better, -1 lower better). %Rank scores are lower-is-better.
MULTIMER_PRESENTATION_BASELINE = {"EL (%Rank score)": -1}
MULTIMER_PRESENTATION_FEATURES = {
    "EL (%Rank score)": -1,
    "Binding affinity (%Rank score)": -1,
    "Proteasomal processing score": +1,   # orthogonal processing signal absent from Müller
}
MULTIMER_RECOGNITION_FEATURES = {
    "Foreignness score": +1,               # orthogonal foreignness signal absent from Müller
    "RNA expression (TPM)": +1,
    "Agretopicity": -1,
}
EXCLUDED_FEATURES = {"RF classifier score": "paper's own label-trained model, not a raw feature"}


def sha256_file(path: Path) -> str:
    digest = sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class MultimerCorpus:
    frame: pd.DataFrame
    reconciliation: dict
    provenance: dict


def load_cd8_multimer(
    path: Path = DEFAULT_PATH, *, strict: bool = True, allow_checksum_mismatch: bool = False
) -> MultimerCorpus:
    digest = sha256_file(path)
    if digest != EXPECTED_SHA256 and not allow_checksum_mismatch:
        raise RuntimeError(
            f"cd8_multimer checksum mismatch (got {digest}, expected {EXPECTED_SHA256}); refusing "
            "to load. Pass allow_checksum_mismatch=True to override intentionally."
        )
    raw = pd.read_excel(path)
    frame = raw.rename(columns=COLUMN_MAP).copy()
    response = frame["Response"].astype(str).str.strip().str.upper()
    unexpected = set(response.unique()) - {"YES", "NO"}
    if unexpected:
        raise ValueError(f"unexpected multimer Response values {unexpected}; verify tested status")
    frame["y"] = response.eq("YES").astype(int)
    # Source-verified: all candidates were multimer-tested (see _TESTED_TOTAL), so NO is a genuine
    # TESTED_NEGATIVE. The reconciliation below fails closed if the totals ever drift.
    tested_verified = len(frame) == _TESTED_TOTAL and int(frame["y"].sum()) == _RESPONSES_TOTAL
    frame["label"] = np.where(
        frame["y"] == 1, "POSITIVE", "TESTED_NEGATIVE" if tested_verified else "UNTESTED"
    )
    frame["mutant_peptide"] = frame["mutant_peptide"].astype(str)
    frame["hla_allele"] = frame["hla_allele"].astype(str)
    frame["candidate_id"] = (
        frame["patient_id"].astype(str) + "|" + frame["mutant_peptide"] + "|" + frame["hla_allele"]
    )

    counts = {
        "rows": int(len(frame)),
        "patients": int(frame["patient_id"].nunique()),
        "positives": int(frame["y"].sum()),
    }
    reconciles = counts == EXPECTED
    reconciliation = {
        "reconciles": reconciles,
        "observed": counts,
        "expected": dict(EXPECTED),
        "source_sha256": digest,
        "patients_with_positive": int(frame.loc[frame.y == 1, "patient_id"].nunique()),
        "all_candidates_multimer_tested": tested_verified,
        "label_state_counts": frame["label"].value_counts().to_dict(),
        "tested_negative_available": tested_verified,
    }
    if strict and not reconciles:
        raise ValueError(
            f"cd8_multimer reconciliation FAILED CLOSED: observed {counts} != expected {EXPECTED}"
        )

    provenance = {
        "source_id": "cd8_multimer_2025",
        "citation": "CD8 pMHC-multimer neoantigen cohort 2025 (PMC12318345 mmc2).",
        "path": str(path),
        "sha256": digest,
        "role": "DECISION_PROBLEM_CORPUS (independent external cohort; pMHC-multimer assay)",
        "excluded_features": EXCLUDED_FEATURES,
        "tested_negative_available": tested_verified,
        "label_semantics": (
            "All candidates multimer-tested (mmc1 reconciliation): YES=POSITIVE, NO=TESTED_NEGATIVE. "
            "Supports genuine tested-pos-vs-tested-neg AUROC/AP and supervised negatives."
        ),
        "honesty": "Independent of Müller NCI (16/7035 peptide overlap); multimer assay; RF score excluded.",
    }
    return MultimerCorpus(frame, reconciliation, provenance)
