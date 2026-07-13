"""Gartner NCI candidate universe WITH genuine per-peptide ascertainment (three-state labels).

This is the SAME NCI TIL cohort that Müller harmonized — but the pre-harmonization Gartner Nmers
export keeps the `Screening Status` column the "min" file dropped, so it distinguishes:
    screened-positive  (Screening Status '1' / 'CD8')  -> POSITIVE
    screened-negative  (Screening Status '0' / '-')    -> TESTED_NEGATIVE   (genuine)
    unscreened                                          -> UNTESTED
That makes it the corrected substrate for the metrics the Müller min file cannot support:
non-circular AUROC/AP and supervised training on genuinely tested positives vs tested negatives.

Testing set: 26 patients / 8,782 mutation-level candidates (25mers), with NetMHCpan EL/BA,
MixMHCpred, MHCflurry, HLAthena %rank predictors + expression/VAF deciles. ~57% are UNTESTED.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_PATH = Path("data/raw/gartner_nci/NmersTestingSet.txt")
EXPECTED_SHA256 = "3630d9209fb26d86e8e1f03dbce0a5811c15ae08a9fb83f3581f1abd1a0aeca6"
EXPECTED = {"rows": 8782, "patients": 26, "positives": 46, "tested_negatives": 3722, "untested": 5014}

COLUMN_MAP = {
    "ID": "patient_id",
    "Mut Epitope": "mutant_peptide",
    "Top netMHCpan4.0 EL ranked minimal": "netmhcpan_el_rank",
    "Top netMHCpan4.0 BA ranked minimal": "netmhcpan_ba_rank",
    "Top mixMHCPred ranked minimal": "mixmhcpred_rank",
    "Top MHCflurry ranked minimal": "mhcflurry_rank",
    "Top HLAthena ranked minimal": "hlathena_rank",
    "Gene Expression Decile for this sample(1=lowest expression-10=highest expression)": "expr_decile",
    "Exome VAF Decile": "vaf_decile",
}

# %rank predictors are lower-is-better; deciles are higher-is-better. Orientation is fixed domain
# knowledge, not tuned. Presentation baseline is NetMHCpan EL %rank (oriented).
PRESENTATION_BASELINE = {"netmhcpan_el_rank": -1}
PRESENTATION_FEATURES = {
    "netmhcpan_el_rank": -1, "netmhcpan_ba_rank": -1, "mixmhcpred_rank": -1,
    "mhcflurry_rank": -1, "hlathena_rank": -1,
}
# Gartner Testing carries only expression/VAF as orthogonal (recognition-adjacent) features; it has
# no agretopicity/foreignness. Its value is the genuine tested-negative denominator, not new features.
RECOGNITION_FEATURES = {"expr_decile": +1, "vaf_decile": +1}

POSITIVE_STATUS = {"1", "CD8"}
TESTED_NEGATIVE_STATUS = {"0", "-"}
UNTESTED_STATUS = {"unscreened"}


def sha256_file(path: Path) -> str:
    digest = sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class GartnerCorpus:
    frame: pd.DataFrame
    reconciliation: dict
    provenance: dict


def _label_from_status(status: pd.Series) -> pd.Series:
    s = status.astype(str).str.strip()
    label = pd.Series("UNTESTED", index=s.index)
    label[s.isin(POSITIVE_STATUS)] = "POSITIVE"
    label[s.isin(TESTED_NEGATIVE_STATUS)] = "TESTED_NEGATIVE"
    unknown = ~s.isin(POSITIVE_STATUS | TESTED_NEGATIVE_STATUS | UNTESTED_STATUS)
    if unknown.any():
        raise ValueError(f"unexpected Screening Status values: {sorted(s[unknown].unique())}")
    return label


def load_gartner_nci(
    path: Path = DEFAULT_PATH, *, strict: bool = True, allow_checksum_mismatch: bool = False
) -> GartnerCorpus:
    digest = sha256_file(path)
    if digest != EXPECTED_SHA256 and not allow_checksum_mismatch:
        raise RuntimeError(
            f"gartner NCI checksum mismatch (got {digest}, expected {EXPECTED_SHA256}); refusing "
            "to load. Pass allow_checksum_mismatch=True to override intentionally."
        )
    raw = pd.read_csv(path, sep="\t", dtype=str)
    frame = raw.rename(columns=COLUMN_MAP).copy()
    frame["label"] = _label_from_status(raw["Screening Status"])
    frame["y"] = (frame["label"] == "POSITIVE").astype(int)
    for col in [*PRESENTATION_FEATURES, *RECOGNITION_FEATURES]:
        frame[col] = pd.to_numeric(frame[col], errors="coerce")
    frame["mutant_peptide"] = frame["mutant_peptide"].astype(str)
    frame["hla_allele"] = "NA"  # Gartner Nmers is mutation-level (best minimal); no per-allele row
    frame["candidate_id"] = frame["patient_id"].astype(str) + "|" + frame["mutant_peptide"]
    # A handful of exact (patient, 25mer) duplicates: keep the most-informative screening state.
    priority = {"POSITIVE": 0, "TESTED_NEGATIVE": 1, "UNTESTED": 2}
    frame = (
        frame.assign(_p=frame["label"].map(priority))
        .sort_values(["candidate_id", "_p"])
        .drop_duplicates("candidate_id", keep="first")
        .drop(columns="_p")
        .reset_index(drop=True)
    )

    counts = {
        "rows": int(len(raw)),
        "patients": int(frame["patient_id"].nunique()),
        "positives": int((frame["label"] == "POSITIVE").sum()),
        "tested_negatives": int((frame["label"] == "TESTED_NEGATIVE").sum()),
        "untested": int((frame["label"] == "UNTESTED").sum()),
    }
    # Reconcile pre-dedup screening-state counts against the pinned three-state totals.
    reconciles = all(
        _label_from_status(raw["Screening Status"]).value_counts().get(k, 0) == v
        for k, v in {"POSITIVE": EXPECTED["positives"], "TESTED_NEGATIVE": EXPECTED["tested_negatives"],
                     "UNTESTED": EXPECTED["untested"]}.items()
    ) and counts["rows"] == EXPECTED["rows"]
    reconciliation = {
        "reconciles": reconciles,
        "observed": counts,
        "expected": dict(EXPECTED),
        "source_sha256": digest,
        "label_state_counts": frame["label"].value_counts().to_dict(),
        "tested_negative_available": True,
        "patients_with_positive": int(frame.loc[frame.y == 1, "patient_id"].nunique()),
    }
    if strict and not reconciles:
        raise ValueError(
            f"gartner NCI reconciliation FAILED CLOSED: {counts} != pinned {EXPECTED}"
        )

    provenance = {
        "source_id": "gartner_nci_testing",
        "citation": "Gartner et al. NCI TIL neoantigen screening (Nmers testing export).",
        "path": str(path),
        "sha256": digest,
        "role": "DECISION_PROBLEM_CORPUS with genuine three-state ascertainment (Screening Status)",
        "tested_negative_available": True,
        "label_semantics": (
            "Screening Status 1/CD8 -> POSITIVE; 0/- -> TESTED_NEGATIVE; unscreened -> UNTESTED. "
            "Supports full-universe positive-unlabeled retrieval AND non-circular tested-subset "
            "AUROC/AP + supervised negatives (tested-only)."
        ),
        "honesty": (
            "Same NCI cohort as Müller; ~57% of candidates were never screened (UNTESTED). "
            "Mutation-level (25mer) candidates; NetMHCpan/MixMHCpred/MHCflurry/HLAthena %rank + "
            "expression/VAF deciles. NOT a PRIME head-to-head (no genuine PRIME scores)."
        ),
    }
    return GartnerCorpus(frame, reconciliation, provenance)
