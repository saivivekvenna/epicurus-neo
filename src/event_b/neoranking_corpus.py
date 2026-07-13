"""Müller NCI neoantigen-ranking corpus — a real within-patient DECISION-PROBLEM substrate.

Unlike the Zhao/CEDAR recognition *label* corpora (tested-subset discrimination), this is a
complete per-patient candidate universe: every predicted neo-peptide for each patient, with a
validated-immunogenicity label. That is exactly what the north star needs to measure within-patient
top-k retrieval — "did the ranker put the experimentally recognized neoantigen in the patient's
top 20 out of their whole candidate list?" — rather than only ranking an already-tested subset.

Source: Müller et al. 2023 harmonized NCI training split (the "neoranking" benchmark), mirrored
locally at data/raw/neoranking_mirror/MullerNCItrain.train-data.min.tsv (checksum-pinned).
    292,495 candidates / 56 patients / 82 validated positives (43 patients with >=1).

HONESTY / SCOPE (corrected after the M7 label/ascertainment audit):
    * These are NCI TIL neo-peptides that overlap IEDB/CEDAR — i.e. PRIME's own training universe.
      There are NO genuine PRIME scores here, so this is NOT a clean "beat PRIME" head-to-head.
    * LABEL STATE: VALIDATED=1 is POSITIVE (immunogenic). VALIDATED=0 is **UNTESTED**, NOT
      tested-negative: the harmonized "min" file carries NO per-peptide assay indicator, and the
      same NCI cohort's Gartner raw (`Screening Status`) shows ~57-69% of candidates were never
      screened. So this corpus supports full-universe positive-UNLABELED top-k retrieval, but NOT
      a tested-pos-vs-tested-neg AUROC/AP or supervised negative training (those are BLOCKED here;
      use the Gartner three-state or multimer corpus for genuine tested negatives).
    * ln_NumTested is a PER-PATIENT covariate: ln of the patient's candidate-universe size (56
      unique values = patient count; exp(ln_NumTested) ~ candidate count, corr 0.9999). It is NOT
      "how many times a peptide was assayed" and is never a model feature.
"""

from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_PATH = Path("data/raw/neoranking_mirror/MullerNCItrain.train-data.min.tsv")
EXPECTED_SHA256 = "d0ea2457513d0178fb323c593168265d0f040ec7092bc9e55bd0d25a78cf3dc7"
EXPECTED = {"rows": 292495, "patients": 56, "positives": 82}

# Raw source column -> canonical name.
COLUMN_MAP = {
    "PatientID": "patient_id",
    "HLA_type_x": "hla_allele",
    "MT_pep_x": "mutant_peptide",
    "VALIDATED": "y",
}

# Feature orientation: +1 higher-is-better, -1 lower-is-better. Domain knowledge (EL/stability/
# expression higher is better; nM affinity and agretopicity ratio lower is better), each confirmed
# against the global label (see progress.md); orientation is a fixed property of the score, NOT a
# tuned hyperparameter.
PRESENTATION_FEATURES = {"Score_EL": +1, "MT_BindAff": -1, "BindStab": +1}
RECOGNITION_FEATURES = {"Agretopicity": -1, "Quantification": +1}
EXCLUDED_FEATURES = {
    "ln_NumTested": "PER-PATIENT covariate = ln(patient candidate-universe size); not a per-peptide "
    "assay count and never a model feature."
}

# The presentation-only incumbent (what a presentation-anchored model like PRIME leans on).
PRESENTATION_BASELINE = "Score_EL"  # higher is better


def sha256_file(path: Path) -> str:
    digest = sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class NeorankingCorpus:
    frame: pd.DataFrame          # canonical candidate frame (one row per candidate)
    reconciliation: dict         # pinned-count fail-closed record
    provenance: dict


def _oriented(values: pd.Series, sign: int) -> np.ndarray:
    """Return the feature so that HIGHER is always better, robust to skew/heavy tails.

    nM affinity spans orders of magnitude, so it is log-compressed; agretopicity is a positive
    ratio and is likewise log-compressed; bounded scores (EL, stability, expression) are used as-is.
    """
    v = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    return float(sign) * v


def load_neoranking_nci(
    path: Path = DEFAULT_PATH, *, strict: bool = True, allow_checksum_mismatch: bool = False
) -> NeorankingCorpus:
    """Load, checksum-verify, and fail-closed reconcile the Müller NCI decision-problem corpus."""
    digest = sha256_file(path)
    if digest != EXPECTED_SHA256 and not allow_checksum_mismatch:
        raise RuntimeError(
            f"neoranking NCI checksum mismatch (got {digest}, expected {EXPECTED_SHA256}); "
            "refusing to load. Pass allow_checksum_mismatch=True to override intentionally."
        )
    raw = pd.read_csv(path, sep="\t")
    frame = raw.rename(columns=COLUMN_MAP).copy()
    frame["y"] = frame["y"].astype(int)
    # Three-state, source-honest: VALIDATED=1 is POSITIVE; VALIDATED=0 is UNTESTED (the min file
    # has no per-peptide assay indicator). NEVER coerce these zeros to TESTED_NEGATIVE.
    frame["label"] = np.where(frame["y"] == 1, "POSITIVE", "UNTESTED")
    frame["candidate_id"] = (
        frame["patient_id"].astype(str) + "|" + frame["mutant_peptide"].astype(str)
        + "|" + frame["hla_allele"].astype(str)
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
        "label_state_counts": frame["label"].value_counts().to_dict(),
        "tested_negative_available": False,
    }
    if strict and not reconciles:
        raise ValueError(
            f"neoranking NCI reconciliation FAILED CLOSED: observed {counts} != expected {EXPECTED}"
        )

    provenance = {
        "source_id": "neoranking_nci_muller",
        "citation": "Müller et al., harmonized NCI neoantigen-ranking training split (neoranking).",
        "path": str(path),
        "sha256": digest,
        "role": "DECISION_PROBLEM_CORPUS (complete within-patient candidate universes)",
        "presentation_baseline": PRESENTATION_BASELINE,
        "excluded_features": EXCLUDED_FEATURES,
        "tested_negative_available": False,
        "label_semantics": (
            "POSITIVE = VALIDATED 1; VALIDATED 0 = UNTESTED (no per-peptide assay indicator). "
            "Full-universe top-k is positive-UNLABELED retrieval; AUROC/AP/supervised-negative "
            "arms are BLOCKED (no tested negatives)."
        ),
        "honesty": (
            "NCI TIL peptides overlap PRIME/IEDB training; no genuine PRIME score present. Not a "
            "PRIME head-to-head. Global AUROC over the full universe is positives-vs-unlabeled and "
            "ascertainment-biased; use Gartner three-state / multimer for tested-negative metrics."
        ),
    }
    return NeorankingCorpus(frame, reconciliation, provenance)


def oriented_feature_matrix(frame: pd.DataFrame, feature_signs: dict[str, int]) -> np.ndarray:
    """Stack the requested features, each oriented so higher is better (raw scale; scaling is the
    model's fold-local job). Missing values are left as NaN for the caller to impute in-fold."""
    cols = [_oriented(frame[name], sign) for name, sign in feature_signs.items()]
    return np.vstack(cols).T if cols else np.empty((len(frame), 0))


def presentation_ceiling_analysis(frame: pd.DataFrame, baseline: str = PRESENTATION_BASELINE) -> dict:
    """Localize the north-star ceiling: where do validated positives fall in each patient's
    presentation ranking, and are the presentation-MISSED ones distinguishable from background?

    This is training-free (pure ranking + describe). It answers the decisive question: is a
    presentation-anchored ranker (like PRIME) structurally leaving recoverable immunogenic
    peptides on the table, and does any available feature distinguish them?
    """
    f = frame.copy()
    f["pres_rank"] = f.groupby("patient_id")[baseline].rank(ascending=False, method="first")
    pos = f[f["y"] == 1]
    neg = f[f["y"] == 0]
    bands = {"rank_1_5": (1, 5), "rank_6_20": (6, 20), "rank_21_100": (21, 100),
             "rank_over_100": (101, 10**12)}
    band_counts = {name: int(pos["pres_rank"].between(lo, hi).sum()) for name, (lo, hi) in bands.items()}
    missed = pos[pos["pres_rank"] > 100]
    caught = pos[pos["pres_rank"] <= 20]
    feat_cols = ["Score_EL", "MT_BindAff", "BindStab", "Agretopicity", "Quantification"]
    return {
        "total_positives": int(len(pos)),
        "positives_by_presentation_rank": band_counts,
        "top20_caught": int(pos["pres_rank"].le(20).sum()),
        "unreachable_over_rank_100": int(len(missed)),
        "caught_medians": {c: float(caught[c].median()) for c in feat_cols} if len(caught) else {},
        "missed_medians": {c: float(missed[c].median()) for c in feat_cols} if len(missed) else {},
        "background_unlabeled_medians": {
            "Quantification": float(neg["Quantification"].median()),
            "Agretopicity": float(neg["Agretopicity"].median()),
        },
        "label_state_note": (
            "The comparison background here is UNLABELED (VALIDATED=0 = untested), not tested-"
            "negative. This is a positive-rank / positive-unlabeled diagnostic, not a discrimination "
            "metric against confirmed negatives."
        ),
        "interpretation": (
            "Presentation-missed positives are weak binders (high nM) but carry elevated expression "
            "and foreignness vs the unlabeled background. This localizes where positives fall in the "
            "presentation ranking; it does NOT establish that presentation is 'near-solved', which "
            "would require tested negatives (see Gartner three-state / multimer)."
        ),
    }


def shared_peptide_diagnostics(frame: pd.DataFrame) -> dict:
    """Report peptides shared across patients (a secondary leakage consideration for patient-CV:
    an identical hotspot peptide can appear for >1 patient)."""
    per_pep_patients = frame.groupby("mutant_peptide")["patient_id"].nunique()
    pos = frame[frame.y == 1]
    pos_pep_patients = pos.groupby("mutant_peptide")["patient_id"].nunique()
    return {
        "unique_peptides": int(frame["mutant_peptide"].nunique()),
        "peptides_in_multiple_patients": int((per_pep_patients > 1).sum()),
        "positive_peptides": int(pos["mutant_peptide"].nunique()),
        "positive_peptides_in_multiple_patients": int((pos_pep_patients > 1).sum()),
    }
