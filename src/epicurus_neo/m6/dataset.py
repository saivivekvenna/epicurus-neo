from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from event_b.models import BiologicalEvent, ResponseLabel

CORPUS_DIR = Path("outputs/event_b_backbone/combined")
CANDIDATE_RESOLVED_STUDIES = (
    "braun_rcc_2025",
    "hu_neovax_2021",
    "mkras_vax_2026",
    "pdac_neovax_2023",
)
# ``gene``/``protein_change`` ride along for the pdac presentation antigen-join
# (Task 5); they are strings and are excluded from every feature set by the
# banned-column list, so they never leak into a model.
_FEATURE_COLUMNS = [
    "candidate_id",
    "patient_id",
    "study_id",
    "mutant_peptide",
    "wildtype_peptide",
    "peptide_length",
    "hla_alleles",
    "mhc_class",
    "gene",
    "protein_change",
]


def parse_alleles(value: object) -> list[str]:
    """Normalize the polymorphic ``hla_alleles`` field to a clean list of strings."""
    if isinstance(value, (list, tuple, np.ndarray)):
        return [str(item).strip() for item in value if str(item).strip()]
    if pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    if text.startswith("["):
        try:
            return [str(item).strip() for item in json.loads(text) if str(item).strip()]
        except json.JSONDecodeError:
            return []
    return [item.strip() for item in text.replace(";", ",").split(",") if item.strip()]


def _first_allele(value: object) -> str:
    alleles = parse_alleles(value)
    return sorted(alleles)[0] if alleles else ""


def load_label_frame(corpus_dir: str | Path = CORPUS_DIR) -> pd.DataFrame:
    """Load the candidate-resolved Event-B label frame (one binary label per candidate)."""
    corpus_dir = Path(corpus_dir)
    candidates = pd.read_parquet(corpus_dir / "candidates.parquet")
    assays = pd.read_parquet(corpus_dir / "assays.parquet")
    primary = assays[
        assays.event_type.astype(str).eq(BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value)
        & assays.candidate_id.notna()
    ]
    frame = primary[["candidate_id", "response_label"]].merge(
        candidates[_FEATURE_COLUMNS], on="candidate_id", how="left", validate="one_to_one"
    )
    keep = [ResponseLabel.POSITIVE.value, ResponseLabel.TESTED_NEGATIVE.value]
    frame = frame[frame.response_label.isin(keep)].copy()
    frame["label"] = (frame.response_label == ResponseLabel.POSITIVE.value).astype(int)
    frame["hla_allele"] = frame.hla_alleles.map(_first_allele)
    return frame.drop(columns=["response_label"]).reset_index(drop=True)
