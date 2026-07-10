"""Event-A (pre-existing reactivity) label frame for the M6B transfer teacher.

Event-A is loaded from the IMPROVE corpus and is *never* merged with Event-B labels:
the teacher trains on these rows only, then scores Event-B candidates. The schema is
identical to the Event-B backbone, so the M6A feature builders apply unchanged.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from event_b.models import BiologicalEvent, ResponseLabel

from epicurus_neo.m6.dataset import _FEATURE_COLUMNS, _first_allele

EVENT_A_CORPUS_DIR = Path("outputs/event_b_corpus_combined")
EVENT_A_STUDY = "improve"


def load_event_a_frame(corpus_dir: str | Path = EVENT_A_CORPUS_DIR) -> pd.DataFrame:
    """Load the IMPROVE Event-A frame (one binary recognition label per assay row)."""
    corpus_dir = Path(corpus_dir)
    candidates = pd.read_parquet(corpus_dir / "candidates.parquet")
    assays = pd.read_parquet(corpus_dir / "assays.parquet")
    primary = assays[
        assays.event_type.astype(str).eq(BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value)
        & assays.candidate_id.notna()
        & assays.study_id.astype(str).eq(EVENT_A_STUDY)
    ]
    frame = primary[["candidate_id", "response_label"]].merge(
        candidates[_FEATURE_COLUMNS], on="candidate_id", how="left", validate="many_to_one"
    )
    keep = [ResponseLabel.POSITIVE.value, ResponseLabel.TESTED_NEGATIVE.value]
    frame = frame[frame.response_label.isin(keep)].copy()
    frame["label"] = (frame.response_label == ResponseLabel.POSITIVE.value).astype(int)
    frame["hla_allele"] = frame.hla_alleles.map(_first_allele)
    return frame.drop(columns=["response_label"]).reset_index(drop=True)
