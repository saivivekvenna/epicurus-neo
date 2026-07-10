"""Deterministic, versioned recognition-corpus split manifests."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from hashlib import sha256
import json
from pathlib import Path

import pandas as pd


SPLIT_VERSION = "event-b-splits-1.0.0"


class SplitType(str, Enum):
    PATIENT_HOLDOUT = "PATIENT_HOLDOUT"
    STUDY_HOLDOUT = "STUDY_HOLDOUT"
    HLA_HOLDOUT = "HLA_HOLDOUT"
    PEPTIDE_CLUSTER_HOLDOUT = "PEPTIDE_CLUSTER_HOLDOUT"
    CANCER_TYPE_HOLDOUT = "CANCER_TYPE_HOLDOUT"
    TEMPORAL_HOLDOUT = "TEMPORAL_HOLDOUT"


@dataclass(frozen=True)
class SplitManifest:
    split_type: str
    split_version: str
    seed: int
    evaluation_fraction: float
    peptide_similarity_threshold: float | None
    temporal_cutoff: str | None
    assignments: tuple[dict, ...]

    def write(self, path: str | Path) -> None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(json.dumps(asdict(self), indent=2, sort_keys=True) + "\n")


class _UnionFind:
    def __init__(self, size: int):
        self.parent = list(range(size))

    def find(self, value: int) -> int:
        while self.parent[value] != value:
            self.parent[value] = self.parent[self.parent[value]]
            value = self.parent[value]
        return value

    def union(self, left: int, right: int) -> None:
        left_root, right_root = self.find(left), self.find(right)
        if left_root != right_root:
            self.parent[max(left_root, right_root)] = min(left_root, right_root)


def peptide_similarity(left: str, right: str) -> float:
    """Normalized Levenshtein similarity for leakage clustering."""
    left, right = str(left).upper(), str(right).upper()
    if not left and not right:
        return 1.0
    previous = list(range(len(right) + 1))
    for i, char_left in enumerate(left, start=1):
        current = [i]
        for j, char_right in enumerate(right, start=1):
            current.append(
                min(current[-1] + 1, previous[j] + 1, previous[j - 1] + (char_left != char_right))
            )
        previous = current
    return 1.0 - previous[-1] / max(len(left), len(right), 1)


def peptide_cluster_ids(peptides: pd.Series, *, threshold: float = 0.8) -> list[str]:
    if not 0 <= threshold <= 1:
        raise ValueError("peptide similarity threshold must be between zero and one")
    values = peptides.fillna("").astype(str).str.upper().tolist()
    unique = sorted(set(values))
    union = _UnionFind(len(unique))
    for left in range(len(unique)):
        for right in range(left + 1, len(unique)):
            if peptide_similarity(unique[left], unique[right]) >= threshold:
                union.union(left, right)
    roots = {value: union.find(index) for index, value in enumerate(unique)}
    labels = {
        root: "pepcluster:"
        + sha256("|".join(v for v in unique if roots[v] == root).encode()).hexdigest()[:16]
        for root in sorted(set(roots.values()))
    }
    return [labels[roots[value]] for value in values]


def _union_by_values(union: _UnionFind, values: pd.Series) -> None:
    first: dict[str, int] = {}
    for index, value in enumerate(values.fillna("<MISSING>").astype(str)):
        if value in first:
            union.union(first[value], index)
        else:
            first[value] = index


def _multivalues(value) -> tuple[str, ...]:
    text = str(value).strip()
    if text.startswith("["):
        try:
            return tuple(sorted(str(item).strip().upper() for item in json.loads(text)))
        except json.JSONDecodeError:
            pass
    return tuple(
        sorted(item.strip().upper() for item in text.replace(";", ",").split(",") if item.strip())
    )


def _union_by_overlapping_sets(union: _UnionFind, values: pd.Series) -> pd.Series:
    first: dict[str, int] = {}
    canonical = []
    for index, value in enumerate(values):
        items = _multivalues(value)
        canonical.append("|".join(items))
        for item in items:
            if item in first:
                union.union(first[item], index)
            else:
                first[item] = index
    return pd.Series(canonical)


def generate_split_manifest(
    frame: pd.DataFrame,
    split_type: SplitType | str,
    *,
    seed: int = 0,
    evaluation_fraction: float = 0.2,
    peptide_similarity_threshold: float = 0.8,
    temporal_cutoff: str | None = None,
) -> SplitManifest:
    split_type = SplitType(split_type)
    required = {"candidate_id", "patient_id", "study_id", "mutant_peptide"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Split input missing columns: {sorted(missing)}")
    if not 0 < evaluation_fraction < 1:
        raise ValueError("evaluation_fraction must be between zero and one")
    if frame.candidate_id.duplicated().any():
        raise ValueError("candidate_id must be unique in a split manifest")
    work = frame.reset_index(drop=True).copy()
    union = _UnionFind(len(work))

    if split_type is SplitType.PATIENT_HOLDOUT:
        target = work.patient_id.astype(str)
    elif split_type is SplitType.STUDY_HOLDOUT:
        target = work.study_id.astype(str)
    elif split_type is SplitType.HLA_HOLDOUT:
        if "hla_alleles" not in work:
            raise ValueError("HLA_HOLDOUT requires hla_alleles")
        target = _union_by_overlapping_sets(union, work.hla_alleles)
        _union_by_values(union, work.patient_id)
    elif split_type is SplitType.PEPTIDE_CLUSTER_HOLDOUT:
        target = pd.Series(
            peptide_cluster_ids(work.mutant_peptide, threshold=peptide_similarity_threshold)
        )
        _union_by_values(union, work.patient_id)
    elif split_type is SplitType.CANCER_TYPE_HOLDOUT:
        if "cancer_type" not in work:
            raise ValueError("CANCER_TYPE_HOLDOUT requires cancer_type")
        target = work.cancer_type.astype(str)
        _union_by_values(union, work.patient_id)
    else:
        if temporal_cutoff is None or "sample_date" not in work:
            raise ValueError("TEMPORAL_HOLDOUT requires temporal_cutoff and sample_date")
        target = work.patient_id.astype(str)
    if split_type is not SplitType.HLA_HOLDOUT:
        _union_by_values(union, target)

    components: dict[int, list[int]] = {}
    for index in range(len(work)):
        components.setdefault(union.find(index), []).append(index)
    if len(components) < 2:
        raise ValueError(f"{split_type.value} has fewer than two leakage-safe components")

    if split_type is SplitType.TEMPORAL_HOLDOUT:
        cutoff = pd.to_datetime(temporal_cutoff, errors="raise", utc=True)
        eval_roots = {
            root
            for root, indices in components.items()
            if pd.to_datetime(work.loc[indices, "sample_date"], errors="coerce", utc=True).max()
            >= cutoff
        }
        if not eval_roots or len(eval_roots) == len(components):
            raise ValueError(
                "temporal cutoff does not produce non-empty train and evaluation splits"
            )
    else:
        ordered = sorted(
            components,
            key=lambda root: sha256(f"{seed}|{split_type.value}|{root}".encode()).hexdigest(),
        )
        eval_count = min(len(ordered) - 1, max(1, round(len(ordered) * evaluation_fraction)))
        eval_roots = set(ordered[:eval_count])

    assignments = []
    for index, row in work.iterrows():
        root = union.find(index)
        assignments.append(
            {
                "candidate_id": str(row.candidate_id),
                "patient_id": str(row.patient_id),
                "study_id": str(row.study_id),
                "group_id": str(target.iloc[index]),
                "component_id": f"component:{root}",
                "split": "evaluation" if root in eval_roots else "train",
            }
        )
    manifest = SplitManifest(
        split_type.value,
        SPLIT_VERSION,
        seed,
        evaluation_fraction,
        peptide_similarity_threshold if split_type is SplitType.PEPTIDE_CLUSTER_HOLDOUT else None,
        temporal_cutoff,
        tuple(sorted(assignments, key=lambda row: row["candidate_id"])),
    )
    assert_split_leakage_safe(manifest)
    return manifest


def assert_split_leakage_safe(manifest: SplitManifest) -> None:
    frame = pd.DataFrame(manifest.assignments)
    for column in ("candidate_id", "patient_id"):
        cross = frame.groupby(column).split.nunique()
        if (cross > 1).any():
            raise AssertionError(f"{column} leaks across split through alternate rows")
    if manifest.split_type == SplitType.STUDY_HOLDOUT.value:
        cross = frame.groupby("study_id").split.nunique()
        if (cross > 1).any():
            raise AssertionError("study leaks across STUDY_HOLDOUT")
    if set(frame.split) != {"train", "evaluation"}:
        raise AssertionError(
            "split manifest must contain non-empty train and evaluation partitions"
        )
