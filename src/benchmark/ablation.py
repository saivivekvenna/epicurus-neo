"""Generate blind, balanced question sets for the §8.2 masking replication."""

from __future__ import annotations

import json
from hashlib import sha256
from pathlib import Path

import numpy as np
import pandas as pd


def _candidate_id(row: pd.Series) -> str:
    identity = f"{row['mutant_peptide']}|{row['hla_allele']}|{row['patient_id']}"
    return sha256(identity.encode()).hexdigest()[:16]


def _amino_acid_changes(wildtype: str, mutant: str) -> list[str]:
    changes = [
        f"{old}>{new}" for old, new in zip(str(wildtype), str(mutant), strict=False) if old != new
    ]
    if len(str(wildtype)) != len(str(mutant)):
        changes.append(f"length:{len(str(wildtype))}>{len(str(mutant))}")
    return changes


def _balanced_indices(frame: pd.DataFrame, rng: np.random.Generator, n: int) -> np.ndarray:
    if n % 2:
        raise ValueError("balanced question-set size must be even")
    labels = pd.to_numeric(frame["label"], errors="raise").to_numpy()
    positives = np.flatnonzero(labels == 1)
    negatives = np.flatnonzero(labels == 0)
    half = n // 2
    if len(positives) < half or len(negatives) < half:
        raise ValueError("not enough positive and negative rows for balanced sampling")
    picked = np.concatenate(
        [rng.choice(positives, half, replace=False), rng.choice(negatives, half, replace=False)]
    )
    rng.shuffle(picked)
    return picked


def generate_masking_sets(
    frame: pd.DataFrame,
    *,
    seeds: range = range(10),
    n_per_condition: int = 50,
) -> dict[int, dict[str, list[dict[str, object]]]]:
    """Create ten disjoint-within-seed A/B sets with labels kept in answer keys."""
    required = {
        "patient_id",
        "mutant_peptide",
        "wildtype_peptide",
        "hla_allele",
        "label",
        "Gene_Symbol",
        "Mutation_Consequence",
        "Protein_position",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Ablation source is missing columns: {sorted(missing)}")

    result: dict[int, dict[str, list[dict[str, object]]]] = {}
    for seed in seeds:
        rng = np.random.default_rng(seed)
        a_indices = _balanced_indices(frame, rng, n_per_condition)
        remaining = frame.drop(frame.index[a_indices]).reset_index(drop=True)
        b_local_indices = _balanced_indices(remaining, rng, n_per_condition)
        condition_a: list[dict[str, object]] = []
        condition_b: list[dict[str, object]] = []
        answers: list[dict[str, object]] = []
        for condition, selected in (
            ("A_UNMASKED", frame.iloc[a_indices]),
            ("B_MASKED", remaining.iloc[b_local_indices]),
        ):
            for _, row in selected.iterrows():
                candidate_id = _candidate_id(row)
                if condition == "A_UNMASKED":
                    question = {
                        "candidate_id": candidate_id,
                        "condition": condition,
                        "gene": row["Gene_Symbol"],
                        "mutation_consequence": row["Mutation_Consequence"],
                        "protein_position": str(row["Protein_position"]),
                        "amino_acid_changes": _amino_acid_changes(
                            row["wildtype_peptide"], row["mutant_peptide"]
                        ),
                        "hla": row["hla_allele"],
                    }
                    condition_a.append(question)
                else:
                    question = {
                        "candidate_id": candidate_id,
                        "condition": condition,
                        "mutant_peptide": row["mutant_peptide"],
                        "wildtype_peptide": row["wildtype_peptide"],
                        "hla": row["hla_allele"],
                    }
                    condition_b.append(question)
                answers.append(
                    {
                        "candidate_id": candidate_id,
                        "condition": condition,
                        "label": int(row["label"]),
                    }
                )
        result[int(seed)] = {
            "condition_a": condition_a,
            "condition_b": condition_b,
            "answer_key": answers,
        }
    return result


def write_masking_sets(
    sets: dict[int, dict[str, list[dict[str, object]]]], output_dir: str | Path
) -> None:
    output = Path(output_dir)
    questions = output / "questions"
    answers = output / "answer_keys"
    questions.mkdir(parents=True, exist_ok=True)
    answers.mkdir(parents=True, exist_ok=True)
    for seed, payload in sets.items():
        question_path = questions / f"seed_{seed:02d}.jsonl"
        answer_path = answers / f"seed_{seed:02d}.jsonl"
        with question_path.open("w") as handle:
            for row in payload["condition_a"] + payload["condition_b"]:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
        with answer_path.open("w") as handle:
            for row in payload["answer_key"]:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
