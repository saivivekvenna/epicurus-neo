"""Score IMPROVE with the frozen NeoPrecis-Immuno model.

IMPROVE publishes NetMHC percentile ranks but NeoPrecis expects raw NetMHC
scores. This adapter uses ``1 - rank / 100`` as an explicit approximation, so
its output is diagnostic and not a fully reproduced NeoPrecis benchmark.
"""

from __future__ import annotations

import argparse
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

from epicurus_neo.metrics import group_metrics, summarize_group_metrics
from epicurus_neo.neoprecis_adapter import (
    normalize_neoprecis_allele,
    wildtype_pseudo_core,
)

DEFAULT_MEMBER = (
    "data/03_data_for_CV/IMPROVE/"
    "03_3_final_peptide_features_Partition_for_CV.txt"
)


def read_improve(path: Path, member: str) -> pd.DataFrame:
    with zipfile.ZipFile(path) as archive, archive.open(member) as source:
        return pd.read_csv(source, sep="\t")


def load_model(repo: Path):
    sys.path.insert(0, str(repo))
    from neoprecis.CRD import PeptCRD

    model_dir = repo / "neoprecis" / "CRD"
    return PeptCRD(
        model_dir / "ref.h5",
        model_dir / "PeptCRD_config.yaml",
        model_dir / "PeptCRD_checkpoint.ckpt",
    )


def score_frame(raw: pd.DataFrame, model, batch_size: int) -> np.ndarray:
    import torch

    alleles = raw["HLA_allele"].map(normalize_neoprecis_allele)
    mutant_cores = raw["Core"].astype(str)
    wildtype_cores = pd.Series(
        [
            wildtype_pseudo_core(str(wildtype), str(mutant), str(core))
            for wildtype, mutant, core in zip(
                raw["Norm_peptide"],
                raw["Mut_peptide"],
                mutant_cores,
                strict=True,
            )
        ],
        index=raw.index,
    )
    approximate_binding_score = (
        1.0 - pd.to_numeric(raw["RankEL_4.1"], errors="coerce").fillna(100.0) / 100.0
    ).clip(0.0, 1.0)
    scores = np.full(len(raw), np.nan, dtype=np.float32)

    valid = (
        wildtype_cores.str.len().eq(9)
        & mutant_cores.str.len().eq(9)
        & alleles.isin(model.allele_dict)
    )
    valid_indices = np.flatnonzero(valid.to_numpy())
    for start in range(0, len(valid_indices), batch_size):
        indices = valid_indices[start : start + batch_size]
        wildtype_tokens = [model._tokenization(wildtype_cores.iloc[index]) for index in indices]
        mutant_tokens = [model._tokenization(mutant_cores.iloc[index]) for index in indices]
        peptides = torch.tensor(
            list(zip(mutant_tokens, wildtype_tokens, strict=True)),
            dtype=torch.long,
        )
        features = torch.tensor(
            [
                [0, model.allele_dict[alleles.iloc[index]], approximate_binding_score.iloc[index]]
                for index in indices
            ],
            dtype=torch.float32,
        )
        with torch.no_grad():
            _, predictions = model.model(peptides, features)
        scores[indices] = predictions.squeeze(1).cpu().numpy()
    return scores


def canonical_scored_frame(raw: pd.DataFrame, scores: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": [f"improve:{index}" for index in range(len(raw))],
            "patient_id": "improve:" + raw["Patient"].astype(str),
            "mutant_peptide": raw["Mut_peptide"].astype(str),
            "wildtype_peptide": raw["Norm_peptide"].astype(str),
            "hla_allele": raw["HLA_allele"].astype(str),
            "label": raw["response"].map({1: "positive", 0: "negative"}),
            "neoprecis_immuno_score": scores,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--neoprecis-repo", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--member", default=DEFAULT_MEMBER)
    parser.add_argument("--batch-size", type=int, default=1024)
    args = parser.parse_args()

    raw = read_improve(args.data, args.member)
    model = load_model(args.neoprecis_repo)
    scores = score_frame(raw, model, args.batch_size)
    scored = canonical_scored_frame(raw, scores)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    scored.to_csv(args.output, index=False)

    evaluated = scored.dropna(subset=["neoprecis_immuno_score"])
    summary = summarize_group_metrics(
        group_metrics(
            evaluated,
            group_col="patient_id",
            score_col="neoprecis_immuno_score",
            k=20,
        )
    )
    report = {
        "score_conversion": "approximate_raw_score = 1 - RankEL_4.1 / 100",
        "evaluated_candidates": len(evaluated),
        "total_candidates": len(scored),
        "summary": summary,
    }
    report_path = args.output.with_suffix(".json")
    report_path.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
