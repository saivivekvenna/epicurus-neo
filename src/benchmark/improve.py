"""Read and verify the official IMPROVE archives shipped in IMPROVE_paper."""

from __future__ import annotations

from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd

from benchmark.headroom import headroom_table, random_expectation, tie_break_canary
from benchmark.metrics import capture_fraction, hits_at_k, p_at_least_one
from benchmark.stats import mde, paired_bootstrap


DATA_MEMBER = "data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt"
RF_MEMBER = "results copy/5_fold_CV/TME_excluded/pred_df_TME_excluded.txt"
RF_WO_PRIME_MEMBER = "results copy/5_fold_CV/TME_excluded/pred_df_wo_primeTME_excluded.txt"


def _canonical(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.rename(
        columns={
            "Patient": "patient_id",
            "Mut_peptide": "mutant_peptide",
            "Norm_peptide": "wildtype_peptide",
            "HLA_allele": "hla_allele",
            "response": "label",
        }
    )


def load_improve_data(repo: str | Path) -> pd.DataFrame:
    archive = Path(repo) / "data.zip"
    if not archive.is_file():
        raise FileNotFoundError(f"Missing {archive}; clone SRHgroup/IMPROVE_paper")
    with ZipFile(archive) as bundle, bundle.open(DATA_MEMBER) as handle:
        return _canonical(pd.read_csv(handle, sep="\t"))


def load_improve_rf(repo: str | Path, *, without_prime: bool = False) -> pd.DataFrame:
    archive = Path(repo) / "results.zip"
    if not archive.is_file():
        raise FileNotFoundError(f"Missing {archive}; clone SRHgroup/IMPROVE_paper")
    member = RF_WO_PRIME_MEMBER if without_prime else RF_MEMBER
    with ZipFile(archive) as bundle, bundle.open(member) as handle:
        return _canonical(pd.read_csv(handle, sep=r"\s+", engine="python"))


def regression_values(repo: str | Path) -> dict[str, object]:
    """Compute every official-CV quantity registered in Milestone 1 §4."""
    data = load_improve_data(repo)
    data["oracle"] = data["label"]
    baseline_specs = {
        "oracle_mean_hits_at_20": ("oracle", False),
        "netmhcpan_rankel_4_1": ("RankEL_4.1", True),
        "prime": ("Prime", False),
        "prioscore": ("PrioScore", False),
        "foreignness": ("Foreigness", False),
        "dai_4_1": ("DAI_4.1", False),
    }
    values: dict[str, object] = {
        "rows": len(data),
        "positives": int(data["label"].sum()),
        "patients": int(data["patient_id"].nunique()),
        "partitions": int(data["Partition"].nunique()),
        "random_expectation": random_expectation(data),
        "positive_count_sd": float(data.groupby("patient_id")["label"].sum().std(ddof=1)),
        "tie_break_canary": tie_break_canary(data),
    }
    for name, (column, ascending) in baseline_specs.items():
        values[name] = float(hits_at_k(data, score_col=column, ascending=ascending).mean())

    rf = load_improve_rf(repo)
    rf_hits = hits_at_k(rf, score_col="prediction_rf")
    prime_hits = hits_at_k(rf, score_col="Prime")
    paired = paired_bootstrap(rf_hits, prime_hits)
    rf_capture = capture_fraction(rf, score_col="prediction_rf")
    prime_capture = capture_fraction(rf, score_col="Prime")
    capture_paired = paired_bootstrap(rf_capture, prime_capture)
    rf_clinical = p_at_least_one(rf, score_col="prediction_rf")
    prime_clinical = p_at_least_one(rf, score_col="Prime")
    clinical_paired = paired_bootstrap(rf_clinical, prime_clinical)
    values.update(
        {
            "rf": float(rf_hits.mean()),
            "rf_without_prime": float(
                hits_at_k(
                    load_improve_rf(repo, without_prime=True), score_col="prediction_rf"
                ).mean()
            ),
            "rf_vs_prime": paired,
            "rf_vs_prime_capture_fraction": capture_paired,
            "rf_p_at_least_one": float(np.nanmean(rf_clinical)),
            "prime_p_at_least_one": float(np.nanmean(prime_clinical)),
            "rf_vs_prime_p_at_least_one": clinical_paired,
            "unreachable_patients": int(np.isnan(rf_clinical).sum()),
            "paired_sd_diff": float(np.std(rf_hits - prime_hits, ddof=1)),
            "paired_mde": float(mde(rf_hits, prime_hits)),
            "headroom": headroom_table(data, base_score_col="Prime"),
        }
    )
    return values
