"""Label-free feature construction for the Zhao PRIME neck-to-neck benchmark.

Every column here is built from peptide sequence, HLA restriction, and the published
per-peptide presentation predictors ONLY. No response label, patient id, study id, or
post-outcome variable is used to construct a feature (the label is attached at the end
purely so the ranking harness can score; it never enters feature construction).

Incumbent (a-priori, not label-selected): ``mixmhcPred-3.0`` — the MixMHCpred binding
score that open-source PRIME is built upon. It is the fair stand-in for PRIME's binding
backbone while a true PRIME run is a documented blocker (see zhao_benchmark).

Feature families (each independently ablatable):
    hydrophobicity          GRAVY, non-anchor (TCR-facing) GRAVY, aromatic fraction
    charge                  net charge, positive/negative residue counts
    mutant_vs_wt            mutation position/anchor, KD & charge deltas at the mutated site
    presentation_predictors orthogonal published predictors (NOT the incumbent), oriented
    hla_calibration         per-HLA one-hot for HLA-specific residual calibration
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from event_b.adapters.zhao_dc import TABLE_FILE, immunogenic_label, read_peptide_sheet, stage_zhao_supplements
from event_b.models import ResponseLabel, stable_candidate_id


# Kyte-Doolittle hydropathy.
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5, "G": -0.4,
    "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6, "S": -0.8,
    "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
AROMATIC = set("FWY")
POSITIVE_AA = set("KR")
NEGATIVE_AA = set("DE")

# The incumbent is a-priori MixMHCpred-3.0 (PRIME's binding backbone).
INCUMBENT_COLUMN = "mixmhcPred-3.0"
INCUMBENT_NAME = "mixmhcpred3_score"

# Published predictor columns and their orientation to "higher = more immunogenic-favourable".
# rank/affinity columns are lower-is-better and are negated in log space; likelihood/score
# columns are higher-is-better and kept. Orientation is by documented predictor semantics.
LOWER_IS_BETTER_AFFINITY = {"NetMHC-3", "NetMHC-4", "MHCflurryBA"}  # rank / nM affinity
HIGHER_IS_BETTER = {
    "BigMHC EL", "NetMHCstab", "NetCTLpan TAP", "NetCTLpan Cleavage",
    "MHCflurryProc", "MHCflurryPres", INCUMBENT_COLUMN,
}
PRESENTATION_PREDICTOR_COLUMNS = sorted(
    (LOWER_IS_BETTER_AFFINITY | HIGHER_IS_BETTER) - {INCUMBENT_COLUMN}
)

# Deterministic presentation floor (guard): a peptide with MHCflurry binding affinity worse
# than this many nM is treated as effectively unpresentable and may never be promoted above a
# presentable candidate by the challenger. The floor is set conservatively at "no meaningful
# binding" (not merely "weak binding") because every Zhao peptide was actually administered and
# assayed, so only a genuine non-binder is deterministically unpresentable here. On this
# vaccinated subset the guard is therefore near-inert by construction (an honest limitation:
# a presentation floor cannot filter an already-vaccinated set); it earns its keep only once
# the challenger is applied to a full somatic candidate universe.
PRESENTATION_AFFINITY_FLOOR_NM = 20000.0


def _oriented(column: str, values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    if column in LOWER_IS_BETTER_AFFINITY:
        # negative log affinity/rank: higher = stronger binder, well-scaled heavy tails.
        return -np.log1p(numeric.clip(lower=0))
    return numeric


def _seq_features(peptide: str) -> dict[str, float]:
    peptide = str(peptide).strip().upper()
    length = len(peptide)
    kd = [KD.get(aa, 0.0) for aa in peptide]
    gravy = float(np.mean(kd)) if kd else 0.0
    # Non-anchor / TCR-facing = exclude position 2 (index 1) and the C-terminus.
    anchor_idx = {1, length - 1}
    nonanchor = [kd[i] for i in range(length) if i not in anchor_idx]
    nonanchor_gravy = float(np.mean(nonanchor)) if nonanchor else gravy
    aromatic_frac = sum(1 for aa in peptide if aa in AROMATIC) / length if length else 0.0
    n_pos = sum(1 for aa in peptide if aa in POSITIVE_AA)
    n_neg = sum(1 for aa in peptide if aa in NEGATIVE_AA)
    return {
        "gravy": gravy,
        "nonanchor_gravy": nonanchor_gravy,
        "aromatic_frac": aromatic_frac,
        "net_charge": float(n_pos - n_neg),
        "n_positive": float(n_pos),
        "n_negative": float(n_neg),
    }


def _mut_features(mut: str, wt: str, position: object, length: int) -> dict[str, float]:
    mut = str(mut).strip().upper()
    wt = str(wt).strip().upper()
    diff = [i for i in range(min(len(mut), len(wt))) if mut[i] != wt[i]]
    if len(diff) == 1:
        idx = diff[0]
    else:
        try:
            idx = int(float(position)) - 1
        except (TypeError, ValueError):
            idx = -1
    if 0 <= idx < len(mut) and idx < len(wt):
        kd_delta = KD.get(mut[idx], 0.0) - KD.get(wt[idx], 0.0)
        charge_delta = (
            (1.0 if mut[idx] in POSITIVE_AA else (-1.0 if mut[idx] in NEGATIVE_AA else 0.0))
            - (1.0 if wt[idx] in POSITIVE_AA else (-1.0 if wt[idx] in NEGATIVE_AA else 0.0))
        )
        pos1 = idx + 1
        is_anchor = 1.0 if (pos1 == 2 or pos1 == length) else 0.0
        pos_frac = pos1 / length if length else 0.0
    else:
        kd_delta = charge_delta = is_anchor = pos_frac = 0.0
    return {
        "mut_kd_delta": kd_delta,
        "mut_charge_delta": charge_delta,
        "mut_is_anchor": is_anchor,
        "mut_position_frac": pos_frac,
    }


FEATURE_FAMILIES: dict[str, list[str]] = {
    "hydrophobicity": ["gravy", "nonanchor_gravy", "aromatic_frac"],
    "charge": ["net_charge", "n_positive", "n_negative"],
    "mutant_vs_wt": ["mut_kd_delta", "mut_charge_delta", "mut_is_anchor", "mut_position_frac"],
    "presentation_predictors": [f"pred_{c}" for c in PRESENTATION_PREDICTOR_COLUMNS],
    # hla_calibration columns are added dynamically (one per observed HLA allele).
    "hla_calibration": [],
}


def build_zhao_feature_frame(raw_dir: str | Path) -> pd.DataFrame:
    """Return one label-free feature row per candidate identity (deduplicated repeats).

    The label column (POSITIVE/TESTED_NEGATIVE and numeric ``y``) is attached for the
    harness only; it is not consumed by any feature above.
    """
    paths = stage_zhao_supplements(raw_dir)
    sheet = read_peptide_sheet(paths[TABLE_FILE])
    sheet = sheet[sheet["Patient ID"].notna()].copy()

    rows: list[dict] = []
    seen: set[str] = set()
    for _, row in sheet.iterrows():
        source_patient = str(int(float(row["Patient ID"])))
        peptide = str(row["peptide(mut)"]).strip().upper()
        wt = str(row["peptide(wt)"]).strip().upper()
        hla = str(row["HLA type"]).strip().upper()
        identity = {
            "study_id": "zhao_dc_2026",
            "patient_id": f"zhao_dc_2026:{source_patient}",
            "sample_id": "",
            "timepoint": "",
            "genomic_variant": "",
            "transcript": "",
            "mutant_peptide": peptide,
            "hla_alleles": hla,
        }
        candidate_id = stable_candidate_id(identity)
        if candidate_id in seen:
            continue  # repeated assay of an identical candidate; scored once for ranking
        seen.add(candidate_id)

        label = immunogenic_label(row["ELSPOT ratio"])[0]
        if label == ResponseLabel.UNTESTED.value:
            continue  # unscorable; never enters the eval universe as a silent negative
        length = len(peptide)
        record = {
            "candidate_id": candidate_id,
            "patient_id": f"zhao_dc_2026:{source_patient}",
            "mutant_peptide": peptide,
            "hla_allele": hla,
            "hla_gene": hla[4:5] if len(hla) > 4 else "",
            "label": label,
            "y": 1 if label == ResponseLabel.POSITIVE.value else 0,
            INCUMBENT_NAME: pd.to_numeric(row[INCUMBENT_COLUMN], errors="coerce"),
            "mhcflurry_affinity_nm": pd.to_numeric(row["MHCflurryBA"], errors="coerce"),
        }
        record.update(_seq_features(peptide))
        record.update(_mut_features(peptide, wt, row.get("position"), length))
        for column in PRESENTATION_PREDICTOR_COLUMNS:
            record[f"pred_{column}"] = float(_oriented(column, pd.Series([row[column]])).iloc[0])
        rows.append(record)

    frame = pd.DataFrame(rows)
    # Presentation guard flag (deterministic): effectively unpresentable weak binders.
    frame["presentation_ok"] = ~(frame["mhcflurry_affinity_nm"] > PRESENTATION_AFFINITY_FLOOR_NM)
    # Per-HLA one-hot for HLA-specific residual calibration.
    hla_dummies = pd.get_dummies(frame["hla_allele"], prefix="hla").astype(float)
    FEATURE_FAMILIES["hla_calibration"] = sorted(hla_dummies.columns)
    frame = pd.concat([frame, hla_dummies], axis=1)
    return frame.sort_values("candidate_id").reset_index(drop=True)


def static_feature_columns(frame: pd.DataFrame, families: list[str] | None = None) -> list[str]:
    families = families or list(FEATURE_FAMILIES)
    cols: list[str] = []
    for family in families:
        cols.extend(c for c in FEATURE_FAMILIES[family] if c in frame.columns)
    return cols


def feature_availability() -> list[dict]:
    """The reviewer's requested correction-feature families vs Zhao recoverability."""
    return [
        {"family": "mutant_vs_wt_contrast", "recoverable": True,
         "basis": "peptide(mut) and peptide(wt) both published per row; SNV single-substitution."},
        {"family": "tcr_facing_hydrophobicity_charge", "recoverable": True,
         "basis": "recomputed from sequence over non-anchor (pos!=2,!=C-term) residues."},
        {"family": "hla_specific_residual_calibration", "recoverable": True,
         "basis": "restricting HLA published per peptide; per-allele one-hot (63 alleles)."},
        {"family": "presentation_predictors", "recoverable": True,
         "basis": "NetMHC-3/4, MHCflurry BA/Proc/Pres, BigMHC EL, NetMHCstab, NetCTLpan published."},
        {"family": "assay_neighbor_retrieval", "recoverable": True,
         "basis": "train-fold-only positive/negative sequence-neighbor retrieval (OOF-safe); "
                  "implemented as an optional fold feature in the harness."},
        {"family": "peptide_novelty_self_similarity", "recoverable": "PARTIAL",
         "basis": "no bundled human proteome reference; a within-corpus novelty proxy is possible "
                  "but a true self-similarity feature needs an external proteome (blocker)."},
        {"family": "contradiction_reliability_weighting", "recoverable": "PARTIAL",
         "basis": "only 2 repeated-assay candidates in Zhao (both concordant); little signal to "
                  "weight until CEDAR/backbone contradictions are pooled."},
        {"family": "assay_and_vaccine_context", "recoverable": False,
         "basis": "single assay (ELISpot) and single vaccine platform (DC) across all rows -> "
                  "no within-Zhao contrast; becomes available only after cross-study pooling."},
        {"family": "prime_incumbent_score", "recoverable": False,
         "basis": "PRIME/MixMHCpred final immunogenicity score not in the supplement and not "
                  "runnable locally; MixMHCpred-3.0 binding score is the a-priori stand-in and a "
                  "reproducible PRIME command + input artifact is emitted."},
    ]
