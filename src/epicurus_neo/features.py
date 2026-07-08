from __future__ import annotations

import math

import numpy as np
import pandas as pd

from epicurus_neo.schema import normalize_peptide


NON_FEATURE_COLUMNS = {
    "candidate_id",
    "source_dataset",
    "study_id",
    "patient_id",
    "hla_allele",
    "hla_allele_norm",
    "mutant_peptide",
    "mutant_peptide_norm",
    "wildtype_peptide",
    "wildtype_peptide_norm",
    "mutant_hla_key",
    "wildtype_hla_key",
    "label",
    "label_weight",
    "assay_type",
    "split",
    "ID",
    "Patient",
    "patient",
    "Response",
    "response",
    "Response Type",
    "response_type",
    "Screening Status",
    "target_value",
    "immunogenicity",
    "Immunogenicity",
    "reactivity",
    "Reactivity",
    "TIL Reactivity",
}


def infer_numeric_feature_columns(frame: pd.DataFrame) -> list[str]:
    """Infer usable numeric feature columns from a canonical candidate table."""
    columns: list[str] = []
    for column in frame.columns:
        if column in NON_FEATURE_COLUMNS:
            continue
        if pd.api.types.is_numeric_dtype(frame[column]):
            columns.append(column)
    return sorted(columns)


def safe_log_inverse(value: object) -> float:
    """Score a positive quantity where smaller is better, such as binding nM."""
    if value is None or pd.isna(value):
        return float("nan")
    value_float = float(value)
    if value_float <= 0:
        return float("nan")
    return -math.log10(value_float)


def zscore(series: pd.Series) -> pd.Series:
    values = pd.to_numeric(series, errors="coerce")
    std = values.std(skipna=True)
    if pd.isna(std) or std == 0:
        return pd.Series(np.zeros(len(series)), index=series.index)
    return (values - values.mean(skipna=True)) / std


HYDROPHOBICITY = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

CHARGE = {
    "D": -1,
    "E": -1,
    "K": 1,
    "R": 1,
    "H": 0.5,
}

AMINO_ACIDS = tuple("ACDEFGHIKLMNPQRSTVWY")
AROMATIC = {"F", "W", "Y"}
POLAR = {"S", "T", "N", "Q", "C", "Y"}
ACIDIC = {"D", "E"}
BASIC = {"K", "R", "H"}


def anchor_positions(length: int) -> set[int]:
    """Return simple class-I anchor positions using zero-based indices."""
    if length < 2:
        return set()
    return {1, length - 1}


def mutation_deltas(mutant: object, wildtype: object) -> dict[str, float]:
    mutant_peptide = normalize_peptide(mutant)
    wildtype_peptide = normalize_peptide(wildtype)
    if not mutant_peptide or not wildtype_peptide or len(mutant_peptide) != len(wildtype_peptide):
        return {
            "mutation_count": float("nan"),
            "mutation_anchor_count": float("nan"),
            "mutation_tcr_face_count": float("nan"),
            "mutation_hydrophobicity_delta": float("nan"),
            "mutation_charge_delta": float("nan"),
        }

    anchors = anchor_positions(len(mutant_peptide))
    changed_positions = [
        idx
        for idx, (mut_aa, wt_aa) in enumerate(zip(mutant_peptide, wildtype_peptide, strict=True))
        if mut_aa != wt_aa
    ]
    hydrophobicity_delta = sum(
        HYDROPHOBICITY.get(mutant_peptide[idx], 0.0) - HYDROPHOBICITY.get(wildtype_peptide[idx], 0.0)
        for idx in changed_positions
    )
    charge_delta = sum(
        CHARGE.get(mutant_peptide[idx], 0.0) - CHARGE.get(wildtype_peptide[idx], 0.0)
        for idx in changed_positions
    )

    return {
        "mutation_count": float(len(changed_positions)),
        "mutation_anchor_count": float(sum(1 for idx in changed_positions if idx in anchors)),
        "mutation_tcr_face_count": float(sum(1 for idx in changed_positions if idx not in anchors)),
        "mutation_hydrophobicity_delta": float(hydrophobicity_delta),
        "mutation_charge_delta": float(charge_delta),
    }


def add_contrastive_features(frame: pd.DataFrame) -> pd.DataFrame:
    """Add mutant-vs-wildtype features tied to recognition plausibility."""
    if "mutant_peptide" not in frame.columns or "wildtype_peptide" not in frame.columns:
        return frame.copy()

    out = frame.copy()
    deltas = [
        mutation_deltas(mutant, wildtype)
        for mutant, wildtype in zip(out["mutant_peptide"], out["wildtype_peptide"], strict=True)
    ]
    delta_frame = pd.DataFrame(deltas, index=out.index)
    for column in delta_frame.columns:
        if column not in out.columns:
            out[column] = delta_frame[column]
    return out


def peptide_sequence_features(peptide: object) -> dict[str, float]:
    sequence = normalize_peptide(peptide)
    if not sequence:
        result = {
            "seq_len": float("nan"),
            "seq_hydrophobicity_mean": float("nan"),
            "seq_charge_sum": float("nan"),
            "seq_aromatic_fraction": float("nan"),
            "seq_polar_fraction": float("nan"),
            "seq_acidic_fraction": float("nan"),
            "seq_basic_fraction": float("nan"),
            "seq_cysteine_fraction": float("nan"),
            "seq_proline_fraction": float("nan"),
            "seq_glycine_fraction": float("nan"),
        }
        result.update({f"seq_aa_frac_{aa}": float("nan") for aa in AMINO_ACIDS})
        return result

    length = len(sequence)
    result = {
        "seq_len": float(length),
        "seq_hydrophobicity_mean": float(
            sum(HYDROPHOBICITY.get(aa, 0.0) for aa in sequence) / length
        ),
        "seq_charge_sum": float(sum(CHARGE.get(aa, 0.0) for aa in sequence)),
        "seq_aromatic_fraction": float(sum(aa in AROMATIC for aa in sequence) / length),
        "seq_polar_fraction": float(sum(aa in POLAR for aa in sequence) / length),
        "seq_acidic_fraction": float(sum(aa in ACIDIC for aa in sequence) / length),
        "seq_basic_fraction": float(sum(aa in BASIC for aa in sequence) / length),
        "seq_cysteine_fraction": float(sequence.count("C") / length),
        "seq_proline_fraction": float(sequence.count("P") / length),
        "seq_glycine_fraction": float(sequence.count("G") / length),
    }
    result.update({f"seq_aa_frac_{aa}": float(sequence.count(aa) / length) for aa in AMINO_ACIDS})
    return result


def add_sequence_features(frame: pd.DataFrame) -> pd.DataFrame:
    if "mutant_peptide" not in frame.columns:
        return frame.copy()
    out = frame.copy()
    features = [peptide_sequence_features(peptide) for peptide in out["mutant_peptide"]]
    feature_frame = pd.DataFrame(features, index=out.index)
    for column in feature_frame.columns:
        if column not in out.columns:
            out[column] = feature_frame[column]
    return out


def add_baseline_scores(frame: pd.DataFrame) -> pd.DataFrame:
    """Add deterministic baseline rank scores when source feature columns exist."""
    out = add_sequence_features(add_contrastive_features(frame))

    if "binding_affinity_nm" in out.columns:
        out["baseline_binding_score"] = out["binding_affinity_nm"].map(safe_log_inverse)
    elif "binding_score" in out.columns:
        out["baseline_binding_score"] = pd.to_numeric(out["binding_score"], errors="coerce")

    presentation_candidates = [
        "presentation_score",
        "bigmhc_el_score",
        "mhcflurry_presentation_score",
    ]
    for column in presentation_candidates:
        if column in out.columns:
            out["baseline_presentation_score"] = pd.to_numeric(out[column], errors="coerce")
            break

    components: list[pd.Series] = []
    if "baseline_binding_score" in out.columns:
        components.append(zscore(out["baseline_binding_score"]))
    if "baseline_presentation_score" in out.columns:
        components.append(zscore(out["baseline_presentation_score"]))
    if "expression_tpm" in out.columns:
        components.append(zscore(np.log1p(pd.to_numeric(out["expression_tpm"], errors="coerce"))))
    if "mutant_wildtype_binding_delta" in out.columns:
        components.append(zscore(out["mutant_wildtype_binding_delta"]))
    if "foreignness_score" in out.columns:
        components.append(zscore(out["foreignness_score"]))

    if components:
        out["baseline_pvac_style_score"] = sum(components) / len(components)

    if "Nmer score" in out.columns:
        out["baseline_gartner_nmer_score"] = pd.to_numeric(out["Nmer score"], errors="coerce")
    if "Top netMHCpan4.0 EL ranked minimal" in out.columns:
        out["baseline_netmhcpan_el_score"] = -pd.to_numeric(
            out["Top netMHCpan4.0 EL ranked minimal"], errors="coerce"
        )
    if "Top MHCflurry ranked minimal" in out.columns:
        out["baseline_mhcflurry_score"] = -pd.to_numeric(
            out["Top MHCflurry ranked minimal"], errors="coerce"
        )

    return out
