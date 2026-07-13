"""Per-patient reconstructability scoring and study-identity confound diagnostics.

This module answers a single operational question: *how good is the data?* It grades every
patient against the "excellent patient" ideal — a fully reconstructable vaccine-target decision
problem (tumor evidence -> full candidate universe -> selection -> tested outcomes) — and it
measures the corpus-level shortcut risk that made past recognition models good at *matching
studies* but poor at *predicting recognition*.

Nothing here fits a model or invents evidence. Every axis scores only source-present columns, so
a corpus that lacks WES/RNA features or a candidate denominator scores low by construction, which
is the finding, not a bug. When a future loader adds those columns the scores rise automatically.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd

from event_b.corpus import EventBCorpus
from event_b.models import BiologicalEvent, MHCClass, ResponseLabel


DATA_QUALITY_VERSION = "event-b-data-quality-1.1.0"

# Columns where a reconstructed WES / RNA evidence loader would deposit numeric features. They do
# not exist in the label-only corpus today; listing them keeps the axes honest and future-proof.
WES_CANDIDATE_COLUMNS = (
    "dna_vaf",
    "tumor_dna_vaf",
    "clonality_ccf",
    "cancer_cell_fraction",
    "tumor_purity",
    "copy_number",
)
RNA_CANDIDATE_COLUMNS = (
    "rna_vaf",
    "expression_tpm",
    "gene_expression_tpm",
    "rna_mutant_reads",
    "rna_depth",
)
# Evidence families that would corroborate a WES / RNA feature at inference time. TUMOR_CLONALITY is
# a WES-derived clonality/CCF call. RNA has no qualifying family: RNA support must come from the
# numeric RNA_CANDIDATE_COLUMNS. IMMUNOPEPTIDOMICS is deliberately absent here — mass-spec peptide
# presentation is direct presentation evidence, not transcription evidence, and under the core
# WES+RNA accessibility contract it must never inflate the rna_features score. It is surfaced as an
# unweighted diagnostic (`diag_immunopeptidomics_present`) instead.
WES_EVIDENCE_FAMILIES = ("TUMOR_CLONALITY",)
RNA_EVIDENCE_FAMILIES: tuple[str, ...] = ()
IMMUNOPEPTIDOMICS_FAMILIES = ("IMMUNOPEPTIDOMICS", "MASS_SPEC", "MHC_LIGANDOMICS")

# Middle reachability-funnel stages. A candidate that dropped here is only informative for a
# *reconstructed pre-selection universe* if the drop was actually assessed (not "not_assessed" /
# "not_applicable"). See ``_complete_reconstructed_universe``.
MIDDLE_FUNNEL_STAGES = ("survives_gating", "presentation_candidate", "ranking_stage", "top_k")


@dataclass(frozen=True)
class QualityAxis:
    key: str
    weight: float
    ideal_item: str
    description: str


# Weights emphasise the north-star-critical, currently-missing evidence: the candidate denominator
# (ideal item #4, "the single most important item") and the WES/RNA features that tell us a target
# is expressed. Presence of a positive label alone (epitope-database tier) is deliberately cheap.
RUBRIC: tuple[QualityAxis, ...] = (
    QualityAxis("candidate_denominator", 0.22, "candidate_universe (#4)",
                "Candidates beyond the vaccinated set: a tested or reconstructed pre-selection universe."),
    QualityAxis("rna_features", 0.14, "tumor_rna (#2)",
                "RNA evidence the target is transcribed: TPM, RNA VAF, mutant RNA reads."),
    QualityAxis("functional_labels", 0.14, "functional_assays / tested_negatives (#6, #7)",
                "Both an Event-B POSITIVE and an explicit TESTED_NEGATIVE — a real ranking problem."),
    QualityAxis("wes_features", 0.12, "tumor_wes (#1)",
                "WES evidence: DNA VAF, cancer-cell fraction/clonality, tumor purity."),
    QualityAxis("prepost_timing", 0.10, "pre_post_assays (#6)",
                "Both PRE- and POST-vaccine assays so Event B is separable from pre-existing Event A."),
    QualityAxis("hla_genotype", 0.08, "hla_genotype (#3)",
                "Patient- and candidate-level HLA calls."),
    QualityAxis("identity_mapping", 0.08, "candidate_identity (#5)",
                "Assay-to-candidate linkage and mutant/wild-type verification."),
    QualityAxis("response_attribution", 0.05, "cd4_cd8_attribution (#6)",
                "Class I / class II attribution on the positive responses."),
    QualityAxis("selection_provenance", 0.04, "manufacturing_selection (#4, #8)",
                "Vaccine-inclusion status, generation provenance, and administration record."),
    QualityAxis("patient_context", 0.03, "patient_context (#9)",
                "Cancer type, stage, treatment context, tumor context."),
)


def _series(frame: pd.DataFrame, column: str) -> pd.Series:
    if column not in frame.columns:
        return pd.Series([], dtype="object")
    return frame[column]


def _present(value: object) -> bool:
    if value is None or value is pd.NA or (isinstance(value, float) and pd.isna(value)):
        return False
    text = str(value).strip()
    return bool(text) and text.lower() not in {"unknown", "nan", "none", "na"}


def _present_fraction(frame: pd.DataFrame, column: str) -> float:
    if column not in frame.columns or frame.empty:
        return 0.0
    return float(frame[column].map(_present).mean())


def _numeric_present_columns(frame: pd.DataFrame, columns: Iterable[str]) -> list[str]:
    hit: list[str] = []
    for column in columns:
        if column not in frame.columns:
            continue
        values = pd.to_numeric(frame[column], errors="coerce")
        if values.notna().any():
            hit.append(column)
    return hit


def _json_alleles(value: object) -> list[str]:
    if not _present(value):
        return []
    text = str(value).strip()
    if text.startswith("["):
        try:
            return [str(item) for item in json.loads(text) if _present(item)]
        except json.JSONDecodeError:
            return []
    return [item.strip() for item in text.replace(";", ",").split(",") if item.strip()]


def _primary_event_b(assays: pd.DataFrame) -> pd.DataFrame:
    if assays.empty or "event_type" not in assays.columns:
        return assays.iloc[0:0]
    mask = assays["event_type"].astype(str).eq(
        BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
    )
    if "candidate_id" in assays.columns:
        mask &= assays["candidate_id"].notna()
    return assays[mask].copy()


def _candidate_linked(assays: pd.DataFrame) -> pd.DataFrame:
    """All candidate-resolved assays, regardless of biological event.

    Pre-vaccine reactivity is recorded as Event A and post-vaccine induction as Event B, so any axis
    that reasons about *timing* must look across both events; filtering to Event B first would make a
    pre-vaccine sample structurally invisible.
    """
    if assays.empty or "candidate_id" not in assays.columns:
        return assays.iloc[0:0]
    return assays[assays["candidate_id"].notna()].copy()


def _stage_informative(value: object) -> bool:
    """True when a funnel stage carries a real reached/dropped call (not 'not_assessed'/'n_a')."""
    norm = str(value).strip().lower()
    if not norm or norm in {"nan", "none", "na"}:
        return False
    return not norm.startswith("not_")


def _has_candidate_alternatives(cands: pd.DataFrame) -> bool:
    """The corpus records candidates the patient could have received beyond the vaccinated set."""
    inclusion = _series(cands, "vaccine_inclusion").astype(str).str.upper()
    return bool((inclusion == "NOT_INCLUDED").any())


def _complete_reconstructed_universe(cands: pd.DataFrame, funnel: pd.DataFrame) -> bool:
    """Strict test for a genuinely reconstructed pre-selection candidate universe.

    A reconstructed denominator is not "some vaccinated candidate reached a middle stage" — the
    included path is always annotated. It requires evidence that the *full* universe was scored and
    that we can see where the non-vaccinated alternatives fell out. We require, conservatively:

      1. candidate alternatives exist (NOT_INCLUDED candidates are present), and
      2. every NOT_INCLUDED candidate carries an assessed middle-stage funnel call (its drop point is
         reconstructed, not merely 'not_assessed'), and
      3. generation provenance is recorded (a candidate-generation run identity proxy).

    Nothing in the current label corpus satisfies (2): NOT_INCLUDED candidates have no assessed
    middle stages, so this returns False for every patient today. A real candidate-generation
    reconstruction that runs the full WES/RNA universe will populate those drops and flip it True.
    """
    if not _has_candidate_alternatives(cands) or funnel is None or funnel.empty:
        return False
    if "generation_provenance" not in cands.columns or not bool(
        _series(cands, "generation_provenance").map(_present).any()
    ):
        return False
    inclusion = _series(cands, "vaccine_inclusion").astype(str).str.upper()
    not_included_ids = set(cands.loc[inclusion == "NOT_INCLUDED", "candidate_id"].astype(str))
    if not not_included_ids or "candidate_id" not in funnel.columns:
        return False
    stages = [stage for stage in MIDDLE_FUNNEL_STAGES if stage in funnel.columns]
    if not stages:
        return False
    drops = funnel[funnel["candidate_id"].astype(str).isin(not_included_ids)]
    if drops.empty:
        return False
    informative = drops[stages].apply(lambda col: col.map(_stage_informative)).any(axis=1)
    covered_ids = set(drops.loc[informative, "candidate_id"].astype(str))
    # Every alternative must have a reconstructed drop point (no unexplained candidate-id gaps).
    return not_included_ids.issubset(covered_ids)


# --- axis scorers -------------------------------------------------------------------------------
# Each scorer receives the patient's candidate rows, primary Event-B assays, patient row, and any
# recognition-evidence rows, and returns a score in [0, 1].


def _score_candidate_denominator(cands, assays, patient, evidence, funnel) -> float:
    # Two distinct concepts, deliberately separated (v1.1): partial credit for *alternatives*
    # (NOT_INCLUDED candidates exist beyond the vaccinated set), full credit only for a *complete
    # reconstructed universe* (the strict test that alone enables full-universe reachability).
    if _complete_reconstructed_universe(cands, funnel):
        return 1.0
    return 0.5 if _has_candidate_alternatives(cands) else 0.0


def _score_rna_features(cands, assays, patient, evidence, funnel) -> float:
    columns = _numeric_present_columns(cands, RNA_CANDIDATE_COLUMNS)
    from_evidence = (
        not evidence.empty
        and "evidence_family" in evidence.columns
        and evidence["evidence_family"].astype(str).isin(RNA_EVIDENCE_FAMILIES).any()
    )
    if not columns and not from_evidence:
        return 0.0
    # Reward breadth of RNA evidence: expression, allele fraction, mutant-read support.
    return float(min(1.0, (len(columns) + int(from_evidence)) / 3.0))


def _score_wes_features(cands, assays, patient, evidence, funnel) -> float:
    columns = _numeric_present_columns(cands, WES_CANDIDATE_COLUMNS)
    from_evidence = (
        not evidence.empty
        and "evidence_family" in evidence.columns
        and evidence["evidence_family"].astype(str).isin(WES_EVIDENCE_FAMILIES).any()
    )
    if not columns and not from_evidence:
        return 0.0
    return float(min(1.0, (len(columns) + int(from_evidence)) / 3.0))


def _score_functional_labels(cands, assays, patient, evidence, funnel) -> float:
    primary = _primary_event_b(assays)
    labels = _series(primary, "response_label").astype(str)
    has_pos = bool((labels == ResponseLabel.POSITIVE.value).any())
    has_neg = bool((labels == ResponseLabel.TESTED_NEGATIVE.value).any())
    return 0.5 * float(has_pos) + 0.5 * float(has_neg)


def _score_prepost_timing(cands, assays, patient, evidence, funnel) -> float:
    # Timing must span events: pre-vaccine reactivity is Event A, post-vaccine induction is Event B.
    # Look at all candidate-linked assays so a paired pre/post design is credited (a prerequisite
    # for isolating vaccine-induced Event B from pre-existing Event A).
    linked = _candidate_linked(assays)
    timing = _series(linked, "relative_to_vaccine").astype(str).str.upper()
    has_pre = bool(timing.str.startswith("PRE").any())
    has_post = bool(timing.str.startswith("POST").any())
    return 0.5 * float(has_pre) + 0.5 * float(has_post)


def _score_hla_genotype(cands, assays, patient, evidence, funnel) -> float:
    patient_hla = 1.0 if _json_alleles(patient.get("hla_alleles")) else 0.0
    if cands.empty:
        candidate_hla = 0.0
    else:
        candidate_hla = float(
            _series(cands, "hla_alleles").map(lambda value: bool(_json_alleles(value))).mean()
        )
    # Patient-level high-resolution genotype is the ideal (item #3); candidate-level HLA still
    # enables peptide->allele presentation mapping, so both are credited.
    return 0.5 * patient_hla + 0.5 * candidate_hla


def _score_identity_mapping(cands, assays, patient, evidence, funnel) -> float:
    primary = _primary_event_b(assays)
    if primary.empty:
        link = 0.0
    else:
        link = float(_series(primary, "candidate_id").map(_present).mean())
    verified = 0.0
    if not cands.empty and "mutant_wildtype_verified" in cands.columns:
        verified = float(
            cands["mutant_wildtype_verified"].astype(str).str.upper().eq("VERIFIED").mean()
        )
    return 0.6 * link + 0.4 * verified


def _score_response_attribution(cands, assays, patient, evidence, funnel) -> float:
    primary = _primary_event_b(assays)
    positives = primary[_series(primary, "response_label").astype(str).eq(ResponseLabel.POSITIVE.value)]
    if positives.empty or cands.empty or "mhc_class" not in cands.columns:
        return 0.0
    class_by_candidate = cands.set_index("candidate_id")["mhc_class"].astype(str)
    known = {MHCClass.CLASS_I.value, MHCClass.CLASS_II.value}
    attributed = positives["candidate_id"].map(class_by_candidate).isin(known)
    return float(attributed.mean())


def _score_selection_provenance(cands, assays, patient, evidence, funnel) -> float:
    if cands.empty:
        return 0.0
    inclusion = float(_series(cands, "vaccine_inclusion").map(_present).mean())
    origin = _present_fraction(cands, "vaccine_inclusion_origin")
    generation = _present_fraction(cands, "generation_provenance")
    return float(np.mean([inclusion, origin, generation]))


def _score_patient_context(cands, assays, patient, evidence, funnel) -> float:
    fields = ("cancer_type", "disease_stage", "treatment_context", "tumor_context")
    return float(np.mean([1.0 if _present(patient.get(field)) else 0.0 for field in fields]))


_SCORERS = {
    "candidate_denominator": _score_candidate_denominator,
    "rna_features": _score_rna_features,
    "functional_labels": _score_functional_labels,
    "wes_features": _score_wes_features,
    "prepost_timing": _score_prepost_timing,
    "hla_genotype": _score_hla_genotype,
    "identity_mapping": _score_identity_mapping,
    "response_attribution": _score_response_attribution,
    "selection_provenance": _score_selection_provenance,
    "patient_context": _score_patient_context,
}


def _tier(scores: dict[str, float], has_pos: bool, has_neg: bool) -> str:
    denom = scores["candidate_denominator"]
    has_features = scores["wes_features"] > 0.0 or scores["rna_features"] > 0.0
    both_labels = has_pos and has_neg
    if not has_pos:
        return "INSUFFICIENT"
    if has_pos and not has_neg:
        return "POSITIVES_ONLY"
    if both_labels and denom >= 1.0 and has_features:
        return "DECISION_PROBLEM_READY"
    if both_labels and denom >= 0.5 and has_features:
        return "PARTIAL_DECISION_PROBLEM"
    if both_labels and has_features:
        return "LABELS_WITH_EVIDENCE"
    return "DISCRIMINATION_ONLY"


def patient_data_quality(corpus: EventBCorpus) -> pd.DataFrame:
    """Score every patient's reconstructability against the ideal-patient rubric."""
    patients = corpus.patients
    if patients.empty:
        return pd.DataFrame(columns=["patient_id", "study_id", "reconstructability_score"])

    candidates = corpus.candidates
    all_assays = corpus.assays
    evidence = corpus.recognition_evidence
    funnel = corpus.candidate_funnel_links

    rows: list[dict] = []
    for _, patient in patients.iterrows():
        patient_id = str(patient["patient_id"])
        cands = candidates[candidates.patient_id.astype(str).eq(patient_id)] if not candidates.empty else candidates
        # Pass every assay for the patient (any biological event); each scorer narrows to the event
        # semantics it needs. Timing spans Event A (pre) and Event B (post); labels use Event B only.
        pat_assays = all_assays[all_assays.patient_id.astype(str).eq(patient_id)] if not all_assays.empty else all_assays
        pat_primary = _primary_event_b(pat_assays)
        pat_evidence = evidence[evidence.patient_id.astype(str).eq(patient_id)] if ("patient_id" in evidence.columns and not evidence.empty) else evidence.iloc[0:0]
        pat_funnel = funnel[funnel.patient_id.astype(str).eq(patient_id)] if ("patient_id" in funnel.columns and not funnel.empty) else funnel.iloc[0:0]

        scores = {
            key: float(np.clip(scorer(cands, pat_assays, patient, pat_evidence, pat_funnel), 0.0, 1.0))
            for key, scorer in _SCORERS.items()
        }
        labels = _series(pat_primary, "response_label").astype(str)
        has_pos = bool((labels == ResponseLabel.POSITIVE.value).any())
        has_neg = bool((labels == ResponseLabel.TESTED_NEGATIVE.value).any())

        weighted = float(sum(scores[axis.key] * axis.weight for axis in RUBRIC))
        tier = _tier(scores, has_pos, has_neg)
        timing = _series(_candidate_linked(pat_assays), "relative_to_vaccine").astype(str).str.upper()
        # Full-universe reachability requires the strict reconstructed universe, not merely denom>=1.0.
        universe_complete = _complete_reconstructed_universe(cands, pat_funnel)
        immunopeptidomics = bool(
            "evidence_family" in pat_evidence.columns
            and pat_evidence["evidence_family"].astype(str).isin(IMMUNOPEPTIDOMICS_FAMILIES).any()
        )
        supports = {
            "supports_vaccinated_subset_discrimination": has_pos and has_neg,
            "supports_full_universe_reachability": universe_complete,
            "supports_feature_conditioned_recognition": (has_pos and has_neg)
            and (scores["wes_features"] > 0.0 or scores["rna_features"] > 0.0),
            "supports_event_b_isolation": bool(timing.str.startswith("PRE").any())
            and bool(timing.str.startswith("POST").any()),
        }
        gap = [
            axis.key
            for axis in sorted(RUBRIC, key=lambda a: a.weight, reverse=True)
            if scores[axis.key] < 0.5
        ]
        rows.append(
            {
                "patient_id": patient_id,
                "study_id": str(patient.get("study_id", "")),
                "cancer_type": str(patient.get("cancer_type", "")),
                "n_candidates": int(len(cands)),
                "n_positive": int((labels == ResponseLabel.POSITIVE.value).sum()),
                "n_tested_negative": int((labels == ResponseLabel.TESTED_NEGATIVE.value).sum()),
                **{f"axis_{key}": value for key, value in scores.items()},
                "reconstructability_score": round(weighted, 4),
                "reconstruction_tier": tier,
                **{key: bool(value) for key, value in supports.items()},
                "diag_immunopeptidomics_present": immunopeptidomics,
                "reconstruction_gap": ", ".join(gap),
            }
        )
    frame = pd.DataFrame(rows)
    return frame.sort_values("reconstructability_score", ascending=False, kind="mergesort").reset_index(drop=True)


def study_data_quality(patient_frame: pd.DataFrame) -> pd.DataFrame:
    """Roll patient scorecards up to per-study readiness."""
    if patient_frame.empty:
        return pd.DataFrame(columns=["study_id", "patient_n", "mean_score"])
    rows = []
    for study_id, group in patient_frame.groupby("study_id", sort=True):
        rows.append(
            {
                "study_id": study_id,
                "patient_n": int(len(group)),
                "mean_score": round(float(group["reconstructability_score"].mean()), 4),
                "max_score": round(float(group["reconstructability_score"].max()), 4),
                "best_tier": group.sort_values("reconstructability_score", ascending=False)
                .iloc[0]["reconstruction_tier"],
                "any_reconstructed_denominator": bool((group["axis_candidate_denominator"] >= 1.0).any()),
                "any_wes_or_rna": bool(
                    ((group["axis_wes_features"] > 0.0) | (group["axis_rna_features"] > 0.0)).any()
                ),
                "discrimination_ready_n": int(group["supports_vaccinated_subset_discrimination"].sum()),
            }
        )
    return pd.DataFrame(rows).sort_values("mean_score", ascending=False, kind="mergesort").reset_index(drop=True)


def _auroc(scores: np.ndarray, labels: np.ndarray) -> float:
    """Rank-based AUROC with tie handling (Mann-Whitney). Chance = 0.5."""
    positives = labels == 1
    negatives = labels == 0
    n_pos = int(positives.sum())
    n_neg = int(negatives.sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = pd.Series(scores).rank(method="average").to_numpy()
    rank_sum_pos = float(order[positives].sum())
    u = rank_sum_pos - n_pos * (n_pos + 1) / 2.0
    return float(u / (n_pos * n_neg))


def _cross_fitted_study_scores(primary: pd.DataFrame) -> np.ndarray:
    """Leave-one-patient-out study base rate for each row.

    Each candidate is scored with its study's positive rate computed *without that candidate's
    patient*, so the score has not seen the labels it is later evaluated against. Rows whose study
    has no other patient (the rate would be undefined) are returned as NaN and excluded downstream.
    """
    study = primary.study_id.astype(str)
    patient = primary.patient_id.astype(str)
    y = primary.y.to_numpy(dtype=float)
    study_pos = study.map(primary.groupby(study).y.sum())
    study_n = study.map(study.value_counts())
    pat_key = study + "\x1f" + patient
    pat_pos = pat_key.map(primary.assign(_k=pat_key.to_numpy()).groupby("_k").y.sum())
    pat_n = pat_key.map(pat_key.value_counts())
    num = study_pos.to_numpy(dtype=float) - pat_pos.to_numpy(dtype=float)
    den = study_n.to_numpy(dtype=float) - pat_n.to_numpy(dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        scores = np.where(den > 0, num / den, np.nan)
    return scores


def study_identity_confound(corpus: EventBCorpus) -> dict:
    """Quantify how much of the label is explained by study identity alone.

    ``study_shortcut_auroc`` (the apparent, in-sample form) is the numeric shape of "good at
    matching studies, bad at predicting recognition": score each candidate by its study's base
    positive rate and you separate positives from tested negatives without learning any biology.
    Because that rate is computed from the very labels it is scored against, it is optimistic;
    ``study_shortcut_auroc_cross_fitted`` recomputes each study rate with the scored patient held
    out (leave-one-patient-out), which is the honest confounding index. Both are lower bounds on the
    exploitable shortcut, since finer within-study cues also leak study identity.
    """
    primary = _primary_event_b(corpus.assays)
    primary = primary[
        primary.response_label.astype(str).isin(
            [ResponseLabel.POSITIVE.value, ResponseLabel.TESTED_NEGATIVE.value]
        )
    ].copy()
    result = {
        "version": DATA_QUALITY_VERSION,
        "candidate_resolved_primary_labels": int(len(primary)),
        "n_studies_with_primary_labels": 0,
        "n_studies_with_both_labels": 0,
        "overall_positive_rate": None,
        "positive_rate_spread": None,
        "study_shortcut_auroc": None,
        "study_shortcut_auroc_cross_fitted": None,
        "cross_fit_coverage": None,
        "per_study_positive_rate": {},
        "repeated_antigen_fraction": None,
        "interpretation": "no candidate-resolved primary labels",
    }
    if primary.empty:
        return result

    primary["y"] = primary.response_label.astype(str).eq(ResponseLabel.POSITIVE.value).astype(int)
    per_study = primary.groupby(primary.study_id.astype(str)).y
    rate = per_study.mean()
    both = per_study.agg(lambda s: s.nunique() == 2)
    scores = primary.study_id.astype(str).map(rate).to_numpy(dtype=float)
    auroc = _auroc(scores, primary.y.to_numpy())

    cf_scores = _cross_fitted_study_scores(primary)
    y = primary.y.to_numpy()
    valid = ~np.isnan(cf_scores)
    cf_auroc = _auroc(cf_scores[valid], y[valid]) if valid.any() else float("nan")
    cf_coverage = float(valid.mean())

    peptides = _series(corpus.candidates, "mutant_peptide").astype(str)
    per_patient = corpus.candidates.groupby(peptides).patient_id.nunique() if not corpus.candidates.empty else pd.Series(dtype=int)
    repeated_fraction = (
        float((per_patient > 1).reindex(peptides).fillna(False).mean()) if len(peptides) else 0.0
    )

    # The cross-fitted index is the defensible one; fall back to apparent when it is undefined.
    headline = cf_auroc if not np.isnan(cf_auroc) else auroc
    result.update(
        {
            "n_studies_with_primary_labels": int(rate.shape[0]),
            "n_studies_with_both_labels": int(both.sum()),
            "overall_positive_rate": round(float(primary.y.mean()), 4),
            "positive_rate_spread": round(float(rate.max() - rate.min()), 4),
            "study_shortcut_auroc": None if np.isnan(auroc) else round(auroc, 4),
            "study_shortcut_auroc_cross_fitted": None if np.isnan(cf_auroc) else round(cf_auroc, 4),
            "cross_fit_coverage": round(cf_coverage, 4),
            "per_study_positive_rate": {study: round(float(value), 4) for study, value in rate.items()},
            "repeated_antigen_fraction": round(repeated_fraction, 4),
            "interpretation": _confound_interpretation(headline),
        }
    )
    return result


def _confound_interpretation(auroc: float) -> str:
    if np.isnan(auroc):
        return "only one label present; study shortcut undefined"
    if auroc >= 0.75 or auroc <= 0.25:
        return "STRONG study-identity shortcut: study base rate alone separates the labels"
    if auroc >= 0.65 or auroc <= 0.35:
        return "MODERATE study-identity shortcut present"
    return "WEAK study-identity shortcut: study base rate is close to chance"


def corpus_data_quality_summary(corpus: EventBCorpus, patient_frame: pd.DataFrame | None = None) -> dict:
    if patient_frame is None:
        patient_frame = patient_data_quality(corpus)
    confound = study_identity_confound(corpus)
    tier_counts = (
        patient_frame["reconstruction_tier"].value_counts().to_dict() if not patient_frame.empty else {}
    )
    excellent = (
        not patient_frame.empty
        and (patient_frame["reconstruction_tier"] == "DECISION_PROBLEM_READY").any()
    )
    return {
        "version": DATA_QUALITY_VERSION,
        "patient_n": int(len(patient_frame)),
        "mean_reconstructability_score": round(float(patient_frame["reconstructability_score"].mean()), 4)
        if not patient_frame.empty
        else 0.0,
        "any_excellent_patient": bool(excellent),
        "patients_with_reconstructed_denominator": int(
            (patient_frame.get("axis_candidate_denominator", pd.Series(dtype=float)) >= 1.0).sum()
        ),
        "patients_with_wes_or_rna_features": int(
            (
                (patient_frame.get("axis_wes_features", pd.Series(dtype=float)) > 0.0)
                | (patient_frame.get("axis_rna_features", pd.Series(dtype=float)) > 0.0)
            ).sum()
        ),
        "patients_discrimination_ready": int(
            patient_frame.get("supports_vaccinated_subset_discrimination", pd.Series(dtype=bool)).sum()
        ),
        "tier_counts": tier_counts,
        "study_identity_confound": confound,
    }


def render_data_quality_markdown(
    summary: dict, patient_frame: pd.DataFrame, study_frame: pd.DataFrame
) -> str:
    confound = summary["study_identity_confound"]
    lines = [
        "# Data-quality and reconstructability audit",
        "",
        f"Version `{summary['version']}`. No model is fitted; axes score only source-present evidence.",
        "",
        "## Headline",
        "",
        f"- Patients scored: **{summary['patient_n']}**",
        f"- Mean reconstructability score: **{summary['mean_reconstructability_score']:.3f}** (0-1)",
        f"- Any 'excellent' (decision-problem-ready) patient: **{summary['any_excellent_patient']}**",
        f"- Patients with a reconstructed candidate denominator: **{summary['patients_with_reconstructed_denominator']}**",
        f"- Patients with any WES/RNA feature: **{summary['patients_with_wes_or_rna_features']}**",
        f"- Patients usable for vaccinated-subset discrimination: **{summary['patients_discrimination_ready']}**",
        "",
        "## Reconstruction tiers",
        "",
    ]
    for tier, count in sorted(summary["tier_counts"].items(), key=lambda kv: -kv[1]):
        lines.append(f"- `{tier}`: {count}")
    lines += [
        "",
        "## Study-identity confound (good-at-matching / bad-at-predicting risk)",
        "",
        f"- `study_shortcut_auroc_cross_fitted`: **{confound['study_shortcut_auroc_cross_fitted']}** "
        f"(leave-one-patient-out; the honest confounding index)",
        f"- `study_shortcut_auroc` (apparent, in-sample): {confound['study_shortcut_auroc']} "
        f"(0.5 = chance; higher = study identity alone predicts the label)",
        f"- `cross_fit_coverage`: {confound['cross_fit_coverage']}",
        f"- `positive_rate_spread`: {confound['positive_rate_spread']}",
        f"- `overall_positive_rate`: {confound['overall_positive_rate']}",
        f"- `n_studies_with_both_labels`: {confound['n_studies_with_both_labels']}",
        f"- `repeated_antigen_fraction`: {confound['repeated_antigen_fraction']}",
        f"- Interpretation: {confound['interpretation']}",
        "",
        "### Per-study positive rate",
        "",
    ]
    for study, value in confound["per_study_positive_rate"].items():
        lines.append(f"- `{study}`: {value}")
    lines += ["", "## Per-study reconstructability", "",
              "| Study | Patients | Mean | Max | Best tier | WES/RNA | Denominator |",
              "|---|---:|---:|---:|---|---|---|"]
    for _, row in study_frame.iterrows():
        lines.append(
            f"| {row['study_id']} | {row['patient_n']} | {row['mean_score']:.3f} | "
            f"{row['max_score']:.3f} | {row['best_tier']} | "
            f"{'yes' if row['any_wes_or_rna'] else 'no'} | "
            f"{'yes' if row['any_reconstructed_denominator'] else 'no'} |"
        )
    lines += ["", "## Top patients by reconstructability", "",
              "| Patient | Study | Score | Tier | Gap (highest-value missing) |",
              "|---|---|---:|---|---|"]
    for _, row in patient_frame.head(15).iterrows():
        lines.append(
            f"| {row['patient_id']} | {row['study_id']} | {row['reconstructability_score']:.3f} | "
            f"{row['reconstruction_tier']} | {row['reconstruction_gap']} |"
        )
    lines.append("")
    return "\n".join(lines)
