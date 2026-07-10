"""Conservative corpus sufficiency and leakage-safe split feasibility audits."""

from __future__ import annotations

import json
from typing import Iterable

import pandas as pd

from event_b.corpus import EventBCorpus
from event_b.models import BiologicalEvent, MHCClass, ResponseLabel
from event_b.registry import StudyRegistry, StudyStatus
from event_b.review import ReviewIssue
from event_b.splits import SplitType, generate_split_manifest


VERDICTS = {
    "INSUFFICIENT_PUBLIC_EVENT_B_DATA",
    "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA",
    "PUBLIC_EVENT_B_MINIMUM_MET_BUT_STUDY_HOLDOUT_NOT_VIABLE",
    "PUBLIC_EVENT_B_BACKBONE_VALIDATED_READY_FOR_BASELINE_EXPERIMENTS",
    "PUBLIC_EVENT_B_DATA_HETEROGENEOUS_REQUIRES_STRATIFIED_MODELLING",
    "PUBLIC_EVENT_B_CORPUS_BLOCKED_BY_LABEL_COMPARABILITY",
}


def _json_values(value: object) -> list[str]:
    if value is None or value is pd.NA or (isinstance(value, float) and pd.isna(value)):
        return []
    text = str(value).strip()
    if not text:
        return []
    if text.startswith("["):
        try:
            return [str(item) for item in json.loads(text)]
        except json.JSONDecodeError:
            return []
    return [item.strip() for item in text.replace(";", ",").split(",") if item.strip()]


def _primary_labels(corpus: EventBCorpus) -> pd.DataFrame:
    assays = corpus.assays
    primary = assays[
        assays.event_type.astype(str).eq(
            BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
        )
        & assays.candidate_id.notna()
    ].copy()
    if primary.empty:
        return primary
    candidate_columns = [
        "candidate_id",
        "mutant_peptide",
        "hla_alleles",
        "mhc_class",
    ]
    patient_columns = ["patient_id", "cancer_type"]
    merged = primary.merge(
        corpus.candidates[candidate_columns], on="candidate_id", how="left", validate="one_to_one"
    ).merge(
        corpus.patients[patient_columns], on="patient_id", how="left", validate="many_to_one"
    )
    if "sample_date" not in merged:
        merged["sample_date"] = pd.NA
    return merged


def _fraction(values: pd.Series, denominator: int) -> float:
    return float(values.max() / denominator) if denominator and len(values) else 0.0


def _split_has_both_labels(assignments: pd.DataFrame, primary: pd.DataFrame) -> tuple[bool, dict]:
    joined = assignments[["candidate_id", "split"]].merge(
        primary[["candidate_id", "patient_id", "response_label"]],
        on="candidate_id",
        how="inner",
        validate="one_to_one",
    )
    counts = (
        joined.groupby(["split", "response_label"]).size().unstack(fill_value=0).to_dict("index")
    )
    patient_counts = (
        joined[joined.response_label.eq(ResponseLabel.POSITIVE.value)]
        .groupby("split")
        .patient_id.nunique()
        .to_dict()
    )
    usable = True
    for split in ("train", "evaluation"):
        split_counts = counts.get(split, {})
        usable &= split_counts.get(ResponseLabel.POSITIVE.value, 0) > 0
        usable &= split_counts.get(ResponseLabel.TESTED_NEGATIVE.value, 0) > 0
        usable &= patient_counts.get(split, 0) > 0
    return bool(usable), {"label_counts": counts, "positive_patient_counts": patient_counts}


def split_feasibility(corpus: EventBCorpus, registry: StudyRegistry) -> dict:
    primary = _primary_labels(corpus)
    eligible = primary[
        primary.response_label.isin(
            [ResponseLabel.POSITIVE.value, ResponseLabel.TESTED_NEGATIVE.value]
        )
    ].copy()
    base = eligible[
        [
            "candidate_id",
            "patient_id",
            "study_id",
            "mutant_peptide",
            "hla_alleles",
            "cancer_type",
            "sample_date",
        ]
    ].drop_duplicates("candidate_id")
    report = {}

    def search(split_type: SplitType, *, precondition=True, reason="insufficient metadata"):
        if not precondition:
            report[split_type.value] = {"feasible": False, "reason": reason}
            return
        last_error = None
        for seed in range(64):
            try:
                manifest = generate_split_manifest(base, split_type, seed=seed)
                usable, details = _split_has_both_labels(pd.DataFrame(manifest.assignments), eligible)
                if usable:
                    report[split_type.value] = {
                        "feasible": True,
                        "seed": seed,
                        **details,
                    }
                    return
            except ValueError as error:
                last_error = str(error)
        report[split_type.value] = {
            "feasible": False,
            "reason": last_error or "no deterministic seed leaves both labels on both sides",
        }

    search(SplitType.PATIENT_HOLDOUT, precondition=eligible.patient_id.nunique() >= 2)
    search(SplitType.STUDY_HOLDOUT, precondition=eligible.study_id.nunique() >= 2)
    known_hla = base.hla_alleles.map(_json_values)
    hla_groups = {allele for values in known_hla for allele in values}
    search(
        SplitType.HLA_HOLDOUT,
        precondition=len(hla_groups) >= 2,
        reason="fewer than two source-resolved HLA groups",
    )
    search(SplitType.PEPTIDE_CLUSTER_HOLDOUT, precondition=base.mutant_peptide.nunique() >= 2)
    search(SplitType.CANCER_TYPE_HOLDOUT, precondition=base.cancer_type.nunique() >= 2)
    valid_dates = pd.to_datetime(base.sample_date, errors="coerce", utc=True).dropna()
    if valid_dates.nunique() >= 2:
        cutoff = valid_dates.sort_values().iloc[len(valid_dates) // 2].isoformat()
        try:
            manifest = generate_split_manifest(
                base,
                SplitType.TEMPORAL_HOLDOUT,
                temporal_cutoff=cutoff,
            )
            usable, details = _split_has_both_labels(pd.DataFrame(manifest.assignments), eligible)
            report[SplitType.TEMPORAL_HOLDOUT.value] = {
                "feasible": usable,
                "temporal_cutoff": cutoff,
                **details,
            }
        except ValueError as error:
            report[SplitType.TEMPORAL_HOLDOUT.value] = {
                "feasible": False,
                "reason": str(error),
            }
    else:
        report[SplitType.TEMPORAL_HOLDOUT.value] = {
            "feasible": False,
            "reason": "fewer than two source-resolved candidate dates",
        }

    design = {entry.canonical_study_id: entry.antigen_design for entry in registry.studies}
    shared = eligible[eligible.study_id.map(design).eq("SHARED")].copy()
    shared_result = {"feasible": False, "reason": "no shared-antigen candidate labels"}
    if not shared.empty:
        for peptide in sorted(shared.mutant_peptide.unique()):
            assignments = eligible[["candidate_id"]].copy()
            evaluation_ids = set(shared.loc[shared.mutant_peptide.eq(peptide), "candidate_id"])
            assignments["split"] = assignments.candidate_id.map(
                lambda value: "evaluation" if value in evaluation_ids else "train"
            )
            usable, details = _split_has_both_labels(assignments, eligible)
            if usable:
                shared_result = {
                    "feasible": True,
                    "evaluation_group": peptide,
                    **details,
                }
                break
        if not shared_result["feasible"]:
            shared_result["reason"] = "no shared-antigen group leaves both labels on both sides"
    report["SHARED_ANTIGEN_GROUP_HOLDOUT"] = shared_result
    return report


def sufficiency_audit(
    corpus: EventBCorpus,
    registry: StudyRegistry,
    issues: Iterable[ReviewIssue] = (),
) -> dict:
    assays = corpus.assays.copy()
    primary = _primary_labels(corpus)
    event_b = assays.event_type.astype(str).eq(
        BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
    )
    positive = assays.response_label.astype(str).eq(ResponseLabel.POSITIVE.value)
    negative = assays.response_label.astype(str).eq(ResponseLabel.TESTED_NEGATIVE.value)
    event_a = assays.event_type.astype(str).eq(
        BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value
    )
    spreading = assays.event_type.astype(str).eq(BiologicalEvent.EPITOPE_SPREADING.value)

    primary_positive = primary.response_label.astype(str).eq(ResponseLabel.POSITIVE.value)
    primary_negative = primary.response_label.astype(str).eq(ResponseLabel.TESTED_NEGATIVE.value)
    primary_untested = primary.response_label.astype(str).eq(ResponseLabel.UNTESTED.value)
    event_b_patients = set(assays.loc[event_b & (positive | negative), "patient_id"].astype(str))
    positive_patients = set(assays.loc[event_b & positive, "patient_id"].astype(str))
    negative_patients = set(assays.loc[event_b & negative, "patient_id"].astype(str))
    event_b_studies = set(assays.loc[event_b & (positive | negative), "study_id"].astype(str))
    # Peptide-ranking evidence tier: only Event-B assays mapped to a specific candidate can
    # train "which peptide works". Patient-level-only cohorts (e.g. no per-pool identity) are
    # reported for eligibility/abstention questions but never counted toward this tier.
    candidate_resolved_patients = set(primary.patient_id.astype(str))
    candidate_resolved_positive_patients = set(
        primary.loc[primary_positive, "patient_id"].astype(str)
    )
    candidate_resolved_studies = set(primary.study_id.astype(str))
    patient_level_only_patients = event_b_patients - candidate_resolved_patients

    publication_ids = set()
    for value in corpus.studies.publication_ids:
        publication_ids.update(_json_values(value))

    design = {entry.canonical_study_id: entry.antigen_design for entry in registry.studies}
    status = {entry.canonical_study_id: entry.ingestion_status.value for entry in registry.studies}
    primary_design = primary.study_id.map(design).fillna("UNKNOWN")

    per_study = []
    all_study_ids = sorted(
        set(corpus.studies.study_id.astype(str))
        | {entry.canonical_study_id for entry in registry.studies}
    )
    for study_id in all_study_ids:
        study_assays = assays[assays.study_id.astype(str).eq(study_id)]
        study_event_b = study_assays.event_type.astype(str).eq(
            BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value
        )
        study_primary = primary[primary.study_id.astype(str).eq(study_id)]
        per_study.append(
            {
                "study_id": study_id,
                "ingestion_status": status.get(study_id, "ACCEPTED_LEGACY"),
                "patient_n": int(
                    corpus.patients.loc[
                        corpus.patients.study_id.astype(str).eq(study_id), "patient_id"
                    ].nunique()
                ),
                "event_b_positive_patient_n": int(
                    study_assays.loc[
                        study_event_b
                        & study_assays.response_label.astype(str).eq(ResponseLabel.POSITIVE.value),
                        "patient_id",
                    ].nunique()
                ),
                "patients_with_tested_negative_n": int(
                    study_assays.loc[
                        study_event_b
                        & study_assays.response_label.astype(str).eq(
                            ResponseLabel.TESTED_NEGATIVE.value
                        ),
                        "patient_id",
                    ].nunique()
                ),
                "primary_label_n": int(len(study_primary)),
                "primary_positive_n": int(
                    study_primary.response_label.astype(str).eq(ResponseLabel.POSITIVE.value).sum()
                ),
                "primary_tested_negative_n": int(
                    study_primary.response_label.astype(str)
                    .eq(ResponseLabel.TESTED_NEGATIVE.value)
                    .sum()
                ),
                "primary_untested_n": int(
                    study_primary.response_label.astype(str).eq(ResponseLabel.UNTESTED.value).sum()
                ),
                "assay_observation_n": int(len(study_assays)),
                "patient_level_only_event_b_n": int(
                    (study_event_b & study_assays.candidate_id.isna()).sum()
                ),
            }
        )

    study_frame = pd.DataFrame(per_study)
    primary_n = len(primary)
    event_b_observation_n = int(event_b.sum())
    positive_n = int(primary_positive.sum())
    negative_n = int(primary_negative.sum())
    dominance = {
        "largest_study_patient_fraction": _fraction(
            study_frame.set_index("study_id").patient_n, int(corpus.patients.patient_id.nunique())
        ),
        "largest_study_observation_fraction": _fraction(
            study_frame.set_index("study_id").assay_observation_n, len(assays)
        ),
        "largest_study_primary_label_fraction": _fraction(
            study_frame.set_index("study_id").primary_label_n, primary_n
        ),
        "largest_study_positive_fraction": _fraction(
            study_frame.set_index("study_id").primary_positive_n, positive_n
        ),
        "largest_study_negative_fraction": _fraction(
            study_frame.set_index("study_id").primary_tested_negative_n, negative_n
        ),
    }

    repeated_peptides = primary.groupby("mutant_peptide").candidate_id.nunique()
    split_report = split_feasibility(corpus, registry)
    # The registered gate measures peptide-ranking readiness, so it counts candidate-resolved
    # patients/studies/positives only. A patient-level-only cohort can lift the headline patient
    # count without adding a single peptide label, so it must not be able to satisfy this gate.
    minimum_met = (
        len(candidate_resolved_patients) >= 100
        and len(candidate_resolved_studies) >= 2
        and len(candidate_resolved_positive_patients) >= 30
    )
    study_holdout = split_report[SplitType.STUDY_HOLDOUT.value]["feasible"]
    no_overwhelming_dominance = dominance["largest_study_primary_label_fraction"] < 0.7
    class_counts = primary.mhc_class.astype(str).value_counts().to_dict()
    class_coverage = bool(
        class_counts.get(MHCClass.CLASS_I.value, 0)
        and class_counts.get(MHCClass.CLASS_II.value, 0)
    )
    explicit_negatives = negative_n > 0
    label_comparable_primary = bool(primary_n and primary.candidate_id.notna().all())

    if not minimum_met:
        verdict = "INSUFFICIENT_CANDIDATE_RESOLVED_PUBLIC_EVENT_B_DATA"
    elif not study_holdout:
        verdict = "PUBLIC_EVENT_B_MINIMUM_MET_BUT_STUDY_HOLDOUT_NOT_VIABLE"
    elif not label_comparable_primary:
        verdict = "PUBLIC_EVENT_B_CORPUS_BLOCKED_BY_LABEL_COMPARABILITY"
    elif not no_overwhelming_dominance or not class_coverage:
        verdict = "PUBLIC_EVENT_B_DATA_HETEROGENEOUS_REQUIRES_STRATIFIED_MODELLING"
    else:
        verdict = "PUBLIC_EVENT_B_BACKBONE_VALIDATED_READY_FOR_BASELINE_EXPERIMENTS"
    assert verdict in VERDICTS

    return {
        "verdict": verdict,
        "required_counts": {
            "total_studies": int(corpus.studies.study_id.nunique()),
            "event_b_studies": len(event_b_studies),
            "publications": len(publication_ids),
            "unique_patients": int(corpus.patients.patient_id.nunique()),
            "event_b_patients": len(event_b_patients),
            "event_b_positive_patients": len(positive_patients),
            "candidate_resolved_patient_n": len(candidate_resolved_patients),
            "candidate_resolved_positive_patient_n": len(candidate_resolved_positive_patients),
            "candidate_resolved_study_n": len(candidate_resolved_studies),
            "patient_level_only_patient_n": len(patient_level_only_patients),
            "patients_with_explicit_tested_negatives": len(negative_patients),
            "vaccine_components": int(corpus.candidates.candidate_id.nunique()),
            "unique_patient_candidate_pairs": int(
                corpus.candidates[["patient_id", "candidate_id"]].drop_duplicates().shape[0]
            ),
            "primary_candidate_labels": primary_n,
            "assay_observations": int(len(assays)),
            "event_b_assay_observations": event_b_observation_n,
            "event_b_positives": positive_n,
            "event_b_tested_negatives": negative_n,
            "event_b_untested_candidates": int(primary_untested.sum()),
            "event_a_observations": int(event_a.sum()),
            "epitope_spreading_observations": int(spreading.sum()),
            "class_i_observations": int(class_counts.get(MHCClass.CLASS_I.value, 0)),
            "class_ii_observations": int(class_counts.get(MHCClass.CLASS_II.value, 0)),
            "unknown_class_observations": int(
                class_counts.get(MHCClass.UNKNOWN.value, 0)
                + class_counts.get(MHCClass.BOTH.value, 0)
            ),
            "cd4_responses": int(
                (primary_positive & primary.mhc_class.astype(str).eq(MHCClass.CLASS_II.value)).sum()
            ),
            "cd8_responses": int(
                (primary_positive & primary.mhc_class.astype(str).eq(MHCClass.CLASS_I.value)).sum()
            ),
            "shared_antigen_observations": int(primary_design.eq("SHARED").sum()),
            "personalized_antigen_observations": int(primary_design.eq("PERSONALIZED").sum()),
            "review_queue_size": len(tuple(issues)),
        },
        "per_study": per_study,
        "independence_and_dominance": {
            "patient_n": int(corpus.patients.patient_id.nunique()),
            "study_n": int(corpus.studies.study_id.nunique()),
            "peptide_n": int(corpus.candidates.mutant_peptide.nunique()),
            "primary_label_n": primary_n,
            "repeated_antigen_primary_label_fraction": float(
                repeated_peptides[repeated_peptides > 1].sum() / primary_n
            )
            if primary_n
            else 0.0,
            "publication_overlap_entries": [
                {
                    "study_id": entry.canonical_study_id,
                    "overlap": list(entry.known_publication_overlap),
                }
                for entry in registry.studies
                if entry.known_publication_overlap
            ],
            "independent_cohort_n": sum(
                entry.ingestion_status is StudyStatus.ACCEPTED
                for entry in registry.studies
            ),
            **dominance,
        },
        "coverage": {
            "cancer_type": corpus.patients.cancer_type.astype(str).value_counts().to_dict(),
            "vaccine_platform": corpus.studies.vaccine_platform.astype(str).value_counts().to_dict(),
            "assay_type": assays.assay_type.astype(str).value_counts().to_dict(),
            "hla_allele": pd.Series(
                [allele for value in primary.hla_alleles for allele in _json_values(value)],
                dtype="object",
            ).value_counts().to_dict(),
            "mhc_class": class_counts,
            "antigen_design": primary_design.value_counts().to_dict(),
            "baseline_verified_studies": ["braun_rcc_2025", "mkras_vax_2026", "pdac_neovax_2023"],
            "baseline_author_asserted_studies": ["hu_neovax_2021"],
            "baseline_unresolved_studies": ["nous_209_2025"],
            "raw_measurement_label_studies": ["braun_rcc_2025"],
            "author_call_label_studies": [
                "hu_neovax_2021",
                "mkras_vax_2026",
                "pdac_neovax_2023",
                "nous_209_2025",
            ],
        },
        "label_quality": {
            "explicit_candidate_negatives": negative_n,
            "inferred_negatives": 0,
            "candidate_resolved_primary_assays": primary_n,
            "pool_or_patient_level_only_event_b_assays": int(
                (event_b & assays.candidate_id.isna()).sum()
            ),
            "patient_level_only_responses": int(
                (
                    event_b
                    & assays.candidate_id.isna()
                    & assays.source_interpretation.astype(str).str.contains("patient-level", na=False)
                ).sum()
            ),
            "ambiguous_candidate_labels": int(primary_untested.sum()),
            "contradictions": sum(issue.code == "CONTRADICTORY_LABELS" for issue in issues),
            "missing_provenance": sum(issue.code == "MISSING_PROVENANCE" for issue in issues),
            "unresolved_patient_mapping": sum(
                entry.ingestion_status is StudyStatus.BLOCKED_PATIENT_MAPPING
                or "patient map" in (entry.current_blocker or "").lower()
                for entry in registry.studies
            ),
        },
        "split_feasibility": split_report,
        "registered_minimum": {
            "thresholds": {
                "candidate_resolved_patients": 100,
                "candidate_resolved_studies": 2,
                "candidate_resolved_positive_patients": 30,
            },
            "met": minimum_met,
            "study_holdout_feasible": study_holdout,
            "no_overwhelming_primary_label_dominance": no_overwhelming_dominance,
            "explicit_tested_negatives_available": explicit_negatives,
            "class_i_and_class_ii_coverage": class_coverage,
            "label_comparability_for_primary_candidate_set": label_comparable_primary,
        },
    }


def render_sufficiency_markdown(audit: dict) -> str:
    counts = audit["required_counts"]
    dominance = audit["independence_and_dominance"]
    minimum = audit["registered_minimum"]
    lines = [
        "# Public Event-B data sufficiency audit",
        "",
        f"**Verdict: {audit['verdict']}**",
        "",
        "No recognition model was fitted. The registered gate is evaluated on candidate-resolved "
        "(peptide-level) patients only; patient-level-only cohorts are reported but never counted "
        "toward peptide-ranking readiness. CD4/CD8 counts below use source-resolved candidate class, "
        "not inferred cellular phenotype.",
        "",
        "## Evidence tiers (peptide-ranking sample vs. total)",
        "",
        f"- `total_event_b_patient_n`: {counts['event_b_patients']}",
        f"- `candidate_resolved_patient_n`: {counts['candidate_resolved_patient_n']}",
        f"- `patient_level_only_patient_n`: {counts['patient_level_only_patient_n']}",
        f"- `candidate_resolved_positive_patient_n`: {counts['candidate_resolved_positive_patient_n']}",
        f"- `candidate_resolved_study_n`: {counts['candidate_resolved_study_n']}",
        f"- `candidate_level_primary_label_n`: {counts['primary_candidate_labels']}",
        "",
        "## Global counts",
        "",
    ]
    for key, value in counts.items():
        lines.append(f"- `{key}`: {value}")
    lines += ["", "## Registered minimum", ""]
    for key, value in minimum.items():
        lines.append(f"- `{key}`: {value}")
    lines += ["", "## Dominance", ""]
    for key, value in dominance.items():
        if key.startswith("largest_study"):
            lines.append(f"- `{key}`: {value:.4f}")
    lines += ["", "## Per-study status", ""]
    for row in audit["per_study"]:
        lines.append(
            f"- `{row['study_id']}`: {row['ingestion_status']}; patients={row['patient_n']}; "
            f"primary={row['primary_label_n']} (+{row['primary_positive_n']}/"
            f"-{row['primary_tested_negative_n']}/untested={row['primary_untested_n']}); "
            f"patient-level-only={row['patient_level_only_event_b_n']}"
        )
    lines += ["", "## Split feasibility", ""]
    for split, row in audit["split_feasibility"].items():
        lines.append(f"- `{split}`: {'viable' if row['feasible'] else 'not viable'}")
    lines.append("")
    return "\n".join(lines)
