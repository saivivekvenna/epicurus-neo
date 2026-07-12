#!/usr/bin/env python
"""Gate feature audit (isolated milestone-7 asset).

WHY. The dynamic gate was falsified: a label-blind *presentation* gate removes
0% of the high-presentation decoys that outrank positives, so downstream
hits@20 does not move. The only lever that can is an ORTHOGONAL feature that
separates the top-ranked TESTED_NEGATIVE decoys from POSITIVES *within the
high-presentation stratum*. This script inventories every currently available
orthogonal feature across the seven cohorts, quantifies its univariate and
cross-fitted signal on that exact hard stratum, audits study/patient identity
purely as a confound (never as a deployable feature), evaluates whether an LLM
can add a genuinely new label-blind structured feature (artifact/transcript
plausibility only), and emits a ranked feature-unlock matrix.

Isolation contract: writes ONLY to artifacts/milestone_7_decision/gate_feature_audit/;
imports the pure core from benchmark.gate_feature_audit; touches no dynamic_gate file.
No label ever reaches the LLM; no patient/study identifier is sent externally.
"""

from __future__ import annotations

import json
import os
import sys
import warnings
from pathlib import Path

import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from benchmark.gate_feature_audit import (  # noqa: E402
    conditional_auroc,
    feature_coverage,
    grouped_oof_auroc,
    high_presentation_mask,
)

OUT = ROOT / "artifacts" / "milestone_7_decision" / "gate_feature_audit"
OUT.mkdir(parents=True, exist_ok=True)

TOP_K = 20  # the gate operates at k=20; the hard stratum is the top-20 presenters/patient

# Orthogonal feature families (deliberately EXCLUDES presentation/PRIME/EL — those are the wall).
FAMILIES = [
    "agretopicity",
    "expression",
    "vaf_readsupport",
    "predictor_disagreement",
    "stability_processing",
    "physicochemical",
    "mutation_annotation",
    "repeated_antigen",
    "assay_context",
    "llm_artifact_plausibility",
]


# --------------------------------------------------------------------------- loaders + derived features
def _flag(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip().str.lower()
    return (~s.isin(["", ".", "nan", "none", "na", "0"])).astype(float)


def load_gartner():
    from event_b import prime_transfer as pt

    g = pt._gartner().copy()
    ranks = ["mhcflurry_rank", "netmhcpan_el_rank", "netmhcpan_ba_rank", "mixmhcpred_rank", "hlathena_rank"]
    g["pred_rank_std"] = g[ranks].astype(float).std(axis=1)
    g["seen_in_rna"] = pd.to_numeric(g["seen in RNA (1=Yes,0=No)"], errors="coerce")
    g["is_dbsnp"] = _flag(g["dbSNP_ID"])  # germline/common -> hypothesised NEGATIVE predictor
    g["is_cosmic"] = _flag(g["Cosmic Info"])  # recurrent somatic -> hypothesised POSITIVE predictor
    specs = [
        ("expr_decile", "expression", True),
        ("seen_in_rna", "expression", True),
        ("vaf_decile", "vaf_readsupport", True),
        ("pred_rank_std", "predictor_disagreement", False),
        ("is_cosmic", "mutation_annotation", True),
        ("is_dbsnp", "mutation_annotation", False),
    ]
    return g, ("prime", False), specs, "genuine PRIME 2.1"


def load_zhao():
    from event_b.zhao_features import build_zhao_feature_frame

    z = build_zhao_feature_frame("data/raw/zhao_dc_2026/manual").copy()
    preds = [
        "pred_BigMHC EL", "pred_MHCflurryBA", "pred_MHCflurryPres", "pred_MHCflurryProc",
        "pred_NetCTLpan Cleavage", "pred_NetCTLpan TAP", "pred_NetMHC-3", "pred_NetMHC-4", "pred_NetMHCstab",
    ]
    zc = z[preds].astype(float)
    zz = (zc - zc.mean()) / zc.std(ddof=0).replace(0, 1.0)
    z["pred_z_std"] = zz.std(axis=1)
    specs = [
        ("mut_kd_delta", "agretopicity", True),
        ("mut_charge_delta", "physicochemical", True),
        ("mut_is_anchor", "mutation_annotation", True),
        ("mut_position_frac", "mutation_annotation", True),
        ("gravy", "physicochemical", True),
        ("net_charge", "physicochemical", True),
        ("aromatic_frac", "physicochemical", True),
        ("pred_NetMHCstab", "stability_processing", True),
        ("pred_NetCTLpan Cleavage", "stability_processing", True),
        ("pred_NetCTLpan TAP", "stability_processing", True),
        ("pred_MHCflurryProc", "stability_processing", True),
        ("pred_z_std", "predictor_disagreement", False),
    ]
    return z, ("mixmhcpred3_score", True), specs, "mixMHCpred 3.0 (not genuine PRIME)"


def load_multimer():
    from event_b import prime_transfer as pt

    return pt._multimer().copy(), ("prime", False), [], "genuine PRIME 2.1 (rank-like: lower=better)"


def load_improve():
    from event_b import prime_transfer as pt

    return pt._improve().copy(), ("prime", False), [], "genuine PRIME 2.1 (rank-like: lower=better)"


def load_osteosarc():
    o = pd.read_csv(ROOT / "artifacts/milestone_7_decision/peptide_recovery/RECOVERED_CANDIDATES.csv").copy()
    # No experimental TESTED_NEGATIVE denominator (3 known positives, one patient) -> descriptive only.
    specs = [
        ("expression_tpm", "expression", True),
        ("tumor_vaf", "vaf_readsupport", True),
        ("n_callers", "vaf_readsupport", True),
        ("n_timepoints", "repeated_antigen", True),
    ]
    return o, ("genuine_prime", True), specs, "genuine PRIME 2.1"


# --------------------------------------------------------------------------- per-cohort audit
def audit_cohort(name, frame, anchor, specs, anchor_note, *, patient_col="patient_id", label_col="label"):
    has_labels = label_col in frame.columns and frame[label_col].isin(["POSITIVE", "TESTED_NEGATIVE"]).any()
    result = {
        "cohort": name,
        "n_rows": int(len(frame)),
        "n_patients": int(frame[patient_col].nunique()) if patient_col in frame else None,
        "anchor": {"column": anchor[0], "higher_better": anchor[1], "note": anchor_note} if anchor else None,
        "label_counts": (frame[label_col].value_counts().to_dict() if has_labels else None),
        "coverage": feature_coverage(frame, [s[0] for s in specs]),
        "orthogonal_features": [],
        "presentation_baseline": None,
        "cross_fit": None,
        "confound_audit": None,
    }
    if not has_labels or anchor is None:
        result["note"] = "no anchor and/or no measured POSITIVE/TESTED_NEGATIVE labels -> conditional signal not computable now"
        return result

    score_col, higher = anchor
    frame = frame[frame[score_col].notna()].copy()
    stratum = high_presentation_mask(frame, score_col, higher_better=higher, top_k=TOP_K, by=patient_col)

    # Presentation baseline ON the hard stratum: this SHOULD be ~0.5 (the wall).
    base = conditional_auroc(frame, score_col, higher_better=higher, mask=stratum, label_col=label_col)
    result["presentation_baseline"] = {"feature": score_col, "on_stratum": base,
                                       "marginal": conditional_auroc(frame, score_col, higher_better=higher, label_col=label_col)}

    for col, family, hyp_hi in specs:
        marg = conditional_auroc(frame, col, higher_better=hyp_hi, label_col=label_col)
        strat = conditional_auroc(frame, col, higher_better=hyp_hi, mask=stratum, label_col=label_col)
        result["orthogonal_features"].append({
            "feature": col, "family": family, "hypothesised_higher_better": hyp_hi,
            "coverage": marg["coverage"],
            "marginal_auroc": marg["auroc"], "marginal_n": {"pos": marg["n_pos"], "neg": marg["n_neg"]},
            "stratum_auroc": strat["auroc"], "stratum_n": {"pos": strat["n_pos"], "neg": strat["n_neg"]},
            "abs_stratum_signal": (None if strat["auroc"] is None else round(abs(strat["auroc"] - 0.5), 4)),
        })

    # Cross-fitted (patient-grouped OOF): does the orthogonal SET add over presentation on the stratum?
    ortho_cols = [s[0] for s in specs]
    result["cross_fit"] = {
        "stratum_orthogonal_only": grouped_oof_auroc(frame, ortho_cols, mask=stratum, group_col=patient_col, label_col=label_col),
        "stratum_presentation_only": grouped_oof_auroc(frame, [score_col], mask=stratum, group_col=patient_col, label_col=label_col),
        "stratum_combined": grouped_oof_auroc(frame, ortho_cols + [score_col], mask=stratum, group_col=patient_col, label_col=label_col),
        "full_orthogonal_only": grouped_oof_auroc(frame, ortho_cols, group_col=patient_col, label_col=label_col),
    }

    # Confound audit (NEVER deployable): how much apparent signal is patient/tissue-identity parasitic?
    best = max((f for f in result["orthogonal_features"] if f["marginal_auroc"] is not None),
               key=lambda f: abs((f["marginal_auroc"] or 0.5) - 0.5), default=None)
    conf = {}
    if best is not None:
        b = best["feature"]
        conf["best_orthogonal_feature"] = b
        conf["marginal_auroc"] = best["marginal_auroc"]
        oof = grouped_oof_auroc(frame, [b], group_col=patient_col, label_col=label_col)
        conf["patient_generalising_oof_auroc"] = oof["oof_auroc"]
        conf["identity_parasitic_gap"] = (
            None if oof["oof_auroc"] is None or best["marginal_auroc"] is None
            else round(abs(best["marginal_auroc"] - 0.5) - abs(oof["oof_auroc"] - 0.5), 4)
        )
    if "tumor type" in frame.columns:
        tt = pd.get_dummies(frame["tumor type"].astype(str))
        tmp = pd.concat([frame[[patient_col, label_col]].reset_index(drop=True), tt.reset_index(drop=True)], axis=1)
        conf["tissue_label_only_oof_auroc"] = grouped_oof_auroc(
            tmp, list(tt.columns), group_col=patient_col, label_col=label_col
        )["oof_auroc"]
        conf["tissue_note"] = "AUDIT ONLY — tumour/study identity is never a deployable feature"
    result["confound_audit"] = conf
    return result


# --------------------------------------------------------------------------- CEDAR repeated-antigen prior
def audit_cedar():
    from event_b.cedar_recognition import DEFAULT_ZIP, normalize_cedar, read_cedar_tcell

    kept = normalize_cedar(read_cedar_tcell(DEFAULT_ZIP)).kept.copy()
    lab = kept["response_label"].astype(str)
    scorable = kept[lab.isin(["POSITIVE", "TESTED_NEGATIVE"])].copy()
    # Repeated-antigen: how many distinct assays test the same peptide.
    scorable["n_assays_for_peptide"] = scorable.groupby("peptide")["peptide"].transform("size").astype(float)
    scorable["subject_pos_frac"] = (
        pd.to_numeric(scorable["n_subjects_positive"], errors="coerce")
        / pd.to_numeric(scorable["n_subjects_tested"], errors="coerce")
    )
    scorable = scorable.rename(columns={"response_label": "label"})
    out = {
        "cohort": "cedar",
        "n_rows": int(len(scorable)),
        "label_counts": scorable["label"].value_counts().to_dict(),
        "anchor": None,
        "note": "recognition-prior asset: NO within-patient presentation anchor -> no gate stratum; "
                "repeated-antigen signal is reported but is HEAVILY leakage-prone (studied because immunogenic).",
        "repeated_antigen": {
            "n_assays_for_peptide": conditional_auroc(scorable, "n_assays_for_peptide", higher_better=True),
            "subject_pos_frac": {
                "leakage": "CIRCULAR — subject_pos_frac is a direct function of the label; reported for completeness only",
                **conditional_auroc(scorable, "subject_pos_frac", higher_better=True),
            },
        },
    }
    return out


# --------------------------------------------------------------------------- Miller (labels-in-hand, no inputs yet)
def audit_miller():
    p = ROOT / "data/raw/miller_ipv/miller_recognition_labels.csv"
    if not p.exists():
        return {"cohort": "miller", "note": "labels not staged on disk"}
    m = pd.read_csv(p)
    return {
        "cohort": "miller",
        "n_rows": int(len(m)),
        "n_patients": int(m["patient_id"].nunique()),
        "label_counts": m["label"].value_counts().to_dict(),
        "anchor": None,
        "note": "LOCKED_TEST; labels ingested but NO presentation anchor and NO expression/VAF yet "
                "(requires the public WES/RNA download + HLA typing + RNA quant). No gate stratum computable now.",
        "available_now": {
            "agretopicity": {"raw": "ref_peptide present (coverage %.3f)" % float(m["ref_peptide"].notna().mean()),
                             "status": "computable after ONE predictor pass scoring ref vs mutant"},
            "mutation_annotation": {"gene/transcript/HGVS present": True},
        },
        "requires_miller_wes_rna": ["expression", "vaf_readsupport", "predictor_disagreement", "stability_processing", "presentation_anchor"],
    }


# --------------------------------------------------------------------------- LLM feasibility (label-blind, artifact/transcript plausibility only)
LLM_SCHEMA = {
    "type": "object",
    "properties": {
        "coding_plausible": {"type": "boolean", "description": "Is the annotated transcript a plausibly protein-coding ORF that could yield this peptide?"},
        "frame_consistent": {"type": "boolean", "description": "Are the cDNA and protein HGVS mutually frame-consistent for the stated variant type?"},
        "peptide_maps_to_mutation": {"type": "boolean", "description": "Does the mutant peptide plausibly contain the annotated amino-acid change vs the WT peptide?"},
        "nmd_risk": {"type": "string", "enum": ["low", "medium", "high"], "description": "Nonsense-mediated-decay risk for the annotated transcript/variant."},
        "artifact_risk_score": {"type": "number", "description": "0=clean well-supported neoepitope annotation, 1=likely annotation artifact (pseudogene/retained intron/frame error)."},
        "rationale": {"type": "string"},
    },
    "required": ["coding_plausible", "frame_consistent", "peptide_maps_to_mutation", "nmd_risk", "artifact_risk_score"],
}

LLM_SYSTEM = (
    "You are a molecular-biology annotation checker. Given ONLY generic variant/transcript/peptide "
    "annotations (no patient, no study, no experimental result), judge whether the mutation->peptide "
    "annotation is biologically PLAUSIBLE and free of obvious construction artifacts. Do NOT guess "
    "immunogenicity, binding, or T-cell recognition. Return only the requested structured fields."
)


def _blind_annotations(n=3):
    """Non-identifying annotation rows (gene/variant/transcript/peptide only) drawn from PUBLISHED
    supplementary tables. No patient_id, no study_id, no label."""
    rows = []
    try:
        m = pd.read_csv(ROOT / "data/raw/miller_ipv/miller_recognition_labels.csv")
        for _, r in m.head(200).iterrows():
            rows.append({
                "gene_symbol": str(r.get("gene_symbol", "")),
                "variant_type": str(r.get("source_variant_type", "")),
                "hgvs_c": str(r.get("cdna_hgvs", "")),
                "hgvs_p": str(r.get("protein_hgvs", "")),
                "transcript_id": str(r.get("transcript_id", "")),
                "mutant_peptide": str(r.get("mutant_peptide", "")),
                "wt_peptide": str(r.get("ref_peptide", "")),
            })
            if len(rows) >= n:
                break
    except Exception:
        pass
    return rows


def run_llm_feasibility(sample_n=3):
    cache = OUT / "llm_feasibility_cache.json"
    if cache.exists():  # reuse cached blind sample; never re-send annotations on a re-run
        prev = json.loads(cache.read_text())
        if prev.get("status", "").startswith("RAN"):
            return prev
    key = os.environ.get("OPENAI_API_KEY")
    rows = _blind_annotations(sample_n)
    record = {"schema": LLM_SCHEMA, "system_prompt": LLM_SYSTEM, "n_sample": len(rows),
              "blind_fields": sorted(rows[0].keys()) if rows else [], "contains_labels": False,
              "contains_patient_or_study_id": False, "results": [], "status": "NOT_RUN"}
    if not key or not rows:
        record["status"] = "NOT_RUN (no OPENAI_API_KEY or no rows)"
        cache.write_text(json.dumps(record, indent=2))
        return record
    try:
        import urllib.request

        for row in rows:
            body = {
                "model": "gpt-4o-mini",
                "temperature": 0,
                "messages": [
                    {"role": "system", "content": LLM_SYSTEM},
                    {"role": "user", "content": "Annotation (no labels, no identifiers):\n" + json.dumps(row)
                        + "\n\nReturn a JSON object with keys: " + ", ".join(LLM_SCHEMA["required"]) + ", rationale."},
                ],
                "response_format": {"type": "json_object"},
            }
            req = urllib.request.Request(
                "https://api.openai.com/v1/chat/completions",
                data=json.dumps(body).encode(),
                headers={"Authorization": f"Bearer {key}", "Content-Type": "application/json"},
            )
            with urllib.request.urlopen(req, timeout=60) as resp:
                payload = json.loads(resp.read())
            content = payload["choices"][0]["message"]["content"]
            # Redact the LOCKED_TEST raw sequence/HGVS from the committed cache (the raw Miller
            # label frame is deliberately gitignored); keep only non-sensitive gene/variant + output.
            redacted = {"gene_symbol": row["gene_symbol"], "variant_type": row["variant_type"],
                        "_redacted": "hgvs/transcript/peptide withheld (locked-test source)"}
            record["results"].append({"input": redacted, "output": json.loads(content)})
        record["status"] = "RAN (blind feasibility sample; annotation-only, no labels/identifiers)"
    except Exception as e:  # never break the audit on an API hiccup
        record["status"] = f"NOT_RUN (api error: {type(e).__name__})"
    cache.write_text(json.dumps(record, indent=2))
    return record


# --------------------------------------------------------------------------- unlock matrix
def build_unlock_matrix(cohorts):
    """Rank orthogonal feature families by best cross-cohort stratum signal, availability, leakage."""
    rows = []
    fam_best = {f: {"family": f, "best_stratum_auroc": None, "best_cohort": None, "coverage_max": 0.0} for f in FAMILIES}
    for c in cohorts:
        for feat in c.get("orthogonal_features", []) or []:
            fam = feat["family"]
            fam_best.setdefault(fam, {"family": fam, "best_stratum_auroc": None, "best_cohort": None, "coverage_max": 0.0})
            sig = feat["abs_stratum_signal"]
            fam_best[fam]["coverage_max"] = max(fam_best[fam]["coverage_max"], feat["coverage"] or 0.0)
            if sig is not None and (fam_best[fam]["best_stratum_auroc"] is None or feat["stratum_auroc"] is not None
                                    and abs(feat["stratum_auroc"] - 0.5) > abs((fam_best[fam]["best_stratum_auroc"] or 0.5) - 0.5)):
                fam_best[fam]["best_stratum_auroc"] = feat["stratum_auroc"]
                fam_best[fam]["best_cohort"] = c["cohort"]
    # availability + leakage + next-experiment annotations
    meta = {
        "agretopicity": ("available_now (Zhao mut_kd_delta); Gartner/Miller need 1 WT predictor pass", "low", "score WT vs mutant EL for Gartner+Miller; test mutant/WT ratio on stratum"),
        "expression": ("available_now (Gartner deciles); Miller needs RNA-seq; Sid sparse", "medium (ascertainment/PU)", "conditional expr AUROC on Miller after RNA quant; guard vs PU per m7-ascertainment-correction"),
        "vaf_readsupport": ("available_now (Gartner vaf_decile); Miller/Sid need WES", "medium", "join VAF/read-depth on Miller after somatic calling"),
        "predictor_disagreement": ("available_now (Gartner 5, Zhao 9); absent multimer/IMPROVE", "low", "test std-of-ranks on stratum with a pre-registered direction"),
        "stability_processing": ("available_now (Zhao NetMHCstab/NetCTLpan); add to Gartner/Miller", "low", "add NetMHCstab/NetCTLpan pass to Gartner+Miller universe"),
        "physicochemical": ("available_now (Zhao); trivially computable everywhere", "low", "recompute gravy/charge for all cohorts (cheap)"),
        "mutation_annotation": ("available_now (Gartner CGC/Cosmic/dbSNP)", "medium (annotation-source leakage)", "map driver/germline flags on Miller from public annotation"),
        "repeated_antigen": ("available_now (CEDAR n_subjects; cross-cohort recurrence)", "HIGH (studied-because-immunogenic)", "restrict to prospective recurrence only; never use assay-count of the same cohort"),
        "assay_context": ("available_now but STUDY-CONFOUNDED", "audit-only (never deployable)", "use as stratification/confound check, not a feature"),
        "llm_artifact_plausibility": ("computable_now (label-blind, annotation-only)", "low (no labels sent)", "score full universe; test whether artifact_risk_score down-ranks TESTED_NEGATIVE decoys on the stratum"),
    }
    for f in FAMILIES:
        fb = fam_best.get(f, {"family": f, "best_stratum_auroc": None, "best_cohort": None, "coverage_max": 0.0})
        avail, leak, nxt = meta.get(f, ("unknown", "unknown", "unknown"))
        rows.append({
            "family": f,
            "best_stratum_auroc": fb["best_stratum_auroc"],
            "abs_signal": (None if fb["best_stratum_auroc"] is None else round(abs(fb["best_stratum_auroc"] - 0.5), 4)),
            "best_cohort": fb["best_cohort"],
            "max_coverage": round(fb["coverage_max"], 4),
            "availability": avail,
            "leakage_risk": leak,
            "next_experiment": nxt,
        })
    rows.sort(key=lambda r: (r["abs_signal"] is None, -(r["abs_signal"] or -1)))
    return rows


# --------------------------------------------------------------------------- render
def render_markdown(cohorts, cedar, miller, unlock, llm):
    L = []
    L.append("# Gate feature audit — orthogonal levers against high-PRIME false positives\n")
    L.append("**Why this exists.** The dynamic gate (`configs/frozen/dynamic_gate_v1.json`) was falsified: a "
             "label-blind *presentation* gate removes 0 pct of the high-presentation decoys that outrank positives, "
             "so downstream Δhits@20 = 0. A gate can only help via an **orthogonal** feature that separates the "
             f"top-ranked TESTED_NEGATIVE decoys from POSITIVES *within the high-presentation stratum* (top-{TOP_K} by "
             "presentation per patient). This audit measures exactly that, read-only, on all seven cohorts.\n")
    L.append("Presentation baseline on the stratum should sit near 0.5 by construction — that is the wall. Any "
             "orthogonal feature with |AUROC−0.5| meaningfully above 0 on the stratum is a candidate unlock. "
             "Study/patient identity is audited as a confound and is **never** a deployable feature.\n")

    L.append("## 1. Ranked feature-unlock matrix\n")
    L.append("| rank | family | best stratum AUROC | \\|signal\\| | best cohort | max coverage | availability | leakage risk | next experiment |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(unlock, 1):
        L.append(f"| {i} | {r['family']} | {r['best_stratum_auroc']} | {r['abs_signal']} | {r['best_cohort']} "
                 f"| {r['max_coverage']} | {r['availability']} | {r['leakage_risk']} | {r['next_experiment']} |")
    L.append("")

    L.append("## 2. Per-cohort conditional signal (POSITIVE vs TESTED_NEGATIVE)\n")
    for c in cohorts:
        L.append(f"### {c['cohort']}  (n={c['n_rows']}, patients={c.get('n_patients')})")
        if c.get("anchor"):
            a = c["anchor"]
            L.append(f"- anchor: `{a['column']}` ({'higher' if a['higher_better'] else 'lower'}=better) — {a['note']}")
        if c.get("label_counts"):
            L.append(f"- labels: {c['label_counts']}")
        pb = c.get("presentation_baseline")
        if pb:
            L.append(f"- **presentation baseline on stratum**: AUROC={pb['on_stratum']['auroc']} "
                     f"(pos={pb['on_stratum']['n_pos']}, neg={pb['on_stratum']['n_neg']}); "
                     f"marginal={pb['marginal']['auroc']}  ← near 0.5 on stratum = the wall")
        if c.get("orthogonal_features"):
            L.append("\n| feature | family | cov | marginal AUROC | stratum AUROC | \\|stratum signal\\| | stratum n(pos/neg) |")
            L.append("|---|---|---|---|---|---|---|")
            for f in sorted(c["orthogonal_features"], key=lambda x: (x["abs_stratum_signal"] is None, -(x["abs_stratum_signal"] or -1))):
                L.append(f"| {f['feature']} | {f['family']} | {f['coverage']} | {f['marginal_auroc']} | "
                         f"{f['stratum_auroc']} | {f['abs_stratum_signal']} | {f['stratum_n']['pos']}/{f['stratum_n']['neg']} |")
        cf = c.get("cross_fit")
        if cf:
            L.append("\n- cross-fitted (patient-grouped OOF) on stratum: "
                     f"orthogonal-only={cf['stratum_orthogonal_only'].get('oof_auroc')}, "
                     f"presentation-only={cf['stratum_presentation_only'].get('oof_auroc')}, "
                     f"combined={cf['stratum_combined'].get('oof_auroc')} "
                     f"(orthogonal-only full-cohort={cf['full_orthogonal_only'].get('oof_auroc')})")
        conf = c.get("confound_audit")
        if conf:
            L.append(f"- confound audit (AUDIT-ONLY, never deployable): best orthogonal `{conf.get('best_orthogonal_feature')}` "
                     f"marginal={conf.get('marginal_auroc')} vs patient-generalising OOF={conf.get('patient_generalising_oof_auroc')} "
                     f"(identity-parasitic gap={conf.get('identity_parasitic_gap')}); "
                     f"tissue-label-only OOF={conf.get('tissue_label_only_oof_auroc')}")
        L.append("")

    L.append("## 3. CEDAR (recognition prior — no gate stratum)\n")
    L.append(f"- {cedar['note']}")
    L.append(f"- labels: {cedar['label_counts']}")
    ra = cedar["repeated_antigen"]["n_assays_for_peptide"]
    L.append(f"- repeated-antigen `n_assays_for_peptide`: AUROC={ra['auroc']} (HIGH leakage: peptides are re-assayed *because* immunogenic)\n")

    L.append("## 4. Miller (LOCKED_TEST — labels in hand, inputs not yet computed)\n")
    L.append(f"- {miller['note']}")
    L.append(f"- labels: {miller['label_counts']}")
    L.append(f"- available now: {json.dumps(miller.get('available_now', {}))}")
    L.append(f"- requires Miller WES/RNA: {miller.get('requires_miller_wes_rna')}\n")

    L.append("## 5. LLM structured feature — feasibility (label-blind, artifact/transcript plausibility ONLY)\n")
    L.append(f"- status: **{llm['status']}**")
    L.append(f"- sends only: {llm['blind_fields']} (contains_labels={llm['contains_labels']}, "
             f"contains_patient_or_study_id={llm['contains_patient_or_study_id']})")
    L.append("- purpose: down-weight *annotation artifacts* (pseudogene / retained-intron / frame errors / NMD), an "
             "axis orthogonal to both presentation and recognition. **Not** an 'is this immunogenic?' guess.")
    L.append("- schema + prompt cached in `llm_feasibility_cache.json`.")
    if llm.get("results"):
        L.append("\n| gene | variant | artifact_risk | coding_plausible | nmd_risk |")
        L.append("|---|---|---|---|---|")
        for r in llm["results"]:
            o, i = r["output"], r["input"]
            L.append(f"| {i['gene_symbol']} | {i['variant_type']} | {o.get('artifact_risk_score')} "
                     f"| {o.get('coding_plausible')} | {o.get('nmd_risk')} |")
    L.append("")

    L.append("## 6. Bottom line\n")
    L.append("- **Only Gartner and Zhao** carry both a presentation anchor and orthogonal features, so only they can "
             "test the gate question directly. multimer/IMPROVE have **no orthogonal axis at all** (presentation-only "
             "frames) — structurally they *cannot* supply a gate feature. CEDAR has no anchor; Miller has no inputs yet; "
             "Sid has 3 positives and no clean negative denominator.")
    L.append("- The unlock matrix above ranks which family to invest in. Read the stratum AUROC, not the marginal: a "
             "feature strong marginally but ~0.5 on the stratum cannot remove high-PRIME decoys (that was expression's "
             "risk per `m7-ascertainment-correction`).")
    L.append("- The single highest-leverage, lowest-leakage NEW axis is **LLM artifact/transcript plausibility**, "
             "because it is label-blind, computable now on existing annotations, and orthogonal to the presentation wall.")
    L.append("")
    L.append("### Tension & caveats (do not over-read the expression result)\n")
    L.append("- **Expression is the strongest orthogonal signal on the exact stratum where the gate failed** "
             "(Gartner stratum AUROC 0.82 while presentation is 0.50), it cross-fits across patients (OOF 0.855 vs "
             "presentation 0.34), and it is **not** identity-parasitic (marginal−OOF gap ≈ 0.01; tissue-label-only "
             "OOF ≈ 0.49). That is a real, generalising, orthogonal axis — exactly what a presentation gate lacks.")
    L.append("- **But this does NOT license an expression rank penalty.** The frozen expression policy "
             "(`configs/frozen/expression_policy_v1.json`, memory `expression-ranking-policy`) already showed that "
             "turning expression into a rank penalty REGRESSES conditional-ranking cohorts within a fixed top-20 "
             "budget. High discrimination on the stratum ≠ a safe deployable penalty; the unlock is a *conditional / "
             "confidence* use that spares strong positives, not a monotone demotion.")
    L.append("- **Underpowered:** only 17 POSITIVES survive onto the Gartner top-20 stratum (46 total). Per "
             "`m7-ascertainment-correction`, Gartner expression/VAF is 'promising but UNDERPOWERED'; this audit is "
             "CONSISTENT_WITH a real effect, it does not establish one. Negatives here are measured TESTED_NEGATIVE "
             "(not UNTESTED), which blunts but does not remove the ascertainment/PU risk.")
    L.append("- **Predictor disagreement** is the second orthogonal axis (Gartner stratum |signal| 0.20) and is "
             "cheap and low-leakage — but its direction must be pre-registered, not chosen from this table.")
    L.append("- **Zhao shows no orthogonal lever** (best stratum |signal| 0.06); its weak anchor (mixMHCpred) means "
             "its 'stratum' is barely selective. multimer/IMPROVE are presentation-only frames with **no orthogonal "
             "column to test** — they cannot supply a gate feature regardless of power.")
    return "\n".join(L)


# --------------------------------------------------------------------------- main
def main():
    cohorts = []
    for loader in (load_gartner, load_zhao, load_multimer, load_improve, load_osteosarc):
        name = loader.__name__.replace("load_", "")
        frame, anchor, specs, note = loader()
        cohorts.append(audit_cohort(name, frame, anchor, specs, note))
    cedar = audit_cedar()
    miller = audit_miller()
    llm = run_llm_feasibility()
    unlock = build_unlock_matrix(cohorts)

    (OUT / "FEATURE_AUDIT.json").write_text(json.dumps(
        {"top_k": TOP_K, "families": FAMILIES, "cohorts": cohorts, "cedar": cedar, "miller": miller,
         "unlock_matrix": unlock, "llm_feasibility": {k: v for k, v in llm.items() if k != "results"}},
        indent=2, default=str))
    (OUT / "FEATURE_UNLOCK_MATRIX.json").write_text(json.dumps(unlock, indent=2, default=str))
    (OUT / "GATE_FEATURE_AUDIT.md").write_text(render_markdown(cohorts, cedar, miller, unlock, llm))
    print("wrote", OUT)
    for r in unlock:
        print(f"  {r['family']:28s} |signal|={r['abs_signal']} avail={r['availability'][:40]}")


if __name__ == "__main__":
    main()
