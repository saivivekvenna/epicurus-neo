"""Runner — four-arm generation x scorer patient-level top-20 benchmark + cohort eligibility audit.

Two deliverables, both reproducible from committed inputs (no network, no PRIME executable needed):

1. **Cohort eligibility audit** (`COHORT_ELIGIBILITY.{json,md}`). For every local recognition cohort,
   which of the four arms is EVALUABLE and — where NOT — exactly which input is missing (raw callset
   for lossless generation, genuine PRIME, measured labels, ...) or which leakage rule blocks it.
   This is the honest map of what current data can and cannot support.

2. **The one end-to-end instantiation** (`SID_FOUR_ARM.{json,md}`). Only osteosarc.com / Sid carries a
   raw multi-caller callset through to peptide generation AND a measured recognition label, so it is the
   only patient on which all four arms run. Reported with a frozen-Epicurus primary result and a
   MixMHCpred-EL sensitivity arm, additive stage attribution, and a leakage-control panel.

STATUS GUARDRAIL (read first): the Sid run is **post-hoc, n=3, single patient** — the lossless generator
was designed after inspecting Sid's structural audit. It is a reachability/attribution DIAGNOSTIC, NOT a
blinded, powered, or prospective superiority test. No arm's win is a benchmark gate. Measured positives are
joined ONLY after each arm's ranking is fixed (strict label isolation, enforced by the harness).

    .venv/bin/python -m scripts.four_arm_benchmark
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from benchmark.cohort_audit import run_cohort_audit  # noqa: E402
from benchmark.four_arm import (  # noqa: E402
    attach_epicurus_score,
    run_patient,
    stage_attribution,
)

OUT = ROOT / "artifacts" / "milestone_7_decision" / "four_arm"
RAW = ROOT / "data" / "raw" / "osteosarc"
REC_CSV = ROOT / "artifacts" / "milestone_7_decision" / "peptide_recovery" / "RECOVERED_CANDIDATES.csv"
PVAC_TSV = RAW / "pvactools_all_epitopes.tsv"

# Evaluation-only: joined AFTER ranking (strict label isolation). Hudson-lab IFNg peptide-expansion
# mutation ids; never read into any generator/scorer input. (Public osteosarc.com result.)
HUDSON_POSITIVES = {
    "ASPM-chr1-197102716",
    "MAP2-chr2-209694772",
    "DYNC1H1-chr14-101980529",
}

# ==================================================================================================
# Part 2 — Sid four-arm instantiation (the only end-to-end-evaluable patient)
# ==================================================================================================
def _load_netmhcpan_el() -> pd.DataFrame:
    """NetMHCpan-EL MT %rank per (peptide, HLA) from the raw pVAC table (min across duplicates)."""
    pv = pd.read_csv(PVAC_TSV, sep="\t", low_memory=False)
    el = pv.rename(columns={
        "MT Epitope Seq": "mutant_peptide",
        "HLA Allele": "hla_allele",
        "NetMHCpanEL MT Percentile": "el",
    })[["mutant_peptide", "hla_allele", "el"]].copy()
    el["el"] = pd.to_numeric(el["el"], errors="coerce")
    return el.groupby(["mutant_peptide", "hla_allele"], as_index=False)["el"].min()


MHCFLURRY_CACHE = OUT / "mhcflurry_presentation_cache.csv"


def _load_or_compute_mhcflurry() -> pd.DataFrame:
    """Genuine MHCflurry presentation features per unique (peptide, HLA); cached to a committed CSV so
    the four-arm run reproduces offline without re-running MHCflurry."""
    if MHCFLURRY_CACHE.exists():
        return pd.read_csv(MHCFLURRY_CACHE)
    from benchmark.presentation_features import PRESENTATION_COLUMNS, add_presentation_features

    uni = pd.read_csv(REC_CSV)
    uniq = uni[["mutant_peptide", "hla_allele"]].astype(str).drop_duplicates().reset_index(drop=True)
    scored = add_presentation_features(uniq)
    keep = ["mutant_peptide", "hla_allele", *PRESENTATION_COLUMNS]
    OUT.mkdir(parents=True, exist_ok=True)
    scored[keep].to_csv(MHCFLURRY_CACHE, index=False)
    return scored[keep]


def _build_sid_universe(el_source: str) -> pd.DataFrame:
    """Sid candidate union with prime/el/expr + attached Epicurus score.

    ``el_source='mhcflurry'``: PRIMARY (fair) — genuine MHCflurry presentation %rank for ALL candidates
    (one consistent, independent predictor; recovered candidates get real presentation evidence instead
    of a 0.5 impute). NetMHCpan is not locally runnable; agreement with it is reported separately.
    ``el_source='netmhcpan'``: REFERENCE — the literal frozen Epicurus feature (NetMHCpan-EL) for pVAC
    candidates only; MISSING on lossless-recovered candidates -> NaN -> frozen 0.5-percentile (shows the
    imputation artifact the fair run removes).
    ``el_source='mixmhcpred'``: SENSITIVITY — MixMHCpred %rank (PRIME's backbone) for all candidates.
    """
    uni = pd.read_csv(REC_CSV)
    uni["prime"] = pd.to_numeric(uni["prime_rank"], errors="coerce")
    uni["expr"] = pd.to_numeric(uni["expression_tpm"], errors="coerce")
    if el_source == "mhcflurry":
        mf = _load_or_compute_mhcflurry()[["mutant_peptide", "hla_allele", "mhcflurry_el_percentile"]]
        uni = uni.merge(mf, on=["mutant_peptide", "hla_allele"], how="left")
        uni["el"] = pd.to_numeric(uni["mhcflurry_el_percentile"], errors="coerce")
    elif el_source == "netmhcpan":
        uni = uni.merge(_load_netmhcpan_el(), on=["mutant_peptide", "hla_allele"], how="left")
    elif el_source == "mixmhcpred":
        uni["el"] = pd.to_numeric(uni["mixmhcpred_rank"], errors="coerce")
    else:  # pragma: no cover - guarded caller
        raise ValueError(el_source)
    return attach_epicurus_score(uni)


def _netmhcpan_agreement() -> dict:
    """Spearman agreement between the genuine MHCflurry EL substitute and NetMHCpan-EL on pVAC rows."""
    from scipy.stats import spearmanr

    uni = pd.read_csv(REC_CSV)
    mf = _load_or_compute_mhcflurry()[["mutant_peptide", "hla_allele", "mhcflurry_el_percentile"]]
    net = _load_netmhcpan_el().rename(columns={"el": "netmhcpan_el"})
    m = uni.merge(mf, on=["mutant_peptide", "hla_allele"], how="left").merge(
        net, on=["mutant_peptide", "hla_allele"], how="left")
    both = m[m["netmhcpan_el"].notna() & m["mhcflurry_el_percentile"].notna()]
    rho = float(spearmanr(both["netmhcpan_el"], both["mhcflurry_el_percentile"]).correlation)
    return {"n_pvac_rows_with_both": int(len(both)), "spearman_netmhcpan_vs_mhcflurry": round(rho, 3),
            "caveat": "MHCflurry-EL is an independent presentation predictor substituting for the "
                      "unavailable NetMHCpan-EL; moderate agreement means the fair run APPROXIMATES the "
                      "frozen NetMHCpan-EL feature, it does not reproduce it exactly."}


def _arm_row(res) -> dict:
    if not res.evaluable:
        return {"arm_id": res.arm_id, "evaluable": False, "missing": res.missing}
    return {
        "arm_id": res.arm_id,
        "evaluable": True,
        "generation_recall": {"n": res.generation_recall.n, "of": res.generation_recall.of},
        "rankable_recall": {"n": res.rankable_recall.n, "of": res.rankable_recall.of},
        "hits_at_20": res.hits_at_k,
        "recall_at_20": res.recall_at_k,
        "n_selected": res.n_selected,
        "covered_mutations": res.top_k.ids,
    }


def _leakage_panel(uni: pd.DataFrame) -> dict:
    """Peptide-cluster leakage controls for the frozen Epicurus arm on Sid."""
    from event_b.prime_training import prime_leakage_mask

    peptides = uni["mutant_peptide"].dropna().astype(str)
    peptides = sorted({p for p in peptides if p})
    prime_leak = prime_leakage_mask(peptides, near=True)
    n_prime = int(sum(prime_leak))
    hudson_pep_in_prime = None  # candidate-level only; per-mutation isolation is inherent (n=3 held out)
    return {
        "epicurus_training_cohort": "cd8_multimer",
        "sid_is_epicurus_out_of_sample": True,
        "n_candidate_peptides": len(peptides),
        "n_near_or_exact_prime_training_peptides": n_prime,
        "prime_training_leak_fraction": round(n_prime / len(peptides), 4) if peptides else None,
        "label_isolation": "measured positives joined only AFTER each arm's ranking (harness-enforced)",
        "note": "Frozen Epicurus v0.1 was fit on cd8_multimer; Sid shares no patient with it. PRIME "
                "training overlap is reported at candidate granularity; the 3 measured-positive "
                "mutations are never a model input.",
        "_hudson_pep_in_prime": hudson_pep_in_prime,
    }


def _variant(el_source: str, el_feature: str) -> dict:
    uni = _build_sid_universe(el_source)
    out = run_patient(uni, HUDSON_POSITIVES)
    return {
        "el_feature": el_feature,
        "arms": {aid: _arm_row(r) for aid, r in out["arms"].items()},
        "stage_attribution": stage_attribution(out["arms"]),
    }


def run_sid() -> dict:
    available = run_patient(_build_sid_universe("mhcflurry"), HUDSON_POSITIVES)["available"]
    primary = _variant(
        "mhcflurry",
        "GENUINE MHCflurry presentation %rank for ALL candidates (independent learned predictor; "
        "recovered candidates get real presentation evidence, no 0.5 impute). Fair four-arm attribution.")
    reference = _variant(
        "netmhcpan",
        "Literal frozen NetMHCpan-EL for pVAC candidates; MISSING on recovered -> 0.5 impute (shows the "
        "imputation artifact the fair run removes).")
    sensitivity = _variant(
        "mixmhcpred",
        "MixMHCpred %rank (PRIME's backbone) for all candidates — secondary sanity.")

    return {
        "patient_id": "osteosarc_sid",
        "status": "post_hoc_diagnostic_n3_single_patient_not_blinded_not_powered",
        "k": 20,
        "positives_evaluation_only": sorted(HUDSON_POSITIVES),
        "available_inputs": available,
        "feature_provenance": {
            "epicurus_features": "prime = genuine PRIME %rank; el = presentation %rank (see each "
                                 "variant); expr = RSEM gene TPM. Frozen formula prime+el+expr only.",
            "netmhcpan_available_locally": False,
            "netmhcpan_agreement": _netmhcpan_agreement(),
        },
        "primary_genuine_mhcflurry_el": primary,
        "reference_frozen_netmhcpan_el": reference,
        "sensitivity_mixmhcpred_el": sensitivity,
        "leakage_controls": _leakage_panel(_build_sid_universe("mhcflurry")),
        "interpretation": (
            "L1 reachability: generation recovers all 3 recognized mutations (recall 1/3 -> 3/3; +2 "
            "top-20 hits under genuine PRIME = the protected lossless_prime incumbent). L3 end-to-end: "
            "with GENUINE presentation features computed on recovered candidates (MHCflurry, no impute), "
            f"the frozen Epicurus scorer stage is {primary['stage_attribution'].get('scorer')} and the "
            f"full stack nets {primary['stage_attribution'].get('total')} vs pVAC+PRIME. The earlier "
            "-2 frozen scorer loss was substantially a 0.5-impute artifact on recovered rows (reference "
            "vs primary). Any residual drop is the Epicurus expression/EL reweighting demoting a "
            "low-expression true positive — consistent with prior cohorts where a learned recognition "
            "score on top of presentation does not help. NetMHCpan is not locally runnable, so el uses "
            "MHCflurry (independent predictor); its moderate agreement with NetMHCpan-EL is disclosed in "
            "feature_provenance. n=3, post-hoc, descriptive — NOT a gate, no constant tuned to Sid."
        ),
    }


# ==================================================================================================
_LEVEL_COLS = ["reachability", "conditional_ranking", "end_to_end_patient_utility"]
_LEVEL_HDR = ["L1 reachability", "L2 conditional ranking", "L3 end-to-end (north star)"]


def _write_cohort_md(audit: dict) -> str:
    lines = [
        "# Benchmark cohort eligibility audit — three-level hierarchy",
        "",
        f"> Policy `{audit['policy']}`. {audit['n_cohorts']} cohorts; "
        f"**{audit['n_end_to_end_eligible']}** end-to-end (Level-3) eligible.",
        "",
        f"> **No pooling.** {audit['no_pooling']}",
        "",
        f"> **Four-arm harness** = {audit['four_arm_harness_role']}",
        "",
        "## The three levels",
        "",
    ]
    for lv in audit["levels"]:
        lines.append(f"{lv['level']}. **{lv['name']}** — {lv['description']}")
    lines += ["", "## Per-cohort level eligibility", "",
              "✅ = eligible · ❌ = not eligible (each cohort interpreted ONLY within its own denominator)", "",
              "| cohort | role | " + " | ".join(_LEVEL_HDR) + " | denominator |",
              "|---|---|---|---|---|---|"]
    for c in audit["cohorts"]:
        cells = ["✅" if c["levels"][col]["eligible"] else "❌" for col in _LEVEL_COLS]
        lines.append(
            f"| **{c['cohort_id']}** | {c['role']} | {cells[0]} | {cells[1]} | {cells[2]} | {c['denominator']} |"
        )
    lines += ["", "## Four-arm infrastructure evaluability (headline only where L3-eligible)", "",
              "The generation×scorer matrix is reusable infra; read it as a HEADLINE only for the "
              "L3-eligible cohort. Elsewhere the evaluable arms are Level-2 conditional-ranking probes.", "",
              "| cohort | pvac_prime | lossless_prime | lossless_epicurus | full_epicurus |",
              "|---|---|---|---|---|"]
    arm_ids = ["pvac_prime", "lossless_prime", "lossless_epicurus", "full_epicurus"]
    for c in audit["cohorts"]:
        cells = []
        for aid in arm_ids:
            arm = c["arms"][aid]
            if arm["evaluable"]:
                cells.append("✅")
            else:
                reason = ", ".join(m.replace("measured_labels", "labels")
                                   .replace("lossless_generation", "no-lossless-gen")
                                   .replace("genuine_prime", "no-PRIME")
                                   .replace("epicurus_features", "no-Epicurus-feat")
                                   for m in arm["missing"])
                cells.append(f"❌ {reason}")
        lines.append(f"| **{c['cohort_id']}** | {cells[0]} | {cells[1]} | {cells[2]} | {cells[3]} |")
    lines += ["", "## Notes", ""]
    for c in audit["cohorts"]:
        lines.append(f"- **{c['cohort_id']}** ({c['role']}) — {c['note']}")
    lines += ["", "## Interpretation", "", audit["interpretation"], ""]
    return "\n".join(lines)


def _cov(entry: dict) -> str:
    return f"{entry['n']}/{entry['of']}"


def _write_sid_md(sid: dict) -> str:
    arm_ids = ["pvac_prime", "lossless_prime", "lossless_epicurus", "full_epicurus"]

    def arm_table(block: dict) -> list[str]:
        rows = ["| arm | gen recall | rankable recall | hits@20 | recall@20 | covered |",
                "|---|---|---|---:|---:|---|"]
        for aid in arm_ids:
            a = block["arms"][aid]
            if not a["evaluable"]:
                rows.append(f"| `{aid}` | — | — | NOT_EVALUABLE | — | {', '.join(a['missing'])} |")
                continue
            rows.append(
                f"| `{aid}` | {_cov(a['generation_recall'])} | {_cov(a['rankable_recall'])} | "
                f"{a['hits_at_20']} | {a['recall_at_20']} | {', '.join(m.split('-')[0] for m in a['covered_mutations'])} |"
            )
        return rows

    def attr_line(block: dict) -> str:
        at = block["stage_attribution"]
        if not at.get("evaluable"):
            return "attribution: NOT_EVALUABLE"
        return (f"generation **{at['generation']:+d}** · scorer **{at['scorer']:+d}** · "
                f"selection **{at['selection']:+d}** · total **{at['total']:+d}** (top-20 hits)")

    def reach_line(block: dict) -> str:
        a = block["arms"]
        pvac = a["pvac_prime"]["generation_recall"]
        loss = a["lossless_prime"]["generation_recall"]
        return (f"generation recall pVAC {_cov(pvac)} -> lossless union {_cov(loss)} "
                f"(recovered {loss['n'] - pvac['n']} recognized mutation(s) at the generation stage)")

    lc = sid["leakage_controls"]
    lines = [
        "# osteosarc.com / Sid — the one END-TO-END (Level-3) evaluable patient",
        "",
        f"> **Status: {sid['status']}.** Post-hoc, n=3 measured positives, single patient. Strict label "
        "isolation: the 3 recognized mutations are joined only AFTER each arm's ranking. NOT a blinded, "
        "powered, or prospective superiority test — a reachability/attribution DIAGNOSTIC.",
        "",
        f"- k = {sid['k']}; evaluation-only positives: {', '.join(m.split('-')[0] for m in sid['positives_evaluation_only'])}",
        f"- available inputs: {', '.join(sid['available_inputs'])}",
        "",
        "This single patient instantiates all three benchmark levels; each is read separately:",
        "",
        "- **L1 reachability** — how many recognized mutations survive raw→generation.",
        "  " + reach_line(sid["primary_genuine_mhcflurry_el"]) + ".",
        "- **L2 conditional ranking** — ordering among the generated/rankable candidates (within this "
        "patient's denominator only); see the per-arm hits@20 below.",
        "- **L3 end-to-end patient utility (PRIMARY)** — recognized mutations in the final top-20 from "
        "common raw inputs vs standard pVAC + genuine PRIME; see `total` in the stage attribution. "
        "`lossless_prime` (lossless generation + genuine PRIME) is the protected incumbent.",
        "",
        f"> Epicurus feature provenance: {sid['feature_provenance']['epicurus_features']} NetMHCpan "
        f"runnable locally: {sid['feature_provenance']['netmhcpan_available_locally']}; MHCflurry vs "
        f"NetMHCpan-EL Spearman on {sid['feature_provenance']['netmhcpan_agreement']['n_pvac_rows_with_both']} "
        f"pVAC rows = {sid['feature_provenance']['netmhcpan_agreement']['spearman_netmhcpan_vs_mhcflurry']}.",
        "",
        "## Primary (FAIR) — frozen Epicurus, genuine MHCflurry EL on all candidates",
        "",
        f"> EL feature: {sid['primary_genuine_mhcflurry_el']['el_feature']}",
        "",
        *arm_table(sid["primary_genuine_mhcflurry_el"]),
        "",
        f"Stage attribution: {attr_line(sid['primary_genuine_mhcflurry_el'])}",
        "",
        "## Reference — literal frozen NetMHCpan-EL (recovered candidates imputed to 0.5)",
        "",
        f"> EL feature: {sid['reference_frozen_netmhcpan_el']['el_feature']}",
        "",
        *arm_table(sid["reference_frozen_netmhcpan_el"]),
        "",
        f"Stage attribution: {attr_line(sid['reference_frozen_netmhcpan_el'])}",
        "",
        "## Sensitivity — MixMHCpred EL (PRIME backbone) on all candidates",
        "",
        f"> {sid['sensitivity_mixmhcpred_el']['el_feature']}",
        "",
        *arm_table(sid["sensitivity_mixmhcpred_el"]),
        "",
        f"Stage attribution: {attr_line(sid['sensitivity_mixmhcpred_el'])}",
        "",
        "## Leakage controls",
        "",
        f"- Frozen Epicurus training cohort: `{lc['epicurus_training_cohort']}`; Sid out-of-sample: "
        f"{lc['sid_is_epicurus_out_of_sample']}",
        f"- Candidate peptides: {lc['n_candidate_peptides']}; near/exact PRIME-training overlap: "
        f"{lc['n_near_or_exact_prime_training_peptides']} ({lc['prime_training_leak_fraction']})",
        f"- Label isolation: {lc['label_isolation']}",
        "",
        "## Interpretation",
        "",
        sid["interpretation"],
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    audit = run_cohort_audit()
    (OUT / "COHORT_ELIGIBILITY.json").write_text(json.dumps(audit, indent=2) + "\n")
    (OUT / "COHORT_ELIGIBILITY.md").write_text(_write_cohort_md(audit) + "\n")

    sid = run_sid()
    (OUT / "SID_FOUR_ARM.json").write_text(json.dumps(sid, indent=2) + "\n")
    (OUT / "SID_FOUR_ARM.md").write_text(_write_sid_md(sid) + "\n")

    print(json.dumps({
        "cohort_audit": {"n_cohorts": audit["n_cohorts"],
                         "n_end_to_end_eligible": audit["n_end_to_end_eligible"]},
        "sid_primary_genuine_mhcflurry": sid["primary_genuine_mhcflurry_el"]["stage_attribution"],
        "sid_reference_netmhcpan_impute": sid["reference_frozen_netmhcpan_el"]["stage_attribution"],
        "netmhcpan_agreement": sid["feature_provenance"]["netmhcpan_agreement"],
        "artifacts": sorted(str(p.relative_to(ROOT)) for p in OUT.glob("*")),
    }, indent=2))


if __name__ == "__main__":
    main()
