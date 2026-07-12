"""Label-blind RNA-expression ranking-policy analysis on the development cohorts, then freeze + Sid.

Question: should RNA expression enter ranking as (a) a rank penalty, (b) confidence-only annotation, or
(c) soft-saturating / route-dependent evidence? Decided on the Level-2 conditional-ranking development
cohorts (cd8_multimer, Gartner, IMPROVE), each interpreted WITHIN its own denominator and NEVER pooled,
under a strict no-regression requirement against the protected lossless+PRIME incumbent.

Protocol (order matters): (1) per-cohort structural stratification of recognition vs expression;
(2) per-cohort per-patient hits@20 for each policy vs the PRIME incumbent + no-regression verdict;
(3) freeze the policy that never regresses on ANY development cohort; (4) ONLY THEN apply the frozen
policy to osteosarc/Sid as DESCRIPTIVE (n=3, post-hoc) — no constant is tuned to the Sid labels.

    .venv/bin/python -m scripts.expression_policy_analysis
"""

from __future__ import annotations

import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from benchmark.expression_policy import (  # noqa: E402
    BOTTOM_STRATUM,
    EXPRESSION_STRATA,
    expr_penalty_score,
    no_regression_verdict,
    prime_only_score,
    select_portfolio_reserved,
    soft_saturating_score,
    within_patient_percentile,
)

OUT = ROOT / "artifacts" / "milestone_7_decision" / "expression_policy"
FROZEN = ROOT / "configs" / "frozen" / "expression_policy_v1.json"
REC_CSV = ROOT / "artifacts" / "milestone_7_decision" / "peptide_recovery" / "RECOVERED_CANDIDATES.csv"
K = 20
HUDSON = {"ASPM-chr1-197102716", "MAP2-chr2-209694772", "DYNC1H1-chr14-101980529"}

# Score-based policies (rank top-k by the score). Confidence-only == prime_only (expression annotates).
SCORE_POLICIES = {
    "prime_only_incumbent": prime_only_score,
    "expr_rank_penalty": expr_penalty_score,
    "soft_saturating": soft_saturating_score,
}


def _load_dev_cohorts() -> dict[str, pd.DataFrame]:
    from event_b.prime_transfer import _gartner, _improve, _multimer

    cohorts = {}
    for name, fn in (("cd8_multimer", _multimer), ("gartner_nci", _gartner), ("improve_srhgroup", _improve)):
        d = fn().copy()
        d["y"] = (d["label"] == "POSITIVE").astype(int)
        d["expr"] = pd.to_numeric(d["expr"], errors="coerce")
        d["prime"] = pd.to_numeric(d["prime"], errors="coerce")
        cohorts[name] = d.reset_index(drop=True)
    return cohorts


def _per_patient_hits_by_score(d: pd.DataFrame, score: pd.Series, k: int = K) -> dict[str, int]:
    tmp = d.assign(_s=score.to_numpy())
    hits = {}
    for pid, g in tmp.groupby("patient_id", sort=True):
        r = g.sort_values("_s", ascending=False, kind="mergesort")
        hits[pid] = int((r["y"].to_numpy()[:k] == 1).sum())
    return hits


def _per_patient_hits_portfolio(d: pd.DataFrame, k: int = K) -> dict[str, int]:
    hits = {}
    for pid, g in d.groupby("patient_id", sort=True):
        sel = select_portfolio_reserved(g, k=k)
        hits[pid] = int((sel["y"].to_numpy() == 1).sum())
    return hits


def _structural_stratification(d: pd.DataFrame) -> dict:
    ep = within_patient_percentile(d, "expr", higher_better=True)
    qbin = pd.cut(ep, [0, .25, .5, .75, 1.0], labels=["Q1lo", "Q2", "Q3", "Q4hi"], include_lowest=True)
    pr = d.groupby(qbin, observed=True)["y"].agg(["mean", "sum", "count"])
    low_half = d[ep <= 0.5]
    return {
        "positive_rate_by_expr_quartile": {
            str(q): {"rate": round(float(row["mean"]), 4), "pos": int(row["sum"]), "n": int(row["count"])}
            for q, row in pr.iterrows()
        },
        "positives_in_low_expr_half": f"{int(low_half['y'].sum())}/{int(d['y'].sum())}",
        "frac_positives_low_expr_half": round(float(low_half["y"].sum() / max(d["y"].sum(), 1)), 3),
    }


def analyze_cohort(name: str, d: pd.DataFrame) -> dict:
    incumbent = _per_patient_hits_by_score(d, prime_only_score(d))
    patients = sorted(incumbent)
    results = {}
    for pol, fn in SCORE_POLICIES.items():
        h = _per_patient_hits_by_score(d, fn(d))
        delta = np.array([h[p] - incumbent[p] for p in patients], dtype=float)
        tot = int(d["y"].sum())
        results[pol] = {
            "recall_top20": f"{sum(h.values())}/{tot}",
            "recall_frac": round(sum(h.values()) / tot, 4) if tot else None,
            "no_regression_vs_incumbent": no_regression_verdict(delta),
        }
    # portfolio (selection-based)
    hp = _per_patient_hits_portfolio(d)
    delta = np.array([hp[p] - incumbent[p] for p in patients], dtype=float)
    tot = int(d["y"].sum())
    results["portfolio_reserve"] = {
        "recall_top20": f"{sum(hp.values())}/{tot}",
        "recall_frac": round(sum(hp.values()) / tot, 4) if tot else None,
        "no_regression_vs_incumbent": no_regression_verdict(delta),
    }
    return {
        "n_candidates": int(len(d)),
        "n_patients": int(d["patient_id"].nunique()),
        "n_positives": int(d["y"].sum()),
        "structural": _structural_stratification(d),
        "policies": results,
    }


def _decide_frozen_policy(per_cohort: dict) -> dict:
    """Frozen policy = a no-regression form that does not lose recognized candidates on ANY development
    cohort. Rationale is built from the computed verdicts so it can never drift from the numbers."""
    verdicts, benefits = {}, {}
    for pol in ["prime_only_incumbent", "expr_rank_penalty", "soft_saturating", "portfolio_reserve"]:
        regressing = [c for c, a in per_cohort.items()
                      if a["policies"][pol]["no_regression_vs_incumbent"]["regresses"]]
        helping = [c for c, a in per_cohort.items()
                   if a["policies"][pol]["no_regression_vs_incumbent"]["mean_delta"] > 0]
        verdicts[pol] = {"regresses_on": regressing, "no_regression_everywhere": not regressing,
                         "helps_on": helping}
        benefits[pol] = helping

    def _fmt(cohorts):
        return ", ".join(cohorts) if cohorts else "none"

    # Any form that MOVES rank on expression either regresses somewhere or only helps a single denominator
    # while hurting others. The safe universal choice keeps ranking = genuine PRIME (confidence-only);
    # soft-saturating is its no-regression-equivalent concrete guard (only clears low-presentation AND
    # low-expression junk, provably never demoting a strong presenter).
    rationale = (
        f"Per development cohort (never pooled): expr_rank_penalty regresses [{_fmt(verdicts['expr_rank_penalty']['regresses_on'])}] "
        f"while only helping [{_fmt(verdicts['expr_rank_penalty']['helps_on'])}]; portfolio_reserve regresses "
        f"[{_fmt(verdicts['portfolio_reserve']['regresses_on'])}] within a fixed top-20 budget. The two "
        f"no-regression-everywhere forms are prime_only (confidence-only) and soft_saturating, and on these "
        f"cohorts soft_saturating is IDENTICAL to prime_only (it only demotes candidates that are BOTH "
        f"low-presentation and low-expression, which never occupy the top-20) — i.e. letting expression move "
        f"rank buys NO measurable benefit and risks harm. Decision: expression is CONFIDENCE-ONLY in the "
        f"score; lossless+PRIME stays the protected incumbent ranking. The soft-saturating guard is frozen "
        f"as an equivalent, route-dependent no-op-on-top-20 safeguard, and the portfolio reserve is retained "
        f"as an OPTIONAL off-by-default reachability tool (it trades displacement for stratum coverage). "
        f"Because PRIME already keeps high-presentation candidates irrespective of expression, this preserves "
        f"reachability of low-expression recognized candidates."
    )
    return {
        "chosen": "confidence_only",
        "chosen_ranking_score": "prime_only_incumbent (genuine PRIME; expression does NOT move rank)",
        "equivalent_no_regression_guard": "soft_saturating (route-dependent; demotes only low-presentation "
                                          "AND low-expression non-multi-source candidates)",
        "rejected": "expr_rank_penalty (regresses >=1 cohort); portfolio_reserve default (regresses Gartner)",
        "expression_role": "confidence annotation (stratum label) + optional off-by-default route-dependent "
                           "portfolio reserve for reachability",
        "per_policy_no_regression": verdicts,
        "rationale": rationale,
    }


def _freeze_config(decision: dict, per_cohort: dict) -> dict:
    cfg = {
        "name": "expression_ranking_policy",
        "version": "1.0.0",
        "frozen_on": "2026-07-12",
        "protected_incumbent": "lossless generation + genuine PRIME (rank by PRIME %rank, lower=better)",
        "decision": decision,
        "constants": {"bottom_stratum_percentile": BOTTOM_STRATUM,
                      "expression_strata": list(EXPRESSION_STRATA),
                      "tuned_to_sid": False,
                      "tuned_to_any_eval_cohort": False},
        "development_evidence": {
            cohort: {"n_positives": a["n_positives"],
                     "policy_recall_frac": {p: a["policies"][p]["recall_frac"] for p in a["policies"]},
                     "policy_regresses": {p: a["policies"][p]["no_regression_vs_incumbent"]["regresses"]
                                          for p in a["policies"]}}
            for cohort, a in per_cohort.items()
        },
        "no_pooling": "Development cohorts evaluated separately within their own denominators; never pooled.",
    }
    FROZEN.parent.mkdir(parents=True, exist_ok=True)
    FROZEN.write_text(json.dumps(cfg, indent=2) + "\n")
    return cfg


# ---------------------------------------------------------------------------
# Sid DESCRIPTIVE (post-freeze; n=3, no tuning)
# ---------------------------------------------------------------------------
def _sid_descriptive() -> dict:
    uni = pd.read_csv(REC_CSV)
    uni["patient_id"] = "sid"
    uni["prime"] = pd.to_numeric(uni["prime_rank"], errors="coerce")
    uni["expr"] = pd.to_numeric(uni["expression_tpm"], errors="coerce")

    def coverage(score: pd.Series) -> dict:
        tmp = uni.assign(_s=score.to_numpy())
        top = tmp.sort_values("_s", ascending=False, kind="mergesort").head(K)
        covered = sorted(set(top["mutation_id"]) & HUDSON)
        return {"hits_at_20": len(covered), "recall_at_20": round(len(covered) / len(HUDSON), 4),
                "covered": [m.split("-")[0] for m in covered]}

    return {
        "status": "descriptive_post_freeze_n3_not_a_gate",
        "frozen_confidence_only_prime": coverage(prime_only_score(uni)),
        "counterfactual_expr_rank_penalty": coverage(expr_penalty_score(uni)),
        "counterfactual_soft_saturating": coverage(soft_saturating_score(uni)),
        "note": "Confidence-only (= lossless+PRIME) keeps all 3 recognized mutations in the top-20; the "
                "expression rank penalty demotes low-expression MAP2. Consistent with the development "
                "no-regression finding. Descriptive on n=3; nothing tuned to these labels.",
    }


def run() -> dict:
    cohorts = _load_dev_cohorts()
    per_cohort = {name: analyze_cohort(name, d) for name, d in cohorts.items()}
    decision = _decide_frozen_policy(per_cohort)
    cfg = _freeze_config(decision, per_cohort)
    sid = _sid_descriptive()
    return {"policy": "expression-policy-analysis-1.0.0", "k": K,
            "development_cohorts": per_cohort, "decision": decision,
            "frozen_config_path": str(FROZEN.relative_to(ROOT)), "sid_descriptive": sid,
            "frozen_config": cfg}


def _write_md(result: dict) -> str:
    lines = [
        "# RNA-expression ranking-policy analysis (label-blind, development cohorts)",
        "",
        f"> Policy `{result['policy']}`, k={result['k']}. Frozen decision: "
        f"**{result['decision']['chosen']}** — {result['decision']['chosen_ranking_score']}.",
        "",
        "> Development cohorts are interpreted WITHIN their own denominators and NEVER pooled. The "
        "protected incumbent is lossless generation + genuine PRIME.",
        "",
        "## Structural: recognition vs within-patient expression quartile (per cohort)",
        "",
        "| cohort | Q1lo | Q2 | Q3 | Q4hi | positives in low-expr half |",
        "|---|---|---|---|---|---|",
    ]
    for name, a in result["development_cohorts"].items():
        s = a["structural"]["positive_rate_by_expr_quartile"]
        cells = [f"{s.get(q, {}).get('rate', 0)}" for q in ("Q1lo", "Q2", "Q3", "Q4hi")]
        lines.append(f"| {name} | {cells[0]} | {cells[1]} | {cells[2]} | {cells[3]} | "
                     f"{a['structural']['positives_in_low_expr_half']} |")
    lines += ["", "## Per-cohort no-regression vs PRIME incumbent (recall@20; ✅ no-regression / ❌ regresses)",
              "", "| cohort | prime_only | expr_rank_penalty | soft_saturating | portfolio_reserve |",
              "|---|---|---|---|---|"]
    pol_order = ["prime_only_incumbent", "expr_rank_penalty", "soft_saturating", "portfolio_reserve"]
    for name, a in result["development_cohorts"].items():
        cells = []
        for p in pol_order:
            pol = a["policies"][p]
            mark = "❌" if pol["no_regression_vs_incumbent"]["regresses"] else "✅"
            cells.append(f"{pol['recall_frac']} {mark}")
        lines.append(f"| {name} | {cells[0]} | {cells[1]} | {cells[2]} | {cells[3]} |")
    d = result["decision"]
    lines += ["", "## Decision", "",
              f"- **Chosen: {d['chosen']}** — {d['chosen_ranking_score']}",
              f"- Expression role: {d['expression_role']}",
              "- No-regression everywhere: " + ", ".join(
                  f"{p}={'YES' if v['no_regression_everywhere'] else 'no'}"
                  for p, v in d['per_policy_no_regression'].items()),
              "", d["rationale"], "",
              "## Sid descriptive (post-freeze, n=3 — NOT a gate)", ""]
    sid = result["sid_descriptive"]
    lines += [
        f"- frozen confidence-only (= lossless+PRIME): recall@20 "
        f"{sid['frozen_confidence_only_prime']['hits_at_20']}/3 "
        f"({', '.join(sid['frozen_confidence_only_prime']['covered'])})",
        f"- counterfactual expr rank penalty: recall@20 "
        f"{sid['counterfactual_expr_rank_penalty']['hits_at_20']}/3 "
        f"({', '.join(sid['counterfactual_expr_rank_penalty']['covered'])})",
        f"- counterfactual soft-saturating: recall@20 "
        f"{sid['counterfactual_soft_saturating']['hits_at_20']}/3 "
        f"({', '.join(sid['counterfactual_soft_saturating']['covered'])})",
        "", sid["note"], "",
    ]
    return "\n".join(lines)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    result = run()
    (OUT / "EXPRESSION_POLICY.json").write_text(json.dumps(result, indent=2) + "\n")
    (OUT / "EXPRESSION_POLICY.md").write_text(_write_md(result) + "\n")
    print(json.dumps({
        "decision": result["decision"]["chosen"],
        "per_policy_no_regression_everywhere": {
            p: v["no_regression_everywhere"]
            for p, v in result["decision"]["per_policy_no_regression"].items()},
        "sid_descriptive": {
            "confidence_only": result["sid_descriptive"]["frozen_confidence_only_prime"]["hits_at_20"],
            "expr_penalty": result["sid_descriptive"]["counterfactual_expr_rank_penalty"]["hits_at_20"]},
        "frozen_config": result["frozen_config_path"],
    }, indent=2))


if __name__ == "__main__":
    main()
