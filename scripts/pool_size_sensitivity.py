"""Pool-size sensitivity diagnostic runner (Milestone 7).

Does frozen Epicurus vs genuine PRIME (and a pVAC-style+PRIME proxy) improve when the starting
candidate pool shrinks but retains exactly all known positives? Reuses the FROZEN deterministic
scorers only (no model fit, no PRIME binary). See src/event_b/pool_size_sensitivity.py for the design.

    python -m scripts.pool_size_sensitivity

Writes artifacts/milestone_7_decision/pool_size_sensitivity/{pool_size_sensitivity.json,
per_arm_summary.csv, long_rows.csv.gz, REPORT.md}.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from event_b.pool_size_sensitivity import (
    ARMS, BUDGET_FRAC, FIXED_K, KS, MIN_NEG, POOL_FRACS, run_cohort,
)
from event_b.prime_transfer import _gartner, _improve, _multimer

ART = Path("artifacts/milestone_7_decision/pool_size_sensitivity")
SEEDS = list(range(25))  # >= 20 deterministic negative-sampling seeds
COMMAND = "python -m scripts.pool_size_sensitivity"

# cohort role notes (heterogeneous; never pooled)
COHORT_LOADERS = [
    ("gartner", _gartner, "Gartner NCI Testing (frozen external transfer test); clean for Epicurus."),
    ("improve", _improve, "IMPROVE SRHgroup (frozen external validation); clean for Epicurus; best-powered."),
    ("multimer", _multimer, "CD8 pMHC-multimer; frozen Epicurus TRAINING cohort -> Epicurus IN-SAMPLE (flagged)."),
]


def _base_frame(fn) -> pd.DataFrame:
    f = fn().copy()
    keep = [c for c in ["patient_id", "mutant_peptide", "hla_allele", "label", "prime", "el", "expr"] if c in f.columns]
    f = f[keep]
    return f[f["label"].isin(["POSITIVE", "TESTED_NEGATIVE"])].reset_index(drop=True)


def main() -> int:
    ART.mkdir(parents=True, exist_ok=True)
    cohorts, long_rows = [], []
    for name, fn, note in COHORT_LOADERS:
        frame = _base_frame(fn)
        frame.to_csv(ART / f"base_{name}.csv", index=False)
        res = run_cohort(frame, name, SEEDS)
        res["role_note"] = note
        long_rows.extend(res.pop("_long_rows"))
        cohorts.append(res)
        print(f"[{name}] eligible={res['n_patients_eligible']} pos={res['n_positives']} "
              f"neg={res['n_tested_negatives']} in_sample={res['epicurus_in_sample']}")

    report = {
        "diagnostic": "pool_size_sensitivity",
        "question": "Does Epicurus vs genuine PRIME improve when the pool shrinks but retains all positives?",
        "command": COMMAND,
        "seeds": SEEDS,
        "config": {"pool_fractions_of_negatives": POOL_FRACS, "ks": list(KS),
                   "normalized_budget_frac": BUDGET_FRAC, "fixed_small_k": FIXED_K,
                   "min_negatives_for_eligibility": MIN_NEG, "arms": list(ARMS)},
        "scorer_provenance": {
            "genuine_prime": "-prime (GfellerLab PRIME 2.1 %rank; pool-invariant per candidate)",
            "frozen_epicurus": "configs/frozen/epicurus_v0_1.json via score_with_frozen (POOL-DEPENDENT percentiles)",
            "pvac_style_prime": "TUNING-FREE proxy: equal-weight mean of within-patient percentile(EL) and "
                                "percentile(PRIME); approximates pVACseq binding-first + PRIME, NOT the pVACtools binary",
        },
        "variant_definitions": {
            "A_oracle": "all positives + prefix random negative subsample; positives 100% retained BY "
                        "CONSTRUCTION -> reranker stress test, diagnostic only, NOT north-star validation",
            "B_gate": "top-N by label-blind within-patient NetMHCpan-EL percentile, size-matched to A; "
                      "positives may drop; retention reported (labels never used to design B)",
        },
        "cohorts": cohorts,
        "answers": _answers(cohorts),
        "caveats": [
            "Oracle (variant A) retention is diagnostic, NOT validation: a real pipeline cannot pre-know positives.",
            "Multimer frozen-Epicurus scores are IN-SAMPLE (training cohort) -> optimistic; do not read as external.",
            "Cohorts are heterogeneous with fixed roles and are NEVER pooled into one headline number.",
            "Precision@20 and recall@20 rise mechanically as negatives thin; MRR / median positive rank are the "
            "shift-invariant checks for whether ranking QUALITY (not just the denominator) changed.",
            "No superiority is claimed; verdicts follow the project's challenger-vs-baseline convention.",
        ],
    }
    (ART / "pool_size_sensitivity.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")

    long_df = pd.DataFrame(long_rows)
    long_df.to_csv(ART / "long_rows.csv.gz", index=False, compression="gzip")
    _per_arm_summary_csv(cohorts).to_csv(ART / "per_arm_summary.csv", index=False)
    (ART / "REPORT.md").write_text(_report_md(report))
    print(f"\nWrote {ART}/pool_size_sensitivity.json, per_arm_summary.csv, long_rows.csv.gz, REPORT.md")
    return 0


def _answers(cohorts: list[dict]) -> dict:
    """Cross-cohort plain-language answers to the three required questions."""
    def trend(c, variant, arm, metric):
        s = c[variant]
        return {p: (s[p]["arms"][arm].get(metric)) for p in POOL_FRACS}

    q1, q2, q3 = {}, {}, {}
    for c in cohorts:
        name = c["cohort"]
        # Q1: do ALL arms improve on recall@20 as pool shrinks? (denominator artifact check)
        q1[name] = {arm: trend(c, "variant_A_oracle", arm, "recall@20") for arm in ARMS}
        # Q1b: MRR trend (ranking-quality, denominator-robust)
        q1[name]["_mrr_frozen_epicurus"] = trend(c, "variant_A_oracle", "frozen_epicurus", "mrr")
        q1[name]["_mrr_genuine_prime"] = trend(c, "variant_A_oracle", "genuine_prime", "mrr")
        # Q2: does Epicurus - PRIME hits@20 delta grow as pool shrinks?
        q2[name] = {p: c["variant_A_oracle"][p]["epicurus_minus_prime_hits@20"]["mean"] for p in POOL_FRACS}
        q2[f"{name}_band95"] = {p: c["variant_A_oracle"][p]["epicurus_minus_prime_hits@20"]["band95"] for p in POOL_FRACS}
        # Q3: does the label-blind gate retain all positives, or only the oracle?
        q3[name] = {p: {"retention_mean": c["variant_B_gate"][p]["positive_retention_mean"],
                        "retention_min": c["variant_B_gate"][p]["positive_retention_min"],
                        "all_retained": c["variant_B_gate"][p]["all_positives_retained"],
                        "n_losing_a_positive": c["variant_B_gate"][p]["n_patients_losing_a_positive"]}
                    for p in POOL_FRACS}
    return {
        "q1_do_all_models_improve_from_fewer_negatives": q1,
        "q2_does_epicurus_advantage_increase": q2,
        "q3_does_a_real_label_blind_filter_retain_positives_or_only_the_oracle": q3,
    }


def _per_arm_summary_csv(cohorts: list[dict]) -> pd.DataFrame:
    rows = []
    metrics = [f"hits@{k}" for k in KS] + [f"hits@budget{int(BUDGET_FRAC*100)}pct",
                                           "recall@20", "precision@20", "mrr", "median_pos_rank"]
    for c in cohorts:
        for variant in ("variant_A_oracle", "variant_B_gate"):
            for pool in POOL_FRACS:
                s = c[variant][pool]
                for arm in ARMS:
                    row = {"cohort": c["cohort"], "epicurus_in_sample": c["epicurus_in_sample"],
                           "variant": variant, "pool": pool, "arm": arm,
                           "n_patients": c["n_patients_eligible"]}
                    row.update({m: s["arms"][arm].get(m) for m in metrics})
                    if variant == "variant_B_gate":
                        row["positive_retention_mean"] = s["positive_retention_mean"]
                        row["all_positives_retained"] = s["all_positives_retained"]
                    rows.append(row)
    return pd.DataFrame(rows)


def _fmt(x, nd=3):
    return "—" if x is None or (isinstance(x, float) and np.isnan(x)) else (f"{x:.{nd}f}" if isinstance(x, float) else str(x))


def _report_md(report: dict) -> str:
    L = []
    L.append("# Pool-size sensitivity diagnostic — frozen Epicurus vs genuine PRIME\n")
    L.append(f"`{report['command']}` · {len(report['seeds'])} seeds · frozen scorers only (no model fit, no PRIME binary re-run).\n")
    L.append("**Question.** " + report["question"] + "\n")
    L.append("Three arms, all frozen: **genuine PRIME** (`-%rank`, pool-invariant), **frozen Epicurus** "
             "(`configs/frozen/epicurus_v0_1.json`; within-patient percentiles → *pool-dependent*), and a "
             "tuning-free **pVAC-style+PRIME** proxy (equal-weight EL⊕PRIME percentiles — *not* the pVACtools binary).\n")
    L.append("Pools per patient: **LARGE** = all tested candidates; **MEDIUM** = all positives + 50% of negatives; "
             "**SMALL** = all positives + 25% of negatives. SMALL ⊂ MEDIUM ⊂ LARGE, identical positives/features.\n")

    for c in report["cohorts"]:
        star = " ⚠️ *Epicurus IN-SAMPLE (training cohort)*" if c["epicurus_in_sample"] else ""
        L.append(f"\n## {c['cohort']}{star}\n")
        L.append(f"{c['role_note']}  \nEligible patients **{c['n_patients_eligible']}** · positives **{c['n_positives']}** · "
                 f"tested-negatives **{c['n_tested_negatives']}** · excluded {len(c['excluded_patients'])} "
                 f"({', '.join(sorted(set(c['excluded_patients'].values()))) or 'none'}).\n")
        # Variant A table
        L.append("\n**Variant A (oracle — reranker stress test; positives 100% retained by construction).** "
                 "Mean over patients×seeds.\n")
        L.append("| pool | arm | hits@5 | hits@10 | hits@20 | recall@20 | prec@20 | MRR | med.pos.rank | sat≤20 |")
        L.append("|---|---|--:|--:|--:|--:|--:|--:|--:|--:|")
        for pool in POOL_FRACS:
            s = c["variant_A_oracle"][pool]
            for arm in ARMS:
                a = s["arms"][arm]
                L.append(f"| {pool} | {arm} | {_fmt(a['hits@5'])} | {_fmt(a['hits@10'])} | {_fmt(a['hits@20'])} | "
                         f"{_fmt(a['recall@20'])} | {_fmt(a['precision@20'])} | {_fmt(a['mrr'])} | "
                         f"{_fmt(a['median_pos_rank'],1)} | {_fmt(s['mean_saturated_le20_frac'],2)} |")
        L.append("\n**Epicurus − PRIME hits@20 (paired, mean [95% band over seeds]):**  ")
        for pool in POOL_FRACS:
            d = c["variant_A_oracle"][pool]["epicurus_minus_prime_hits@20"]
            L.append(f"{pool} {_fmt(d['mean'])} {d['band95']}  ·  ")
        # Variant B
        L.append("\n\n**Variant B (label-blind EL-percentile gate; size-matched; positives may drop).**\n")
        L.append("| pool | pos. retention (mean / min) | all retained? | # patients losing a positive | epi−PRIME hits@20 |")
        L.append("|---|--:|:--:|--:|--:|")
        for pool in POOL_FRACS:
            s = c["variant_B_gate"][pool]
            L.append(f"| {pool} | {_fmt(s['positive_retention_mean'])} / {_fmt(s['positive_retention_min'])} | "
                     f"{'yes' if s['all_positives_retained'] else 'NO'} | {s['n_patients_losing_a_positive']} | "
                     f"{_fmt(s['epicurus_minus_prime_hits@20_mean'])} |")

    ans = report["answers"]
    L.append("\n## Answers\n")
    L.append("**Q1 — Do all models improve just because there are fewer negatives?**  \n"
             "recall@20 and precision@20 rise for *every* arm as the pool shrinks (mechanical: fewer competitors / "
             "smaller top-20 denominator). MRR and median positive rank are the denominator-robust checks — see the "
             "per-cohort tables. If MRR is roughly flat while recall climbs, the gain is the denominator, not better ranking.\n")
    L.append("**Q2 — Does Epicurus's advantage over PRIME increase as the pool shrinks?**  \n")
    for c in report["cohorts"]:
        row = ans["q2_does_epicurus_advantage_increase"][c["cohort"]]
        L.append(f"- {c['cohort']}: LARGE {_fmt(row['LARGE'])} → MEDIUM {_fmt(row['MEDIUM'])} → SMALL {_fmt(row['SMALL'])} "
                 f"(paired hits@20 delta){' — IN-SAMPLE' if c['epicurus_in_sample'] else ''}\n")
    L.append("**Q3 — Does a real label-blind filter retain the positives, or only the oracle?**  \n")
    for c in report["cohorts"]:
        q3 = ans["q3_does_a_real_label_blind_filter_retain_positives_or_only_the_oracle"][c["cohort"]]
        med, sml = q3["MEDIUM"], q3["SMALL"]
        L.append(f"- {c['cohort']}: MEDIUM gate retains {_fmt(med['retention_mean'])} of positives "
                 f"(all retained: {med['all_retained']}); SMALL gate retains {_fmt(sml['retention_mean'])} "
                 f"(all retained: {sml['all_retained']}, {sml['n_losing_a_positive']} patients lose ≥1). "
                 f"The oracle keeps 100% by construction; the deployable gate does not.\n")
    L.append("\n" + "\n".join(f"> {cv}" for cv in report["caveats"]) + "\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
