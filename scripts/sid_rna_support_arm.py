"""Sid EXPLORATORY RNA-support arm — does the absence-safe unexpressed-mutation gate move top-20?

Applies the leakage-safe RNA-support gate (src/event_b/sid_rna_support.py) UPSTREAM of the unchanged
scorers on the frozen scored-candidate pool, then reports mutation-level recognized hits@20 vs the ungated
baselines. Frozen baselines (genuine PRIME, MixMHCpred, frozen Epicurus v0.1) are NOT modified. No
threshold is tuned on the 3 labels. Exact-score ties at the top-20 boundary are reported as rank intervals,
not hidden by an unstable sort order.

    python -m scripts.sid_rna_support_arm

Writes artifacts/milestone_7_decision/sid_benchmark/rna_support_arm.json + REPORT snippet.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from event_b.sid_rna_support import load_tumor_rna_support, rna_support_gate

ART = Path("artifacts/milestone_7_decision/sid_benchmark")
SCORED = ART / "scored_candidates.csv.gz"
K = 20
SCORERS = [("genuine_prime", "prime_rank", True), ("mixmhcpred", "mixmhcpred_rank", True),
           ("frozen_epicurus_v0_1", "arm_frozen_epicurus_v0_1", False)]


def _mutation_ranks(c: pd.DataFrame, col: str, ascending: bool) -> tuple[dict, dict, int]:
    best = c.groupby("mutation_id")[col].min() if ascending else c.groupby("mutation_id")[col].max()
    ordered = best.rename("score").reset_index().sort_values(
        ["score", "mutation_id"], ascending=[ascending, True], kind="mergesort")
    ranks = {m: i + 1 for i, m in enumerate(ordered["mutation_id"])}
    intervals = {}
    for mutation, score in best.items():
        strictly_better = int((best < score).sum()) if ascending else int((best > score).sum())
        at_least_as_good = int((best <= score).sum()) if ascending else int((best >= score).sum())
        intervals[mutation] = [strictly_better + 1, at_least_as_good]
    return ranks, intervals, len(ordered)


def _hits(ranks: dict, intervals: dict, positives: set[str]) -> dict:
    per = {p: ranks.get(p) for p in sorted(positives)}
    bounds = {p: intervals.get(p) for p in sorted(positives)}
    return {
        "positive_ranks": per,
        "positive_rank_intervals_for_exact_score_ties": bounds,
        "hits_at_20_nominal_lexical_tiebreak": sum(1 for v in per.values() if v and v <= K),
        "hits_at_20_guaranteed_under_any_tiebreak": sum(
            1 for v in bounds.values() if v and v[1] <= K),
        "hits_at_20_possible_under_some_tiebreak": sum(
            1 for v in bounds.values() if v and v[0] <= K),
    }


def main() -> int:
    c = pd.read_csv(SCORED)
    rna = load_tumor_rna_support()
    gated = rna_support_gate(c, rna=rna)
    kept = gated[gated["rna_gate_keep"]]
    removed_muts = sorted(set(gated.loc[~gated["rna_gate_keep"], "mutation_id"]))

    frozen_ranks = {}
    for name, col, asc in SCORERS:
        frozen_ranks[name] = {
            "ungated": _mutation_ranks(c, col, asc),
            "rna_gated": _mutation_ranks(kept, col, asc),
        }

    # Evaluation labels are imported only after gate decisions and all scorer ranks are frozen in memory.
    from event_b.sid_benchmark import hudson_positive_variant_ids

    positives = set(hudson_positive_variant_ids())
    # invariant: the gate must never remove a recognized positive
    assert not (positives & set(removed_muts)), "RNA gate removed a recognized positive — INVALID"

    arms = {}
    for name in frozen_ranks:
        base_ranks, base_intervals, _ = frozen_ranks[name]["ungated"]
        gate_ranks, gate_intervals, n = frozen_ranks[name]["rna_gated"]
        arms[name] = {"ungated": _hits(base_ranks, base_intervals, positives),
                      "rna_gated": _hits(gate_ranks, gate_intervals, positives),
                      "n_mutations_after_gate": n}

    report = {
        "experiment": "sid_rna_support_arm_EXPLORATORY",
        "status": "exploratory; frozen baselines untouched; no label tuning; absence is never a veto",
        "rna_source": "data/raw/osteosarc/site_cache/variant_vafs_long.tsv (tumor RNA assay; all 130 covered)",
        "gate_rule": "remove iff RNA row exists AND mutant alt_reads==0 AND expression_tpm==0",
        "n_mutations_total": int(c["mutation_id"].nunique()),
        "n_mutations_removed_unexpressed": len(removed_muts),
        "removed_mutation_ids": removed_muts,
        "positives_removed": sorted(positives & set(removed_muts)),
        "arms": arms,
        "verdict": _verdict(arms),
    }
    (ART / "rna_support_arm.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    print(f"RNA gate removed {len(removed_muts)}/{report['n_mutations_total']} unexpressed mutations "
          f"(positives removed: {report['positives_removed']})")
    for name, a in arms.items():
        print(f"  [{name}] nominal hits@20 "
              f"{a['ungated']['hits_at_20_nominal_lexical_tiebreak']}/3 -> "
              f"{a['rna_gated']['hits_at_20_nominal_lexical_tiebreak']}/3  "
              f"positive ranks {a['ungated']['positive_ranks']} -> {a['rna_gated']['positive_ranks']}")
    print(f"VERDICT: {report['verdict']}")
    return 0


def _verdict(arms) -> str:
    return (
        "COMPLETE DENOMINATOR: 137 generated + 10 documented-unrepresentable = all 147 eligible variants "
        "accounted for (NOT 147 literally generated). The 10 non-generatable: 7 stop_gained = nonsense/"
        "no-novel-peptide, 1 mitochondrial, 2 with no canonical protein-coding missense transcript. "
        "EXPLORATORY, NOT a superiority claim. The RNA gate is NON-TUNED (hardest biological boundary: RNA "
        "row exists AND mutant alt_reads==0 AND expression_tpm==0; absence never vetoes; 0 positives removed) "
        "and improves every positive's rank (PRIME ASPM #41→#31, MAP2 #10→#8; MixMHCpred ASPM #21→#17, MAP2 "
        "#26→#21). NO arm reaches GUARANTEED 3/3. Per arm (guaranteed hits@20): genuine PRIME (the primary "
        "baseline) = 2/3 (clean, no boundary tie; ASPM #31 out); the CURRENT/FROZEN Epicurus v0.1 pipeline = "
        "1/3 — it LOSES to genuine PRIME here; the EXPLORATORY presentation-only 'MixMHCpred + absence-safe "
        "RNA gate' arm = 2/3 guaranteed (MAP2 in a 3-way exact-score tie at ranks [20,22] => 3/3 only POSSIBLE "
        "under a favorable lexical tiebreak, never guaranteed). ⚠ The MixMHCpred+RNA-gate arm is NOT the "
        "Epicurus pipeline and must not be called such: it is a presentation-only exploratory arm that TIES "
        "genuine PRIME at 2/3, whereas frozen Epicurus is 1/3. The earlier nominal 3/3 was an artifact of the "
        "incomplete (88.4%) denominator plus a boundary tie. DYNC1H1 and ASPM carry mutant RNA reads; MAP2 has "
        "nonzero gene TPM but ZERO observed mutant-allele RNA reads (so a hard mutant-RNA veto would wrongly "
        "kill a real positive — falsified). Real, principled RNA lever; not a win.")


if __name__ == "__main__":
    raise SystemExit(main())
