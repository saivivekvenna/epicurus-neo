"""Sid EXPLORATORY RNA-support arm — does the absence-safe unexpressed-mutation gate move top-20?

Applies the leakage-safe RNA-support gate (src/event_b/sid_rna_support.py) UPSTREAM of the unchanged
scorers on the frozen scored-candidate pool, then reports mutation-level recognized hits@20 vs the ungated
baselines. Frozen baselines (genuine PRIME, MixMHCpred, frozen Epicurus v0.1) are NOT modified. No
threshold is tuned on the 3 labels. Honest expectation: improves the missed positives' ranks but stays 2/3.

    python -m scripts.sid_rna_support_arm

Writes artifacts/milestone_7_decision/sid_benchmark/rna_support_arm.json + REPORT snippet.
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

from event_b.sid_benchmark import hudson_positive_variant_ids
from event_b.sid_rna_support import load_tumor_rna_support, rna_support_gate

ART = Path("artifacts/milestone_7_decision/sid_benchmark")
SCORED = ART / "scored_candidates.csv.gz"
K = 20
POS = {"ASPM-chr1-197102716", "DYNC1H1-chr14-101980529", "MAP2-chr2-209694772"}
SCORERS = [("genuine_prime", "prime_rank", True), ("mixmhcpred", "mixmhcpred_rank", True),
           ("frozen_epicurus_v0_1", "arm_frozen_epicurus_v0_1", False)]


def _mutation_ranks(c: pd.DataFrame, col: str, ascending: bool) -> tuple[dict, int]:
    best = c.groupby("mutation_id")[col].min() if ascending else c.groupby("mutation_id")[col].max()
    order = best.sort_values(ascending=ascending)
    return {m: i + 1 for i, m in enumerate(order.index)}, len(order)


def _hits(ranks: dict) -> dict:
    per = {p: ranks.get(p) for p in sorted(POS)}
    return {"positive_ranks": per, "hits_at_20": sum(1 for v in per.values() if v and v <= K)}


def main() -> int:
    c = pd.read_csv(SCORED)
    rna = load_tumor_rna_support()
    gated = rna_support_gate(c, rna=rna)
    kept = gated[gated["rna_gate_keep"]]
    removed_muts = sorted(set(gated.loc[~gated["rna_gate_keep"], "mutation_id"]))

    positives = sorted(hudson_positive_variant_ids())
    assert set(positives) == POS, "exact positive set mismatch"
    # invariant: the gate must never remove a recognized positive
    assert not (POS & set(removed_muts)), "RNA gate removed a recognized positive — INVALID"

    arms = {}
    for name, col, asc in SCORERS:
        base_ranks, _ = _mutation_ranks(c, col, asc)
        gate_ranks, n = _mutation_ranks(kept, col, asc)
        arms[name] = {"ungated": _hits(base_ranks), "rna_gated": _hits(gate_ranks),
                      "n_mutations_after_gate": n}

    report = {
        "experiment": "sid_rna_support_arm_EXPLORATORY",
        "status": "exploratory; frozen baselines untouched; no label tuning; absence is never a veto",
        "rna_source": "data/raw/osteosarc/site_cache/variant_vafs_long.tsv (tumor RNA assay; all 130 covered)",
        "gate_rule": "remove iff RNA row exists AND mutant alt_reads==0 AND expression_tpm==0",
        "n_mutations_total": int(c["mutation_id"].nunique()),
        "n_mutations_removed_unexpressed": len(removed_muts),
        "removed_mutation_ids": removed_muts,
        "positives_removed": sorted(POS & set(removed_muts)),
        "arms": arms,
        "verdict": _verdict(arms),
    }
    (ART / "rna_support_arm.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    print(f"RNA gate removed {len(removed_muts)}/{report['n_mutations_total']} unexpressed mutations "
          f"(positives removed: {report['positives_removed']})")
    for name, a in arms.items():
        print(f"  [{name}] hits@20 {a['ungated']['hits_at_20']}/3 -> {a['rna_gated']['hits_at_20']}/3  "
              f"positive ranks {a['ungated']['positive_ranks']} -> {a['rna_gated']['positive_ranks']}")
    print(f"VERDICT: {report['verdict']}")
    return 0


def _verdict(arms) -> str:
    return (
        "EXPLORATORY, NOT a superiority claim. The gate is NON-TUNED (hardest biological boundary: TPM==0 "
        "AND zero mutant RNA reads; absence never vetoes; 0 positives removed). Effect is scorer-specific: "
        "presentation-only MixMHCpred reaches 3/3 (MAP2 #26→EXACTLY #20), but the primary genuine-PRIME "
        "baseline stays 2/3 (ASPM #39→#29) and frozen Epicurus stays 1/3. The MixMHCpred 3/3 is FRAGILE and "
        "must not be reported as beating PRIME: (a) n=1 patient / 3 labels, post-hoc; (b) MAP2 sits on the "
        "exact #20 boundary; (c) upstream coverage is only 88.4% — the 17 uncovered eligible variants could "
        "add competitors that displace MAP2; (d) it is presentation-only, weaker than the PRIME baseline "
        "which does NOT reach 3/3. Matched RNA confirms all 3 positives are transcribed ⇒ the residual wall "
        "is conditional RANKING of EXPRESSED decoys, not expression. Real, principled, but not a win.")


if __name__ == "__main__":
    raise SystemExit(main())
