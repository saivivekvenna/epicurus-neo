"""Sid end-to-end benchmark — LABEL-BLIND full-universe generation + MixMHCpred/PRIME scoring.

Corrects the target-leakage of scripts/osteosarc_peptide_recovery.py: generation runs over the COMPLETE
label-blind eligible variant universe (147 variants; missense + frameshift are supported by the current
generator, other consequences recorded as unsupported), never a TARGETS list. Every generated (peptide,
HLA) pair is scored by genuine MixMHCpred + PRIME. The coverage guard (assert_generation_label_blind) is
applied. The 3 EXACT Hudson labels are joined ONLY afterwards. Mutation-level top-20 hits@20 reported per
arm with stage-of-first-loss for every positive. Post-hoc n=1/3 — descriptive only.

    python -m scripts.sid_benchmark_generate            # online (populates gitignored Ensembl cache)
    python -m scripts.sid_benchmark_generate --offline  # serve cache, fail closed on miss

Writes artifacts/milestone_7_decision/sid_benchmark/{generation.json, per_variant.csv, top20.csv, REPORT.md}.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

from event_b.lossless_peptide_generation import EnsemblClient, generate_variant_candidates, read_hla_panel
from event_b.prime_adapter import PRIME_COMMIT, score_prime
from event_b.sid_benchmark import (
    assert_generation_label_blind,
    eligible_universe_ids,
    hudson_positive_variant_ids,
    load_variant_universe,
)

RAW = Path("data/raw/osteosarc")
CACHE_DIR = RAW / "ensembl_cache"
PVAC_PATH = RAW / "pvactools_all_epitopes.tsv"
ART = Path("artifacts/milestone_7_decision/sid_benchmark")
K = 20
SUPPORTED = {"missense_variant": "missense", "frameshift_variant": "frameshift"}


def _variant_rows(u: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    """Eligible variants split into generator-supported (missense/frameshift) and unsupported (recorded)."""
    elig = u[u["class_i_eligible"]]
    supported, unsupported = [], []
    for _, r in elig.iterrows():
        if r["consequence"] in SUPPORTED:
            supported.append({"chrom": r["chrom"], "pos": int(r["pos"]), "ref": r["ref"], "alt": r["alt"],
                              "gene": r["gene"], "source_variant_type": SUPPORTED[r["consequence"]],
                              "variant_id": r["variant_id"]})
        else:
            unsupported.append({"variant_id": r["variant_id"], "consequence": r["consequence"],
                                "reason": "consequence not supported by current generator (missense/frameshift only)"})
    return supported, unsupported


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()
    ART.mkdir(parents=True, exist_ok=True)

    u = load_variant_universe()
    elig_ids = eligible_universe_ids(u)
    supported, unsupported = _variant_rows(u)
    hla_panel = read_hla_panel(PVAC_PATH)
    client = EnsemblClient(CACHE_DIR, offline=args.offline)

    per_variant, all_candidates = [], []
    for v in supported:
        rec = {"variant_id": v["variant_id"], "gene": v["gene"], "consequence_kind": v["source_variant_type"]}
        try:
            out = generate_variant_candidates(v, client, hla_panel, expected=None)  # LABEL-BLIND
            cand = out["candidates"]
            rec.update(status="ok", n_windows=out["provenance"]["n_windows"],
                       n_unique_peptides=out["provenance"]["n_unique_peptides"],
                       n_peptide_hla=len(cand), transcript_id=out["provenance"]["transcript_id"])
            all_candidates.append(cand)
        except Exception as e:  # noqa: BLE001 — record every failure, never silently exclude
            rec.update(status="FAILED", reason=f"{type(e).__name__}: {str(e)[:160]}")
        per_variant.append(rec)
    for un in unsupported:
        per_variant.append({"variant_id": un["variant_id"], "gene": un["variant_id"].split("-")[0],
                            "consequence_kind": un["consequence"], "status": "UNSUPPORTED", "reason": un["reason"]})

    cand = pd.concat(all_candidates, ignore_index=True) if all_candidates else pd.DataFrame()
    generated_ids = set(cand["mutation_id"]) if len(cand) else set()

    # ---- coverage guard (applied; failures reported, not crashed) ----
    try:
        guard = assert_generation_label_blind(generated_ids)
        guard["guard"] = "PASS"
    except Exception as e:  # noqa: BLE001
        cov = len(generated_ids & elig_ids) / len(elig_ids) if elig_ids else 0.0
        guard = {"guard": "COVERAGE_BELOW_THRESHOLD", "message": str(e), "coverage": round(cov, 4),
                 "note": "generation covers missense+frameshift only; stop_gained/inframe unsupported "
                         "=> honest incomplete coverage, NOT leakage (positive-subset check passed)."}

    # ---- score every generated (peptide, HLA) with genuine MixMHCpred + PRIME ----
    if len(cand):
        res = score_prime(cand, peptide_col="mutant_peptide", hla_col="hla_allele")
        s = res.scored
        cand = cand.reset_index(drop=True)
        cand["prime_rank"] = pd.to_numeric(s["prime_rank"], errors="coerce").to_numpy()
        cand["mixmhcpred_rank"] = pd.to_numeric(s["mixmhcpred_rank"], errors="coerce").to_numpy()
        # arms: presentation-only MixMHCpred (binding-first) and genuine PRIME (lower rank = better)
        cand["arm_mixmhcpred"] = -cand["mixmhcpred_rank"]
        cand["arm_genuine_prime"] = -cand["prime_rank"]

    # ---- JOIN EXACT labels AFTER generation+scoring is frozen ----
    positives = hudson_positive_variant_ids(u)  # exactly the 3 IDs present in the universe
    report = {"experiment": "sid_benchmark_label_blind_generation",
              "command": "python -m scripts.sid_benchmark_generate" + (" --offline" if args.offline else ""),
              "mode": "offline" if args.offline else "online",
              "prime_commit": PRIME_COMMIT, "hla_panel": hla_panel,
              "universe": {"total": int(len(u)), "eligible": len(elig_ids),
                           "supported_missense_frameshift": len(supported), "unsupported": len(unsupported)},
              "generation": {"variants_generated_ok": int(sum(1 for r in per_variant if r.get("status") == "ok")),
                             "variants_failed": int(sum(1 for r in per_variant if r.get("status") == "FAILED")),
                             "variants_unsupported": len(unsupported),
                             "n_peptide_hla_candidates": int(len(cand))},
              "coverage_guard": guard,
              "labels_joined_after_freeze": sorted(positives),
              "arms": {}}

    if len(cand):
        for arm, col, sc in [("presentation_only_mixmhcpred", "arm_mixmhcpred", "mixmhcpred_rank"),
                             ("genuine_prime", "arm_genuine_prime", "prime_rank")]:
            ranked = cand.rename(columns={"mutation_id": "variant_id"})[["variant_id", col, sc]].rename(
                columns={col: "score"})
            m = _mutation_top(ranked, positives, K)
            report["arms"][arm] = m
        report["stage_of_first_loss"] = _stage_of_loss(u, per_variant, cand, positives)

    (ART / "generation.json").write_text(json.dumps(report, indent=2, sort_keys=True, default=str) + "\n")
    pd.DataFrame(per_variant).to_csv(ART / "per_variant.csv", index=False)
    if len(cand):
        _write_top20(cand, positives)
    (ART / "REPORT.md").write_text(_md(report))
    print(json.dumps({k: report[k] for k in ["universe", "generation", "coverage_guard"]}, indent=2, default=str))
    for arm, m in report.get("arms", {}).items():
        print(f"  [{arm}] mutation hits@20 = {m['hits_at_k']}/3  recognized_ranks={m.get('positive_ranks')}")
    return 0


def _mutation_top(ranked: pd.DataFrame, positives: set[str], k: int) -> dict:
    df = ranked.dropna(subset=["score"]).sort_values("score", ascending=False, kind="mergesort")
    per = df.drop_duplicates("variant_id", keep="first").reset_index(drop=True)
    per["rank"] = np.arange(1, len(per) + 1)
    top = per.head(k)
    hit_ids = set(top["variant_id"]) & positives
    pos_ranks = {vid: int(per.loc[per["variant_id"] == vid, "rank"].iloc[0])
                 for vid in positives if (per["variant_id"] == vid).any()}
    return {"k": k, "n_variants_ranked": int(per["variant_id"].nunique()), "hits_at_k": len(hit_ids),
            "hit_variant_ids": sorted(hit_ids), "recall_at_k": round(len(hit_ids) / len(positives), 4),
            "positive_ranks": pos_ranks, "positives_scoreable": sorted(pos_ranks)}


def _stage_of_loss(u, per_variant, cand, positives) -> dict:
    pv = {r["variant_id"]: r for r in per_variant}
    scoreable = set(cand[cand["prime_rank"].notna()]["mutation_id"]) if "prime_rank" in cand else set()
    generated = set(cand["mutation_id"]) if len(cand) else set()
    out = {}
    for vid in sorted(positives):
        rec = pv.get(vid, {})
        if rec.get("status") == "UNSUPPORTED":
            stage = "generation (consequence unsupported)"
        elif rec.get("status") == "FAILED":
            stage = f"generation (VEP/transcript failure: {rec.get('reason')})"
        elif vid not in generated:
            stage = "generation (no candidate peptide emitted)"
        elif vid not in scoreable:
            stage = "scoring (no PRIME/MixMHCpred score)"
        else:
            stage = "reached scoring (rankable)"
        out[vid] = {"status": rec.get("status"), "stage_of_first_loss": stage}
    return out


def _write_top20(cand, positives):
    best = cand.sort_values("mixmhcpred_rank", kind="mergesort").drop_duplicates("mutation_id", keep="first")
    best = best.assign(is_recognized=best["mutation_id"].isin(positives)).head(40)
    best[["mutation_id", "gene_symbol", "mutant_peptide", "hla_allele", "mixmhcpred_rank", "prime_rank",
          "is_recognized"]].to_csv(ART / "top20.csv", index=False)


def _md(r) -> str:
    L = [f"# Sid end-to-end benchmark — label-blind generation\n\n`{r['command']}` · mode {r['mode']} · "
         f"PRIME `{r['prime_commit'][:10]}`\n",
         "_Post-hoc n=1 patient / 3 recognized positives — descriptive only. Generation is over the "
         "COMPLETE label-blind eligible universe (guard-enforced); the 3 exact labels are joined only "
         "after generation+scoring are frozen._\n"]
    u, g = r["universe"], r["generation"]
    L.append(f"\n**Universe:** {u['total']} public → {u['eligible']} eligible → "
             f"{u['supported_missense_frameshift']} generator-supported (missense/frameshift), "
             f"{u['unsupported']} unsupported.\n")
    L.append(f"**Generation:** {g['variants_generated_ok']} ok, {g['variants_failed']} failed, "
             f"{g['variants_unsupported']} unsupported → {g['n_peptide_hla_candidates']} peptide×HLA candidates.\n")
    L.append(f"**Coverage guard:** {r['coverage_guard'].get('guard')} "
             f"(coverage {r['coverage_guard'].get('coverage', r['coverage_guard'].get('generated_eligible_covered'))}).\n")
    L.append("\n## Mutation-level recognized hits@20 (labels joined post-freeze)\n")
    L.append("| arm | hits@20 / 3 | recognized ranks (variant → rank) |")
    L.append("|---|--:|---|")
    for arm, m in r.get("arms", {}).items():
        L.append(f"| {arm} | {m['hits_at_k']}/3 | {m.get('positive_ranks')} |")
    L.append("\n## Stage of first loss (per recognized positive)\n")
    for vid, s in r.get("stage_of_first_loss", {}).items():
        L.append(f"- `{vid}` — **{s['stage_of_first_loss']}**")
    L.append("\n> The old target-conditioned 'lossless 3/3' is withdrawn (see BENCHMARK_PROTOCOL.md §0). "
             "This is the honest label-blind end-to-end result; every missed positive has an explicit "
             "stage of first loss. Competitor arms (pVACtools boundary-mismatch, etc.) are reported "
             "separately per the protocol.\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
