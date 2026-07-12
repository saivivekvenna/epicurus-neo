"""Runner — input-only lossless peptide recovery for osteosarc.com / Sid (EXPLORATORY).

Frozen exploratory protocol:
``docs/superpowers/specs/2026-07-12-osteosarc-lossless-peptide-recovery-exploratory-protocol.md``
(paired copy ``artifacts/milestone_7_decision/peptide_recovery/EXPLORATORY_PROTOCOL.md``).

POST-HOC STATUS (read first): Sid's structural audit motivated this generator and its feasibility
scores were already inspected, so this is a reproducibility / reachability diagnostic on the one
patient that motivated it — NOT preregistered, blind, or independent validation. This run does NOT
show Epicurus beats PRIME. The selection score is genuine PRIME itself (``genuine_prime = -PRIME
%rank``); any top-20 gain is because better candidate GENERATION lets genuine PRIME score targets it
previously never received.

Pipeline (frozen): raw GRCh38 allele -> Ensembl VEP MANE/canonical -> Ensembl CDS/protein -> all
standard-AA 8-14 mutation/novel-frame windows -> patient HLA panel (read from pVAC) -> genuine PRIME
-> union with the original pVAC candidate set (PRIME from ``_cache_prime.tsv``) on a stable genomic
candidate identity -> frozen evidence router + route-aware top-20 (``genuine_prime``). The recognition
outcome labels are joined ONLY after ranking, isolated as evaluation.

Run online (populates the gitignored Ensembl cache) then ``--offline`` (serves the cache, fails
closed on a miss); both produce identical recovered candidates.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from event_b.lossless_peptide_generation import (
    POLICY_ID,
    EnsemblClient,
    generate_variant_candidates,
    read_hla_panel,
    union_candidates,
)
from event_b.prime_adapter import PRIME_COMMIT, score_prime
from epicurus_neo.evidence_router import (
    DEFAULT_ROUTER_POLICY,
    route_candidates,
    select_route_aware_topk,
)

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "osteosarc"
CACHE_DIR = RAW / "ensembl_cache"
OUT = ROOT / "artifacts" / "milestone_7_decision" / "peptide_recovery"

PVAC_PATH = RAW / "pvactools_all_epitopes.tsv"
CACHE_PRIME_PATH = RAW / "_cache_prime.tsv"
VAF_PATH = RAW / "site_cache" / "variant_vafs_long.tsv"
RSEM_PATH = RAW / "rsem.2025.01.genes.results"

# Programmatic, asserted (not private label knowledge) — read back from the pVAC table too.
EXPECTED_HLA_PANEL = [
    "HLA-A*01:01", "HLA-B*08:01", "HLA-B*27:05", "HLA-C*01:02", "HLA-C*07:01",
]

# Frozen recovery targets (INPUTS only): gene + genomic locus. The raw ref/alt/VAF/provenance are
# read from the public VAF table; nothing recognition-derived is used to pick or shape a peptide.
# DYNC1H1 is the input-only positive control (already in pVAC; the generator must reproduce it).
TARGETS = [
    {"gene": "ASPM", "chrom": "chr1", "pos": 197102716, "kind": "missense", "ensg": "ENSG00000066279",
     "expected": {"transcript_id": "ENST00000367409", "mane_refseq": "NM_018136.5",
                  "protein_start": 2179, "amino_acids": "G/R", "hgvsc": "c.6535G>A"},
     "expected_windows": 77, "expected_peptide_present": "RRVRVRRTLR"},
    {"gene": "MAP2", "chrom": "chr2", "pos": 209694772, "kind": "frameshift", "ensg": "ENSG00000078018",
     "expected": {"transcript_id": "ENST00000682079", "mane_refseq": "NM_001375505.1",
                  "protein_start": 868, "amino_acids": "GYCVFNKYTV/X", "hgvsc": "c.2603_2630del"},
     "expected_windows": 259, "expected_peptide_present": "RVVPFTKAL"},
    {"gene": "DYNC1H1", "chrom": "chr14", "pos": 101980529, "kind": "missense", "ensg": "ENSG00000197102",
     "expected": {"transcript_id": "ENST00000360184", "mane_refseq": "NM_001376.5",
                  "protein_start": 314, "amino_acids": "V/I", "hgvsc": "c.940G>A"},
     "expected_windows": None, "expected_peptide_present": "KRFHATISF", "positive_control": True},
]

# ================= EVALUATION-ONLY (joined AFTER ranking, never a generator input) =================
# The three IFNgamma/TCR-expanded osteosarc.com/Sid target mutation IDs. Held out of generation and
# ranking entirely; used only to measure candidate/rankable/top-20 coverage after the ranking is
# fixed. (Public Hudson-lab result; see the reconstruction. Not read from any assay/vaccine table.)
HUDSON_EVAL_MUTATION_IDS = {
    "ASPM-chr1-197102716",
    "MAP2-chr2-209694772",
    "DYNC1H1-chr14-101980529",
}
# ===================================================================================================

POLICY = DEFAULT_ROUTER_POLICY


def _mutation_id(gene: str, chrom: str, pos: int) -> str:
    return f"{gene}-{chrom}-{int(pos)}"


# ---------------------------------------------------------------------------
# Input provenance (raw allele, VAF, multi-caller support, gene TPM)
# ---------------------------------------------------------------------------
def load_target_provenance() -> dict:
    vafs = pd.read_csv(VAF_PATH, sep="\t", low_memory=False)
    rsem = pd.read_csv(RSEM_PATH, sep="\t")
    rsem["ensg"] = rsem["gene_id"].astype(str).str.split(".").str[0]
    tpm_by_ensg = dict(zip(rsem["ensg"], pd.to_numeric(rsem["TPM"], errors="coerce")))

    provenance: dict = {}
    for target in TARGETS:
        rows = vafs[(vafs["gene"].astype(str) == target["gene"])
                    & (vafs["pos"].astype(str) == str(target["pos"]))]
        if rows.empty:
            raise ValueError(f"no VAF rows for {target['gene']} @ {target['pos']}")
        ref = str(rows["ref"].iloc[0])
        alt = str(rows["alt"].iloc[0])
        tumor = rows[rows["tissue"].astype(str).str.lower().eq("tumor")
                     & (pd.to_numeric(rows["alt_reads"], errors="coerce").fillna(0) > 0)]
        n_callers = int(tumor["pipeline"].astype(str).str.strip().replace("", pd.NA).dropna().nunique())
        n_timepoints = int(tumor["timepoint"].astype(str).str.strip().replace("", pd.NA).dropna().nunique())
        tumor_vaf = float(pd.to_numeric(tumor["vaf"], errors="coerce").max()) if len(tumor) else float("nan")
        provenance[target["gene"]] = {
            "chrom": target["chrom"], "pos": target["pos"], "ref": ref, "alt": alt,
            "n_callers": n_callers, "n_timepoints": n_timepoints, "tumor_vaf": tumor_vaf,
            "callers": sorted(tumor["pipeline"].dropna().astype(str).unique().tolist()),
            "expression_tpm": float(tpm_by_ensg.get(target["ensg"], float("nan"))),
        }
    return provenance


# ---------------------------------------------------------------------------
# Incumbent pVAC candidate set (+ PRIME from the frozen cache)
# ---------------------------------------------------------------------------
def build_pvac_frame() -> pd.DataFrame:
    pvac = pd.read_csv(PVAC_PATH, sep="\t", low_memory=False)
    cache = pd.read_csv(CACHE_PRIME_PATH, sep="\t", low_memory=False)

    frame = pd.DataFrame({
        "patient_id": "sid",
        "gene_symbol": pvac["Gene Name"].astype(str),
        "chrom": pvac["Chromosome"].astype(str),
        "pos": pd.to_numeric(pvac["Stop"], errors="coerce").astype("Int64"),  # 1-based end == SNV pos
        "ref": pvac["Reference"].astype(str),
        "alt": pvac["Variant"].astype(str),
        "mutant_peptide": pvac["MT Epitope Seq"].astype(str),
        "hla_allele": pvac["HLA Allele"].astype(str),
        "mhc_class": "I",
        "expression_tpm": pd.to_numeric(pvac["Gene Expression"], errors="coerce"),
        "candidate_source": "pvactools_2025_01",
    })
    frame["source_variant_type"] = pvac["Variant Type"].astype(str).map(
        lambda v: "SNV" if v == "missense" else v
    )
    frame["mutation_id"] = [
        _mutation_id(g, c, p) for g, c, p in zip(frame["gene_symbol"], frame["chrom"], frame["pos"])
    ]

    merged = frame.merge(
        cache[["mutant_peptide", "hla_allele", "prime_rank", "prime_score", "mixmhcpred_rank"]],
        on=["mutant_peptide", "hla_allele"], how="left",
    )
    merged["binding_percentile_rank"] = merged["mixmhcpred_rank"]
    merged["genuine_prime"] = -pd.to_numeric(merged["prime_rank"], errors="coerce")
    return merged


# ---------------------------------------------------------------------------
# Recovered candidate generation + genuine PRIME scoring
# ---------------------------------------------------------------------------
def generate_recovered(client: EnsemblClient, hla_panel: list[str], prov: dict) -> tuple[pd.DataFrame, dict]:
    frames, gen_provenance, per_variant = [], {}, {}
    for target in TARGETS:
        variant = {
            "gene": target["gene"], "chrom": target["chrom"], "pos": target["pos"],
            "ref": prov[target["gene"]]["ref"], "alt": prov[target["gene"]]["alt"],
            "source_variant_type": target["kind"],
        }
        result = generate_variant_candidates(variant, client, hla_panel, expected=target["expected"])
        windows = result["windows"]
        # Fail closed on the frozen window count and the frozen must-be-present peptide.
        if target["expected_windows"] is not None and len(windows) != target["expected_windows"]:
            raise ValueError(
                f"{target['gene']}: expected {target['expected_windows']} windows, got {len(windows)}"
            )
        if target["expected_peptide_present"] not in set(windows):
            raise ValueError(
                f"{target['gene']}: frozen peptide {target['expected_peptide_present']} not generated"
            )
        candidates = result["candidates"]
        # Attach input-only provenance features used by the router (no recognition data).
        candidates["expression_tpm"] = prov[target["gene"]]["expression_tpm"]
        candidates["n_callers"] = prov[target["gene"]]["n_callers"]
        candidates["n_timepoints"] = prov[target["gene"]]["n_timepoints"]
        candidates["tumor_vaf"] = prov[target["gene"]]["tumor_vaf"]
        frames.append(candidates)
        gen_provenance[target["gene"]] = result["provenance"]
        per_variant[target["gene"]] = {
            "n_windows": len(windows),
            "n_unique_peptides": len(set(windows)),
            "n_peptide_hla_pairs": int(len(candidates)),
            "positive_control": bool(target.get("positive_control", False)),
        }
    recovered = pd.concat(frames, ignore_index=True)
    return recovered, {"per_variant": per_variant, "ensembl": gen_provenance}


def score_recovered_prime(recovered: pd.DataFrame) -> pd.DataFrame:
    result = score_prime(recovered, peptide_col="mutant_peptide", hla_col="hla_allele")
    scored = result.scored.copy()
    scored["binding_percentile_rank"] = scored["mixmhcpred_rank"]
    scored["genuine_prime"] = -pd.to_numeric(scored["prime_rank"], errors="coerce")
    return scored


# ---------------------------------------------------------------------------
# Ranking (frozen router + route-aware top-k) and coverage metrics
# ---------------------------------------------------------------------------
def _rankable_pool(frame: pd.DataFrame) -> pd.DataFrame:
    routed = route_candidates(frame)
    return routed[routed["router_eligible"].astype(bool) & routed["rankable"].astype(bool)]


def _pure_prime_top_mutations(frame: pd.DataFrame, k: int) -> tuple[set, list]:
    pool = _rankable_pool(frame).copy()
    pool = pool[pd.to_numeric(pool["genuine_prime"], errors="coerce").notna()]
    pool["_tie_key"] = [
        hashlib.md5(f"{p}|{h}".encode()).hexdigest()
        for p, h in zip(pool["mutant_peptide"].astype(str), pool["hla_allele"].astype(str))
    ]
    ordered = pool.sort_values(["genuine_prime", "_tie_key"], ascending=[False, True], kind="mergesort")
    top = ordered.head(k)
    return set(top["mutation_id"]), top[["mutation_id", "gene_symbol", "mutant_peptide", "hla_allele",
                                         "prime_rank", "genuine_prime"]].to_dict("records")


def _route_aware_top_mutations(frame: pd.DataFrame) -> tuple[set, pd.DataFrame]:
    routed = route_candidates(frame)
    selected = select_route_aware_topk(routed, score_column="genuine_prime")
    chosen = selected[selected["route_selected"].astype(bool)].copy()
    return set(chosen["mutation_id"]), chosen


def _coverage(mutation_ids: set) -> dict:
    hit = sorted(HUDSON_EVAL_MUTATION_IDS & set(mutation_ids))
    return {"n": len(hit), "of": len(HUDSON_EVAL_MUTATION_IDS), "mutations": hit}


def _candidate_recall(frame: pd.DataFrame, require_rankable: bool) -> dict:
    if require_rankable:
        present = set(_rankable_pool(frame)["mutation_id"])
    else:
        present = set(frame["mutation_id"])
    return _coverage(present)


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------
def run(offline: bool) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    hla_panel = read_hla_panel(PVAC_PATH)
    if hla_panel != EXPECTED_HLA_PANEL:
        raise ValueError(f"pVAC HLA panel {hla_panel} != frozen expected {EXPECTED_HLA_PANEL}")

    prov = load_target_provenance()
    client = EnsemblClient(CACHE_DIR, offline=offline)
    recovered, gen_meta = generate_recovered(client, hla_panel, prov)
    recovered = score_recovered_prime(recovered)

    pvac = build_pvac_frame()
    keep = ["patient_id", "mutation_id", "gene_symbol", "chrom", "pos", "ref", "alt",
            "source_variant_type", "mhc_class", "mutant_peptide", "hla_allele", "expression_tpm",
            "binding_percentile_rank", "prime_rank", "prime_score", "mixmhcpred_rank",
            "genuine_prime", "candidate_source"]
    for col in ("n_callers", "n_timepoints", "tumor_vaf"):
        if col not in pvac:
            pvac[col] = pd.NA
    for col in keep + ["n_callers", "n_timepoints", "tumor_vaf"]:
        if col not in recovered:
            recovered[col] = pd.NA
    cols = keep + ["n_callers", "n_timepoints", "tumor_vaf"]
    union = union_candidates([pvac[cols], recovered[cols]])

    # ---- coverage: pVAC-only vs augmented (union), at four stages ----
    k = POLICY.k
    pvac_pure_muts, _ = _pure_prime_top_mutations(pvac, k)
    aug_pure_muts, aug_pure_rows = _pure_prime_top_mutations(union, k)
    pvac_ra_muts, _ = _route_aware_top_mutations(pvac)
    aug_ra_muts, aug_ra_chosen = _route_aware_top_mutations(union)

    per_variant = gen_meta["per_variant"]
    best_rank = {}
    for gene in per_variant:
        sub = recovered[recovered["gene_symbol"] == gene]
        ranks = pd.to_numeric(sub["prime_rank"], errors="coerce").dropna()
        best_rank[gene] = float(ranks.min()) if len(ranks) else None

    # Content hash (mode-invariant): recovered candidate identity + Ensembl SHAs + PRIME best ranks.
    recovered_identity = sorted(
        f"{r.mutation_id}|{r.mutant_peptide}|{r.hla_allele}|{r.prime_rank}"
        for r in recovered.itertuples()
    )
    ensembl_shas = {
        gene: {kind: rec["sha256"] for kind, rec in meta["ensembl"].items()}
        for gene, meta in gen_meta["ensembl"].items()
    }
    content_hash = hashlib.sha256(
        json.dumps({"identity": recovered_identity, "ensembl": ensembl_shas,
                    "best_rank": best_rank}, sort_keys=True).encode()
    ).hexdigest()

    result = {
        "policy_id": POLICY_ID,
        "router_policy_id": POLICY.policy_id,
        "mode": "offline" if offline else "online",
        "status": "post_hoc_reachability_diagnostic_not_independent_validation",
        "prime_commit": PRIME_COMMIT,
        "hla_panel": hla_panel,
        "generation": {
            "per_variant": {
                g: {**per_variant[g], "best_prime_rank": best_rank[g]} for g in per_variant
            },
            "n_recovered_pairs": int(len(recovered)),
            "n_recovered_scored": int(pd.to_numeric(recovered["prime_rank"], errors="coerce").notna().sum()),
        },
        "candidate_universe": {
            "n_pvac_candidates": int(len(pvac)),
            "n_union_candidates": int(len(union)),
        },
        "recall": {
            "candidate_generation": {
                "pvac_only": _candidate_recall(pvac, require_rankable=False),
                "augmented": _candidate_recall(union, require_rankable=False),
            },
            "rankable": {
                "pvac_only": _candidate_recall(pvac, require_rankable=True),
                "augmented": _candidate_recall(union, require_rankable=True),
            },
        },
        "top20_coverage": {
            "pure_prime": {"pvac_only": _coverage(pvac_pure_muts), "augmented": _coverage(aug_pure_muts)},
            "route_aware": {"pvac_only": _coverage(pvac_ra_muts), "augmented": _coverage(aug_ra_muts)},
        },
        "content_hash": content_hash,
        "ensembl_provenance": ensembl_shas,
        "interpretation": (
            "Post-hoc reachability fix on the motivating patient. NOT a claim that Epicurus beats "
            "PRIME: the selection score IS genuine PRIME (genuine_prime = -PRIME %rank); better "
            "candidate GENERATION lets genuine PRIME score targets it previously never received. "
            "A benefit claim needs a future untouched cohort."
        ),
    }

    _write_artifacts(union, aug_ra_chosen, aug_pure_rows, recovered, result, gen_meta, prov, offline)
    return result


def _write_artifacts(union, aug_ra_chosen, aug_pure_rows, recovered, result, gen_meta, prov, offline):
    OUT.mkdir(parents=True, exist_ok=True)

    union.to_csv(OUT / "RECOVERED_CANDIDATES.csv", index=False)

    top = aug_ra_chosen.copy()
    top = top.sort_values("route_rank", kind="mergesort") if "route_rank" in top else top
    hudson = HUDSON_EVAL_MUTATION_IDS  # evaluation-only join, AFTER ranking
    top_cols = ["route_rank", "route_selection_kind", "primary_route", "mutation_id", "gene_symbol",
                "mutant_peptide", "hla_allele", "prime_rank", "genuine_prime", "candidate_source"]
    top_cols = [c for c in top_cols if c in top.columns]
    top_out = top[top_cols].copy()
    top_out["hudson_eval_target"] = top_out["mutation_id"].isin(hudson)  # evaluation only
    top_out.to_csv(OUT / "TOP20.csv", index=False)

    (OUT / "RESULT.json").write_text(json.dumps(result, indent=2) + "\n")

    provenance = {
        "policy_id": POLICY_ID,
        "router_policy_id": DEFAULT_ROUTER_POLICY.policy_id,
        "mode": result["mode"],
        "content_hash": result["content_hash"],
        "prime_commit": PRIME_COMMIT,
        "input_provenance": prov,
        "ensembl": gen_meta["ensembl"],
        "cache_dir": str(CACHE_DIR.relative_to(ROOT)),
        "note": "Full third-party Ensembl sequences live only in the gitignored cache; this file "
                "carries URLs + SHA-256 + short junction context so an offline rerun reproduces the run.",
    }
    (OUT / "PROVENANCE.json").write_text(json.dumps(provenance, indent=2, default=str) + "\n")

    _write_report(result)


def _write_report(result: dict) -> None:
    r = result["recall"]
    t = result["top20_coverage"]
    gen = result["generation"]["per_variant"]

    def cov(entry):
        return f"{entry['n']}/{entry['of']}" + (f" ({', '.join(entry['mutations'])})" if entry["mutations"] else "")

    lines = [
        "# Lossless peptide recovery — osteosarc.com / Sid (EXPLORATORY, post-hoc)",
        "",
        f"> Generator policy `{result['policy_id']}` composed with router policy "
        f"`{result['router_policy_id']}`. Mode: **{result['mode']}**. Genuine PRIME commit "
        f"`{result['prime_commit'][:10]}`.",
        "",
        "> **Status:** post-hoc reachability diagnostic on the patient that motivated it — NOT "
        "preregistered / blind / independent. This does **not** show Epicurus beats PRIME: the "
        "selection score IS genuine PRIME (`genuine_prime = -PRIME %rank`); better candidate "
        "GENERATION lets genuine PRIME score targets it previously never received.",
        "",
        "## Generation (input-only; no assay/vaccine/label input)",
        "",
        "| Variant | Windows | Unique peptides | Peptide×HLA pairs | Best genuine-PRIME %rank | Role |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for gene, g in gen.items():
        role = "positive control" if g["positive_control"] else "recovered"
        best = f"{g['best_prime_rank']:.4f}" if g["best_prime_rank"] is not None else "n/a"
        lines.append(
            f"| {gene} | {g['n_windows']} | {g['n_unique_peptides']} | {g['n_peptide_hla_pairs']} | {best} | {role} |"
        )
    lines += [
        "",
        f"- pVAC candidate rows: **{result['candidate_universe']['n_pvac_candidates']}**; "
        f"union rows: **{result['candidate_universe']['n_union_candidates']}**; "
        f"recovered pairs scored by genuine PRIME: "
        f"**{result['generation']['n_recovered_scored']}/{result['generation']['n_recovered_pairs']}**.",
        "",
        "## Coverage of the 3 Hudson-expanded targets (labels joined AFTER ranking, evaluation only)",
        "",
        "| Stage | pVAC-only | Augmented (pVAC + lossless recovery) |",
        "|---|---|---|",
        f"| candidate generation recall | {cov(r['candidate_generation']['pvac_only'])} | {cov(r['candidate_generation']['augmented'])} |",
        f"| rankable recall (peptide+HLA) | {cov(r['rankable']['pvac_only'])} | {cov(r['rankable']['augmented'])} |",
        f"| pure genuine-PRIME top-20 | {cov(t['pure_prime']['pvac_only'])} | {cov(t['pure_prime']['augmented'])} |",
        f"| route-aware top-20 | {cov(t['route_aware']['pvac_only'])} | {cov(t['route_aware']['augmented'])} |",
        "",
        f"- Content hash (mode-invariant): `{result['content_hash']}`.",
        "",
        "## Interpretation guardrail",
        "",
        result["interpretation"],
    ]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--offline", action="store_true",
                        help="Serve Ensembl responses from the gitignored cache; fail closed on a miss.")
    args = parser.parse_args()
    result = run(offline=args.offline)
    print(json.dumps({
        "mode": result["mode"],
        "content_hash": result["content_hash"],
        "generation": {g: {"windows": v["n_windows"], "best_prime_rank": v["best_prime_rank"]}
                       for g, v in result["generation"]["per_variant"].items()},
        "recall": result["recall"],
        "top20_coverage": result["top20_coverage"],
    }, indent=2))


if __name__ == "__main__":
    main()
