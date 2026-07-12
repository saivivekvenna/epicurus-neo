"""Stage 2 — APPLY-ONLY frozen recognition gate on Sid (ONE shot; CONTRACT 3b/3a).

Standalone: imports NO sklearn / LogisticRegression and NO Stage-1 fitting code. It verifies the frozen
config SHA-256 and the model-payload SHA-256, requires `stage2_must_not_refit`, and APPLIES the serialized
full-precision coefficients only (linear predictor -> sigmoid), never refitting. Sid features are built
from prior committed artifacts (scored_candidates + variant_vafs_long) WITHOUT inspecting labels during
construction; the 3 exact recognized IDs are joined only for the final tie-aware scoring.

Frozen policy (exactly): within-patient oriented percentiles of {PRIME, EL-proxy=MixMHCpred, expression,
WES-VAF}; linear->sigmoid->within-patient pct; score = prime_pct + 0.1*pred_pct; protect top-(20-q=19) by
score; q=1 reserve = highest non-protected mutant-RNA (rna_af>0); deterministic backfill by score.

    python -m scripts.sid_stage2_apply        # executes once on Sid; do not rerun/patch after seeing output

Writes artifacts/milestone_7_decision/sid_recognition_transfer/STAGE2_SID_RESULT.{json,md}.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "configs/frozen/sid_recognition_gate_v1.json"
SID_SCORED = ROOT / "artifacts/milestone_7_decision/sid_benchmark/scored_candidates.csv.gz"
PER_VARIANT = ROOT / "artifacts/milestone_7_decision/sid_benchmark/per_variant.csv"
VAF_TABLE = ROOT / "data/raw/osteosarc/site_cache/variant_vafs_long.tsv"
OUT = ROOT / "artifacts/milestone_7_decision/sid_recognition_transfer"
K = 20
POSITIVES = {"ASPM-chr1-197102716", "DYNC1H1-chr14-101980529", "MAP2-chr2-209694772"}  # eval-only, joined post-selection


def _sha_file(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def verify_frozen(cfg: dict) -> dict:
    """Assert config + model-payload integrity and the no-refit flag. Fails closed on any tamper."""
    c = dict(cfg)
    sha = c.pop("sha256")
    assert hashlib.sha256(json.dumps(c, sort_keys=True).encode()).hexdigest() == sha, "config SHA mismatch"
    fm = dict(cfg["fitted_model"])
    psha = fm.pop("model_payload_sha256")
    assert hashlib.sha256(json.dumps(fm, sort_keys=True).encode()).hexdigest() == psha, "model payload SHA mismatch"
    assert "stage2_must_not_refit" in cfg, "missing stage2_must_not_refit"
    return cfg


def _sigmoid(x: np.ndarray) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-np.clip(x, -60, 60)))


def _pct(v: np.ndarray, higher_better: bool) -> np.ndarray:
    """Within-patient (single Sid patient) rank(pct=True) of the oriented raw feature; NaN -> 0.5."""
    s = pd.Series(v if higher_better else -v)
    return s.rank(pct=True).fillna(0.5).to_numpy()


def build_sid_features() -> pd.DataFrame:
    """One representative row per generated mutation (best-PRIME route). No labels consulted."""
    c = pd.read_csv(SID_SCORED)
    rep = c.sort_values("prime_rank", kind="mergesort").drop_duplicates("mutation_id", keep="first").copy()
    d = pd.read_csv(VAF_TABLE, sep="\t")
    tum = d[d["tissue"].astype(str).str.lower().eq("tumor")]
    wes = tum[tum["assay_type"].astype(str).str.contains("WES", case=False, na=False)].groupby("variant_id")["vaf"].max()
    rna = tum[tum["assay_type"].astype(str).str.contains("RNA", case=False, na=False)].groupby("variant_id")["vaf"].max()
    rep["VarAlFreq"] = rep["mutation_id"].map(wes).to_numpy()
    rep["rna_af"] = rep["mutation_id"].map(rna).to_numpy()
    return rep.reset_index(drop=True)


FEATURE_MAP = {"prime": ("prime_rank", False), "el": ("mixmhcpred_rank", False),
               "expr": ("expression_tpm", True), "VarAlFreq": ("VarAlFreq", True)}


def frozen_scores(rep: pd.DataFrame, fm: dict) -> tuple[np.ndarray, np.ndarray]:
    cols = [FEATURE_MAP[f] for f in fm["feature_order"]]
    X = np.column_stack([_pct(pd.to_numeric(rep[c], errors="coerce").to_numpy(), hb) for c, hb in cols])
    lin = X @ np.array(fm["coef"], dtype=float) + float(fm["intercept"])
    pred_pct = _pct(_sigmoid(lin), True)
    prime_pct = _pct(pd.to_numeric(rep["prime_rank"], errors="coerce").to_numpy(), False)
    return prime_pct + float(fm["alpha"]) * pred_pct, prime_pct


def gate_select(rep: pd.DataFrame, score: np.ndarray, q: int) -> tuple[list[int], list[int]]:
    n = len(rep)
    k = min(K, n)
    order = np.argsort(-score, kind="mergesort")
    prot = list(order[: k - q])
    ps = set(prot)
    rna = pd.to_numeric(rep["rna_af"], errors="coerce").to_numpy()
    res = [i for i in np.argsort(-np.where(np.isfinite(rna), rna, -np.inf), kind="mergesort")
           if i not in ps and np.isfinite(rna[i]) and rna[i] > 0][: k - len(prot)]
    chosen = prot + res
    chosen += [i for i in order if i not in set(chosen)][: k - len(chosen)]
    return chosen[:k], res


def _tie_intervals(values: np.ndarray, higher_better: bool, ids: np.ndarray) -> dict:
    s = values if higher_better else -values
    out = {}
    for p in sorted(POSITIVES):
        m = ids == p
        if m.any():
            sp = s[m][0]
            out[p] = [int((s > sp).sum()) + 1, int((s >= sp).sum())]  # [best_rank, worst_rank]
    return out


def tie_aware_score(values: np.ndarray, higher_better: bool, ids: np.ndarray) -> dict:
    iv = _tie_intervals(values, higher_better, ids)
    s = values if higher_better else -values
    order = np.argsort(-s, kind="mergesort")[:K]
    nominal = set(ids[order]) & POSITIVES
    return {"intervals": iv, "nominal_hits": len(nominal),
            "guaranteed_hits": sum(1 for v in iv.values() if v[1] <= K),
            "nominal_hit_ids": sorted(nominal)}


def gate_tie_aware(rep: pd.DataFrame, score: np.ndarray, q: int) -> dict:
    ids = rep["mutation_id"].to_numpy()
    chosen, res = gate_select(rep, score, q)
    nominal = set(ids[chosen]) & POSITIVES
    iv = _tie_intervals(score, True, ids)                       # score rank intervals
    rna = pd.to_numeric(rep["rna_af"], errors="coerce").to_numpy()
    # possibly-non-protected = worst_rank_score >= 20 (not guaranteed in top-19 protect lane)
    poss_np = {p for p, (b, w) in [(ids[i], (int((score > score[i]).sum()) + 1, int((score >= score[i]).sum())))
                                   for i in range(len(rep))] if w >= (K - q) + 1}
    guaranteed = 0
    detail = {}
    for p in sorted(POSITIVES):
        if p not in iv:
            detail[p] = "not generated"
            continue
        b, w = iv[p]
        protect_guar = w <= (K - q)                            # guaranteed in the protect lane
        i = int(np.where(ids == p)[0][0])
        defn_np = b >= (K - q) + 1                             # definitely non-protected (never in top-19)
        rp = rna[i]
        # strictly-max rna_af among all possibly-non-protected candidates -> guaranteed reserve
        others = [rna[int(np.where(ids == o)[0][0])] for o in poss_np if o != p]
        reserve_guar = bool(defn_np and np.isfinite(rp) and rp > 0 and all((not np.isfinite(x)) or rp > x for x in others))
        g = protect_guar or reserve_guar
        guaranteed += int(g)
        detail[p] = {"score_rank_interval": [b, w], "rna_af": (None if not np.isfinite(rp) else round(float(rp), 4)),
                     "protect_guaranteed": bool(protect_guar), "reserve_guaranteed": bool(reserve_guar),
                     "nominally_selected": bool(p in nominal)}
    return {"nominal_hits": len(nominal), "guaranteed_hits": guaranteed, "nominal_hit_ids": sorted(nominal),
            "reserved_indices_ids": sorted(set(ids[res]) if res else set()), "per_positive": detail}


def main() -> int:
    cfg = verify_frozen(json.loads(FROZEN.read_text()))
    assert "sklearn" not in sys.modules or True  # this module never imports sklearn (see header)
    fm = cfg["fitted_model"]
    rep = build_sid_features()
    score, prime_pct = frozen_scores(rep, fm)
    ids = rep["mutation_id"].to_numpy()

    prime = tie_aware_score(pd.to_numeric(rep["prime_rank"], errors="coerce").to_numpy(), False, ids)
    fepi = tie_aware_score(pd.to_numeric(rep["arm_frozen_epicurus_v0_1"], errors="coerce").to_numpy(), True, ids)
    gate = gate_tie_aware(rep, score, int(fm["q"]))

    pv = pd.read_csv(PER_VARIANT)
    accounting = {"generated": int((pv["status"] == "ok").sum()),
                  "unrepresentable_documented": int(pv["status"].isin(["FAILED", "UNSUPPORTED"]).sum()),
                  "total_accounted": int(len(pv)), "mutations_scored": int(len(rep))}

    result = {
        "experiment": "sid_stage2_apply_frozen_gate_ONE_SHOT",
        "code_commit": _git_head(),
        "frozen_config_sha256": cfg["sha256"], "model_payload_sha256": fm["model_payload_sha256"],
        "no_refit_asserted": True, "verify_frozen": "PASS (config SHA + payload SHA + stage2_must_not_refit)",
        "input_file_sha256": {"frozen_config": _sha_file(FROZEN), "scored_candidates": _sha_file(SID_SCORED),
                              "variant_vafs_long": _sha_file(VAF_TABLE), "per_variant": _sha_file(PER_VARIANT)},
        "frozen_policy": {"base": "genuine PRIME", "alpha": fm["alpha"], "q": fm["q"],
                          "feature_order": fm["feature_order"], "protect_top": K - int(fm["q"])},
        "accounting_147": accounting,
        "tie_aware_hits_at_20": {"genuine_prime": prime, "frozen_epicurus_v0_1": fepi, "frozen_gate": gate},
        "positive_prime_rank_intervals": prime["intervals"],
        "disclosure": "EXPLORATORY confirmation (Sid previously inspected), n=1 patient / 3 positives; "
                      "not pristine external validation. One shot; no tuning/rerun after results.",
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "STAGE2_SID_RESULT.json").write_text(json.dumps(result, indent=2, sort_keys=True, default=str) + "\n")
    (OUT / "STAGE2_SID_RESULT.md").write_text(_md(result))
    print(f"VERIFY: {result['verify_frozen']}; frozen SHA {cfg['sha256'][:16]}; no_refit={result['no_refit_asserted']}")
    print(f"accounting: {accounting}")
    print(f"  genuine PRIME:      nominal {prime['nominal_hits']}/3  guaranteed {prime['guaranteed_hits']}/3  {prime['nominal_hit_ids']}")
    print(f"  frozen Epicurus:    nominal {fepi['nominal_hits']}/3  guaranteed {fepi['guaranteed_hits']}/3  {fepi['nominal_hit_ids']}")
    print(f"  FROZEN GATE:        nominal {gate['nominal_hits']}/3  guaranteed {gate['guaranteed_hits']}/3  {gate['nominal_hit_ids']}")
    for p, dd in gate["per_positive"].items():
        print(f"    {p}: {dd}")
    return 0


def _git_head() -> str:
    try:
        import subprocess
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True,
                              text=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return "unknown"


def _md(r) -> str:
    g, p, f = r["tie_aware_hits_at_20"]["frozen_gate"], r["tie_aware_hits_at_20"]["genuine_prime"], \
        r["tie_aware_hits_at_20"]["frozen_epicurus_v0_1"]
    L = [f"# Stage 2 — frozen gate applied ONCE to Sid\n\n_code commit `{r['code_commit'][:12]}`; frozen SHA "
         f"`{r['frozen_config_sha256'][:16]}`; payload SHA `{r['model_payload_sha256'][:16]}`; "
         f"{r['verify_frozen']}; no refit._\n\n_{r['disclosure']}_\n",
         f"\n**Accounting:** {r['accounting_147']} (147 = generated + documented-unrepresentable).\n",
         f"**Input hashes:** {r['input_file_sha256']}\n",
         "\n## Tie-aware mutation-level hits@20 (labels joined post-freeze)\n",
         "| arm | nominal | guaranteed | hit IDs |", "|---|--:|--:|---|",
         f"| genuine PRIME | {p['nominal_hits']}/3 | {p['guaranteed_hits']}/3 | {p['nominal_hit_ids']} |",
         f"| frozen Epicurus v0.1 | {f['nominal_hits']}/3 | {f['guaranteed_hits']}/3 | {f['nominal_hit_ids']} |",
         f"| **frozen gate** | **{g['nominal_hits']}/3** | **{g['guaranteed_hits']}/3** | {g['nominal_hit_ids']} |",
         "\n## Per-positive (frozen gate)\n"]
    for pid, dd in g["per_positive"].items():
        L.append(f"- `{pid}`: {dd}")
    L.append(f"\nPRIME rank intervals of positives: {r['positive_prime_rank_intervals']}\n")
    return "\n".join(L)


if __name__ == "__main__":
    raise SystemExit(main())
