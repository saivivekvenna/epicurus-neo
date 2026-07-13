"""Miller Hu_287 — candidate-universe FREEZE (label-free) then UNSEAL (single label join).

Two commands enforce the pre-registration split (no threshold or arm can be chosen after seeing outcomes):
  freeze  — build the class-I universe from reconstructed inputs (PASS Mutect2 variants + OptiType HLA +
            salmon TPM), VEP-classify each variant's consequence (only protein-altering ones enumerated;
            synonymous/stop/splice/ambiguous => NOT_ENUMERABLE), score genuine PRIME + frozen Epicurus,
            compute each arm's ORDERED top-20 selection, and write universe + per-arm selections + SHA-256s
            + n_selected/saturation + a LOCK marker. Opens NO label file.
  unseal  — verify the frozen SHA-256s, then open the sealed labels ONCE to compute (a) HLA-agnostic
            mutation reachability and (b) the class-I four-arm hits@20 (unique mutations).

Epicurus v0.1 expects NetMHCpan-EL as `el`; MixMHCpred %rank is NOT a valid substitute, so `el` is set NaN
(frozen 0.5-percentile fallback) unless genuine NetMHCpan-EL is available — disclosed in provenance. pVAC
candidates come from a GENUINE pVACseq run if available; lossless-recovered rows are never relabeled pVAC.

    PYTHONPATH=src python -m scripts.miller_hu287_universe freeze
    PYTHONPATH=src python -m scripts.miller_hu287_universe unseal
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HU = ROOT / "data/raw/miller_ipv/hu_287"
PASS_VCF = HU / "somatic/Hu_287.somatic.pass.vcf.gz"
QUANT = HU / "salmon_quant/quant.sf"
HLA_JSON = ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction/HLA_PROVENANCE.json"
PVAC_CANDIDATES = HU / "pvac/pvac_candidates.csv"       # genuine pVACseq output if the env ran; else absent
LABELS = ROOT / "data/raw/miller_ipv/miller_recognition_labels.csv"          # SEALED until unseal()
ENS_CACHE = HU / "ensembl_cache"
ART = ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction"
FREEZE_DIR = HU / "freeze"

MIN_TUMOR_VAF, MAX_NORMAL_VAF, MIN_DEPTH = 0.05, 0.05, 10
K = 20
_ENUMERABLE = {"missense_variant": "missense", "inframe_insertion": "inframe",
               "inframe_deletion": "inframe", "frameshift_variant": "frameshift"}


def norm_chrom(c) -> str:
    return str(c).replace("chr", "").replace("Chr", "").upper()


def variant_key(chrom, pos, ref, alt) -> str:
    return f"{norm_chrom(chrom)}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def passes_base_filters(tvaf: float, nvaf: float, tdp: int, ndp: int) -> bool:
    return bool(tvaf >= MIN_TUMOR_VAF and nvaf <= MAX_NORMAL_VAF and tdp >= MIN_DEPTH and ndp >= MIN_DEPTH
                and (nvaf == 0 or tvaf / max(nvaf, 1e-9) >= 1.0))


def _vaf_depth(sample) -> tuple[float, int]:
    ad = sample.get("AD")
    if ad and len(ad) >= 2 and sum(ad) > 0:
        return ad[1] / sum(ad), int(sum(ad))
    dp = sample.get("DP") or 0
    af = sample.get("AF")
    af = af[0] if isinstance(af, (tuple, list)) else af
    return (float(af) if af is not None else 0.0), int(dp)


def sha256_file(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def load_filtered_variants(vcf_path: Path) -> pd.DataFrame:
    """PASS somatic variants + frozen §3 base filters (label-free). Consequence assigned later by VEP."""
    import pysam
    rows = []
    for rec in pysam.VariantFile(str(vcf_path)):
        if rec.alts is None:
            continue
        alt = rec.alts[0]
        try:
            tvaf, tdp = _vaf_depth(rec.samples["Hu_287_T"])
            nvaf, ndp = _vaf_depth(rec.samples["Hu_287_N"])
        except KeyError:
            continue
        rows.append({"key": variant_key(rec.chrom, rec.pos, rec.ref, alt), "chrom": norm_chrom(rec.chrom),
                     "pos": int(rec.pos), "ref": rec.ref, "alt": alt,
                     "tumor_vaf": round(tvaf, 4), "normal_vaf": round(nvaf, 4), "tumor_dp": tdp, "normal_dp": ndp,
                     "pass_filters": passes_base_filters(tvaf, nvaf, tdp, ndp)})
    return pd.DataFrame(rows)


def gene_tpm_by_ensg(quant_path: Path) -> dict:
    d = pd.read_csv(quant_path, sep="\t")
    ensg = d["Name"].astype(str).str.split("|").str[1].str.split(".").str[0]
    return d.assign(ensg=ensg).groupby("ensg")["TPM"].sum().to_dict()


def classify_consequence(client, chrom, pos, ref, alt) -> tuple[str | None, str]:
    """VEP most-severe consequence -> (enumerable_kind or None, raw_term). Only protein-altering variants
    are enumerable; synonymous/stop/splice/UTR/intron/ambiguous -> NOT_ENUMERABLE (honest)."""
    from event_b.lossless_peptide_generation import genomic_hgvs
    vep = client.vep_hgvs(genomic_hgvs(chrom, pos, ref, alt))
    payload = vep["json"]
    term = ""
    if isinstance(payload, list) and payload:
        term = str(payload[0].get("most_severe_consequence", ""))
    return _ENUMERABLE.get(term), term


def load_pvac_candidates() -> pd.DataFrame:
    """Genuine pVACseq class-I candidates if a real run produced them; else empty (pvac arm NOT_EVALUABLE).
    NEVER synthesized from lossless rows."""
    if not PVAC_CANDIDATES.exists():
        return pd.DataFrame()
    d = pd.read_csv(PVAC_CANDIDATES)
    d["candidate_source"] = "pvac"
    return d


def build_universe(variants: pd.DataFrame, hla_panel: list[str], tpm_by_ensg: dict):
    """Union genuine pVAC rows (first) with the label-blind lossless class-I universe."""
    from event_b.lossless_peptide_generation import EnsemblClient, generate_variant_candidates
    from event_b.lossless_peptide_generation import union_candidates
    client = EnsemblClient(ENS_CACHE)
    frames, notes = [], []
    for v in variants[variants["pass_filters"]].itertuples():
        kind, term = classify_consequence(client, v.chrom, v.pos, v.ref, v.alt)   # Fix 4: VEP first
        if kind is None:
            notes.append({"key": v.key, "status": "NOT_ENUMERABLE", "consequence": term})
            continue
        variant = {"chrom": v.chrom, "pos": v.pos, "ref": v.ref, "alt": v.alt, "gene": "",
                   "source_variant_type": kind}
        try:
            res = generate_variant_candidates(variant, client, hla_panel, expected=None,
                                              require_mane_refseq=False)
        except Exception as e:
            notes.append({"key": v.key, "status": "NOT_ENUMERABLE", "consequence": term, "reason": str(e)[:180]})
            continue
        cand = res["candidates"].copy()
        cand["patient_id"] = "Hu_287"
        cand["mutation_id"] = v.key
        cand["candidate_source"] = "lossless_recovery"
        ensg = str(res["provenance"].get("gene_id", "")).split(".")[0]
        cand["expr"] = float(tpm_by_ensg.get(ensg, 0.0))
        cand["expression_tpm"] = cand["expr"]
        cand["vep_consequence"] = term
        frames.append(cand)
    lossless = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    pvac = load_pvac_candidates()
    # genuine pVAC rows FIRST so the incumbent represents any shared (peptide,HLA) route; lossless = union
    uni = union_candidates([pvac, lossless]) if len(pvac) else lossless
    return uni, notes, bool(len(pvac))


def score_universe(uni: pd.DataFrame) -> pd.DataFrame:
    """genuine PRIME (%rank + MixMHCpred) + frozen Epicurus v0.1. Fix 2: el=NaN (NetMHCpan-EL unavailable;
    MixMHCpred is NOT NetMHCpan-EL), so Epicurus's el term takes the frozen 0.5-percentile fallback."""
    from event_b.prime_adapter import score_prime
    from event_b.prime_transfer import score_with_frozen
    res = score_prime(uni[["mutant_peptide", "hla_allele"]].rename(columns={"mutant_peptide": "peptide"}))
    sc = res.scored.rename(columns={"peptide": "mutant_peptide"})
    uni = uni.merge(sc[["mutant_peptide", "hla_allele", "prime_rank", "mixmhcpred_rank"]].drop_duplicates(),
                    on=["mutant_peptide", "hla_allele"], how="left")
    uni["genuine_prime"] = -pd.to_numeric(uni["prime_rank"], errors="coerce")
    uni["el"] = np.nan                                  # NetMHCpan-EL not available -> neutral (disclosed)
    ep = uni.rename(columns={"prime_rank": "prime"})[["mutant_peptide", "hla_allele", "prime", "el", "expr"]].copy()
    ep["patient_id"] = "Hu_287"
    uni["epicurus"] = score_with_frozen(ep[["patient_id", "prime", "el", "expr"]])
    return uni


def arm_selection(uni: pd.DataFrame, arm) -> pd.DataFrame:
    """One arm's ORDERED top-K selection, LABEL-FREE (positives never consulted)."""
    from benchmark.four_arm import (DEFAULT_ROUTER_POLICY, _generation_rows, _plain_topk, _route_aware_topk)
    gen = _generation_rows(uni, arm.generation)
    if len(gen) == 0:
        return gen.assign(_rank=[]) if "mutation_id" in gen.columns else pd.DataFrame(
            columns=["mutation_id", "mutant_peptide", "hla_allele"])
    sel = (_route_aware_topk(gen, arm.scorer, DEFAULT_ROUTER_POLICY) if arm.selection == "route_aware"
           else _plain_topk(gen, arm.scorer, K, DEFAULT_ROUTER_POLICY))
    cols = [c for c in ("mutation_id", "mutant_peptide", "hla_allele", arm.scorer, "candidate_source") if c in sel.columns]
    return sel[cols].reset_index(drop=True)


def freeze() -> dict:
    """Build + score + select, then LOCK — with NO label access."""
    from benchmark.four_arm import FOUR_ARMS, detect_available, evaluate_eligibility
    for f in (PASS_VCF, QUANT, HLA_JSON):
        if not f.exists():
            return {"status": "NOT_EVALUABLE", "missing": str(f)}
    hla_panel = json.loads(HLA_JSON.read_text())["class_i_alleles"]
    variants = load_filtered_variants(PASS_VCF)
    uni, notes, has_pvac = build_universe(variants, hla_panel, gene_tpm_by_ensg(QUANT))
    if len(uni):
        uni = score_universe(uni)
    FREEZE_DIR.mkdir(parents=True, exist_ok=True)
    # freeze the filtered variant keys (for HLA-agnostic reachability at unseal) + the scored universe
    variants.to_csv(FREEZE_DIR / "variants.csv", index=False)
    (uni if len(uni) else pd.DataFrame()).to_csv(FREEZE_DIR / "universe.csv", index=False)
    available = detect_available(uni, {"__any__"}) if len(uni) else set()
    elig = evaluate_eligibility(available)
    arms_meta, hashes = {}, {"variants.csv": sha256_file(FREEZE_DIR / "variants.csv"),
                             "universe.csv": sha256_file(FREEZE_DIR / "universe.csv")}
    for arm in FOUR_ARMS:
        sel = arm_selection(uni, arm) if len(uni) else pd.DataFrame(columns=["mutation_id"])
        fn = f"select_{arm.arm_id}.csv"
        sel.to_csv(FREEZE_DIR / fn, index=False)
        hashes[fn] = sha256_file(FREEZE_DIR / fn)
        arms_meta[arm.arm_id] = {"evaluable": bool(elig[arm.arm_id].evaluable),
                                 "missing": elig[arm.arm_id].missing,
                                 "n_selected": int(len(sel)), "saturated": bool(len(sel) >= K),
                                 "selection_file": fn, "n_unique_mutations": int(sel["mutation_id"].nunique()) if len(sel) else 0}
    manifest = {"patient_id": "Hu_287", "labels_opened": False, "LOCK": "FROZEN_NO_LABELS",
                "hla_panel": hla_panel, "genuine_pvac_lane": has_pvac,
                "n_variants_pass": int(variants["pass_filters"].sum()), "n_universe_rows": int(len(uni)),
                "not_enumerable": notes, "arms": arms_meta, "sha256": hashes,
                "el_feature": "NaN (NetMHCpan-EL unavailable; MixMHCpred is NOT a valid el substitute) -> frozen 0.5 fallback",
                "pvac_note": ("genuine pVACseq candidates unioned first" if has_pvac
                              else "no genuine pVAC run -> pvac_prime NOT_EVALUABLE (lossless NEVER relabeled pVAC)")}
    (FREEZE_DIR / "FREEZE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    ART.mkdir(parents=True, exist_ok=True)
    (ART / "FREEZE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    return manifest


def verify_frozen_hashes(man: dict, freeze_dir: Path = FREEZE_DIR) -> tuple[bool, str | None]:
    """Recompute every frozen file's SHA-256 and compare to the manifest. Returns (ok, first_bad_file)."""
    for fn, h in man["sha256"].items():
        if sha256_file(freeze_dir / fn) != h:
            return False, fn
    return True, None


def unseal() -> dict:
    """Verify frozen hashes, then open labels ONCE for the two endpoints."""
    man = json.loads((FREEZE_DIR / "FREEZE_MANIFEST.json").read_text())
    ok, bad = verify_frozen_hashes(man)                  # integrity gate BEFORE any label read
    if not ok:
        return {"status": "HASH_MISMATCH", "file": bad}
    labels = pd.read_csv(LABELS)                          # <-- the ONLY label read
    lab = labels[labels["patient_id"] == "Hu_287"].copy()
    lab["key"] = [variant_key(c, p, r, a) for c, p, r, a in zip(lab["chrom"], lab["pos"], lab["ref"], lab["alt"])]
    recognized = set(lab.loc[lab["label"] == "POSITIVE", "key"])
    variants = pd.read_csv(FREEZE_DIR / "variants.csv")
    called = set(variants["key"])
    called_pass = set(variants.loc[variants["pass_filters"], "key"])
    reach = {"n_recognized_mutations": len(recognized),
             "reachability_called": len(recognized & called),
             "reachability_called_and_passed": len(recognized & called_pass),
             "recognized_called_and_passed_keys": sorted(recognized & called_pass),
             "n_tested_mutations": lab["key"].nunique()}
    arms = {}
    for arm_id, meta in man["arms"].items():
        sel = pd.read_csv(FREEZE_DIR / meta["selection_file"])
        hit_muts = sorted(set(sel["mutation_id"]) & recognized) if len(sel) else []
        arms[arm_id] = {"evaluable": meta["evaluable"], "n_selected": meta["n_selected"],
                        "saturated": meta["saturated"], "hits_at_20_unique_mutations": len(hit_muts),
                        "hit_mutation_keys": hit_muts}
    out = {"patient_id": "Hu_287", "endpoint_a_hla_agnostic_reachability": reach,
           "endpoint_b_class_i_four_arm": arms,
           "limitation": "Miller IFN-g 20mers are NOT HLA-I restricted (CD4/class-II possible); endpoint (b) "
                         "is a CD8/class-I mechanistic view and undercounts full biological recall.",
           "el_disclosure": man["el_feature"], "pvac_disclosure": man["pvac_note"]}
    (ART / "FOUR_ARM_RESULT.json").write_text(json.dumps(out, indent=2, default=str) + "\n")
    return out


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    cmd = argv[0] if argv else "freeze"
    if cmd == "freeze":
        m = freeze()
        print("FREEZE:", m.get("LOCK", m.get("status")), "| universe rows", m.get("n_universe_rows"),
              "| pvac", m.get("genuine_pvac_lane"), "| arms",
              {a: (v["n_selected"], v["evaluable"]) for a, v in m.get("arms", {}).items()})
        return 0
    if cmd == "unseal":
        o = unseal()
        if o.get("status"):
            print("UNSEAL:", o["status"])
            return 1
        r = o["endpoint_a_hla_agnostic_reachability"]
        print(f"(a) reachability: {r['reachability_called_and_passed']}/{r['n_recognized_mutations']}")
        print("(b) class-I hits@20: " + json.dumps({a: v["hits_at_20_unique_mutations"] for a, v in o["endpoint_b_class_i_four_arm"].items()}))
        return 0
    print("usage: miller_hu287_universe.py [freeze|unseal]")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
