"""Miller Hu_287 — candidate-universe assembly + four-arm scoring (RECONSTRUCTION_METHOD_PREREG.md).

Builds the label-blind class-I candidate universe from the reconstructed inputs, scores genuine PRIME +
frozen Epicurus, runs the four-arm harness, and ONLY THEN joins the sealed recognition labels to compute
the two frozen endpoints:
  (a) HLA-agnostic mutation reachability (called+filtered somatic variants ∩ recognized mutations), and
  (b) the class-I top-20 mechanistic four-arm hits@20.
Label ISOLATION: the label CSV is opened solely inside `score_against_labels()`; every upstream function is
label-free. Join key is position-based (chrom-normalized, pos, ref, alt) — annotation-independent.

    PYTHONPATH=src python -m scripts.miller_hu287_universe
"""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HU = ROOT / "data/raw/miller_ipv/hu_287"
PASS_VCF = HU / "somatic/Hu_287.somatic.pass.vcf.gz"
QUANT = HU / "salmon_quant/quant.sf"
HLA_JSON = ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction/HLA_PROVENANCE.json"
LABELS = ROOT / "data/raw/miller_ipv/miller_recognition_labels.csv"          # SEALED until score_against_labels
ENS_CACHE = HU / "ensembl_cache"
ART = ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction"

# frozen §3 base-variant filters
MIN_TUMOR_VAF, MAX_NORMAL_VAF, MIN_DEPTH = 0.05, 0.05, 10


def norm_chrom(c: str) -> str:
    return str(c).replace("chr", "").replace("Chr", "").upper()


def variant_key(chrom, pos, ref, alt) -> str:
    return f"{norm_chrom(chrom)}:{int(pos)}:{str(ref).upper()}:{str(alt).upper()}"


def passes_base_filters(tvaf: float, nvaf: float, tdp: int, ndp: int) -> bool:
    """Frozen §3 base-variant filters (identical for all four arms)."""
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


def load_filtered_variants(vcf_path: Path) -> pd.DataFrame:
    """PASS somatic variants passing the frozen §3 base filters (label-free)."""
    import pysam
    rows = []
    vf = pysam.VariantFile(str(vcf_path))
    for rec in vf:
        if rec.alts is None:
            continue
        alt = rec.alts[0]
        try:
            tvaf, tdp = _vaf_depth(rec.samples["Hu_287_T"])
            nvaf, ndp = _vaf_depth(rec.samples["Hu_287_N"])
        except KeyError:
            continue
        keep = passes_base_filters(tvaf, nvaf, tdp, ndp)
        r, a = rec.ref, alt
        if len(r) == 1 and len(a) == 1:
            vtype = "missense"                          # coding SNV (CDS-restricted calling); VEP confirms
        elif (len(a) - len(r)) % 3 != 0:
            vtype = "frameshift"
        else:
            vtype = "inframe"
        rows.append({"key": variant_key(rec.chrom, rec.pos, r, a), "chrom": norm_chrom(rec.chrom),
                     "pos": int(rec.pos), "ref": r, "alt": a, "source_variant_type": vtype,
                     "tumor_vaf": round(tvaf, 4), "normal_vaf": round(nvaf, 4), "tumor_dp": tdp, "normal_dp": ndp,
                     "pass_filters": bool(keep)})
    return pd.DataFrame(rows)


def gene_tpm_by_ensg(quant_path: Path) -> dict:
    d = pd.read_csv(quant_path, sep="\t")
    ensg = d["Name"].astype(str).str.split("|").str[1].str.split(".").str[0]   # ENSG (strip version)
    return d.assign(ensg=ensg).groupby("ensg")["TPM"].sum().to_dict()


def build_universe(variants: pd.DataFrame, hla_panel: list[str], tpm_by_ensg: dict) -> pd.DataFrame:
    """Label-blind lossless class-I universe for the filtered variants (Ensembl REST VEP enumeration)."""
    from event_b.lossless_peptide_generation import EnsemblClient, generate_variant_candidates
    client = EnsemblClient(ENS_CACHE)
    frames, notes = [], []
    for v in variants[variants["pass_filters"]].itertuples():
        variant = {"chrom": v.chrom, "pos": v.pos, "ref": v.ref, "alt": v.alt, "gene": "",
                   "source_variant_type": v.source_variant_type}
        try:
            res = generate_variant_candidates(variant, client, hla_panel, expected=None,
                                              require_mane_refseq=False)
        except Exception as e:                          # NOT_ENUMERABLE (phase/transcript) — recorded, not faked
            notes.append({"key": v.key, "status": "NOT_ENUMERABLE", "reason": str(e)[:200]})
            continue
        cand = res["candidates"].copy()
        cand["patient_id"] = "Hu_287"
        cand["mutation_id"] = v.key                     # position-based key (join-safe, annotation-independent)
        ensg = str(res["provenance"].get("gene_id", "")).split(".")[0]
        cand["expr"] = float(tpm_by_ensg.get(ensg, 0.0))
        cand["expression_tpm"] = cand["expr"]
        frames.append(cand)
    uni = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    return uni, notes


def score_universe(uni: pd.DataFrame) -> pd.DataFrame:
    """Attach genuine PRIME (%rank + MixMHCpred) and frozen Epicurus v0.1. No labels."""
    from event_b.prime_adapter import score_prime
    from event_b.prime_transfer import score_with_frozen
    res = score_prime(uni[["mutant_peptide", "hla_allele"]].rename(columns={"mutant_peptide": "peptide"}))
    sc = res.scored.rename(columns={"peptide": "mutant_peptide"})
    uni = uni.merge(sc[["mutant_peptide", "hla_allele", "prime_rank", "mixmhcpred_rank"]],
                    on=["mutant_peptide", "hla_allele"], how="left")
    uni["genuine_prime"] = -pd.to_numeric(uni["prime_rank"], errors="coerce")     # higher=better for top-k
    ep = uni.rename(columns={"prime_rank": "prime", "mixmhcpred_rank": "el"})
    ep["patient_id"] = "Hu_287"
    uni["epicurus"] = score_with_frozen(ep[["patient_id", "prime", "el", "expr"]])
    return uni


def score_against_labels(uni: pd.DataFrame, variants: pd.DataFrame) -> dict:
    """FINAL step — the ONLY place the sealed labels are read. Computes (a) HLA-agnostic reachability and
    (b) the class-I four-arm hits@20."""
    from benchmark.four_arm import run_patient, stage_attribution
    from benchmark.miller_ingest import mutation_recognition, parse_sra_runinfo  # noqa: F401
    labels = pd.read_csv(LABELS)
    lab = labels[labels["patient_id"] == "Hu_287"].copy()
    lab["key"] = [variant_key(c, p, r, a) for c, p, r, a in zip(lab["chrom"], lab["pos"], lab["ref"], lab["alt"])]
    recognized = sorted(lab.loc[lab["label"] == "POSITIVE", "key"].unique())
    tested = sorted(lab["key"].unique())

    called = set(variants["key"])
    called_pass = set(variants.loc[variants["pass_filters"], "key"])
    reach = {
        "n_recognized_mutations": len(recognized),
        "recognized_called_any": sorted(set(recognized) & called),
        "recognized_called_and_passed": sorted(set(recognized) & called_pass),
        "reachability_called": len(set(recognized) & called),
        "reachability_passed": len(set(recognized) & called_pass),
        "n_tested_mutations": len(tested),
    }
    positives = set(recognized)
    four = run_patient(uni, positives, k=20) if len(uni) else {"available": [], "arms": {}}
    arms = {a: {"evaluable": r.evaluable, "missing": r.missing,
                "hits_at_20": r.hits_at_k, "n_selected": r.n_selected,
                "generation_recall": (r.generation_recall.n if r.generation_recall else None)}
            for a, r in four.get("arms", {}).items()} if four.get("arms") else {}
    attrib = stage_attribution(four["arms"]) if four.get("arms") else {"evaluable": False}
    return {"reachability_hla_agnostic": reach, "class_i_four_arm": arms, "attribution": attrib,
            "recognized_keys": recognized}


def main() -> int:
    for f in (PASS_VCF, QUANT, HLA_JSON):
        if not f.exists():
            print(f"NOT_EVALUABLE: missing {f}")
            return 1
    hla_panel = json.loads(HLA_JSON.read_text())["class_i_alleles"]
    variants = load_filtered_variants(PASS_VCF)
    tpm = gene_tpm_by_ensg(QUANT)
    print(f"variants: {len(variants)} PASS, {int(variants['pass_filters'].sum())} after §3 filters; "
          f"HLA={hla_panel}")
    uni, notes = build_universe(variants, hla_panel, tpm)
    print(f"universe rows: {len(uni)}; NOT_ENUMERABLE: {len(notes)}")
    if len(uni):
        uni = score_universe(uni)
    result = score_against_labels(uni, variants)          # <-- labels read here only
    ART.mkdir(parents=True, exist_ok=True)
    (ART / "FOUR_ARM_RESULT.json").write_text(json.dumps(
        {"patient_id": "Hu_287", "hla": hla_panel, "n_variants_pass": int(variants["pass_filters"].sum()),
         "n_universe_rows": int(len(uni)), "not_enumerable": notes, **result}, indent=2, default=str) + "\n")
    if len(uni):
        uni.to_csv(HU / "universe_scored.csv.gz", index=False, compression="gzip")
    r = result["reachability_hla_agnostic"]
    print(f"(a) HLA-agnostic reachability: {r['reachability_passed']}/{r['n_recognized_mutations']} recognized "
          f"mutations called+passed")
    print(f"(b) class-I four-arm: {json.dumps(result['class_i_four_arm'])}")
    print("wrote", ART / "FOUR_ARM_RESULT.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
