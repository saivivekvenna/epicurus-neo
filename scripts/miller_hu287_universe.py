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
import re
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
HU = ROOT / "data/raw/miller_ipv/hu_287"
PASS_VCF = HU / "somatic/Hu_287.somatic.pass.vcf.gz"
NORM_VCF = HU / "somatic/Hu_287.somatic.pass.norm.vcf.gz"
REF = ROOT / "data/raw/refs/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
QUANT = HU / "salmon_quant/quant.sf"
HLA_JSON = ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction/HLA_PROVENANCE.json"
PVAC_CANDIDATES = HU / "pvac/pvac_candidates.csv"       # genuine pVACseq output if the env ran; else absent
PVAC_PROVENANCE = HU / "pvac/PVAC_PROVENANCE.json"      # required proof of a genuine pVACseq run
_KEY_RE = re.compile(r"^[0-9XYMT]+:\d+:[ACGTN]+:[ACGTN]+$")
LABELS = ROOT / "data/raw/miller_ipv/miller_recognition_labels.csv"          # SEALED until unseal()
ENS_CACHE = HU / "ensembl_cache"
ART = ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction"
FREEZE_DIR = HU / "freeze"

MIN_TUMOR_VAF, MAX_NORMAL_VAF, MIN_DEPTH = 0.05, 0.05, 10
K = 20
CLASS_I_MIN, CLASS_I_MAX = 8, 11                        # prereg §3 class-I lengths (generator/PRIME allow 8-14)
_ENUMERABLE = {"missense_variant": "missense", "inframe_insertion": "inframe",
               "inframe_deletion": "inframe", "frameshift_variant": "frameshift"}
# router-compatible truthful source_variant_type (MISSENSE is 'typical'; inframe/frameshift are atypical)
_ROUTER_TYPE = {"missense": "MISSENSE", "inframe": "INFRAME", "frameshift": "FRAMESHIFT"}
RNA_BAM = HU / "rna/Hu_287_tumor_rna.sorted.bam"
FROZEN_EPICURUS = ROOT / "configs/frozen/epicurus_v0_1.json"
FROZEN_EXPR_POLICY = ROOT / "configs/frozen/expression_policy_v1.json"
FROZEN_ROUTER = ROOT / "configs/frozen/evidence_router_v1.json"


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


def normalize_hla(a) -> str:
    """Normalize an HLA allele to the 'HLA-A*02:01' form used by the panel/PRIME adapter."""
    s = str(a).strip().upper().replace("HLA-", "")
    return f"HLA-{s}" if s else ""


def load_pvac_candidates() -> pd.DataFrame:
    """Genuine pVACseq class-I candidates ONLY. A CSV alone is insufficient — it must be accompanied by a
    provenance file proving a real pvacseq/pvactools run (tool+version), carry the required schema, use
    position-based mutation_ids (chrom:pos:ref:alt), and normalized HLA. Anything else is REJECTED (empty,
    => pvac arm NOT_EVALUABLE). Lossless rows are NEVER relabeled pVAC."""
    if not (PVAC_CANDIDATES.exists() and PVAC_PROVENANCE.exists()):
        return pd.DataFrame()
    try:
        prov = json.loads(PVAC_PROVENANCE.read_text())
    except Exception:
        return pd.DataFrame()
    if str(prov.get("tool", "")).lower() not in {"pvacseq", "pvactools"} or not prov.get("version"):
        return pd.DataFrame()                            # no genuine provenance -> reject
    d = pd.read_csv(PVAC_CANDIDATES)
    if not {"mutation_id", "mutant_peptide", "hla_allele"}.issubset(d.columns):
        return pd.DataFrame()                            # wrong schema -> reject
    if not d["mutation_id"].astype(str).map(lambda x: bool(_KEY_RE.match(x))).all():
        return pd.DataFrame()                            # non-position-based mutation_id -> reject
    d = d.copy()
    d["hla_allele"] = d["hla_allele"].map(normalize_hla)
    d["candidate_source"] = "pvac"
    return d


def normalize_pass_vcf(pass_vcf: Path, ref: Path, out: Path) -> Path:
    """Left-align + split multiallelics against the EXACT reference (bcftools norm) so indel keys are
    canonical before enumeration/freeze. SNVs are unchanged. Fail-closed on error."""
    subprocess.run(["bcftools", "norm", "-f", str(ref), "-m-any", "-Oz", "-o", str(out), str(pass_vcf)],
                   check=True, capture_output=True)
    subprocess.run(["bcftools", "index", "-t", str(out)], check=True, capture_output=True)
    return out


def rna_alt_evidence(variants: pd.DataFrame, rna_bam: Path = RNA_BAM) -> tuple[dict, str]:
    """Label-blind allele-specific RNA evidence at each PASS site from the tumor-RNA BAM (evidence-only,
    NEVER a hard filter). For SNVs, count the query base at the site. For INDELS, exact allele matching from
    a pileup is unreliable, so RNA-alt is NOT_ASSESSED (never fabricated). `rna_mutant_reads` mirrors
    `rna_alt_obs` (the canonical field the evidence router reads). Returns ({key: {...}}, status)."""
    if not rna_bam.exists():
        return {}, "NOT_ASSESSED (tumor-RNA HISAT2 BAM not available)"
    import pysam
    out = {}
    with pysam.AlignmentFile(str(rna_bam)) as bam:
        for v in variants[variants["pass_filters"]].itertuples():
            ref, alt = str(v.ref).upper(), str(v.alt).upper()
            if not (len(ref) == 1 and len(alt) == 1):        # indel -> NOT_ASSESSED (do not fabricate)
                out[v.key] = {"rna_alt_obs": None, "rna_mutant_reads": None, "rna_depth": None,
                              "rna_vaf": None, "indel_rna": "NOT_ASSESSED"}
                continue
            depth = altobs = 0
            try:
                for col in bam.pileup(str(v.chrom), int(v.pos) - 1, int(v.pos), truncate=True, min_base_quality=0):
                    for pr in col.pileups:
                        if pr.is_refskip or pr.query_position is None:
                            continue
                        depth += 1
                        if pr.alignment.query_sequence[pr.query_position].upper() == alt:
                            altobs += 1
            except (ValueError, KeyError):
                continue
            out[v.key] = {"rna_alt_obs": altobs, "rna_mutant_reads": altobs, "rna_depth": depth,
                          "rna_vaf": round(altobs / depth, 4) if depth else 0.0}
    return out, "COMPUTED (SNV allele-specific; indels NOT_ASSESSED)"


def build_universe(variants: pd.DataFrame, hla_panel: list[str], tpm_by_ensg: dict, rna_by_key: dict):
    """Union genuine pVAC rows (first) with the label-blind lossless class-I universe. Attaches true
    gene_symbol, router-truthful source_variant_type, WES + RNA evidence per mutation."""
    from event_b.lossless_peptide_generation import EnsemblClient, generate_variant_candidates
    from event_b.lossless_peptide_generation import union_candidates
    client = EnsemblClient(ENS_CACHE)
    frames, notes = [], []
    for v in variants[variants["pass_filters"]].itertuples():
        kind, term = classify_consequence(client, v.chrom, v.pos, v.ref, v.alt)   # VEP first
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
        cand["gene_symbol"] = res["provenance"].get("gene_symbol", "")          # Fix: true gene (not blank)
        cand["source_variant_type"] = _ROUTER_TYPE[kind]                        # Fix: VEP kind, not "SNV"
        cand["vep_consequence"] = term
        ensg = str(res["provenance"].get("gene_id", "")).split(".")[0]
        cand["expr"] = float(tpm_by_ensg.get(ensg, 0.0))
        cand["expression_tpm"] = cand["expr"]
        for col in ("tumor_vaf", "normal_vaf", "tumor_dp", "normal_dp"):        # Fix: propagate WES evidence
            cand[col] = getattr(v, col)
        rna = rna_by_key.get(v.key, {})
        cand["rna_alt_obs"] = rna.get("rna_alt_obs", np.nan)
        cand["rna_mutant_reads"] = rna.get("rna_mutant_reads", np.nan)   # canonical field the router reads
        cand["rna_depth"] = rna.get("rna_depth", np.nan)
        cand["rna_vaf"] = rna.get("rna_vaf", np.nan)
        frames.append(cand)
    lossless = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    pvac = load_pvac_candidates()
    uni = union_candidates([pvac, lossless]) if len(pvac) else lossless
    return uni, notes, bool(len(pvac))


def filter_class_i_lengths(uni: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """Restrict BOTH pVAC and lossless rows to class-I lengths 8-11 (prereg §3) before scoring/selection.
    The generator/PRIME accept 8-14, so 12-14mers are explicitly dropped here. Returns (filtered, counts)."""
    if not len(uni):
        return uni, {"pre": 0, "post": 0, "dropped_len_12_14": 0}
    ln = uni["mutant_peptide"].astype(str).str.len()
    keep = (ln >= CLASS_I_MIN) & (ln <= CLASS_I_MAX)
    counts = {"pre": int(len(uni)), "post": int(keep.sum()),
              "dropped_len_12_14": int(((ln >= 12) & (ln <= 14)).sum()),
              "policy": f"class-I {CLASS_I_MIN}-{CLASS_I_MAX}mers only (generator/PRIME allow 8-14)"}
    return uni[keep].reset_index(drop=True), counts


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
    # Fix: give the router FACTUAL presentation evidence — MixMHCpred %rank IS a binding percentile rank
    # (label-blind, real predictor). This is NOT NetMHCpan-EL and is not used as Epicurus's `el`.
    uni["binding_percentile_rank"] = pd.to_numeric(uni["mixmhcpred_rank"], errors="coerce")
    uni["binding_rank_provenance"] = "MixMHCpred 3.0 %rank (via genuine PRIME); not NetMHCpan-EL"
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


def _git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def freeze() -> dict:
    """Build + score + select, then LOCK — with NO label access. ONE-SHOT/IMMUTABLE: refuses to overwrite
    an existing valid LOCK."""
    from benchmark.four_arm import DEFAULT_ROUTER_POLICY, FOUR_ARMS, detect_available, evaluate_eligibility
    existing = FREEZE_DIR / "FREEZE_MANIFEST.json"
    if existing.exists():                                # ONE-SHOT: never overwrite an existing manifest
        try:
            prev = json.loads(existing.read_text())
        except Exception:
            return {"status": "FROZEN_CORRUPT", "reason": "existing FREEZE_MANIFEST.json is unreadable; refusing to overwrite"}
        if prev.get("LOCK") != "FROZEN_NO_LABELS":
            return {"status": "FROZEN_INVALID", "reason": "existing manifest lacks a valid LOCK; refusing to overwrite"}
        ok, bad = verify_frozen_hashes(prev, FREEZE_DIR)   # verify derived hashes; fail closed on corruption
        return {"status": "ALREADY_FROZEN" if ok else "FROZEN_HASH_MISMATCH", "bad_file": bad,
                "sha256": prev.get("sha256"), "arms": prev.get("arms")}
    for f in (PASS_VCF, QUANT, HLA_JSON, REF):
        if not f.exists():
            return {"status": "NOT_EVALUABLE", "missing": str(f)}
    hla_panel = json.loads(HLA_JSON.read_text())["class_i_alleles"]
    norm_vcf = normalize_pass_vcf(PASS_VCF, REF, NORM_VCF)     # bcftools norm indels vs exact FASTA
    variants = load_filtered_variants(norm_vcf)
    rna_by_key, rna_status = rna_alt_evidence(variants)       # evidence-only (NOT_ASSESSED if BAM absent)
    uni, notes, has_pvac = build_universe(variants, hla_panel, gene_tpm_by_ensg(QUANT), rna_by_key)
    uni, length_counts = filter_class_i_lengths(uni)          # class-I 8-11 only, BEFORE scoring/selection
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
    def _rel(p: Path) -> str:
        try:
            return str(p.relative_to(ROOT))
        except ValueError:
            return str(p)
    fai = REF.parent / (REF.name + ".fai")
    input_sha = {_rel(p): (sha256_file(p) if p.exists() else "MISSING")
                 for p in (NORM_VCF, HLA_JSON, QUANT, FROZEN_EPICURUS, FROZEN_EXPR_POLICY, FROZEN_ROUTER,
                           fai, REF)}   # fai = contig names/lengths; REF = actual FASTA bases (3GB, hashed)
    manifest = {"patient_id": "Hu_287", "labels_opened": False, "LOCK": "FROZEN_NO_LABELS",
                "hla_panel": hla_panel, "genuine_pvac_lane": has_pvac,
                "indel_normalization": "bcftools norm -f GRCh38.fa -m-any (left-align+split) before enumeration",
                "class_i_length_filter": length_counts, "rna_alt_evidence_status": rna_status,
                "n_variants_pass": int(variants["pass_filters"].sum()), "n_universe_rows": int(len(uni)),
                "not_enumerable": notes, "arms": arms_meta, "sha256": hashes,
                "input_sha256": input_sha, "git_commit": _git_commit(),
                "router_policy_id": DEFAULT_ROUTER_POLICY.policy_id,
                "presentation_evidence": "router binding_percentile_rank = MixMHCpred %rank (real predictor; not NetMHCpan-EL)",
                "el_feature": "NaN (NetMHCpan-EL unavailable; MixMHCpred is NOT a valid el substitute) -> frozen 0.5 fallback",
                "pvac_note": ("genuine pVACseq candidates unioned first" if has_pvac
                              else "no genuine pVAC run -> pvac_prime NOT_EVALUABLE (lossless NEVER relabeled pVAC)")}
    (FREEZE_DIR / "FREEZE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    ART.mkdir(parents=True, exist_ok=True)
    (ART / "FREEZE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    return manifest


def verify_frozen_hashes(man: dict, freeze_dir: Path = FREEZE_DIR) -> tuple[bool, str | None]:
    """Recompute every frozen file's SHA-256 and compare. Fails CLOSED (returns (False, name)) on a missing
    or unreadable file or a malformed/missing sha entry — never raises."""
    sha = man.get("sha256")
    if not isinstance(sha, dict):
        return False, "<malformed sha256 manifest>"
    for fn, h in sha.items():
        p = freeze_dir / fn
        if not p.exists() or not isinstance(h, str):
            return False, fn
        try:
            if sha256_file(p) != h:
                return False, fn
        except OSError:
            return False, fn
    return True, None


def unseal() -> dict:
    """Mechanically ONCE-ONLY: refuse if already unsealed. Verify derived + INPUT hashes, THEN open labels
    once for the two endpoints."""
    marker = FREEZE_DIR / "UNSEALED.json"
    if marker.exists():                                  # once-only: do not re-open LABELS
        return {"status": "ALREADY_UNSEALED", **json.loads(marker.read_text())}
    man = json.loads((FREEZE_DIR / "FREEZE_MANIFEST.json").read_text())
    ok, bad = verify_frozen_hashes(man, FREEZE_DIR)      # derived-file integrity gate BEFORE any label read
    if not ok:
        return {"status": "HASH_MISMATCH", "file": bad}
    for rel, h in man.get("input_sha256", {}).items():   # input integrity gate BEFORE any label read
        cur = sha256_file(ROOT / rel) if (ROOT / rel).exists() else "MISSING"
        if cur != h:
            return {"status": "INPUT_HASH_MISMATCH", "file": rel, "expected": h, "got": cur}
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
        if not meta["evaluable"]:                        # NOT evaluated -> null, never 0 (0 implies evaluated)
            arms[arm_id] = {"evaluable": False, "missing": meta.get("missing"), "n_selected": None,
                            "saturated": None, "hits_at_20_unique_mutations": None, "hit_mutation_keys": None}
            continue
        sel = pd.read_csv(FREEZE_DIR / meta["selection_file"])
        hit_muts = sorted(set(sel["mutation_id"]) & recognized) if len(sel) else []
        arms[arm_id] = {"evaluable": True, "n_selected": meta["n_selected"],
                        "saturated": meta["saturated"], "hits_at_20_unique_mutations": len(hit_muts),
                        "hit_mutation_keys": hit_muts}
    out = {"patient_id": "Hu_287", "endpoint_a_hla_agnostic_reachability": reach,
           "endpoint_b_class_i_four_arm": arms,
           "limitation": "Miller IFN-g 20mers are NOT HLA-I restricted (CD4/class-II possible); endpoint (b) "
                         "is a CD8/class-I mechanistic view and undercounts full biological recall.",
           "el_disclosure": man.get("el_feature"), "pvac_disclosure": man.get("pvac_note")}
    (ART / "FOUR_ARM_RESULT.json").write_text(json.dumps(out, indent=2, default=str) + "\n")
    marker.write_text(json.dumps(                         # marker stores the manifest FILE sha, not the dict
        {"unsealed": True, "manifest_file_sha256": sha256_file(FREEZE_DIR / "FREEZE_MANIFEST.json")}, indent=2) + "\n")
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
