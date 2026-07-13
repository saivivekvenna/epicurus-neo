"""Generic, label-blind candidate-universe FREEZE for an ARBITRARY Miller IPV patient.

This is the parameterized going-forward twin of ``scripts/miller_hu287_universe.py``. That Hu_287 script is
FROZEN: its own file SHA-256 is baked into the Hu_287 ``FREEZE_MANIFEST.json`` (``CODE_FILES`` ->
``_all_inputs`` -> ``verify_input_hashes``), so it must stay byte-identical and is never imported-for-mutation
here. We REUSE every pure/patient-independent helper from it (filters, scoring policy, hash/provenance
verifiers) so the science cannot silently diverge, and reimplement ONLY the three pieces that hard-code
Hu_287: the VCF sample-column names, the ``patient_id`` written into the universe, and the per-patient paths /
code-file set / Ensembl cache.

Scope of THIS module: ``freeze`` ONLY. There is deliberately no ``unseal`` and no recognition-label path here
— unsealing a calibration/final Miller patient is a separate, gated step. ``freeze`` opens NO label file.

    from benchmark.miller_patient import load_patient
    from benchmark.miller_universe_core import UniverseConfig, freeze
    freeze(UniverseConfig.for_patient(load_patient("Hu_315")))
"""

from __future__ import annotations

import importlib
import json
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

# The frozen Hu_287 script is the single source of truth for every pure transform / verifier we reuse.
# Importing it does NOT change its bytes (its SHA, and thus the Hu_287 frozen provenance, is untouched).
u = importlib.import_module("scripts.miller_hu287_universe")

ROOT = u.ROOT
REF = u.REF
K = u.K
FROZEN_EPICURUS = u.FROZEN_EPICURUS
FROZEN_EXPR_POLICY = u.FROZEN_EXPR_POLICY
FROZEN_ROUTER = u.FROZEN_ROUTER

# Reused pure helpers (kept as attributes so tests can monkeypatch them exactly like the Hu_287 suite does).
sha256_file = u.sha256_file
variant_key = u.variant_key
norm_chrom = u.norm_chrom
_rel = u._rel
_is_hex_sha = u._is_hex_sha

# "code that shapes output" for the GENERIC lane: this module + its driver + the shared science libs. It
# MUST also pin (a) scripts/miller_hu287_universe.py — the generic core imports and EXECUTES its pure
# helpers/scorers, so those bytes shape every frozen output — and (b) src/benchmark/miller_patient.py, which
# resolves every per-patient path and sample id fed into the freeze.
CODE_FILES = (
    ROOT / "src/benchmark/miller_universe_core.py",
    ROOT / "scripts/miller_patient_universe.py",
    ROOT / "scripts/miller_hu287_universe.py",
    ROOT / "src/benchmark/miller_patient.py",
    ROOT / "src/benchmark/four_arm.py",
    ROOT / "src/event_b/lossless_peptide_generation.py",
    ROOT / "src/event_b/prime_adapter.py",
    ROOT / "src/event_b/prime_transfer.py",
    ROOT / "src/epicurus_neo/evidence_router.py",
)


@dataclass(frozen=True)
class UniverseConfig:
    """Every per-patient path + identity the freeze needs. Built from a label-blind ``MillerPatient``.

    ``for_patient`` reproduces the exact Hu_287 constants when given Hu_287 (asserted in tests), so the
    generic lane is a faithful parameterization of the frozen pipeline, not a re-derivation of the science.
    """

    patient_id: str
    sample_normal: str
    sample_tumor: str
    raw_dir: Path
    pass_vcf: Path
    norm_vcf: Path
    quant: Path
    hla_json: Path
    rna_bam: Path
    ens_cache: Path
    pvac_candidates: Path
    pvac_provenance: Path
    artifact_dir: Path
    freeze_dir: Path
    ref: Path = REF
    code_files: tuple[Path, ...] = CODE_FILES
    semantic_files: tuple[Path, ...] = field(default=())

    @classmethod
    def for_patient(cls, patient) -> "UniverseConfig":
        raw = patient.raw_dir
        art = patient.artifact_dir
        return cls(
            patient_id=patient.patient_id,
            sample_normal=patient.normal_sample,
            sample_tumor=patient.tumor_sample,
            raw_dir=raw,
            pass_vcf=raw / f"somatic/{patient.patient_id}.somatic.pass.vcf.gz",
            norm_vcf=raw / f"somatic/{patient.patient_id}.somatic.pass.norm.vcf.gz",
            quant=raw / "salmon_quant/quant.sf",
            hla_json=art / "HLA_PROVENANCE.json",
            rna_bam=raw / f"rna/{patient.patient_id}_tumor_rna.sorted.bam",
            ens_cache=raw / "ensembl_cache",
            pvac_candidates=raw / "pvac/pvac_candidates.csv",
            pvac_provenance=raw / "pvac/PVAC_PROVENANCE.json",
            artifact_dir=art,
            freeze_dir=raw / "freeze",
            semantic_files=CODE_FILES + (FROZEN_EPICURUS, FROZEN_EXPR_POLICY, FROZEN_ROUTER),
        )


# ---------------------------------------------------------------------------
# Patient-bound reimplementations (the ONLY places Hu_287 was hard-coded)
# ---------------------------------------------------------------------------
def load_filtered_variants(config: UniverseConfig, vcf_path: Path) -> pd.DataFrame:
    """PASS somatic variants + frozen §3.1 base filters, reading the patient's OWN tumor/normal sample
    columns (the sole generalization of the Hu_287 loader). Label-free; consequence assigned later by VEP."""
    import pysam

    rows = []
    for rec in pysam.VariantFile(str(vcf_path)):
        if rec.alts is None:
            continue
        alt = rec.alts[0]
        try:
            tvaf, tdp, talt = u._vaf_depth(rec.samples[config.sample_tumor])
            nvaf, ndp, nalt = u._vaf_depth(rec.samples[config.sample_normal])
        except KeyError:
            continue
        rows.append({
            "key": u.variant_key(rec.chrom, rec.pos, rec.ref, alt), "chrom": u.norm_chrom(rec.chrom),
            "pos": int(rec.pos), "ref": rec.ref, "alt": alt,
            "tumor_vaf": round(tvaf, 4), "normal_vaf": round(nvaf, 4), "tumor_dp": tdp, "normal_dp": ndp,
            "tumor_alt_reads": talt, "normal_alt_reads": nalt,
            "pass_filters": u.passes_base_filters(nvaf, tdp, ndp, talt),
            "strict5_pass": u.legacy_strict5_filters(tvaf, nvaf, tdp, ndp),
        })
    return pd.DataFrame(rows)


def load_pvac_candidates(config: UniverseConfig) -> pd.DataFrame:
    """Genuine pVACseq class-I candidates for THIS patient only (config paths). Same acceptance contract as
    the Hu_287 loader: a CSV alone is insufficient without a real pvacseq/pvactools provenance file."""
    if not (config.pvac_candidates.exists() and config.pvac_provenance.exists()):
        return pd.DataFrame()
    try:
        prov = json.loads(config.pvac_provenance.read_text())
    except Exception:
        return pd.DataFrame()
    if str(prov.get("tool", "")).lower() not in {"pvacseq", "pvactools"} or not prov.get("version"):
        return pd.DataFrame()
    d = pd.read_csv(config.pvac_candidates)
    if not {"mutation_id", "mutant_peptide", "hla_allele"}.issubset(d.columns):
        return pd.DataFrame()
    if not d["mutation_id"].astype(str).map(lambda x: bool(u._KEY_RE.match(x))).all():
        return pd.DataFrame()
    d = d.copy()
    d["hla_allele"] = d["hla_allele"].map(u.normalize_hla)
    d["candidate_source"] = "pvac"
    return d


def build_universe(config: UniverseConfig, variants: pd.DataFrame, hla_panel: list[str],
                   tpm_by_ensg: dict, rna_by_key: dict):
    """Union genuine pVAC rows (first) with the label-blind lossless class-I universe. Identical to the
    Hu_287 builder except the patient's OWN Ensembl cache + ``patient_id`` are used."""
    import numpy as np
    from event_b.lossless_peptide_generation import (EnsemblClient, generate_variant_candidates,
                                                      union_candidates)

    client = EnsemblClient(config.ens_cache)
    frames, notes = [], []
    for v in variants[variants["pass_filters"]].itertuples():
        kind, term = u.classify_consequence(client, v.chrom, v.pos, v.ref, v.alt)
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
        cand["patient_id"] = config.patient_id
        cand["mutation_id"] = v.key
        cand["candidate_source"] = "lossless_recovery"
        cand["gene_symbol"] = res["provenance"].get("gene_symbol", "")
        cand["source_variant_type"] = u._ROUTER_TYPE[kind]
        cand["vep_consequence"] = term
        ensg = str(res["provenance"].get("gene_id", "")).split(".")[0]
        cand["expr"] = float(tpm_by_ensg.get(ensg, 0.0))
        cand["expression_tpm"] = cand["expr"]
        for col in ("tumor_vaf", "normal_vaf", "tumor_dp", "normal_dp", "tumor_alt_reads"):
            cand[col] = getattr(v, col)
        rna = rna_by_key.get(v.key, {})
        cand["rna_alt_obs"] = rna.get("rna_alt_obs", np.nan)
        cand["rna_mutant_reads"] = rna.get("rna_mutant_reads", np.nan)
        cand["rna_depth"] = rna.get("rna_depth", np.nan)
        cand["rna_vaf"] = rna.get("rna_vaf", np.nan)
        frames.append(cand)
    lossless = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    pvac = load_pvac_candidates(config)
    uni = union_candidates([pvac, lossless]) if len(pvac) else lossless
    used = [{"url": url, "cache_path": u._rel(client.cache_dir / e["file"]), "sha256": e["sha256"]}
            for url, e in client.accessed.items()]
    return uni, notes, bool(len(pvac)), used


# ---------------------------------------------------------------------------
# Provenance gates (parameterized versions of the Hu_287 globals-bound gates)
# ---------------------------------------------------------------------------
def _source_inputs(config: UniverseConfig) -> tuple[Path, ...]:
    """Mandatory inputs that MUST pre-exist before freeze proceeds (RNA BAM+.bai included => freeze is
    NOT_EVALUABLE until the tumor-RNA alignment is complete)."""
    return (config.pass_vcf, config.rna_bam, Path(str(config.rna_bam) + ".bai"), config.hla_json,
            config.quant, config.ref, config.ref.parent / (config.ref.name + ".fai"),
            FROZEN_EPICURUS, FROZEN_EXPR_POLICY, FROZEN_ROUTER, *config.code_files)


def _all_inputs(config: UniverseConfig, has_pvac: bool) -> tuple[Path, ...]:
    extra = (config.pvac_candidates, config.pvac_provenance) if has_pvac else ()
    return (*_source_inputs(config), config.norm_vcf, *extra)


def verify_input_hashes(config: UniverseConfig, man: dict) -> tuple[bool, str | None, str | None]:
    """Fail-CLOSED gate on recorded input hashes: input_sha256 must cover EVERY expected input key for this
    lane with a valid 64-hex sha, each recompute-matching on disk."""
    isha = man.get("input_sha256")
    if not isinstance(isha, dict) or not isha:
        return False, None, "input_sha256 missing/empty/malformed"
    for p in _all_inputs(config, bool(man.get("genuine_pvac_lane", False))):
        key = u._rel(p)
        if key not in isha:
            return False, key, "expected input hash key absent"
        if not u._is_hex_sha(isha[key]):
            return False, key, "recorded hash not a valid sha256 (MISSING/empty/malformed)"
        cur = u.sha256_file(ROOT / key) if (ROOT / key).exists() else None
        if cur != isha[key]:
            return False, key, "input hash mismatch or file absent"
    return True, None, None


def verify_ensembl_used(config: UniverseConfig, records, *, require_nonempty: bool) -> tuple[bool, str | None]:
    """Fail-CLOSED verification of the exact Ensembl responses consumed by generation, resolved against THIS
    patient's cache. Mirrors the Hu_287 verifier (traversal guard, dup-URL reject, recompute-match)."""
    if records is None or not isinstance(records, list):
        return False, "records_missing_or_not_a_list"
    if require_nonempty and not records:
        return False, "empty_used_response_set_with_variants_processed"
    cache_root = config.ens_cache.resolve()
    seen: set[str] = set()
    for r in records:
        if not isinstance(r, dict) or not {"url", "cache_path", "sha256"} <= set(r):
            return False, "malformed_record"
        url, cp, sha = r["url"], r["cache_path"], r["sha256"]
        if not isinstance(url, str) or not url.strip():
            return False, "bad_url"
        if not isinstance(cp, str) or not cp:
            return False, "bad_cache_path"
        if not u._is_hex_sha(sha):
            return False, f"bad_sha:{cp}"
        if url in seen:
            return False, f"duplicate_url:{url}"
        seen.add(url)
        p = (ROOT / cp).resolve()
        try:
            p.relative_to(cache_root)
        except ValueError:
            return False, f"path_escapes_cache:{cp}"
        if not p.is_file():
            return False, f"not_a_regular_file:{cp}"
        try:
            actual = u.sha256_file(p)
        except OSError:
            return False, f"hash_io_error:{cp}"
        if actual != sha:
            return False, f"sha_mismatch:{cp}"
    return True, None


def verify_git_tracked_clean(config: UniverseConfig) -> tuple[bool, dict]:
    """Every repo-local semantic file for the GENERIC lane must be git-tracked and byte-identical to HEAD so
    the exact run is reconstructable by checkout. Returns (ok, {rel: status})."""
    import subprocess

    def _rc(*args) -> int:
        return subprocess.call(["git", "-C", str(ROOT), *args], stdout=subprocess.DEVNULL,
                               stderr=subprocess.DEVNULL)

    info, ok = {}, True
    for f in config.semantic_files:
        rel = u._rel(f)
        if _rc("ls-files", "--error-unmatch", "--", rel) != 0:
            status = "UNTRACKED"
        elif _rc("diff", "--cached", "--quiet", "--", rel) != 0:
            status = "STAGED_MODIFIED"
        elif _rc("diff", "--quiet", "--", rel) != 0:
            status = "UNSTAGED_MODIFIED"
        else:
            status = "CLEAN"
        info[rel] = status
        ok = ok and status == "CLEAN"
    return ok, info


def _git_commit() -> str:
    return u._git_commit()


# ---------------------------------------------------------------------------
# FREEZE (label-free). No unseal here by design.
# ---------------------------------------------------------------------------
def freeze(config: UniverseConfig) -> dict:
    """Build + score + select, then LOCK — with NO label access. ONE-SHOT/IMMUTABLE: refuses to overwrite an
    existing valid LOCK. Hu_287 is refused outright (its dedicated frozen script owns that provenance)."""
    from benchmark.four_arm import (DEFAULT_ROUTER_POLICY, FOUR_ARMS, detect_available, evaluate_eligibility)

    if config.patient_id == "Hu_287":
        return {"status": "REFUSED_HU287", "reason": "Hu_287 is owned by scripts/miller_hu287_universe.py; "
                "the generic lane must not touch its frozen provenance."}

    existing = config.freeze_dir / "FREEZE_MANIFEST.json"
    if existing.exists():                                # ONE-SHOT: never overwrite an existing manifest
        try:
            prev = json.loads(existing.read_text())
        except Exception:
            return {"status": "FROZEN_CORRUPT", "reason": "existing FREEZE_MANIFEST.json is unreadable; refusing to overwrite"}
        if prev.get("LOCK") != "FROZEN_NO_LABELS":
            return {"status": "FROZEN_INVALID", "reason": "existing manifest lacks a valid LOCK; refusing to overwrite"}
        ok, bad = u.verify_frozen_hashes(prev, config.freeze_dir)
        if not ok:
            return {"status": "FROZEN_HASH_MISMATCH", "bad_file": bad}
        iok, ibad, ireason = verify_input_hashes(config, prev)
        if not iok:
            return {"status": "FROZEN_INPUT_HASH_MISMATCH", "bad_input": ibad, "reason": ireason}
        eok, ereason = verify_ensembl_used(config, prev.get("ensembl_used_responses"),
                                           require_nonempty=int(prev.get("n_variants_pass", 0)) > 0)
        if not eok:
            return {"status": "FROZEN_ENSEMBL_USED_MISMATCH", "reason": ereason}
        return {"status": "ALREADY_FROZEN", "sha256": prev.get("sha256"), "arms": prev.get("arms"),
                "git_commit": prev.get("git_commit"), "input_verified": True}

    missing = [u._rel(p) for p in _source_inputs(config) if not Path(p).exists()]
    if missing:
        return {"status": "NOT_EVALUABLE", "missing_inputs": missing}
    tools_ok, tool_commits = u.verify_tool_commits()
    if not tools_ok:
        return {"status": "NOT_EVALUABLE", "tool_commit_issue": tool_commits}
    git_ok, git_tracking = verify_git_tracked_clean(config)
    if not git_ok:
        return {"status": "NOT_EVALUABLE", "git_tracking_issue": git_tracking}
    mod_ok, mod_info = u.verify_frozen_module_integrity()
    if not mod_ok:
        return {"status": "NOT_EVALUABLE", "frozen_module_issue": mod_info}

    hla_panel = json.loads(config.hla_json.read_text())["class_i_alleles"]
    norm_vcf = u.normalize_pass_vcf(config.pass_vcf, config.ref, config.norm_vcf)
    variants = load_filtered_variants(config, norm_vcf)
    for col in ("pass_filters", "strict5_pass"):
        if col not in variants.columns:
            return {"status": "NOT_EVALUABLE", "missing_variant_column": col}
    rna_by_key, rna_status = u.rna_alt_evidence(variants, rna_bam=config.rna_bam)
    uni, notes, has_pvac, ensembl_used = build_universe(config, variants, hla_panel,
                                                        u.gene_tpm_by_ensg(config.quant), rna_by_key)
    n_processed = int(variants["pass_filters"].sum())
    eok, ereason = verify_ensembl_used(config, ensembl_used, require_nonempty=n_processed > 0)
    if not eok:
        return {"status": "NOT_EVALUABLE", "ensembl_used_issue": ereason}
    uni, length_counts = u.filter_class_i_lengths(uni)
    if len(uni):
        uni = u.score_universe(uni)

    input_sha = {}
    for p in _all_inputs(config, has_pvac):
        if not Path(p).exists():
            return {"status": "NOT_EVALUABLE", "missing_input": u._rel(p)}
        s = u.sha256_file(p)
        if not u._is_hex_sha(s):
            return {"status": "NOT_EVALUABLE", "bad_input_hash": u._rel(p)}
        input_sha[u._rel(p)] = s

    config.freeze_dir.mkdir(parents=True, exist_ok=True)
    variants.to_csv(config.freeze_dir / "variants.csv", index=False)
    (uni if len(uni) else pd.DataFrame()).to_csv(config.freeze_dir / "universe.csv", index=False)
    available = detect_available(uni, {"__any__"}) if len(uni) else set()
    elig = evaluate_eligibility(available)
    arms_meta = {}
    hashes = {"variants.csv": u.sha256_file(config.freeze_dir / "variants.csv"),
              "universe.csv": u.sha256_file(config.freeze_dir / "universe.csv")}
    for arm in FOUR_ARMS:
        sel = u.arm_selection(uni, arm) if len(uni) else pd.DataFrame(columns=["mutation_id"])
        fn = f"select_{arm.arm_id}.csv"
        sel.to_csv(config.freeze_dir / fn, index=False)
        hashes[fn] = u.sha256_file(config.freeze_dir / fn)
        arms_meta[arm.arm_id] = {"evaluable": bool(elig[arm.arm_id].evaluable),
                                 "missing": elig[arm.arm_id].missing,
                                 "n_selected": int(len(sel)), "saturated": bool(len(sel) >= K),
                                 "selection_file": fn,
                                 "n_unique_mutations": int(sel["mutation_id"].nunique()) if len(sel) else 0}
    manifest = {"patient_id": config.patient_id, "provenance_lane": "generic-miller-universe-core",
                "labels_opened": False, "LOCK": "FROZEN_NO_LABELS",
                "hla_panel": hla_panel, "genuine_pvac_lane": has_pvac,
                "indel_normalization": "bcftools norm -f GRCh38.fa -m-any (left-align+split) before enumeration",
                "class_i_length_filter": length_counts, "rna_alt_evidence_status": rna_status,
                "sample_names": {"normal": config.sample_normal, "tumor": config.sample_tumor},
                "n_variants_pass": int(variants["pass_filters"].sum()), "n_universe_rows": int(len(uni)),
                "n_variants_pass_strict5": int(variants["strict5_pass"].sum()),
                "base_filter_policy": "prereg §3.1 (2026-07-12): Mutect2 PASS + normal VAF<=0.05 + depth>=10 both "
                                      "+ tumor alt-reads>=3; tumor VAF = continuous annotation (NOT a gate)",
                "not_enumerable": notes, "arms": arms_meta, "sha256": hashes,
                "input_sha256": input_sha, "code_files": [u._rel(c) for c in config.code_files],
                "tool_commits": tool_commits, "git_tracked_clean": git_tracking,
                "frozen_module_integrity": mod_info, "ensembl_used_responses": ensembl_used,
                "git_commit": _git_commit(),
                "router_policy_id": DEFAULT_ROUTER_POLICY.policy_id,
                "presentation_evidence": "router binding_percentile_rank = MixMHCpred %rank (real predictor; not NetMHCpan-EL)",
                "el_feature": "NaN (NetMHCpan-EL unavailable; MixMHCpred is NOT a valid el substitute) -> frozen 0.5 fallback",
                "pvac_note": ("genuine pVACseq candidates unioned first" if has_pvac
                              else "no genuine pVAC run -> pvac_prime NOT_EVALUABLE (lossless NEVER relabeled pVAC)")}
    (config.freeze_dir / "FREEZE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    config.artifact_dir.mkdir(parents=True, exist_ok=True)
    (config.artifact_dir / "FREEZE_MANIFEST.json").write_text(json.dumps(manifest, indent=2, default=str) + "\n")
    return manifest
