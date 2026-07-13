"""Fail-closed post-freeze storage lifecycle for sequential Miller calibration patients.

After a patient's label-blind universe FREEZE is locked, the big raw intermediates (per-run FASTQ,
extracted HLA reads) are pure functions of preserved inputs (the ``.sra`` archives, the aligned BAMs)
and can be reclaimed to make room for the next calibration patient. This module PLANS that reclamation
and, only under an independently re-verified freeze, can execute it.

Hard invariants (all enforced here, all tested):
  * **Dry-run by default.** ``plan_cleanup`` never mutates the filesystem. ``execute_cleanup`` refuses
    unless the patient's ``FROZEN_NO_LABELS`` manifest exists AND independently re-verifies every declared
    output hash, every declared input hash, and the pinned code/config provenance — AND ``confirm=True``.
  * **Preserve, always:** ``.sra`` archives, source download/convert/expression manifests + hashes, the
    entire immutable ``freeze/`` (universe + portfolios + provenance), the Ensembl generation cache, and the
    final required BAM/VCF/quant/HLA artifacts. Anything unrecognized is preserved by fail-closed default.
  * **Reclaim, only:** explicitly regenerable FASTQ and documented disposable intermediates, each carrying
    a concrete regeneration recipe. Nothing else is ever proposed for deletion.
  * **Candidate discovery and deletion are confined to the patient's own raw dir.** Verification
    intentionally reads repo-root code/config/reference inputs (to recompute the freeze's declared input
    hashes and provenance), but the set of files that can ever be classified, proposed, or deleted is drawn
    only from a symlink-free walk of ``raw_dir``. The cohort recognition-label CSV is a sibling of that dir,
    is never opened, and a defensive guard additionally refuses any candidate whose real path escapes the
    patient dir or is the label table.
  * **Symlinks are never followed for deletion** and never dereferenced; they are skipped and reported.
  * **Every manifest-declared path is bounds-checked** before it is hashed: output keys must resolve inside
    ``freeze_dir`` and input keys inside the repo root, with no absolute paths, ``..`` traversal, or symlinks.
"""

from __future__ import annotations

import hashlib
import json
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[2]
COHORT_DIR = ROOT / "data/raw/miller_ipv"
# The cohort recognition-label table. NEVER read or deleted by this module. Named only to guard against it.
LABELS_PATH = COHORT_DIR / "miller_recognition_labels.csv"
_HEX = re.compile(r"^[0-9a-f]{64}$")

# Documented disposable intermediates (regenerable from PRESERVED inputs). Kept as an explicit ALLOW-list:
# a file is reclaimable only if it matches one of these AND matches no preservation rule.
FASTQ_SUFFIXES = (".fastq", ".fq", ".fastq.gz", ".fq.gz")


def _sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _is_hex(value) -> bool:
    return isinstance(value, str) and bool(_HEX.match(value))


def _safe_member(base: Path, key: str) -> Path | None:
    """Resolve a manifest hash KEY under ``base`` only if it is a plain relative member: no absolute path,
    no ``..`` traversal, and NO symlink at any component between ``base`` and the target (a symlinked parent
    whose target happens to stay inside ``base`` is still rejected), with the real path inside ``base``."""
    if not isinstance(key, str) or not key or "\x00" in key:
        return None
    if os.path.isabs(key) or key.startswith(("/", "\\")) or "\\" in key:
        return None
    parts = PurePosixPath(key).parts
    if any(part in ("..", "") for part in parts):
        return None
    cur = base
    for part in parts:                                     # reject a symlink at ANY component (parent or leaf)
        cur = cur / part
        if cur.is_symlink():
            return None
    cand = base / key
    try:
        real = cand.resolve()
        real.relative_to(base.resolve())                   # containment backstop
    except (ValueError, OSError):
        return None
    return cand


def _commit_exists(sha: str, root: Path) -> bool:
    return subprocess.call(["git", "-C", str(root), "cat-file", "-e", f"{sha}^{{commit}}"],
                           stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0


def _git_status(rel: str, root: Path) -> str:
    """CLEAN only if ``rel`` is git-tracked and byte-identical to HEAD (no staged/unstaged diff)."""
    def rc(*args) -> int:
        return subprocess.call(["git", "-C", str(root), *args],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if rc("ls-files", "--error-unmatch", "--", rel) != 0:
        return "UNTRACKED"
    if rc("diff", "--cached", "--quiet", "--", rel) != 0:
        return "STAGED_MODIFIED"
    if rc("diff", "--quiet", "--", rel) != 0:
        return "UNSTAGED_MODIFIED"
    return "CLEAN"


@dataclass
class PatientPaths:
    """The only per-patient identity/paths the lifecycle needs. Decoupled from MillerPatient for fixtures."""

    patient_id: str
    raw_dir: Path
    freeze_dir: Path
    root: Path = ROOT

    @classmethod
    def for_patient(cls, patient) -> "PatientPaths":
        return cls(patient_id=patient.patient_id, raw_dir=patient.raw_dir,
                   freeze_dir=patient.raw_dir / "freeze", root=ROOT)

    @classmethod
    def from_patient_id(cls, patient_id: str) -> "PatientPaths":
        from benchmark.miller_patient import load_patient
        return cls.for_patient(load_patient(patient_id))


# ---------------------------------------------------------------------------
# Independent freeze verification (the destructive-mode gate)
# ---------------------------------------------------------------------------
def verify_frozen_no_labels(freeze_dir: Path, root: Path = ROOT, *, git_status=_git_status,
                            commit_exists=_commit_exists) -> tuple[bool, dict]:
    """Re-derive, from scratch, that ``freeze_dir`` holds a valid FROZEN_NO_LABELS lock. Fail-closed on the
    first problem; ``ok`` is True only if EVERY check passes:

      * LOCK == FROZEN_NO_LABELS and labels_opened is exactly False;
      * every declared OUTPUT hash key is a safe member of freeze_dir and recompute-matches;
      * every declared INPUT hash key is a safe member of repo root and recompute-matches;
      * code_files is nonempty and every entry appears in BOTH git_tracked_clean and input_sha256;
      * every git_tracked_clean value recorded as CLEAN AND is currently git-clean on disk;
      * frozen_module_integrity is a complete {module, module_sha256}, contained, non-symlink, hash-matching;
      * git_commit is a full 40-hex sha that resolves to an existing commit object.

    Records the manifest's own sha256 in the returned detail."""
    detail: dict = {"freeze_dir": _rel(freeze_dir, root), "checks": {}}
    manifest_path = freeze_dir / "FREEZE_MANIFEST.json"
    if manifest_path.is_symlink():
        detail["reason"] = "FREEZE_MANIFEST.json is a symlink; the manifest must be a regular file in freeze_dir"
        return False, detail
    if not manifest_path.is_file():
        detail["reason"] = "no FROZEN_NO_LABELS manifest (freeze absent or incomplete)"
        return False, detail
    detail["manifest_sha256"] = _sha256(manifest_path)
    try:
        man = json.loads(manifest_path.read_text())
    except (OSError, ValueError) as exc:
        detail["reason"] = f"manifest unreadable: {exc}"
        return False, detail

    if man.get("LOCK") != "FROZEN_NO_LABELS":
        detail["reason"] = f"LOCK marker is {man.get('LOCK')!r}, expected 'FROZEN_NO_LABELS'"
        return False, detail
    if man.get("labels_opened") is not False:
        detail["reason"] = f"labels_opened is {man.get('labels_opened')!r}, expected False"
        return False, detail

    # (1) declared OUTPUT hashes: keys must be safe members of freeze_dir and recompute-match.
    outputs = man.get("sha256")
    if not isinstance(outputs, dict) or not outputs:
        detail["reason"] = "sha256 (output hashes) missing/empty/malformed"
        return False, detail
    for name, want in outputs.items():
        target = _safe_member(freeze_dir, name)
        if target is None:
            detail["reason"] = f"unsafe output key (absolute/traversal/symlink/escape): {name!r}"
            return False, detail
        if not _is_hex(want):
            detail["reason"] = f"output hash for {name} is not a valid sha256"
            return False, detail
        if not target.is_file() or _sha256(target) != want:
            detail["reason"] = f"output hash mismatch or file absent: {name}"
            return False, detail
    detail["checks"]["outputs_verified"] = len(outputs)

    # (2) declared INPUT hashes: keys must be safe members of the repo root and recompute-match.
    inputs = man.get("input_sha256")
    if not isinstance(inputs, dict) or not inputs:
        detail["reason"] = "input_sha256 missing/empty/malformed"
        return False, detail
    for rel, want in inputs.items():
        target = _safe_member(root, rel)
        if target is None:
            detail["reason"] = f"unsafe input key (absolute/traversal/symlink/escape): {rel!r}"
            return False, detail
        if not _is_hex(want):
            detail["reason"] = f"input hash for {rel} is not a valid sha256"
            return False, detail
        if not target.is_file() or _sha256(target) != want:
            detail["reason"] = f"input hash mismatch or file absent: {rel}"
            return False, detail
    detail["checks"]["inputs_verified"] = len(inputs)

    # (3) code_files: nonempty, and every entry provenance-pinned in BOTH git_tracked_clean and input_sha256.
    code_files = man.get("code_files")
    tracked = man.get("git_tracked_clean")
    if not isinstance(code_files, list) or not code_files:
        detail["reason"] = "code_files missing/empty"
        return False, detail
    if not isinstance(tracked, dict) or not tracked:
        detail["reason"] = "git_tracked_clean provenance missing/empty/malformed"
        return False, detail
    missing_tracked = [c for c in code_files if c not in tracked]
    missing_input = [c for c in code_files if c not in inputs]
    if missing_tracked or missing_input:
        detail["reason"] = (f"code_files not fully pinned: absent from git_tracked_clean={missing_tracked}, "
                            f"from input_sha256={missing_input}")
        return False, detail

    # (4) pinned code/config provenance: every recorded value must be CLEAN AND still git-clean on disk now.
    recorded_bad = {rel: st for rel, st in tracked.items() if st != "CLEAN"}
    if recorded_bad:
        detail["reason"] = f"git_tracked_clean recorded non-CLEAN statuses: {recorded_bad}"
        return False, detail
    live_bad = {rel: git_status(rel, root) for rel in tracked}
    live_bad = {rel: st for rel, st in live_bad.items() if st != "CLEAN"}
    if live_bad:
        detail["reason"] = f"pinned code/config not clean on disk: {live_bad}"
        return False, detail
    detail["checks"]["provenance_files_clean"] = len(tracked)

    # (5) frozen scoring module: mandatory, complete, contained, non-symlink, hash-matching.
    integ = man.get("frozen_module_integrity")
    if not isinstance(integ, dict) or not integ.get("module") or not integ.get("module_sha256"):
        detail["reason"] = "frozen_module_integrity missing/incomplete (need module + module_sha256)"
        return False, detail
    if not _is_hex(integ["module_sha256"]):
        detail["reason"] = "frozen_module_integrity.module_sha256 is not a valid sha256"
        return False, detail
    mod = _safe_member(root, integ["module"])
    if mod is None:
        detail["reason"] = f"unsafe frozen module path: {integ['module']!r}"
        return False, detail
    if not mod.is_file() or _sha256(mod) != integ["module_sha256"]:
        detail["reason"] = f"frozen module changed on disk or absent: {integ['module']}"
        return False, detail
    detail["checks"]["frozen_module_intact"] = True

    # (6) git commit: full 40-hex sha resolving to an existing commit object.
    commit = man.get("git_commit")
    if not (isinstance(commit, str) and re.fullmatch(r"[0-9a-f]{40}", commit)):
        detail["reason"] = f"git_commit is not a full 40-hex sha: {commit!r}"
        return False, detail
    if not commit_exists(commit, root):
        detail["reason"] = f"git_commit does not resolve to an existing commit object: {commit}"
        return False, detail
    detail["git_commit"] = commit
    return True, detail


# ---------------------------------------------------------------------------
# Classification (pure; no filesystem mutation)
# ---------------------------------------------------------------------------
def _rel(path: Path, base: Path) -> str:
    try:
        return str(path.resolve().relative_to(base.resolve()))
    except ValueError:
        return str(path)


def _preserve_reason(rel: str) -> str | None:
    """Return a preservation reason if ``rel`` (patient-dir-relative, POSIX) must be kept, else None."""
    parts = rel.split("/")
    top = parts[0]
    name = parts[-1]
    if name.endswith(".sra"):
        return "SRA archive (the reclaimable FASTQ regenerate from this)"
    if top == "freeze":
        return "immutable frozen universe/portfolios/provenance"
    if top == "ensembl_cache":
        return "Ensembl generation cache (recorded in freeze provenance)"
    if name.endswith("_MANIFEST.json") or name.endswith("_PROVENANCE.json") or name in {
        "DOWNLOAD_MANIFEST.json", "CONVERT_PROVENANCE.json", "EXPRESSION_QUANT.json", "PATIENT_INPUTS.json",
    }:
        return "source manifest / provenance + hashes"
    if name.endswith((".bam", ".bai", ".crai", ".cram", ".dict", ".fai")):
        return "final required alignment artifact (BAM/index)"
    if name.endswith((".vcf", ".vcf.gz", ".tbi", ".vcf.gz.tbi")):
        return "final required somatic VCF"
    if top == "salmon_quant" or name == "quant.sf":
        return "final required RNA quantification"
    if top == "hla" and (name.endswith("_result.tsv") or name.endswith(".pdf") or name == "HLA_PROVENANCE.json"):
        return "final required HLA typing result"
    return None


def _regular_contained(p: Path, base: Path) -> bool:
    """True iff ``p`` is an existing regular, non-symlink file whose real path stays inside ``base``.
    Uses only stat/resolve — never opens the file (so a label table can never be read here)."""
    try:
        if p.is_symlink() or not p.is_file():
            return False
        p.resolve().relative_to(base.resolve())
        return True
    except (OSError, ValueError):
        return False


def _reclaim_status(rel: str, raw_dir: Path) -> tuple[str | None, str | None]:
    """Classify a documented disposable intermediate. Returns (recipe, blocked_reason):
      * (recipe, None)  -> reclaimable; the preserved regeneration SOURCE was proven to exist;
      * (None, reason)  -> withheld: it *looks* regenerable but its source is missing → preserve fail-closed;
      * (None, None)    -> not a documented disposable intermediate at all.
    Source existence is proven so we never delete a FASTQ/extracted-read whose regeneration input is gone."""
    parts = rel.split("/")
    top, name = parts[0], parts[-1]
    if top == "fastq" and name.endswith(FASTQ_SUFFIXES):
        run = name.split("_")[0].split(".")[0]
        if _regular_contained(raw_dir / f"{run}.sra", raw_dir):
            return f"regenerate: fasterq-dump {run}.sra (preserved SRA archive)", None
        return None, f"withheld: regeneration source {run}.sra is missing or not a regular contained file"
    if top == "hla" and name.endswith(".fq"):
        somatic = raw_dir / "somatic"
        bams = sorted(somatic.glob("*_N.md.bam")) if somatic.is_dir() else []
        src = next((b for b in bams if _regular_contained(b, raw_dir)), None)
        if src is not None:
            return (f"regenerate: re-extract the GRCh38 MHC region from the preserved normal alignment "
                    f"{src.name} via the generic reconstruction HLA stage "
                    "(scripts/miller_patient_reconstruct.py <PATIENT_ID> hla)"), None
        return None, "withheld: no preserved normal alignment (somatic/*_N.md.bam) to regenerate the HLA reads from"
    return None, None


def classify_entries(paths: PatientPaths) -> list[dict]:
    """Walk the patient raw dir (never following symlinked dirs) and classify every regular file.

    Categories: REMOVE (documented regenerable), PRESERVE (required or fail-closed default),
    SKIP_SYMLINK (never followed/deleted), SKIP_ESCAPES (real path outside the patient dir — defensive)."""
    raw = paths.raw_dir
    raw_real = raw.resolve()
    entries: list[dict] = []
    if not raw.is_dir():
        return entries
    for dirpath, dirnames, filenames in os.walk(raw, followlinks=False):
        dirnames.sort()
        for fname in sorted(filenames):
            p = Path(dirpath) / fname
            rel = os.path.relpath(p, raw).replace(os.sep, "/")
            is_link = p.is_symlink()
            # Size without dereferencing symlinks (lstat), so a dangling/hostile link never triggers I/O.
            try:
                size = os.lstat(p).st_size
            except OSError:
                size = 0
            rec = {"rel": rel, "size_bytes": int(size), "is_symlink": bool(is_link)}
            if is_link:
                rec.update(category="SKIP_SYMLINK",
                           reason="symlink — never followed or deleted (label/target safety)")
                entries.append(rec)
                continue
            # Defensive containment: a regular file resolving outside the patient dir, or the label table.
            real = p.resolve()
            if real == LABELS_PATH.resolve():
                rec.update(category="SKIP_ESCAPES", reason="resolves to the recognition-label table — never touched")
                entries.append(rec)
                continue
            try:
                real.relative_to(raw_real)
            except ValueError:
                rec.update(category="SKIP_ESCAPES", reason="real path escapes the patient dir")
                entries.append(rec)
                continue
            preserve = _preserve_reason(rel)
            recipe, blocked = _reclaim_status(rel, raw)
            if preserve is None and recipe is not None:
                rec.update(category="REMOVE", reason="documented regenerable intermediate", regeneration=recipe)
            elif preserve is None and blocked is not None:
                rec.update(category="PRESERVE", reason=blocked)      # regenerable-shaped but source missing
            else:
                rec.update(category="PRESERVE",
                           reason=preserve or "unclassified — preserved by fail-closed default")
            entries.append(rec)
    return entries


def _summary(entries: list[dict]) -> dict:
    out: dict = {}
    for cat in ("PRESERVE", "REMOVE", "SKIP_SYMLINK", "SKIP_ESCAPES"):
        rows = [e for e in entries if e["category"] == cat]
        out[cat] = {"n": len(rows), "bytes": sum(e["size_bytes"] for e in rows)}
    return out


# ---------------------------------------------------------------------------
# Plan (dry-run) and Execute (fail-closed)
# ---------------------------------------------------------------------------
def plan_cleanup(paths: PatientPaths, *, git_status=_git_status, commit_exists=_commit_exists) -> dict:
    """Dry-run plan. NEVER mutates the filesystem. Reclamation is *eligible* only when the freeze
    independently re-verifies; otherwise every reclaim candidate is withheld pending freeze."""
    frozen_ok, verification = verify_frozen_no_labels(paths.freeze_dir, paths.root, git_status=git_status,
                                                      commit_exists=commit_exists)
    entries = classify_entries(paths)
    reclaimable = [e for e in entries if e["category"] == "REMOVE"]
    report = {
        "patient_id": paths.patient_id,
        "raw_dir": _rel(paths.raw_dir, paths.root),
        "mode": "dry_run",
        "frozen_verified": frozen_ok,
        "manifest_sha256": verification.get("manifest_sha256"),
        "verification": verification,
        "labels_guard": ("verification recomputes hashes of the freeze's declared repo-root inputs "
                         "(code/config/reference); candidate discovery and deletion are confined to "
                         f"{_rel(paths.raw_dir, paths.root)}; the recognition-label table is never opened"),
        "summary": _summary(entries),
        "entries": entries,
    }
    if frozen_ok:
        report["reclaimable"] = reclaimable
        report["reclaimable_bytes"] = sum(e["size_bytes"] for e in reclaimable)
        report["withheld_pending_freeze"] = []
    else:
        # Fail-closed lifecycle ordering: nothing is reclaimable until the freeze verifies.
        report["reclaimable"] = []
        report["reclaimable_bytes"] = 0
        report["withheld_pending_freeze"] = reclaimable
    return report


def execute_cleanup(paths: PatientPaths, *, confirm: bool = False, git_status=_git_status,
                    commit_exists=_commit_exists) -> dict:
    """Destructive reclamation — REFUSES unless the freeze independently re-verifies AND confirm=True.

    Even when permitted, each deletion re-checks (exists, regular file, not a symlink, contained in the
    patient dir) immediately before unlinking. Only REMOVE-classified files are ever touched."""
    report = plan_cleanup(paths, git_status=git_status, commit_exists=commit_exists)
    report["mode"] = "execute"
    if not report["frozen_verified"]:
        report["status"] = "REFUSED_UNVERIFIED_FREEZE"
        report["deleted"] = []
        report["deleted_bytes"] = 0
        return report
    if not confirm:
        report["status"] = "REFUSED_NO_CONFIRM"
        report["deleted"] = []
        report["deleted_bytes"] = 0
        return report

    raw_real = paths.raw_dir.resolve()
    deleted: list[dict] = []
    for e in report["reclaimable"]:
        p = paths.raw_dir / e["rel"]
        if not p.exists() or p.is_symlink() or not p.is_file():
            continue
        real = p.resolve()
        if real == LABELS_PATH.resolve():
            continue
        try:
            real.relative_to(raw_real)
        except ValueError:
            continue
        os.remove(p)
        deleted.append({"rel": e["rel"], "size_bytes": e["size_bytes"], "regeneration": e.get("regeneration")})
    report["status"] = "EXECUTED"
    report["deleted"] = deleted
    report["deleted_bytes"] = sum(d["size_bytes"] for d in deleted)
    return report
