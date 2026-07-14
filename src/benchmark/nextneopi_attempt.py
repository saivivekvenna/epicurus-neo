"""Fail-closed finalization of a label-blind nextNEOpi Track-A attempt.

This module deliberately does not know where recognition labels live.  It turns
an execution (including a verified failure) into a checksum-bound attempt that
can be joined to outcomes only by a later evaluator.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable, Mapping


COMPARATOR_ID = "nextneopi_track_a"
POLICY_ID = "nextneopi-track-a-attempt-v1"
SUCCESS = "SUCCEEDED"
NON_SUCCESS = frozenset({"FAILED", "ABSTAINED"})
REQUIRED_AGGREGATE_COLUMNS = frozenset({"ID", "Best Peptide", "Allele", "Tier"})
Normalizer = Callable[[Path, Path, Path], None]


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _file(path: str | Path, *, nonempty: bool = True) -> Path:
    unresolved = Path(path).expanduser()
    if unresolved.is_symlink():
        raise ValueError(f"required evidence file is a symlink: {unresolved}")
    resolved = unresolved.resolve()
    if not resolved.is_file() or resolved.is_symlink():
        raise FileNotFoundError(resolved)
    if nonempty and resolved.stat().st_size == 0:
        raise ValueError(f"required evidence file is empty: {resolved}")
    return resolved


def _attestation(path: Path) -> dict:
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256_file(path)}


def _load_json(path: Path) -> dict:
    try:
        value = json.loads(path.read_text())
    except (OSError, ValueError) as exc:
        raise ValueError(f"invalid JSON evidence: {path}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"JSON evidence must be an object: {path}")
    return value


@dataclass(frozen=True)
class LockedAttempt:
    patient_id: str
    input_manifest: Path
    source: Path
    config: Path
    instrumentation: Path
    execution_log: Path
    execution_trace: Path
    input: dict
    frozen_config: dict


def _lock_common(
    *,
    patient_id: str,
    input_manifest: str | Path,
    source: str | Path,
    config: str | Path,
    instrumentation: str | Path,
    execution_log: str | Path,
    execution_trace: str | Path,
) -> LockedAttempt:
    if not patient_id or any(c in patient_id for c in "\r\n,"):
        raise ValueError("patient_id must be non-empty and manifest-safe")
    manifest_path = _file(input_manifest)
    source_path = _file(source)
    config_path = _file(config)
    patch_path = _file(instrumentation)
    log_path = _file(execution_log)
    trace_path = _file(execution_trace)
    manifest = _load_json(manifest_path)
    frozen = _load_json(config_path)

    if manifest.get("patient_id") != patient_id:
        raise ValueError("INPUT_MANIFEST patient_id does not match attempt")
    if manifest.get("labels_opened") is not False:
        raise ValueError("INPUT_MANIFEST is not label-blind")
    if frozen.get("policy_id") != manifest.get("policy_id"):
        raise ValueError("frozen config policy differs from INPUT_MANIFEST")

    expected_config = (manifest.get("config") or {}).get("sha256")
    expected_source = (manifest.get("upstream") or {}).get("nextneopi_nf_sha256")
    expected_patch = (manifest.get("instrumentation_patch") or {}).get("sha256")
    if not all(isinstance(x, str) and len(x) == 64 for x in
               (expected_config, expected_source, expected_patch)):
        raise ValueError("INPUT_MANIFEST lacks pinned source/config/instrumentation hashes")
    if sha256_file(config_path) != expected_config:
        raise ValueError("frozen config hash differs from INPUT_MANIFEST")
    if sha256_file(source_path) != expected_source:
        raise ValueError("nextNEOpi source hash differs from pinned INPUT_MANIFEST")
    if sha256_file(patch_path) != expected_patch:
        raise ValueError("instrumentation hash differs from pinned INPUT_MANIFEST")
    if (frozen.get("upstream") or {}).get("nextneopi_nf_sha256") != expected_source:
        raise ValueError("frozen config and INPUT_MANIFEST pin different nextNEOpi sources")
    if (frozen.get("runtime") or {}).get("instrumentation_patch_sha256") != expected_patch:
        raise ValueError("frozen config and INPUT_MANIFEST pin different instrumentation")
    return LockedAttempt(
        patient_id, manifest_path, source_path, config_path, patch_path,
        log_path, trace_path, manifest, frozen,
    )


def _open_text(path: Path):
    return gzip.open(path, "rt") if path.name.endswith(".gz") else path.open()


@dataclass(frozen=True)
class VcfAllele:
    chrom: str
    pos: int
    ref: str
    alt: str
    line_number: int

    @property
    def pvac_id(self) -> str:
        # pVACtools uses VEP's zero-based half-open coordinates but retains the
        # original VCF REF/ALT strings in the aggregate ID.  Strip common VCF
        # padding only to derive the affected span (including zero-span insertions).
        prefix = 0
        limit = min(len(self.ref), len(self.alt))
        while prefix < limit and self.ref[prefix] == self.alt[prefix]:
            prefix += 1
        ref_end, alt_end = len(self.ref), len(self.alt)
        while ref_end > prefix and alt_end > prefix and self.ref[ref_end - 1] == self.alt[alt_end - 1]:
            ref_end -= 1
            alt_end -= 1
        start = self.pos - 1 + prefix
        stop = start + (ref_end - prefix)
        return f"{self.chrom}-{start}-{stop}-{self.ref}-{self.alt}"

    @property
    def input_key(self) -> str:
        return f"{self.chrom}:{self.pos}:{self.ref}:{self.alt}"


def _read_vcf_alleles(path: Path) -> tuple[list[str], dict[str, list[VcfAllele]]]:
    headers: list[str] = []
    by_id: dict[str, list[VcfAllele]] = {}
    saw_columns = False
    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith("#"):
                headers.append(raw if raw.endswith("\n") else raw + "\n")
                saw_columns |= raw.startswith("#CHROM\t")
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"malformed pVAC input VCF line {line_number}")
            chrom, pos_text, _identifier, ref, alt_text = fields[:5]
            try:
                pos = int(pos_text)
            except ValueError as exc:
                raise ValueError(f"invalid pVAC input position at line {line_number}") from exc
            if pos < 1 or not chrom or not ref or ref == ".":
                raise ValueError(f"invalid pVAC input allele at line {line_number}")
            for alt in alt_text.split(","):
                if not alt or alt == "." or alt.startswith("<") or "[" in alt or "]" in alt:
                    raise ValueError(f"unsupported pVAC input ALT at line {line_number}: {alt}")
                allele = VcfAllele(chrom, pos, ref, alt, line_number)
                by_id.setdefault(allele.pvac_id, []).append(allele)
    if not saw_columns:
        raise ValueError("pVAC input VCF lacks #CHROM header")
    if not by_id:
        raise ValueError("pVAC input VCF has no variants")
    return headers, by_id


def _read_aggregate(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = REQUIRED_AGGREGATE_COLUMNS.difference(reader.fieldnames or ())
        if missing:
            raise ValueError(f"native aggregate missing columns: {sorted(missing)}")
        rows = list(reader)
        columns = list(reader.fieldnames or ())
    if not rows:
        raise ValueError("native aggregate is empty")
    ids = [row["ID"] for row in rows]
    if any(not item for item in ids) or len(ids) != len(set(ids)):
        raise ValueError("native aggregate IDs must be non-empty and unique")
    return columns, rows


def _write_variant_vcf(path: Path, headers: Iterable[str], mapped: list[tuple[str, VcfAllele]]) -> None:
    meta = [line for line in headers if line.startswith("##")]
    with path.open("w") as handle:
        handle.writelines(meta)
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for token, allele in mapped:
            handle.write(
                f"{allele.chrom}\t{allele.pos}\t{token}\t{allele.ref}\t{allele.alt}"
                "\t.\tPASS\t.\n"
            )


def _bcftools_normalizer(input_vcf: Path, reference: Path, output_vcf: Path) -> None:
    try:
        completed = subprocess.run(
            ["bcftools", "norm", "-f", str(reference), "-m-any", "-Ov", "-o",
             str(output_vcf), str(input_vcf)],
            check=False, capture_output=True, text=True,
        )
    except FileNotFoundError as exc:
        raise RuntimeError("bcftools is required to finalize a Track-A success") from exc
    if completed.returncode != 0:
        raise RuntimeError(f"bcftools norm failed: {completed.stderr.strip()}")


def _parse_normalized(path: Path, expected_tokens: set[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    with _open_text(path) as handle:
        for line_number, raw in enumerate(handle, 1):
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 5:
                raise ValueError(f"malformed normalized VCF line {line_number}")
            chrom, pos, token, ref, alt = fields[:5]
            if token not in expected_tokens:
                raise ValueError(f"normalizer emitted unknown mapping token: {token}")
            if token in result or "," in alt:
                raise ValueError(f"normalizer did not emit one allele for token: {token}")
            result[token] = f"{chrom}:{int(pos)}:{ref}:{alt}"
    missing = expected_tokens.difference(result)
    if missing:
        raise ValueError(f"normalizer dropped mapped variants: {sorted(missing)}")
    if len(set(result.values())) != len(result):
        raise ValueError("distinct aggregate IDs collapse to canonical duplicate mutations")
    return result


def _normalize_set(
    *, name: str, destination: Path, headers: list[str],
    mapped: list[tuple[str, VcfAllele]], reference: Path, normalizer: Normalizer,
) -> tuple[Path, dict[str, str]]:
    source = destination / f"{name}_mapped.vcf"
    normalized = destination / f"{name}_normalized.vcf"
    _write_variant_vcf(source, headers, mapped)
    normalizer(source, reference, normalized)
    normalized = _file(normalized)
    return normalized, _parse_normalized(normalized, {token for token, _ in mapped})


def normalize_variant_rows(
    rows: Iterable[Mapping[str, object]],
    reference_fasta: str | Path,
    normalizer: Normalizer | None = None,
    *,
    row_id_field: str = "row_id",
) -> list[dict]:
    """Canonicalize labeled or unlabeled variant rows without changing their identity/order.

    Each row must contain ``row_id`` (or ``row_id_field``), ``chrom``, ``pos``, ``ref``
    and ``alt``.  The normalizer must return exactly one allele per row.  In particular,
    this contract never decomposes an MNV into independently rescuable components.
    """
    reference = _file(reference_fasta)
    copied = [dict(row) for row in rows]
    if not copied:
        return []
    identifiers: list[str] = []
    mapped: list[tuple[str, VcfAllele]] = []
    for index, row in enumerate(copied, 1):
        required = {row_id_field, "chrom", "pos", "ref", "alt"}
        missing = required.difference(row)
        if missing:
            raise ValueError(f"variant row missing fields: {sorted(missing)}")
        identifier = str(row[row_id_field])
        if not identifier:
            raise ValueError("variant row IDs must be non-empty")
        identifiers.append(identifier)
        try:
            position = int(row["pos"])
        except (TypeError, ValueError) as exc:
            raise ValueError(f"invalid position for row {identifier}") from exc
        chrom, ref, alt = str(row["chrom"]), str(row["ref"]), str(row["alt"])
        if position < 1 or not chrom or not ref or not alt or "," in alt:
            raise ValueError(f"invalid single-allele variant for row {identifier}")
        if alt.startswith("<") or "[" in alt or "]" in alt:
            raise ValueError(f"unsupported ALT for row {identifier}")
        mapped.append((f"ROWMAP{index:08d}", VcfAllele(chrom, position, ref, alt, index)))
    if len(set(identifiers)) != len(identifiers):
        raise ValueError("variant row IDs must be unique")
    normalize = normalizer or _bcftools_normalizer
    with tempfile.TemporaryDirectory(prefix="nextneopi-normalize-") as temp:
        directory = Path(temp)
        # A minimal header is intentional: annotations and genotypes cannot influence identity.
        normalized_path, normalized = _normalize_set(
            name="rows", destination=directory,
            headers=["##fileformat=VCFv4.2\n", "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"],
            mapped=mapped, reference=reference, normalizer=normalize,
        )
        # Force parsing before the temporary evidence disappears.
        if not normalized_path.is_file():  # pragma: no cover - guarded by _normalize_set
            raise ValueError("normalizer output disappeared")
    result = []
    for index, row in enumerate(copied, 1):
        result.append({**row, "canonical_mutation": normalized[f"ROWMAP{index:08d}"]})
    return result


def _canonical_key(value: str) -> bool:
    parts = value.rsplit(":", 3)
    if len(parts) != 4 or not all(parts):
        return False
    try:
        return int(parts[1]) >= 1
    except ValueError:
        return False


def _verify_attestation(value: object, *, generated_root: Path | None = None) -> Path:
    if not isinstance(value, dict):
        raise ValueError("file attestation must be an object")
    if not {"path", "bytes", "sha256"}.issubset(value):
        raise ValueError("file attestation lacks path/bytes/sha256")
    raw = value["path"]
    if not isinstance(raw, str) or not Path(raw).is_absolute() or ".." in Path(raw).parts:
        raise ValueError("attested path must be absolute and traversal-free")
    unresolved = Path(raw)
    if unresolved.is_symlink():
        raise ValueError(f"attested path is a symlink: {raw}")
    path = _file(unresolved)
    if generated_root is not None:
        try:
            path.relative_to(generated_root)
        except ValueError as exc:
            raise ValueError(f"generated artifact escapes attempt directory: {raw}") from exc
    if not isinstance(value["bytes"], int) or value["bytes"] != path.stat().st_size:
        raise ValueError(f"attested size mismatch: {raw}")
    digest = value["sha256"]
    if not isinstance(digest, str) or len(digest) != 64 or sha256_file(path) != digest:
        raise ValueError(f"attested hash mismatch: {raw}")
    return path


def verify_attempt_manifest(
    path: str | Path,
    expected_patient_id: str,
    expected_comparator_id: str = COMPARATOR_ID,
) -> tuple[bool, dict]:
    """Pure, label-free preflight returning an immutable evaluator snapshot.

    No exception crosses this boundary: an unverified attempt is represented as ``ok=False``.
    Successful snapshots expose native-order canonical mutation IDs; verified failures and
    abstentions expose an empty ordered list and remain evaluable.
    """
    try:
        manifest_path = _file(path)
        manifest = _load_json(manifest_path)
        root = manifest_path.parent.resolve()
        if manifest_path.name != "ATTEMPT_MANIFEST.json":
            raise ValueError("unexpected attempt manifest basename")
        if manifest.get("schema_version") != 1 or manifest.get("policy_id") != POLICY_ID:
            raise ValueError("unsupported attempt manifest schema/policy")
        if manifest.get("patient_id") != expected_patient_id:
            raise ValueError("attempt patient_id mismatch")
        if manifest.get("comparator_id") != expected_comparator_id:
            raise ValueError("attempt comparator_id mismatch")
        if manifest.get("scope") != "FULL_PIPELINE_IDENTICAL_RAW_INPUT":
            raise ValueError("attempt scope mismatch")
        if manifest.get("labels_opened") is not False or manifest.get("evaluable") is not True:
            raise ValueError("attempt is not a label-blind evaluable record")
        status = manifest.get("status")
        execution = manifest.get("execution")
        if status not in {SUCCESS, *NON_SUCCESS} or not isinstance(execution, dict):
            raise ValueError("invalid attempt status/execution")
        if execution.get("status") != status:
            raise ValueError("top-level and execution statuses differ")

        input_path = _verify_attestation(manifest.get("input_manifest"))
        pins = manifest.get("pinned_provenance")
        if not isinstance(pins, dict):
            raise ValueError("missing pinned provenance")
        source_path = _verify_attestation(pins.get("source"))
        config_path = _verify_attestation(pins.get("config"))
        patch_path = _verify_attestation(pins.get("instrumentation"))
        _verify_attestation(execution.get("log"))
        _verify_attestation(execution.get("trace"))

        locked_input = _load_json(input_path)
        locked_config = _load_json(config_path)
        if locked_input.get("patient_id") != expected_patient_id or locked_input.get("labels_opened") is not False:
            raise ValueError("attested INPUT_MANIFEST identity/label lock mismatch")
        if (locked_input.get("config") or {}).get("sha256") != sha256_file(config_path):
            raise ValueError("INPUT_MANIFEST does not bind attested config")
        if (locked_input.get("upstream") or {}).get("nextneopi_nf_sha256") != sha256_file(source_path):
            raise ValueError("INPUT_MANIFEST does not bind attested source")
        if (locked_input.get("instrumentation_patch") or {}).get("sha256") != sha256_file(patch_path):
            raise ValueError("INPUT_MANIFEST does not bind attested instrumentation")
        if (locked_config.get("upstream") or {}).get("nextneopi_nf_sha256") != sha256_file(source_path):
            raise ValueError("config does not bind attested source")
        if (locked_config.get("runtime") or {}).get("instrumentation_patch_sha256") != sha256_file(patch_path):
            raise ValueError("config does not bind attested instrumentation")
        if pins.get("upstream_commit") != (locked_config.get("upstream") or {}).get("commit"):
            raise ValueError("attempt and config pin different upstream commits")
        observed = _execution_identity(
            LockedAttempt(
                expected_patient_id, input_path, source_path, config_path, patch_path,
                Path(execution["log"]["path"]), Path(execution["trace"]["path"]),
                locked_input, locked_config,
            ),
            command=execution.get("command"), started_at=execution.get("started_at"),
            finished_at=execution.get("finished_at"), runtime=execution.get("runtime"),
        )
        if execution.get("duration_seconds") != observed["duration_seconds"]:
            raise ValueError("execution duration does not match timestamps")

        mutation_ids: list[str] = []
        if status == SUCCESS:
            if execution.get("exit_code") != 0 or manifest.get("attempt_outcome") != "PORTFOLIO_FROZEN":
                raise ValueError("successful attempt lacks successful execution/freeze")
            aggregate = manifest.get("native_aggregate")
            published = manifest.get("published_pvac_input")
            normalization = manifest.get("normalization")
            portfolio = manifest.get("portfolio")
            if not all(isinstance(x, dict) for x in (aggregate, published, normalization, portfolio)):
                raise ValueError("successful attempt lacks required native evidence")
            aggregate_path = _verify_attestation(aggregate)
            published_vcf_path = _verify_attestation(published.get("vcf"))
            _verify_attestation(published.get("index"))
            _verify_attestation(manifest.get("grch38_reference"))
            all_norm_path = _verify_attestation(
                normalization.get("all_variants_vcf"), generated_root=root
            )
            selected_norm_path = _verify_attestation(
                normalization.get("selected_variants_vcf"), generated_root=root
            )
            map_path = _verify_attestation(normalization.get("map"), generated_root=root)
            portfolio_path = _verify_attestation(portfolio, generated_root=root)
            _aggregate_columns, aggregate_rows = _read_aggregate(aggregate_path)
            _headers, vcf_by_id = _read_vcf_alleles(published_vcf_path)
            if aggregate.get("rows") != len(aggregate_rows):
                raise ValueError("native aggregate row attestation mismatch")
            expected_all_tokens = {f"PVACMAP{i:08d}" for i in range(1, len(aggregate_rows) + 1)}
            all_normalized = _parse_normalized(all_norm_path, expected_all_tokens)
            expected_top_tokens = {
                f"PVACMAP{i:08d}" for i in range(1, min(20, len(aggregate_rows)) + 1)
            }
            selected_normalized = _parse_normalized(selected_norm_path, expected_top_tokens)
            with portfolio_path.open(newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            if int(portfolio.get("size", -1)) != len(rows) or len(rows) > 20:
                raise ValueError("portfolio size does not match ordered rows")
            if [row.get("rank") for row in rows] != [str(i) for i in range(1, len(rows) + 1)]:
                raise ValueError("portfolio ranks are not contiguous native order")
            mutation_ids = [row.get("canonical_mutation", "") for row in rows]
            if any(not _canonical_key(item) for item in mutation_ids):
                raise ValueError("portfolio contains invalid canonical mutation key")
            if len(mutation_ids) != len(set(mutation_ids)):
                raise ValueError("portfolio contains canonical duplicate mutations")
            with map_path.open(newline="") as handle:
                mappings = list(csv.DictReader(handle, delimiter="\t"))
            if len(mappings) != normalization.get("mapped_aggregate_ids"):
                raise ValueError("normalization map row count mismatch")
            if len(mappings) != len(aggregate_rows):
                raise ValueError("normalization map and native aggregate row counts differ")
            for rank, (mapping, aggregate_row) in enumerate(zip(mappings, aggregate_rows), 1):
                token = f"PVACMAP{rank:08d}"
                matches = vcf_by_id.get(aggregate_row["ID"], [])
                if len(matches) != 1:
                    raise ValueError("native aggregate no longer maps one-to-one to pVAC input")
                expected_mapping = {
                    "aggregate_rank": str(rank),
                    "aggregate_id": aggregate_row["ID"],
                    "mapping_token": token,
                    "input_variant": matches[0].input_key,
                    "canonical_variant": all_normalized[token],
                    "selected_top20": str(rank <= 20).lower(),
                }
                if mapping != expected_mapping:
                    raise ValueError(f"normalization map semantic mismatch at rank {rank}")
                if rank <= 20 and selected_normalized[token] != all_normalized[token]:
                    raise ValueError(f"selected/all normalized VCF mismatch at rank {rank}")
            selected = [row.get("canonical_variant") for row in mappings
                        if row.get("selected_top20") == "true"]
            if selected != mutation_ids:
                raise ValueError("normalization map and ordered portfolio differ")
        else:
            exit_code = execution.get("exit_code")
            if status == "FAILED" and (not isinstance(exit_code, int) or exit_code == 0):
                raise ValueError("verified FAILED attempt lacks nonzero exit")
            if status == "ABSTAINED" and exit_code == 0:
                raise ValueError("verified ABSTAINED attempt carries success exit")
            portfolio = manifest.get("portfolio")
            if not isinstance(portfolio, dict) or portfolio.get("size") != 0 or portfolio.get("rows") != []:
                raise ValueError("non-success attempt must have an explicit zero portfolio")
            if any(manifest.get(key) is not None for key in
                   ("native_aggregate", "published_pvac_input", "grch38_reference", "normalization")):
                raise ValueError("non-success attempt claims success-only artifacts")
            if not isinstance(manifest.get("reason"), str) or not manifest["reason"].strip():
                raise ValueError("non-success attempt lacks reason")

        return True, {
            "schema_version": 1,
            "manifest_path": str(manifest_path),
            "manifest_sha256": sha256_file(manifest_path),
            "patient_id": expected_patient_id,
            "comparator_id": expected_comparator_id,
            "status": status,
            "evaluable": True,
            "mutation_ids": mutation_ids,
        }
    except (OSError, ValueError, TypeError) as exc:
        return False, {"reason": str(exc), "evaluable": False, "mutation_ids": []}


def _common_manifest(lock: LockedAttempt, status: str, exit_code: int | None) -> dict:
    return {
        "schema_version": 1,
        "policy_id": POLICY_ID,
        "comparator_id": COMPARATOR_ID,
        "scope": "FULL_PIPELINE_IDENTICAL_RAW_INPUT",
        "patient_id": lock.patient_id,
        "labels_opened": False,
        "evaluable": True,
        "status": status,
        "execution": {
            "status": status,
            "exit_code": exit_code,
            "log": _attestation(lock.execution_log),
            "trace": _attestation(lock.execution_trace),
        },
        "input_manifest": _attestation(lock.input_manifest),
        "pinned_provenance": {
            "source": _attestation(lock.source),
            "config": _attestation(lock.config),
            "instrumentation": _attestation(lock.instrumentation),
            "upstream_commit": (lock.frozen_config.get("upstream") or {}).get("commit"),
        },
    }


def _parse_time(value: str, name: str) -> datetime:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} is required")
    try:
        return datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError(f"{name} must be an ISO-8601 timestamp") from exc


def _execution_identity(
    lock: LockedAttempt, *, command: str, started_at: str, finished_at: str, runtime: str,
) -> dict:
    if command != lock.input.get("command"):
        raise ValueError("executed command differs from locked INPUT_MANIFEST command")
    started = _parse_time(started_at, "started_at")
    finished = _parse_time(finished_at, "finished_at")
    if started.tzinfo is None or finished.tzinfo is None:
        raise ValueError("execution timestamps must include timezone offsets")
    if finished < started:
        raise ValueError("finished_at precedes started_at")
    if not isinstance(runtime, str) or not runtime.strip():
        raise ValueError("runtime identity is required")
    return {
        "command": command,
        "started_at": started_at,
        "finished_at": finished_at,
        "duration_seconds": (finished - started).total_seconds(),
        "runtime": runtime.strip(),
    }


def _empty_destination(output_dir: str | Path) -> Path:
    destination = Path(output_dir).expanduser().resolve()
    if destination.exists():
        if destination.is_symlink() or not destination.is_dir():
            raise ValueError("attempt output must be a real directory")
        if any(destination.iterdir()):
            raise FileExistsError(f"attempt output is nonempty; refusing overwrite: {destination}")
    else:
        destination.mkdir(parents=True)
    return destination


def _write_once_manifest(destination: Path, value: dict) -> Path:
    """Durably publish a manifest without an overwrite window."""
    final = destination / "ATTEMPT_MANIFEST.json"
    if final.exists() or final.is_symlink():
        raise FileExistsError(f"attempt manifest already exists: {final}")
    temporary = destination / f".ATTEMPT_MANIFEST.{os.getpid()}.tmp"
    payload = (json.dumps(value, indent=2, sort_keys=True) + "\n").encode()
    fd = os.open(temporary, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    try:
        with os.fdopen(fd, "wb", closefd=True) as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(temporary, final)
        os.chmod(final, 0o444)
        directory_fd = os.open(destination, os.O_RDONLY)
        try:
            os.fsync(directory_fd)
        finally:
            os.close(directory_fd)
    finally:
        temporary.unlink(missing_ok=True)
    return final


def finalize_success(
    *, patient_id: str, input_manifest: str | Path, source: str | Path,
    config: str | Path, instrumentation: str | Path, aggregate: str | Path,
    pvac_input_vcf: str | Path, pvac_input_vcf_index: str | Path,
    reference: str | Path, execution_log: str | Path, execution_trace: str | Path,
    output_dir: str | Path, execution_status: str = SUCCESS, exit_code: int = 0,
    execution_command: str, started_at: str, finished_at: str, runtime: str,
    normalizer: Normalizer | None = None,
) -> dict:
    """Finalize a successful attempt and freeze its canonical mutation top 20."""
    if execution_status != SUCCESS or exit_code != 0:
        raise ValueError("a successful Track-A attempt requires status SUCCEEDED and exit_code 0")
    lock = _lock_common(
        patient_id=patient_id, input_manifest=input_manifest, source=source, config=config,
        instrumentation=instrumentation, execution_log=execution_log,
        execution_trace=execution_trace,
    )
    execution_identity = _execution_identity(
        lock, command=execution_command, started_at=started_at,
        finished_at=finished_at, runtime=runtime,
    )
    aggregate_path = _file(aggregate)
    vcf_path = _file(pvac_input_vcf)
    vcf_index_path = _file(pvac_input_vcf_index)
    reference_path = _file(reference)
    destination = _empty_destination(output_dir)

    _columns, rows = _read_aggregate(aggregate_path)
    headers, vcf_by_id = _read_vcf_alleles(vcf_path)
    mapped: list[tuple[str, VcfAllele]] = []
    for rank, row in enumerate(rows, 1):
        matches = vcf_by_id.get(row["ID"], [])
        if len(matches) != 1:
            raise ValueError(
                f"aggregate ID must map one-to-one to exact pVAC input VCF allele: "
                f"{row['ID']} ({len(matches)} matches)"
            )
        mapped.append((f"PVACMAP{rank:08d}", matches[0]))

    normalize = normalizer or _bcftools_normalizer
    all_vcf, all_norm = _normalize_set(
        name="all", destination=destination, headers=headers, mapped=mapped,
        reference=reference_path, normalizer=normalize,
    )
    selected_mapped = mapped[:20]
    top_vcf, top_norm = _normalize_set(
        name="top20", destination=destination, headers=headers, mapped=selected_mapped,
        reference=reference_path, normalizer=normalize,
    )
    # Independently running the selected subset must not alter its canonical representation.
    for token, _allele in selected_mapped:
        if top_norm[token] != all_norm[token]:
            raise ValueError(f"selected/all normalization disagrees for {token}")

    mapping_path = destination / "NORMALIZATION_MAP.tsv"
    with mapping_path.open("w", newline="") as handle:
        fields = ("aggregate_rank", "aggregate_id", "mapping_token", "input_variant",
                  "canonical_variant", "selected_top20")
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for rank, ((token, allele), row) in enumerate(zip(mapped, rows), 1):
            writer.writerow({
                "aggregate_rank": rank, "aggregate_id": row["ID"], "mapping_token": token,
                "input_variant": allele.input_key, "canonical_variant": all_norm[token],
                "selected_top20": str(rank <= 20).lower(),
            })

    portfolio_path = destination / "NEXTNEOPI_MUTATION_TOP20.tsv"
    with portfolio_path.open("w", newline="") as handle:
        fields = ("rank", "canonical_mutation", "aggregate_id", "peptide", "hla", "tier")
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for rank, ((token, _allele), row) in enumerate(zip(selected_mapped, rows[:20]), 1):
            writer.writerow({
                "rank": rank, "canonical_mutation": top_norm[token],
                "aggregate_id": row["ID"], "peptide": row["Best Peptide"],
                "hla": row["Allele"], "tier": row["Tier"],
            })

    result = _common_manifest(lock, SUCCESS, 0)
    result["execution"].update(execution_identity)
    result.update({
        "attempt_outcome": "PORTFOLIO_FROZEN",
        "native_aggregate": {**_attestation(aggregate_path), "rows": len(rows)},
        "published_pvac_input": {
            "vcf": _attestation(vcf_path), "index": _attestation(vcf_index_path),
        },
        "grch38_reference": _attestation(reference_path),
        "normalization": {
            "command": "bcftools norm -f <exact-reference> -m-any",
            "all_variants_vcf": _attestation(all_vcf),
            "selected_variants_vcf": _attestation(top_vcf),
            "map": _attestation(mapping_path),
            "mapped_aggregate_ids": len(mapped),
        },
        "portfolio": {
            **_attestation(portfolio_path), "unit": "canonical mutation",
            "size": len(selected_mapped), "ordering": "native aggregate order",
        },
    })
    _write_once_manifest(destination, result)
    return result


def finalize_non_success(
    *, patient_id: str, input_manifest: str | Path, source: str | Path,
    config: str | Path, instrumentation: str | Path, execution_log: str | Path,
    execution_trace: str | Path, output_dir: str | Path, execution_status: str,
    exit_code: int | None, reason: str,
    execution_command: str, started_at: str, finished_at: str, runtime: str,
) -> dict:
    """Record a verified failure/abstention as an evaluable zero-portfolio attempt."""
    if execution_status not in NON_SUCCESS:
        raise ValueError("non-success status must be FAILED or ABSTAINED")
    if execution_status == "FAILED" and (exit_code is None or exit_code == 0):
        raise ValueError("FAILED attempt requires a non-zero exit code")
    if execution_status == "ABSTAINED" and exit_code == 0:
        raise ValueError("ABSTAINED attempt cannot carry successful exit code 0")
    if not reason or not reason.strip():
        raise ValueError("failure/abstention reason is required")
    lock = _lock_common(
        patient_id=patient_id, input_manifest=input_manifest, source=source, config=config,
        instrumentation=instrumentation, execution_log=execution_log,
        execution_trace=execution_trace,
    )
    execution_identity = _execution_identity(
        lock, command=execution_command, started_at=started_at,
        finished_at=finished_at, runtime=runtime,
    )
    destination = _empty_destination(output_dir)
    result = _common_manifest(lock, execution_status, exit_code)
    result["execution"].update(execution_identity)
    result.update({
        "attempt_outcome": execution_status,
        "reason": reason.strip(),
        "native_aggregate": None,
        "published_pvac_input": None,
        "grch38_reference": None,
        "normalization": None,
        "portfolio": {"size": 0, "unit": "canonical mutation", "rows": []},
    })
    _write_once_manifest(destination, result)
    return result
