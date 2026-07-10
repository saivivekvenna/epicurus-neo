"""Checksum-pinned, format-aware source validation without whole-file loading."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
from zipfile import BadZipFile, ZipFile

from event_b.manifest import sha256_file


SUPPORTED_MEDIA_TYPES = {"csv", "tsv", "json", "jsonl", "docx", "pdf", "xlsx"}


@dataclass(frozen=True)
class SourceValidation:
    path: str
    media_type: str
    file_size_bytes: int
    checksum_sha256: str
    valid: bool
    reason: str | None = None


def _looks_like_html(prefix: bytes) -> bool:
    text = prefix.lstrip().lower()
    return text.startswith((b"<!doctype html", b"<html", b"<head", b"<body"))


def _zip_media_type(path: Path) -> str:
    with ZipFile(path) as archive:
        names = set(archive.namelist())
    if "[Content_Types].xml" not in names:
        raise ValueError("ZIP container is not an OOXML document")
    if any(name.startswith("xl/") for name in names):
        return "xlsx"
    if any(name.startswith("word/") for name in names):
        return "docx"
    raise ValueError("OOXML container is neither XLSX nor DOCX")


def detect_media_type(path: str | Path) -> str:
    source = Path(path)
    with source.open("rb") as handle:
        prefix = handle.read(4096)
    if _looks_like_html(prefix):
        raise ValueError("HTML response detected; refusing disguised supplement")
    if prefix.startswith(b"%PDF-"):
        return "pdf"
    if prefix.startswith(b"PK\x03\x04"):
        try:
            return _zip_media_type(source)
        except BadZipFile as error:
            raise ValueError("invalid ZIP/OOXML container") from error
    suffix = source.suffix.lower().lstrip(".")
    if suffix in {"json", "jsonl"}:
        with source.open(encoding="utf-8-sig") as handle:
            if suffix == "json":
                json.load(handle)
            else:
                for line in handle:
                    if line.strip():
                        json.loads(line)
        return suffix
    if suffix in {"csv", "tsv"}:
        prefix.decode("utf-8-sig")
        return suffix
    raise ValueError(f"unsupported or unrecognized source type: {suffix or '<none>'}")


def validate_source_file(
    path: str | Path, *, expected_media_type: str | None = None
) -> SourceValidation:
    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    media_type = detect_media_type(source)
    if expected_media_type and media_type != expected_media_type.lower().lstrip("."):
        raise ValueError(f"expected {expected_media_type}, detected {media_type}")
    return SourceValidation(
        str(source.resolve()),
        media_type,
        source.stat().st_size,
        sha256_file(source),
        True,
    )
