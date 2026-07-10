"""Source manifests and content-addressed raw-document provenance."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from pathlib import Path

from event_b.models import SCHEMA_VERSION


def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


@dataclass(frozen=True)
class SourceDocument:
    document_id: str
    local_path: str | None
    source_url: str | None
    checksum_sha256: str
    media_type: str


@dataclass(frozen=True)
class SourceManifest:
    manifest_id: str
    source_name: str
    source_version: str
    adapter_name: str
    adapter_version: str
    documents: tuple[SourceDocument, ...]
    schema_version: str = SCHEMA_VERSION

    def to_dict(self) -> dict:
        return asdict(self)

    def write(self, path: str | Path) -> None:
        Path(path).parent.mkdir(parents=True, exist_ok=True)
        Path(path).write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n")


def manifest_from_paths(
    source_name: str,
    source_version: str,
    adapter_name: str,
    adapter_version: str,
    paths: list[str | Path],
) -> SourceManifest:
    documents = []
    for raw_path in sorted(map(Path, paths), key=lambda path: str(path)):
        checksum = sha256_file(raw_path)
        documents.append(
            SourceDocument(
                document_id="doc:" + checksum[:20],
                local_path=str(raw_path.resolve()),
                source_url=None,
                checksum_sha256=checksum,
                media_type=raw_path.suffix.lower().lstrip(".") or "binary",
            )
        )
    identity = source_name + source_version + "".join(doc.checksum_sha256 for doc in documents)
    return SourceManifest(
        manifest_id="manifest:" + sha256(identity.encode()).hexdigest()[:20],
        source_name=source_name,
        source_version=source_version,
        adapter_name=adapter_name,
        adapter_version=adapter_version,
        documents=tuple(documents),
    )
