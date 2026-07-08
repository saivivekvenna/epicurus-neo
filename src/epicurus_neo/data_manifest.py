from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class DatasetSource:
    key: str
    name: str
    role: str
    locked_test: bool
    url: str
    direct_download: str | None
    notes: str
    files: dict


def load_dataset_manifest(path: str | Path) -> list[DatasetSource]:
    manifest_path = Path(path)
    payload = yaml.safe_load(manifest_path.read_text()) or {}
    raw_datasets = payload.get("datasets", {})
    sources: list[DatasetSource] = []
    for key, item in sorted(raw_datasets.items()):
        sources.append(
            DatasetSource(
                key=str(key),
                name=str(item["name"]),
                role=str(item["role"]),
                locked_test=bool(item["locked_test"]),
                url=str(item["url"]),
                direct_download=item.get("direct_download"),
                notes=str(item.get("notes", "")).strip(),
                files=dict(item.get("files", {})),
            )
        )
    return sources


def downloadable_sources(sources: list[DatasetSource]) -> list[DatasetSource]:
    return [source for source in sources if source.direct_download]
