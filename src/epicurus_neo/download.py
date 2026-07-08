from __future__ import annotations

import shutil
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from epicurus_neo.data_manifest import load_dataset_manifest


@dataclass(frozen=True)
class DownloadPlan:
    dataset_key: str
    file_key: str
    url: str
    output_path: Path
    notes: str


def dataset_file_plans(
    manifest_path: str | Path,
    *,
    output_dir: str | Path,
    dataset_key: str | None = None,
) -> list[DownloadPlan]:
    output_root = Path(output_dir)
    payload = load_dataset_manifest(manifest_path)
    plans: list[DownloadPlan] = []
    for source in payload:
        if dataset_key is not None and source.key != dataset_key:
            continue
        raw_files = getattr(source, "files", {})
        for file_key, item in raw_files.items():
            plans.append(
                DownloadPlan(
                    dataset_key=source.key,
                    file_key=file_key,
                    url=str(item["url"]),
                    output_path=output_root / source.key / str(item["filename"]),
                    notes=str(item.get("notes", "")),
                )
            )
    return plans


def download_file(url: str, output_path: str | Path, *, overwrite: bool = False) -> Path:
    out = Path(output_path)
    if out.exists() and not overwrite:
        return out
    out.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response, out.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    return out

