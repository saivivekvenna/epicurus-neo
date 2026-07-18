"""Provenance recording: hash inputs and capture which tool produced each artifact."""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any

from epicurus_neo.pipeline.tools import ToolSpec

_CHUNK = 1024 * 1024


def hash_file(path: str | Path) -> str:
    """SHA-256 of a file's bytes; ``"absent"`` if the file does not exist."""
    p = Path(path)
    if not p.exists() or not p.is_file():
        return "absent"
    digest = hashlib.sha256()
    with p.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_CHUNK), b""):
            digest.update(chunk)
    return digest.hexdigest()


def stage_provenance(
    *,
    stage: str,
    tool: ToolSpec | None,
    command: list[str],
    inputs: dict[str, str | Path],
    outputs: dict[str, str | Path],
) -> dict[str, Any]:
    """Build a JSON-serializable provenance record for one stage execution.

    No wall-clock timestamp is embedded so records are reproducible and diffable;
    the calling runner may add one if desired.
    """
    return {
        "stage": stage,
        "tool": None if tool is None else tool.name,
        "tool_binary": None if tool is None else tool.binary,
        "tool_version": None if tool is None else tool.version(),
        "command": list(command),
        "inputs": {name: hash_file(path) for name, path in inputs.items()},
        "outputs": {name: str(path) for name, path in outputs.items()},
    }
