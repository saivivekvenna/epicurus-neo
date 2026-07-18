"""Readiness check: are the external tools and reference data installed?"""

from __future__ import annotations

from pathlib import Path

from epicurus_neo.pipeline.tools import ALL_TOOLS

# Reference artifacts the pipeline expects inside the bundle directory.
REFERENCE_ITEMS: tuple[str, ...] = ("genome.fa", "vep", "salmon_index")


def readiness_report(bundle_dir: str | Path | None = None) -> dict:
    """Report which tools and reference items are present.

    ``ready`` is True only when every tool is on PATH and, when a bundle
    directory is given, every expected reference item exists.
    """
    tools = [
        {
            "name": tool.name,
            "binary": tool.binary,
            "purpose": tool.purpose,
            "available": tool.is_available(),
            "path": tool.resolved_path(),
        }
        for tool in ALL_TOOLS
    ]
    tools_ready = all(entry["available"] for entry in tools)

    references: list[dict] = []
    references_ready = True
    if bundle_dir is not None:
        base = Path(bundle_dir).expanduser()
        for item in REFERENCE_ITEMS:
            present = (base / item).exists()
            references_ready = references_ready and present
            references.append({"item": item, "path": str(base / item), "present": present})

    return {
        "ready": tools_ready and references_ready,
        "tools_ready": tools_ready,
        "references_ready": references_ready if bundle_dir is not None else None,
        "tools": tools,
        "references": references,
        "bundle_dir": str(bundle_dir) if bundle_dir is not None else None,
    }
