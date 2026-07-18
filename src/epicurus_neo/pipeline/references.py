"""Reference-bundle manifest and scaffold.

The pipeline needs a one-time GRCh38 reference bundle (tens-hundreds of GB). Rather
than ship an unverified multi-hundred-GB downloader, Epicurus documents exactly what
the bundle must contain and where each item comes from, and scaffolds the directory
with a REFERENCES.md the user follows on their machine.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class ReferenceItem:
    name: str
    used_by: str
    source: str
    note: str


REFERENCE_MANIFEST: tuple[ReferenceItem, ...] = (
    ReferenceItem(
        "genome.fa",
        "align, call, annotate",
        "GATK GRCh38 resource bundle (Homo_sapiens_assembly38.fasta) or Ensembl GRCh38 primary assembly",
        "Also needs the BWA-MEM2 index (bwa-mem2 index genome.fa) and .fai / .dict.",
    ),
    ReferenceItem(
        "gatk/",
        "call",
        "GATK GRCh38 resource bundle: germline resource (af-only-gnomad) + panel of normals",
        "Used by Mutect2 / FilterMutectCalls.",
    ),
    ReferenceItem(
        "vep/",
        "annotate",
        "Ensembl VEP cache for GRCh38 (vep_install --cache) matching your VEP version",
        "Offline cache directory.",
    ),
    ReferenceItem(
        "salmon_index/",
        "express",
        "Salmon index built from GENCODE GRCh38 transcripts (salmon index)",
        "Match the transcript annotation used downstream.",
    ),
)


def references_manifest() -> list[dict]:
    return [
        {"name": item.name, "used_by": item.used_by, "source": item.source, "note": item.note}
        for item in REFERENCE_MANIFEST
    ]


def _references_markdown() -> str:
    lines = [
        "# Epicurus reference bundle (GRCh38)",
        "",
        "The pipeline expects the following items inside this directory. Each is a one-time",
        "download; sizes range from ~3 GB (genome) to tens of GB (VEP cache, Salmon index).",
        "",
    ]
    for item in REFERENCE_MANIFEST:
        lines += [
            f"## `{item.name}`",
            f"- **Used by stage(s):** {item.used_by}",
            f"- **Source:** {item.source}",
            f"- **Note:** {item.note}",
            "",
        ]
    lines += [
        "After populating this directory, verify readiness with:",
        "",
        "```bash",
        "epicurus doctor --bundle-dir <this directory>",
        "```",
        "",
    ]
    return "\n".join(lines)


def scaffold_references(dest: str | Path) -> dict:
    """Create the bundle directory and write REFERENCES.md; return the manifest.

    Does not download anything: it records exactly what to install and where the
    pipeline looks for it.
    """
    base = Path(dest).expanduser()
    base.mkdir(parents=True, exist_ok=True)
    readme = base / "REFERENCES.md"
    readme.write_text(_references_markdown())
    return {
        "bundle_dir": str(base),
        "instructions": str(readme),
        "items": references_manifest(),
    }
