"""External-tool specifications and availability checks.

Epicurus orchestrates established, validated bioinformatics tools; it does not
reimplement them. Each tool is described by a :class:`ToolSpec` so the pipeline
can (a) check the tool is installed before a long run and (b) record which tool
and version produced each artifact.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass


@dataclass(frozen=True)
class ToolSpec:
    """A single external command-line tool the pipeline depends on."""

    name: str
    binary: str
    version_args: tuple[str, ...] = ("--version",)
    purpose: str = ""

    def is_available(self) -> bool:
        """True if the tool's binary is on PATH."""
        return shutil.which(self.binary) is not None

    def resolved_path(self) -> str | None:
        return shutil.which(self.binary)

    def version(self) -> str | None:
        """Best-effort tool version string; ``None`` if not installed or it errors.

        Never raises: version probing must not be able to crash a run or a
        readiness check.
        """
        path = self.resolved_path()
        if path is None:
            return None
        try:
            completed = subprocess.run(  # noqa: S603 - fixed argv, no shell
                [path, *self.version_args],
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.SubprocessError):
            return None
        output = (completed.stdout or completed.stderr or "").strip()
        return output.splitlines()[0].strip() if output else None


# Canonical tool registry. Binaries are the conventional executable names; a user
# gets them from the provided container or bioconda environment.
BWA_MEM2 = ToolSpec("bwa-mem2", "bwa-mem2", ("version",), "read alignment")
SAMTOOLS = ToolSpec("samtools", "samtools", ("--version",), "BAM utilities")
GATK = ToolSpec("gatk", "gatk", ("--version",), "somatic calling (Mutect2)")
VEP = ToolSpec("vep", "vep", ("--help",), "variant annotation")
SALMON = ToolSpec("salmon", "salmon", ("--version",), "RNA transcript quantification")
OPTITYPE = ToolSpec("optitype", "OptiTypePipeline.py", ("--help",), "class-I HLA typing")
PVACSEQ = ToolSpec("pvacseq", "pvacseq", ("--version",), "neoantigen candidate generation")
MHCFLURRY = ToolSpec("mhcflurry", "mhcflurry-predict", ("--version",), "MHC presentation prediction")

ALL_TOOLS: tuple[ToolSpec, ...] = (
    BWA_MEM2,
    SAMTOOLS,
    GATK,
    VEP,
    SALMON,
    OPTITYPE,
    PVACSEQ,
    MHCFLURRY,
)
