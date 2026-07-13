"""Incumbent-score provenance with a fail-closed guard against mislabeling a proxy as PRIME.

The Zhao supplement ships MixMHCpred-3.0 (the binding backbone that open-source PRIME is
built on) but NOT PRIME's final immunogenicity score, and no PRIME 2.0 executable / training
table is available locally. Any report that calls the incumbent "PRIME" without a genuine
PRIME score is a correctness bug; this module makes that impossible to do silently.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class IncumbentProvenance(str, Enum):
    REAL_PRIME = "REAL_PRIME"
    MIXMHCPRED_PROXY = "MIXMHCPRED_PROXY"


@dataclass(frozen=True)
class IncumbentSpec:
    name: str  # short arm/column name (must NOT claim PRIME unless provenance is REAL_PRIME)
    column: str  # the score column in the feature frame
    provenance: IncumbentProvenance
    description: str
    higher_is_better: bool = True

    def __post_init__(self) -> None:
        claims_prime = "prime" in self.name.lower() or "prime" in self.column.lower()
        if claims_prime and self.provenance is not IncumbentProvenance.REAL_PRIME:
            raise ValueError(
                f"Incumbent {self.name!r}/{self.column!r} names PRIME but provenance is "
                f"{self.provenance.value}; a proxy may never be labeled PRIME."
            )

    @property
    def is_real_prime(self) -> bool:
        return self.provenance is IncumbentProvenance.REAL_PRIME


# The only incumbent currently runnable on Zhao. Explicitly a proxy — never labeled PRIME.
MIXMHCPRED_PROXY = IncumbentSpec(
    name="incumbent_proxy",
    column="mixmhcpred3_score",
    provenance=IncumbentProvenance.MIXMHCPRED_PROXY,
    description=(
        "MixMHCpred-3.0 binding score (PRIME's binding backbone). PROXY incumbent — NOT PRIME. "
        "True PRIME 2.0 requires the GfellerLab executable (+MixMHCpred), absent locally; a "
        "reproducible acquisition adapter + input artifact is emitted instead."
    ),
    higher_is_better=True,
)


def require_real_prime(spec: IncumbentSpec, *, context: str) -> None:
    """Fail closed where a result is only valid with a genuine PRIME score."""
    if not spec.is_real_prime:
        raise RuntimeError(
            f"{context}: requires a REAL_PRIME incumbent, but got provenance "
            f"{spec.provenance.value} ({spec.description}). This arm is BLOCKED until PRIME is "
            "acquired; do not report it as a PRIME result."
        )


def assert_report_labeling(report: dict) -> None:
    """Guard: any report field mentioning 'prime' as a result must carry REAL_PRIME provenance."""
    provenance = report.get("incumbent_provenance")
    for key in report:
        key_lower = key.lower()
        is_blocked_or_acquisition_metadata = any(
            marker in key_lower for marker in ("blocked", "blocker", "acquisition")
        )
        if "prime" in key_lower and not is_blocked_or_acquisition_metadata:
            if provenance != IncumbentProvenance.REAL_PRIME.value:
                raise AssertionError(
                    f"report key {key!r} implies PRIME but incumbent_provenance={provenance!r}"
                )
