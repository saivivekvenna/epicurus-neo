"""Bootstrap and paired-power statistics over patient-level vectors."""

from __future__ import annotations

from dataclasses import dataclass
from statistics import NormalDist

import numpy as np


@dataclass(frozen=True)
class BootstrapCI:
    mean: float
    lo: float
    hi: float
    n: int

    def __iter__(self):
        # Preserve the build-spec's three-value unpacking while exposing retained n.
        yield self.mean
        yield self.lo
        yield self.hi


@dataclass(frozen=True)
class PairedBootstrapCI:
    delta: float
    lo: float
    hi: float
    p_better: float
    n: int

    def __iter__(self):
        # Preserve the build-spec's four-value unpacking while exposing retained n.
        yield self.delta
        yield self.lo
        yield self.hi
        yield self.p_better


class MDE(float):
    """A float carrying the retained paired sample size."""

    n: int

    def __new__(cls, value: float, n: int):
        instance = float.__new__(cls, value)
        instance.n = n
        return instance


def _clean(v: np.ndarray) -> np.ndarray:
    values = np.asarray(v, dtype=float).reshape(-1)
    return values[np.isfinite(values)]


def _paired(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    left = np.asarray(a, dtype=float).reshape(-1)
    right = np.asarray(b, dtype=float).reshape(-1)
    if left.shape != right.shape:
        raise ValueError("paired vectors must have the same shape")
    keep = np.isfinite(left) & np.isfinite(right)
    return left[keep], right[keep]


def bootstrap_ci(v: np.ndarray, n: int = 20_000, seed: int = 0) -> BootstrapCI:
    values = _clean(v)
    if not len(values):
        raise ValueError("bootstrap requires at least one finite value")
    if n <= 0:
        raise ValueError("n must be positive")
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(values), size=(n, len(values)))
    means = values[indices].mean(axis=1)
    lo, hi = np.quantile(means, [0.025, 0.975])
    return BootstrapCI(float(values.mean()), float(lo), float(hi), len(values))


def paired_bootstrap(
    a: np.ndarray,
    b: np.ndarray,
    n: int = 20_000,
    seed: int = 0,
) -> PairedBootstrapCI:
    left, right = _paired(a, b)
    if not len(left):
        raise ValueError("paired bootstrap requires at least one finite pair")
    if n <= 0:
        raise ValueError("n must be positive")
    differences = left - right
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(differences), size=(n, len(differences)))
    means = differences[indices].mean(axis=1)
    lo, hi = np.quantile(means, [0.025, 0.975])
    return PairedBootstrapCI(
        delta=float(differences.mean()),
        lo=float(lo),
        hi=float(hi),
        p_better=float(np.mean(means > 0.0)),
        n=len(differences),
    )


def mde(
    a: np.ndarray,
    b: np.ndarray,
    power: float = 0.80,
    alpha: float = 0.05,
) -> MDE:
    """Return the plan's current-n CI-exclusion threshold using paired SD.

    The plan's registered 0.237 value is z_(1-alpha/2) * SE, not the
    conventional z_alpha + z_power prospective MDE.  ``power`` is retained in
    the public signature for compatibility; prospective planning belongs in
    :func:`n_required` and uses it there.
    """
    del power
    left, right = _paired(a, b)
    if len(left) < 2:
        raise ValueError("mde requires at least two finite pairs")
    if not 0 < alpha < 1:
        raise ValueError("alpha must be between zero and one")
    sd_diff = float(np.std(left - right, ddof=1))
    z_alpha = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    return MDE(z_alpha * sd_diff / np.sqrt(len(left)), len(left))


def n_required(
    sd_diff: float,
    delta: float,
    power: float = 0.80,
    alpha: float = 0.05,
) -> int:
    """Prospective paired sample size using z_alpha + z_power."""
    if sd_diff < 0 or delta <= 0:
        raise ValueError("sd_diff must be non-negative and delta must be positive")
    if not 0 < power < 1 or not 0 < alpha < 1:
        raise ValueError("power and alpha must be between zero and one")
    z_alpha = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    z_power = NormalDist().inv_cdf(power)
    return int(np.ceil(((z_alpha + z_power) * sd_diff / delta) ** 2))
