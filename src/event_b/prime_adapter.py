"""REAL_PRIME scoring adapter — genuine GfellerLab PRIME 2.1 (+ MixMHCpred 3.0).

This runs the actual PRIME executable installed under data/raw/tools (gitignored, non-commercial
academic license; see configs/source_manifests/prime_tool.yml). It is the ONLY source of PRIME
scores in this repo. A MixMHCpred/NetMHCpan-EL baseline is NOT PRIME and must never be labeled so.

Per-(peptide, restricting-HLA) scoring: PRIME is run once per restricting allele, and that allele's
`%Rank_<ALLELE>` (PRIME %rank, lower = better), `Score_<ALLELE>` (PRIME score), and
`%RankBinding_<ALLELE>` (MixMHCpred %rank) are read back. Peptides outside 8-14mers / non-standard
residues, and alleles PRIME cannot score, are marked (never fabricated).
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
import subprocess
import tempfile

import pandas as pd

from event_b.incumbent import IncumbentProvenance, IncumbentSpec

TOOLS = Path("data/raw/tools")
PRIME_DIR = TOOLS / "PRIME"
PRIME_EXEC = PRIME_DIR / "PRIME"
MIX_EXEC = TOOLS / "MixMHCpred" / "MixMHCpred"
TOOLENV_BIN = TOOLS / "toolenv" / "bin"

PRIME_COMMIT = "7b18d4e11042141e7102f7c69be2b0e03d138dab"
MIX_COMMIT = "0a7f9b9e20d1cf02236f4a0a90d16735be879b38"

STD_AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LEN, MAX_LEN = 8, 14
# PRIME.x scales super-linearly (memory) on a single huge peptide file: ~80k peptides score in ~1s,
# but ~300k in one call thrashes for >20 min. Chunk large per-allele batches to stay in the fast
# regime — results are identical (PRIME %rank is per (peptide, allele), independent of batching).
PRIME_MAX_PEPTIDES_PER_CALL = 40000


REAL_PRIME_INCUMBENT = IncumbentSpec(
    name="prime_incumbent",
    column="prime_rank",
    provenance=IncumbentProvenance.REAL_PRIME,
    description=(
        f"Genuine GfellerLab PRIME 2.1 %Rank (commit {PRIME_COMMIT[:10]}, MixMHCpred 3.0 "
        f"{MIX_COMMIT[:10]}); lower %rank = better. Non-commercial academic license."
    ),
    higher_is_better=False,  # PRIME %rank: lower is better
)


def prime_available() -> bool:
    return all(p.exists() for p in (PRIME_EXEC, MIX_EXEC, TOOLENV_BIN / "python3"))


def normalize_allele_for_prime(hla: object) -> str:
    """'HLA-A*02:01' / 'HLA-A02:01' / 'A*02:01' -> 'A0201'; return '' if not a class-I HLA-A/B/C."""
    text = str(hla).strip().upper().replace("HLA-", "").replace("HLA", "").replace(" ", "")
    text = text.replace("*", "").replace(":", "")
    m = re.match(r"^([ABC])(\d{4})$", text)
    return f"{m.group(1)}{m.group(2)}" if m else ""


def _usable_peptide(pep: object) -> bool:
    p = str(pep).strip().upper()
    return MIN_LEN <= len(p) <= MAX_LEN and set(p).issubset(STD_AA)


def _run_prime_for_allele(peptides: list[str], allele: str, workdir: Path) -> pd.DataFrame:
    """Run PRIME for one allele; return per-peptide prime_rank/score/mixmhcpred_rank."""
    pep_file = workdir / f"pep_{allele}.txt"
    out_file = workdir / f"out_{allele}.txt"
    pep_file.write_text("\n".join(peptides) + "\n")
    env = {"PATH": f"{TOOLENV_BIN.resolve()}:{MIX_EXEC.parent.resolve()}:/usr/bin:/bin"}
    proc = subprocess.run(
        [str(PRIME_EXEC.resolve()), "-i", str(pep_file.resolve()), "-o", str(out_file.resolve()),
         "-a", allele, "-mix", str(MIX_EXEC.resolve())],
        capture_output=True, text=True, env=env, cwd=str(PRIME_DIR.resolve()),
    )
    if not out_file.exists() or f"cannot be run for {allele}" in (proc.stdout + proc.stderr):
        return pd.DataFrame(columns=["peptide", "prime_rank", "prime_score", "mixmhcpred_rank"])
    rows = [ln for ln in out_file.read_text().splitlines() if ln and not ln.startswith("#")]
    if len(rows) < 2:
        return pd.DataFrame(columns=["peptide", "prime_rank", "prime_score", "mixmhcpred_rank"])
    frame = pd.read_csv(out_file, sep="\t", comment="#")
    rank_col, score_col, bind_col = f"%Rank_{allele}", f"Score_{allele}", f"%RankBinding_{allele}"
    if rank_col not in frame.columns:
        return pd.DataFrame(columns=["peptide", "prime_rank", "prime_score", "mixmhcpred_rank"])
    return pd.DataFrame({
        "peptide": frame["Peptide"].astype(str),
        "prime_rank": pd.to_numeric(frame[rank_col], errors="coerce"),
        "prime_score": pd.to_numeric(frame[score_col], errors="coerce"),
        "mixmhcpred_rank": pd.to_numeric(frame[bind_col], errors="coerce"),
    })


@dataclass(frozen=True)
class PrimeScoringResult:
    scored: pd.DataFrame       # one row per input (peptide, allele) with prime_rank etc + status
    provenance: dict


def score_prime(pairs: pd.DataFrame, *, peptide_col="peptide", hla_col="hla_allele") -> PrimeScoringResult:
    """Score genuine PRIME for each (peptide, restricting HLA). Never fabricates a missing score."""
    if not prime_available():
        raise RuntimeError(
            "PRIME is not installed at data/raw/tools. Run the acquisition in "
            "configs/source_manifests/prime_tool.yml first."
        )
    work = pairs[[peptide_col, hla_col]].drop_duplicates().copy()
    work["peptide"] = work[peptide_col].astype(str).str.upper()
    work["prime_allele"] = work[hla_col].map(normalize_allele_for_prime)
    work["usable"] = work["peptide"].map(_usable_peptide) & work["prime_allele"].ne("")

    results = []
    with tempfile.TemporaryDirectory() as tmp:
        workdir = Path(tmp)
        for allele, group in work[work["usable"]].groupby("prime_allele"):
            peptides = sorted(group["peptide"].unique())
            for start in range(0, len(peptides), PRIME_MAX_PEPTIDES_PER_CALL):
                chunk = peptides[start:start + PRIME_MAX_PEPTIDES_PER_CALL]
                scored = _run_prime_for_allele(chunk, allele, workdir)
                scored["prime_allele"] = allele
                results.append(scored)
    scored_all = pd.concat(results, ignore_index=True) if results else pd.DataFrame(
        columns=["peptide", "prime_rank", "prime_score", "mixmhcpred_rank", "prime_allele"]
    )

    merged = work.merge(scored_all, on=["peptide", "prime_allele"], how="left")
    merged["prime_status"] = "SCORED"
    merged.loc[~merged["usable"], "prime_status"] = "SKIPPED_UNSCORABLE_INPUT"
    merged.loc[merged["usable"] & merged["prime_rank"].isna(), "prime_status"] = "PRIME_ALLELE_UNSUPPORTED"

    out = pairs.merge(
        merged[[peptide_col, hla_col, "prime_rank", "prime_score", "mixmhcpred_rank", "prime_status"]]
        .rename(columns={peptide_col: peptide_col, hla_col: hla_col}),
        on=[peptide_col, hla_col], how="left",
    )
    provenance = {
        "incumbent_provenance": IncumbentProvenance.REAL_PRIME.value,
        "prime_version": "2.1",
        "prime_commit": PRIME_COMMIT,
        "mixmhcpred_version": "3.0",
        "mixmhcpred_commit": MIX_COMMIT,
        "license": "GfellerLab non-commercial academic; not redistributed.",
        "n_pairs": int(len(out)),
        "n_scored": int((out["prime_status"] == "SCORED").sum()),
        "status_counts": out["prime_status"].value_counts().to_dict(),
    }
    return PrimeScoringResult(out, provenance)
