"""Miller IPV (PRJNA980652, DOI 10.1126/scitranslmed.abj9905) ingestion scaffold.

State of the cohort (2026-07-12): the raw WES/RNA inputs are PUBLIC on SRA (PRJNA980652); the per-peptide
ELISpot label table (data files S1/S2) is paywalled behind science.org (Cloudflare 403 to every public
CLI/API route — Crossref/Unpaywall/PMC/Europe PMC/direct file patterns all fail; the article is not open
access). So this module covers everything that does NOT need the paywalled file:

  * parse the public SRA run metadata -> a per-patient input crosswalk (13 patients x {normal exome,
    tumor exome, tumor RNA}) and ordered download tranches (smallest valid first = one patient trio); and
  * the ingestion CONTRACT the S1/S2 label frame must satisfy once obtained, enforcing the north-star
    invariants: three-state labels (POSITIVE / TESTED_NEGATIVE / UNTESTED), NO collapsing of contradictory
    longitudinal rows, a real conflict surfaced (never dropped), and a deterministic patient crosswalk to
    the SRA BioSamples.

Aggregate labels are known from the paper text (RUNNABLE BUT BLOCKED ON FILE for the per-row table):
754 assayed 20-mers, 199 POSITIVE (26%), 555 TESTED_NEGATIVE, across 13 patients / 349 tested variants.
"""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

_ROOT = Path(__file__).resolve().parents[2]
SRA_RUNINFO_FIXTURE = _ROOT / "artifacts" / "milestone_7_decision" / "external_validation" / "miller_ipv" / "SRA_RUNINFO.csv"

VALID_LABELS = {"POSITIVE", "TESTED_NEGATIVE", "UNTESTED"}
_STD_AA = set("ACDEFGHIKLMNPQRSTVWY")
# one row per this key; a repeated key with a DIFFERENT label is a conflict to surface, never collapse.
LABEL_KEY = ("patient_id", "mutant_peptide", "hla_allele", "assay", "assay_timepoint")


# ---------------------------------------------------------------------------
# Public SRA metadata -> input crosswalk + download tranches
# ---------------------------------------------------------------------------
def _patient_from_sample(name: object) -> str:
    m = re.match(r"(Hu_\d+)", str(name))
    return m.group(1) if m else ""


def parse_sra_runinfo(path: str | Path) -> pd.DataFrame:
    """Normalize the public SRA runinfo CSV into per-run rows with patient id + tumor/normal + exome/RNA."""
    raw = pd.read_csv(path)
    sample = raw["SampleName"].astype(str)
    libname = raw.get("LibraryName", pd.Series("", index=raw.index)).astype(str)
    out = pd.DataFrame({
        "run": raw["Run"].astype(str),
        "biosample": raw["BioSample"].astype(str),
        "sample_name": sample,
        "patient_id": sample.map(_patient_from_sample),
        "tissue": ["tumor" if "tumor" in s.lower() else "normal" if "normal" in s.lower() else ""
                   for s in sample],
        "library_strategy": raw["LibraryStrategy"].astype(str),
        "is_rna": libname.str.contains("RNA", case=False) | raw["LibraryStrategy"].str.contains("RNA", case=False),
        "bases": pd.to_numeric(raw["bases"], errors="coerce"),
        "size_mb": pd.to_numeric(raw["size_MB"], errors="coerce"),
    })
    out["size_gb"] = out["size_mb"] / 1000.0
    out["assay_kind"] = ["tumor_rna" if r.is_rna else f"{r.tissue}_exome" for r in out.itertuples()]
    return out


def patient_input_crosswalk(runs: pd.DataFrame) -> pd.DataFrame:
    """Per-patient input completeness + size (needs normal exome + tumor exome + tumor RNA)."""
    rows = []
    for pid, g in runs.groupby("patient_id", sort=True):
        kinds = set(g["assay_kind"])
        has_ne = "normal_exome" in kinds
        has_te = "tumor_exome" in kinds
        has_rna = "tumor_rna" in kinds
        rows.append({
            "patient_id": pid,
            "has_normal_exome": has_ne,
            "has_tumor_exome": has_te,
            "has_tumor_rna": has_rna,
            "complete": has_ne and has_te and has_rna,
            "n_runs": int(len(g)),
            "size_gb": round(float(g["size_gb"].sum()), 3),
        })
    return pd.DataFrame(rows)


def build_download_tranches(runs: pd.DataFrame) -> list[dict]:
    """Ordered download tranches, one patient per tranche, smallest first.

    The first tranche is the smallest scientifically valid unit — ONE patient's full trio (normal exome +
    tumor exome + tumor RNA) — enough to prove the raw -> generation -> PRIME/Epicurus -> hits@20 loop
    before committing to the full ~0.22 TB.
    """
    xwalk = patient_input_crosswalk(runs).sort_values("size_gb", kind="mergesort")
    tranches = []
    for _, r in xwalk.iterrows():
        g = runs[runs["patient_id"] == r["patient_id"]]
        tranches.append({
            "patient_id": r["patient_id"],
            "n_patients": 1,
            "n_runs": int(len(g)),
            "runs": sorted(g["run"]),
            "library_strategies": sorted(set(g["library_strategy"])),
            "complete_trio": bool(r["complete"]),
            "size_gb": round(float(g["size_gb"].sum()), 4),
        })
    return tranches


# ---------------------------------------------------------------------------
# S1/S2 label ingestion contract
# ---------------------------------------------------------------------------
def _valid_peptide(value: object) -> bool:
    pep = str(value).strip().upper()
    return bool(pep) and not (set(pep) - _STD_AA)


def validate_recognition_labels(frame: pd.DataFrame, *, sra_patients: set[str]) -> dict:
    """Validate an ingested S1/S2 recognition-label frame against the north-star ingestion contract.

    Enforces: required columns; three-state label vocabulary; valid peptides; deterministic crosswalk to
    the SRA patient set; and — critically — it PRESERVES every row (never collapses), reporting the count
    of same-key conflicting labels for human adjudication rather than silently resolving them.
    """
    required = ["patient_id", "mutant_peptide", "label", "assay", "assay_timepoint"]
    missing_cols = [c for c in required if c not in frame.columns]
    if missing_cols:
        return {"ok": False, "error": f"missing columns: {missing_cols}", "n_rows": int(len(frame))}

    labels = frame["label"].astype(str).str.upper()
    invalid_labels = sorted(set(labels) - VALID_LABELS)
    n_invalid_peptides = int((~frame["mutant_peptide"].map(_valid_peptide)).sum())

    patients = set(frame["patient_id"].astype(str))
    patients_without_inputs = sorted(patients - set(sra_patients))

    # conflict detection WITHOUT collapsing: same full key, >1 distinct label
    key_cols = [c for c in LABEL_KEY if c in frame.columns]
    keyed = frame.assign(_label_u=labels)
    conflicts = keyed.groupby(key_cols)["_label_u"].nunique()
    n_conflicting_keys = int((conflicts > 1).sum())

    label_counts = {lab: int((labels == lab).sum()) for lab in sorted(set(labels) & VALID_LABELS)}

    ok = not invalid_labels and n_invalid_peptides == 0 and not missing_cols
    return {
        "ok": bool(ok),
        "n_rows": int(len(frame)),                 # every ingested row is retained
        "label_counts": label_counts,
        "invalid_labels": invalid_labels,
        "n_invalid_peptides": n_invalid_peptides,
        "patients_without_inputs": patients_without_inputs,
        "n_conflicting_keys": n_conflicting_keys,   # surfaced, NOT resolved
        "n_patients": len(patients),
    }
