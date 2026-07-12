"""Sid RNA-support gate (EXPLORATORY) — absence-safe removal of confidently-unexpressed mutations.

Leakage-safe, label-blind, biologically justified: a mutation whose MUTANT allele is not transcribed cannot
present a class-I neoantigen. Matched Sid tumor RNA evidence for all 130 generated mutations is already in
`data/raw/osteosarc/site_cache/variant_vafs_long.tsv` (per-variant tumor RNA assay rows). We remove a
mutation ONLY when RNA data EXISTS and positively shows non-expression (tumor RNA mutant `alt_reads == 0`
AND gene `expression_tpm == 0`). ABSENCE IS NEVER A VETO — a mutation with no RNA row, or any mutant RNA
read, or any expression, is KEPT. No threshold is tuned on the three Hudson labels.

Audit result (see the arm runner): this gate flags 18/130 mutations, never a recognized positive, and
improves the missed positives' ranks (ASPM #39→~29 PRIME / #20→~16 MixMHCpred; MAP2 #26→~21 MixMHCpred) but
does NOT reach 3/3 — the missed positive stays just outside top-20. It is a principled partial, not a win.
The exact per-peptide mutant-read fallback (public T2 RNA BAM via remote pysam) is unnecessary for
variant-level gating and would not change this verdict.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

VAF_TABLE = Path("data/raw/osteosarc/site_cache/variant_vafs_long.tsv")


def load_tumor_rna_support(vaf_table: Path = VAF_TABLE) -> pd.DataFrame:
    """Per variant_id: does a tumor RNA assay row exist, and its max mutant alt_reads / RNA-VAF.
    Label-blind (no recognition label touched)."""
    d = pd.read_csv(vaf_table, sep="\t")
    tumor = d[d["tissue"].astype(str).str.lower().eq("tumor")]
    rna = tumor[tumor["assay_type"].astype(str).str.contains("RNA", case=False, na=False)]
    g = rna.groupby("variant_id").agg(rna_max_alt_reads=("alt_reads", "max"),
                                      rna_max_vaf=("vaf", "max")).reset_index()
    g["rna_assay_present"] = True
    return g


def confidently_unexpressed(mutation_ids: pd.Series, expression_tpm: pd.Series,
                            rna: pd.DataFrame) -> np.ndarray:
    """Absence-safe boolean per row: True ONLY if RNA data exists AND shows no mutant reads AND TPM==0.
    Missing RNA row, any mutant read, or any TPM -> False (KEEP)."""
    r = rna.set_index("variant_id")
    present = mutation_ids.map(lambda m: bool(r["rna_assay_present"].get(m, False))).to_numpy()
    alt = mutation_ids.map(lambda m: float(r["rna_max_alt_reads"].get(m, np.nan))).to_numpy()
    tpm = pd.to_numeric(expression_tpm, errors="coerce").to_numpy()
    with np.errstate(invalid="ignore"):
        return present & (alt == 0) & (tpm == 0)


def rna_support_gate(candidates: pd.DataFrame, *, mutation_col: str = "mutation_id",
                     tpm_col: str = "expression_tpm", rna: pd.DataFrame | None = None) -> pd.DataFrame:
    """Add `rna_unexpressed` (removal flag) and `rna_gate_keep` columns. Never removes a mutation whose
    RNA evidence is absent or positive. Row-level (mutation-consistent) — collapse to mutation for scoring."""
    r = rna if rna is not None else load_tumor_rna_support()
    out = candidates.copy()
    out["rna_unexpressed"] = confidently_unexpressed(out[mutation_col], out[tpm_col], r)
    out["rna_gate_keep"] = ~out["rna_unexpressed"]
    return out
