"""Conservative normalization of the local CEDAR T-cell export into a mutation-derived
cancer recognition-training asset (NOT a decision denominator).

Design goals (reviewer contract):
    * Distinguish TRUE mutation-derived cancer recognition assays from generic (non-mutated)
      tumor antigens, pathogen epitopes, autoimmune/allergic rows, and ambiguous rows.
    * Preserve assay-level positives AND explicit negatives, source PMID + assay IRI (row),
      host / disease context / HLA, and contradictions. NEVER infer a negative.
    * Quantify unique peptides, PMIDs, duplicates, contradictory outcomes, and overlap with
      Zhao and the current Event-B backbone.

Conservative inclusion rule (each row must satisfy all):
    1. mutation-derived: CEDAR ``Related Object | Epitope Relation`` is an in-frame /
       frameshift / fusion / unspecified neo-epitope, OR the ``Epitope | Mutation`` field is
       populated. Explicit ``analog`` and ``mimotope`` relations are EXCLUDED (synthetic /
       non-genetic). Post-translational-modification-only epitopes carry no ``Mutation`` and
       no neo-epitope relation, so they are excluded automatically.
    2. linear peptide epitope with a recoverable 8-25mer amino-acid sequence.
    3. human host (``Host | Name`` starts with "Homo sapiens").
    4. not a pathogen-derived epitope (defense-in-depth on ``Epitope | Source Organism``;
       neo-epitope source organism is typically blank because the mutated protein is self).

Disease context is a TAG, not a hard filter: a bona-fide neo-epitope is cancer-relevant even
when assayed in a healthy donor (priming) or with a blank disease. ``context_class`` records
cancer_patient / healthy_donor / unknown so downstream training can choose.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
import re
import zipfile

import pandas as pd


DEFAULT_ZIP = Path("data/raw/cedar/tcell_full_v3_tsv.zip")
TSV_MEMBER = "tcell_full_v3.tsv"
EXPECTED_SHA256 = "494d5b075198cc95135d0074e1c7d22a2c823a943d0cc664706c526dffeb2101"
CEDAR_EXPORT_URL = "https://cedar.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip"
CEDAR_LICENSE = "Free to use with attribution (CEDAR/IEDB); cite Koşaloğlu-Yalçın et al."

# Flat "Group|Field" source columns -> canonical short names.
COLUMN_MAP = {
    "Assay ID|CEDAR IRI": "cedar_assay_iri",
    "Reference|CEDAR IRI": "cedar_reference_iri",
    "Reference|PMID": "pmid",
    "Epitope|Name": "epitope_name",
    "Epitope|Object Type": "object_type",
    "Epitope|Mutation": "mutation",
    "Epitope|Source Molecule": "source_molecule",
    "Epitope|Source Organism": "source_organism",
    "Related Object|Epitope Relation": "epitope_relation",
    "Host|Name": "host",
    "1st in vivo Process|Disease": "disease",
    "Assay|Method": "assay_method",
    "Assay|Qualitative Measurement": "qualitative_measurement",
    "Assay|Number of Subjects Tested": "n_subjects_tested",
    "Assay|Number of Subjects Positive": "n_subjects_positive",
    "MHC Restriction|Name": "mhc_allele",
    "MHC Restriction|Class": "mhc_class",
}

NEO_RELATIONS = {
    "in-frame neo-epitope",
    "frameshift neo-epitope",
    "fusion neo-epitope",
    "unspecified neo-epitope",
}
EXCLUDE_RELATIONS = {"analog", "mimotope"}

PATHOGEN_REGEX = re.compile(
    r"virus|viral|hepat|papilloma|corona|herpes|influenza|Plasmodium|bacter|Mycobact|"
    r"HIV|retrovir|Epstein|cytomegalo|Listeria|Salmonella|Toxoplasma|Leishmania",
    re.IGNORECASE,
)
CANCER_REGEX = re.compile(
    r"cancer|carcinoma|melanoma|leukemi|lymphoma|myeloma|sarcoma|glioma|glioblast|tumou?r|"
    r"neoplasm|adenocarc|blastoma|mesothelioma|myelodysplas|myelofibros|polycythemia|"
    r"thrombocythemia|astrocytoma|malignan",
    re.IGNORECASE,
)
STD_AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_LEN = 8
MAX_LEN = 25
POSITIVE_MEASURES = {"Positive", "Positive-Low", "Positive-High", "Positive-Intermediate"}
NEGATIVE_MEASURES = {"Negative"}


def sha256_file(path: str | Path) -> str:
    digest = sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _flat_header(zip_path: Path) -> list[str]:
    with zipfile.ZipFile(zip_path) as archive:
        with archive.open(TSV_MEMBER) as handle:
            reader = csv.reader((line.decode() for line in handle), delimiter="\t")
            group = next(reader)
            field = next(reader)
    return [f"{a}|{b}" for a, b in zip(group, field)]


def read_cedar_tcell(zip_path: str | Path = DEFAULT_ZIP, *, nrows: int | None = None) -> pd.DataFrame:
    """Read only the needed columns from the CEDAR T-cell export (direct from the zip)."""
    zip_path = Path(zip_path)
    flat = _flat_header(zip_path)
    missing = [c for c in COLUMN_MAP if c not in flat]
    if missing:
        raise ValueError(f"CEDAR export is missing expected columns: {missing}")
    positions = sorted(flat.index(c) for c in COLUMN_MAP)
    with zipfile.ZipFile(zip_path) as archive, archive.open(TSV_MEMBER) as handle:
        frame = pd.read_csv(
            handle,
            sep="\t",
            header=None,
            names=flat,
            usecols=positions,
            skiprows=2,
            dtype=str,
            nrows=nrows,
            low_memory=False,
        )
    frame = frame[list(COLUMN_MAP)].rename(columns=COLUMN_MAP)
    return frame


def extract_peptide(name: object) -> tuple[str, str]:
    """Return (sequence, reason). A usable sequence is 8-25 standard amino acids."""
    seq = re.sub(r"[^A-Za-z]", "", str(name)).upper()
    if not seq:
        return "", "no_sequence"
    if not set(seq).issubset(STD_AA):
        return "", "nonstandard_residue"
    if len(seq) < MIN_LEN:
        return "", "too_short"
    if len(seq) > MAX_LEN:
        return "", "too_long"
    return seq, "ok"


def response_label(measure: object) -> tuple[str, str] | tuple[None, str]:
    """Map the CEDAR qualitative measure to a three-state label; never infer a negative."""
    text = str(measure).strip()
    if text in POSITIVE_MEASURES:
        return "POSITIVE", text
    if text in NEGATIVE_MEASURES:
        return "TESTED_NEGATIVE", text
    return None, text or "<missing>"


def classify_context(disease: object) -> str:
    text = str(disease).strip()
    if not text or text.lower() in {"nan", "<na>"}:
        return "unknown"
    if text.lower() == "healthy":
        return "healthy_donor"
    if CANCER_REGEX.search(text):
        return "cancer_patient"
    return "other"


@dataclass(frozen=True)
class CedarNormalization:
    kept: pd.DataFrame
    exclusion_ledger: dict


def normalize_cedar(frame: pd.DataFrame) -> CedarNormalization:
    """Apply the conservative filter; return the assay-level asset + an exclusion ledger."""
    work = frame.copy()
    relation = work["epitope_relation"].fillna("")
    has_mutation = work["mutation"].notna() & work["mutation"].astype(str).str.strip().ne("")
    is_neo = relation.isin(NEO_RELATIONS)
    is_excluded_relation = relation.isin(EXCLUDE_RELATIONS)
    mutation_derived = (is_neo | has_mutation) & ~is_excluded_relation

    linear = work["object_type"].fillna("").eq("Linear peptide")
    human_host = work["host"].fillna("").str.startswith("Homo sapiens")
    pathogen = work["source_organism"].fillna("").map(lambda v: bool(PATHOGEN_REGEX.search(v)))

    seq_reason = work["epitope_name"].map(extract_peptide)
    work["peptide"] = [s for s, _ in seq_reason]
    reasons = pd.Series([r for _, r in seq_reason], index=work.index)
    valid_seq = work["peptide"].ne("")

    labels = work["qualitative_measurement"].map(response_label)
    work["response_label"] = [lab for lab, _ in labels]
    work["source_interpretation"] = [raw for _, raw in labels]
    scorable = work["response_label"].notna()

    keep = mutation_derived & linear & human_host & ~pathogen & valid_seq & scorable

    ledger = {
        "total_rows": int(len(work)),
        "excluded_not_mutation_derived": int((~mutation_derived).sum()),
        "excluded_relation_analog_or_mimotope": int(is_excluded_relation.sum()),
        "excluded_not_linear_peptide": int((mutation_derived & ~linear).sum()),
        "excluded_non_human_host": int((mutation_derived & linear & ~human_host).sum()),
        "excluded_pathogen_source": int((mutation_derived & linear & human_host & pathogen).sum()),
        "excluded_unusable_sequence": int(
            (mutation_derived & linear & human_host & ~pathogen & ~valid_seq).sum()
        ),
        "sequence_exclusion_reasons": reasons[
            mutation_derived & linear & human_host & ~pathogen & ~valid_seq
        ].value_counts().to_dict(),
        "excluded_unscorable_measure": int(
            (mutation_derived & linear & human_host & ~pathogen & valid_seq & ~scorable).sum()
        ),
        "kept_rows": int(keep.sum()),
    }

    kept = work.loc[keep].copy()
    kept["peptide_length"] = kept["peptide"].str.len()
    kept["context_class"] = kept["disease"].map(classify_context)
    kept["cancer_context"] = kept["context_class"].eq("cancer_patient")

    # Contradiction + duplicate flags (preserve every row; never collapse).
    kept["_is_positive"] = kept["response_label"].eq("POSITIVE")
    pair = kept.groupby(["peptide", "mhc_allele"], dropna=False)["_is_positive"]
    contradictory_pairs = pair.transform("nunique").gt(1)
    kept["contradiction_flag"] = contradictory_pairs
    kept["is_duplicate_assay"] = kept.duplicated(
        ["peptide", "mhc_allele", "pmid", "qualitative_measurement"], keep="first"
    )
    kept = kept.drop(columns=["_is_positive"])

    columns = [
        "cedar_assay_iri", "cedar_reference_iri", "pmid", "peptide", "peptide_length",
        "mutation", "epitope_relation", "source_molecule", "host", "disease",
        "context_class", "cancer_context", "assay_method", "mhc_allele", "mhc_class",
        "qualitative_measurement", "response_label", "source_interpretation",
        "n_subjects_tested", "n_subjects_positive", "contradiction_flag", "is_duplicate_assay",
    ]
    kept = kept.loc[:, columns].sort_values(["peptide", "mhc_allele", "cedar_assay_iri"]).reset_index(
        drop=True
    )
    return CedarNormalization(kept, ledger)


def cedar_audit(
    kept: pd.DataFrame,
    *,
    zhao_peptides: set[str] | None = None,
    backbone_peptides: set[str] | None = None,
) -> dict:
    """Quantify the asset: uniques, labels, dups, contradictions, and corpus overlap."""
    peptides = set(kept["peptide"])
    pair = kept.groupby(["peptide", "mhc_allele"], dropna=False)["response_label"]
    contradictory_pair_n = int((pair.nunique() > 1).sum())

    audit = {
        "kept_rows": int(len(kept)),
        "unique_peptides": len(peptides),
        "unique_pmids": int(kept["pmid"].nunique()),
        "unique_peptide_hla_pairs": int(kept.groupby(["peptide", "mhc_allele"], dropna=False).ngroups),
        "label_counts": kept["response_label"].value_counts().to_dict(),
        "positive_sublabels": kept.loc[
            kept.response_label.eq("POSITIVE"), "source_interpretation"
        ].value_counts().to_dict(),
        "assay_method_counts": kept["assay_method"].value_counts().head(12).to_dict(),
        "mhc_class_counts": kept["mhc_class"].fillna("<NA>").value_counts().to_dict(),
        "epitope_relation_counts": kept["epitope_relation"].fillna("<NA>").value_counts().to_dict(),
        "context_class_counts": kept["context_class"].value_counts().to_dict(),
        "cancer_context_rows": int(kept["cancer_context"].sum()),
        "duplicate_assay_rows": int(kept["is_duplicate_assay"].sum()),
        "contradictory_peptide_hla_pairs": contradictory_pair_n,
        "contradictory_rows": int(kept["contradiction_flag"].sum()),
        "top_source_molecules": kept["source_molecule"].fillna("<NA>").value_counts().head(12).to_dict(),
    }
    if zhao_peptides is not None:
        overlap = peptides & zhao_peptides
        audit["overlap_zhao"] = {
            "peptides": len(overlap),
            "examples": sorted(overlap)[:20],
            "note": "Shared peptides must be held out together when pooling with Zhao.",
        }
    if backbone_peptides is not None:
        overlap = peptides & backbone_peptides
        audit["overlap_backbone"] = {
            "peptides": len(overlap),
            "examples": sorted(overlap)[:20],
            "note": "Backbone (IMPROVE/direct-screen) is literature-derived; this overlap is a "
            "leakage risk and must be purged or grouped before any pooled evaluation.",
        }
    audit["scope"] = (
        "Recognition-training asset ONLY (assay-level literature aggregation). It is NOT a "
        "per-patient candidate denominator and must not back a top-K decision-problem claim."
    )
    return audit


def render_cedar_markdown(ledger: dict, audit: dict) -> str:
    lines = [
        "# CEDAR mutation-derived cancer recognition asset — normalization audit",
        "",
        f"Source: local CEDAR T-cell export `{TSV_MEMBER}` (sha256 `{EXPECTED_SHA256[:16]}…`). "
        f"{CEDAR_LICENSE}",
        "",
        "## Conservative filter cascade",
        "",
        "| Stage | Rows |",
        "|---|---:|",
        f"| Total assay rows | {ledger['total_rows']} |",
        f"| Excluded: not mutation-derived | {ledger['excluded_not_mutation_derived']} |",
        f"| Excluded: analog / mimotope relation | {ledger['excluded_relation_analog_or_mimotope']} |",
        f"| Excluded: not a linear peptide | {ledger['excluded_not_linear_peptide']} |",
        f"| Excluded: non-human host | {ledger['excluded_non_human_host']} |",
        f"| Excluded: pathogen source | {ledger['excluded_pathogen_source']} |",
        f"| Excluded: unusable sequence | {ledger['excluded_unusable_sequence']} "
        f"({ledger['sequence_exclusion_reasons']}) |",
        f"| Excluded: unscorable measure | {ledger['excluded_unscorable_measure']} |",
        f"| **Kept (recognition asset)** | **{ledger['kept_rows']}** |",
        "",
        "## Asset summary",
        "",
        f"- Unique peptides: **{audit['unique_peptides']}**; unique PMIDs: **{audit['unique_pmids']}**; "
        f"unique (peptide, HLA) pairs: {audit['unique_peptide_hla_pairs']}",
        f"- Labels: `{audit['label_counts']}` (positive sublabels: `{audit['positive_sublabels']}`)",
        f"- MHC class: `{audit['mhc_class_counts']}`; relations: `{audit['epitope_relation_counts']}`",
        f"- Context: `{audit['context_class_counts']}` (cancer-patient rows: {audit['cancer_context_rows']})",
        f"- Duplicate assay rows (preserved): {audit['duplicate_assay_rows']}",
        f"- Contradictory (peptide, HLA) pairs: **{audit['contradictory_peptide_hla_pairs']}** "
        f"({audit['contradictory_rows']} rows; preserved, flagged, never collapsed)",
        f"- Top assay methods: `{audit['assay_method_counts']}`",
        "",
        "## Overlap / leakage",
        "",
    ]
    if "overlap_zhao" in audit:
        lines.append(f"- Overlap with Zhao peptides: **{audit['overlap_zhao']['peptides']}**")
    if "overlap_backbone" in audit:
        lines.append(
            f"- Overlap with Event-B backbone peptides: **{audit['overlap_backbone']['peptides']}** "
            f"— {audit['overlap_backbone']['note']}"
        )
    lines += ["", f"**Scope:** {audit['scope']}", ""]
    return "\n".join(lines)
