"""Label-blind patient configuration for the Miller IPV raw-data pipeline."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from benchmark.miller_ingest import SRA_RUNINFO_FIXTURE, parse_sra_runinfo

ROOT = Path(__file__).resolve().parents[2]


@dataclass(frozen=True)
class MillerPatient:
    patient_id: str
    normal_exome_run: str
    tumor_exome_run: str
    tumor_rna_run: str

    @property
    def slug(self) -> str:
        return self.patient_id.lower()

    @property
    def raw_dir(self) -> Path:
        return ROOT / "data/raw/miller_ipv" / self.slug

    @property
    def artifact_dir(self) -> Path:
        if self.patient_id == "Hu_287":
            return ROOT / "artifacts/milestone_7_decision/external_validation/miller_ipv/hu_287_reconstruction"
        return ROOT / "artifacts/milestone_8_generalization/patients" / self.patient_id

    @property
    def tumor_sample(self) -> str:
        return f"{self.patient_id}_T"

    @property
    def normal_sample(self) -> str:
        return f"{self.patient_id}_N"


def load_patient(patient_id: str, runinfo_path: str | Path = SRA_RUNINFO_FIXTURE) -> MillerPatient:
    """Resolve the three public runs using SRA metadata only; no outcome file is opened."""
    runs = parse_sra_runinfo(runinfo_path)
    group = runs[runs["patient_id"] == patient_id]
    if group.empty:
        raise ValueError(f"unknown Miller patient: {patient_id}")
    by_role: dict[str, list[str]] = {}
    for row in group.itertuples():
        by_role.setdefault(row.assay_kind, []).append(row.run)
    required = ("normal_exome", "tumor_exome", "tumor_rna")
    bad = {role: by_role.get(role, []) for role in required if len(by_role.get(role, [])) != 1}
    if bad:
        raise ValueError(f"patient {patient_id} does not have exactly one run per required role: {bad}")
    return MillerPatient(
        patient_id=patient_id,
        normal_exome_run=by_role["normal_exome"][0],
        tumor_exome_run=by_role["tumor_exome"][0],
        tumor_rna_run=by_role["tumor_rna"][0],
    )

