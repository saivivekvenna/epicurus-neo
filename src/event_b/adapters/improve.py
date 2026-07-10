"""Compatibility adapter: IMPROVE remains Event A, never silently becomes Event B."""

from __future__ import annotations

from hashlib import sha256
import json
from pathlib import Path
from zipfile import ZipFile

import pandas as pd

from event_b.adapters.base import AdapterDeclaration
from event_b.corpus import EventBCorpus
from event_b.manifest import SourceManifest
from event_b.models import (
    AssayType,
    BiologicalEvent,
    MHCClass,
    ResponseLabel,
    ReviewStatus,
    SCHEMAS,
    ValueOrigin,
    VaccineInclusion,
    stable_candidate_id,
)


IMPROVE_MEMBER = "data/03_data_for_CV/IMPROVE/03_3_final_peptide_features_Partition_for_CV.txt"


def _prov(entity: str, entity_id: str, document: str, row: int, fragment: str) -> dict:
    provenance_id = "prov:" + sha256(f"{entity}|{entity_id}|{row}".encode()).hexdigest()[:20]
    return {
        "provenance_id": provenance_id,
        "entity_type": entity,
        "entity_id": entity_id,
        "field_name": "*",
        "source_document": document,
        "row": row,
        "source_fragment": fragment,
        "extraction_method": "deterministic_table_adapter",
        "extraction_confidence": 1.0,
        "value_origin": ValueOrigin.SOURCE_REPORTED.value,
        "review_status": ReviewStatus.ACCEPTED.value,
    }


class ImproveEventAAdapter:
    declaration = AdapterDeclaration(
        "IMPROVE official CV",
        "paper-release",
        "improve_event_a",
        "1.0.0",
        ("studies", "patients", "candidates", "assays", "provenance"),
        (BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value,),
        ("Event-A pre-existing reactivity only; not vaccine-inducible Event B",),
        ("response=0 is an explicitly screened non-reactive candidate",),
        ("vaccine", "clinical_outcome", "vaccination schedule"),
    )

    def __init__(self, data_zip: str | Path):
        self.data_zip = Path(data_zip)

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        with ZipFile(self.data_zip) as archive, archive.open(IMPROVE_MEMBER) as handle:
            return {"raw": pd.read_csv(handle, sep="\t")}

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        raw = pd.DataFrame(extracted["raw"])
        document = f"{self.data_zip.name}:{IMPROVE_MEMBER}"
        studies, patients, candidates, assays, provenance = {}, {}, [], [], []
        for source_row, row in raw.reset_index(drop=True).iterrows():
            study_id = "improve"
            patient_id = f"improve:{row.Patient}"
            study_prov = "prov:" + sha256(f"study|{study_id}".encode()).hexdigest()[:20]
            patient_prov = "prov:" + sha256(f"patient|{patient_id}".encode()).hexdigest()[:20]
            studies[study_id] = {
                "study_id": study_id,
                "title": "IMPROVE neoantigen immunogenicity benchmark",
                "cancer_type": "multi-cancer",
                "source_paths": str(self.data_zip),
                "source_checksums": manifest.documents[0].checksum_sha256,
                "source_manifest_id": manifest.manifest_id,
                "provenance_id": study_prov,
            }
            patients[patient_id] = {
                "patient_id": patient_id,
                "source_patient_id": row.Patient,
                "study_id": study_id,
                "cancer_type": row.cohort,
                "hla_alleles": json.dumps(
                    sorted(raw.loc[raw.Patient.eq(row.Patient), "HLA_allele"].unique())
                ),
                "provenance_id": patient_prov,
            }
            record = {
                "study_id": study_id,
                "patient_id": patient_id,
                "sample_id": str(row.Sample),
                "timepoint": "UNKNOWN",
                "genomic_variant": str(row.Genomic_Position),
                "transcript": str(row.Transcript_ID),
                "mutant_peptide": str(row.Mut_peptide),
                "hla_alleles": [str(row.HLA_allele)],
            }
            candidate_id = stable_candidate_id(record)
            candidate_prov = _prov(
                "candidates", candidate_id, document, source_row + 2, str(row.identity)
            )
            assay_id = "assay:" + sha256(f"improve|{candidate_id}".encode()).hexdigest()[:20]
            assay_prov = _prov(
                "assays", assay_id, document, source_row + 2, f"response={row.response}"
            )
            candidates.append(
                {
                    "candidate_id": candidate_id,
                    "patient_id": patient_id,
                    "study_id": study_id,
                    "sample_id": str(row.Sample),
                    "genomic_variant": str(row.Genomic_Position),
                    "gene": row.Gene_Symbol,
                    "transcript": str(row.Transcript_ID),
                    "protein_change": str(row.Protein_position),
                    "mutant_peptide": str(row.Mut_peptide),
                    "wildtype_peptide": str(row.Norm_peptide),
                    "peptide_length": len(str(row.Mut_peptide)),
                    "hla_alleles": json.dumps([str(row.HLA_allele)]),
                    "mhc_class": MHCClass.CLASS_I.value,
                    "candidate_source": "IMPROVE official candidate table",
                    "vaccine_inclusion": VaccineInclusion.UNKNOWN.value,
                    "vaccine_inclusion_origin": ValueOrigin.UNKNOWN.value,
                    "mutant_wildtype_verified": True,
                    "provenance_id": candidate_prov["provenance_id"],
                }
            )
            assays.append(
                {
                    "assay_id": assay_id,
                    "patient_id": patient_id,
                    "study_id": study_id,
                    "candidate_id": candidate_id,
                    "assay_type": AssayType.OTHER.value,
                    "timepoint": "UNKNOWN",
                    "relative_to_vaccine": "UNKNOWN",
                    "source_interpretation": int(row.response),
                    "event_type": BiologicalEvent.EVENT_A_PREEXISTING_REACTIVITY.value,
                    "response_label": (
                        ResponseLabel.POSITIVE.value
                        if int(row.response) == 1
                        else ResponseLabel.TESTED_NEGATIVE.value
                    ),
                    "explicit_assay_inclusion": True,
                    "review_status": ReviewStatus.ACCEPTED.value,
                    "provenance_id": assay_prov["provenance_id"],
                }
            )
            provenance.extend([candidate_prov, assay_prov])
        for study_id, study in studies.items():
            provenance.append(_prov("studies", study_id, document, 1, "cohort-level metadata"))
            study["provenance_id"] = provenance[-1]["provenance_id"]
        for patient_id, patient in patients.items():
            provenance.append(_prov("patients", patient_id, document, 1, "patient/HLA metadata"))
            patient["provenance_id"] = provenance[-1]["provenance_id"]
        return EventBCorpus(
            studies=SCHEMAS["studies"].normalize(pd.DataFrame(studies.values())),
            patients=SCHEMAS["patients"].normalize(pd.DataFrame(patients.values())),
            candidates=SCHEMAS["candidates"].normalize(pd.DataFrame(candidates)),
            assays=SCHEMAS["assays"].normalize(pd.DataFrame(assays)),
            provenance=SCHEMAS["provenance"].normalize(pd.DataFrame(provenance)),
        )
