"""Nous-209 Lynch-syndrome vaccine adapter with strict pool non-decomposition."""

from __future__ import annotations

import json
from pathlib import Path
import shutil
import urllib.request

import pandas as pd

from event_b.adapters.base import AdapterDeclaration
from event_b.adapters.mkras_vax import _frame, _id, _prov
from event_b.corpus import EventBCorpus
from event_b.manifest import SourceManifest, manifest_from_paths, sha256_file
from event_b.models import (
    AssayType,
    BiologicalEvent,
    MHCClass,
    ResponseLabel,
    ReviewStatus,
    SCHEMAS,
    ValueOrigin,
)


STUDY_ID = "nous_209_2025"
COHORT_ID = "NCT05078866_initial_vaccination"
DOI = "10.1038/s41591-025-04182-9"
NCT = "NCT05078866"
SOURCE_FILES = {
    "41591_2025_4182_MOESM3_ESM.xlsx": "4b0ce79f8e6021dbc3ecd8c04c952b05fbdd7173d595c6fb3a1fc512ddbf8dc5",
    "41591_2025_4182_MOESM6_ESM.xlsx": "62481d675055259c55a4000bb7e2e22183396c8ac23422d4df098f2487a5c038",
}
BASE_URL = (
    "https://static-content.springer.com/esm/"
    "art%3A10.1038%2Fs41591-025-04182-9/MediaObjects/"
)


def stage_sources(raw_dir: str | Path) -> list[Path]:
    root = Path(raw_dir)
    root.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, expected in SOURCE_FILES.items():
        path = root / name
        if not path.exists() or sha256_file(path) != expected:
            url = BASE_URL + name
            request = urllib.request.Request(url, headers={"User-Agent": "epicurus-neo/1.0"})
            try:
                with urllib.request.urlopen(request, timeout=180) as response, path.open("wb") as out:
                    shutil.copyfileobj(response, out)
            except Exception as error:  # noqa: BLE001 - actionable manual-source contract
                raise RuntimeError(
                    f"Download {url} manually to {path}; expected sha256={expected}"
                ) from error
        observed = sha256_file(path)
        if observed != expected:
            raise RuntimeError(f"checksum mismatch for {path}: observed {observed}")
        paths.append(path)
    return paths


def source_manifest(raw_dir: str | Path) -> SourceManifest:
    adapter = Nous209Adapter(raw_dir)
    return manifest_from_paths(
        adapter.declaration.source_name,
        adapter.declaration.source_version,
        adapter.declaration.adapter_name,
        adapter.declaration.adapter_version,
        stage_sources(raw_dir),
    )


def _patient_number(value: object) -> str:
    return str(value).replace("Pt", "").replace("\xa0", " ").strip()


class Nous209Adapter:
    declaration = AdapterDeclaration(
        "Nous-209 Lynch syndrome cancer-interception cohort",
        "Nature-Medicine-2026",
        "nous_209_event_b",
        "1.0.0",
        ("antigens", "studies", "patients", "vaccines", "assays", "entity_relationships", "provenance"),
        (BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,),
        (
            "The patient source reports number of reactive pools, not patient-by-pool identities.",
            "The 115 deconvolved immunogenic FSP identifiers have no public patient mapping here.",
            "FSP sequences and restricting HLA alleles are not in the ingested source tables.",
            "Four participants had some baseline spontaneous pool reactivity, but their identities "
            "are not resolved by the public source table used for patient outcomes.",
        ),
        (
            "Each of 37 evaluable participants has a source-reported positive response at week 9.",
            "Reactive-pool counts are patient-level quantitative observations, not 209 positive labels.",
            "Global FSP-to-pool membership is preserved independently of patient observations.",
        ),
        ("patient_pool_identity", "fsp_sequence", "patient_hla_genotype", "tumor_fsp_confirmation"),
        canonical_study_id=STUDY_ID,
        cohort_id=COHORT_ID,
        source_files=tuple(SOURCE_FILES),
        supported_timepoints=("WEEK_9_PEAK",),
        positivity_rules=("predefined protocol ELISpot positivity; author call at peak",),
        baseline_semantics="Baseline spontaneous responses exist in 4/37 but are not patient-mapped.",
        vaccine_component_structure="GAd prime/MVA boost encoding 209 shared FSPs.",
        assay_target_structure="Sixteen multi-FSP peptide pools; identities absent per patient.",
        candidate_identity_completeness="PATIENT_LEVEL_ONLY_NO_CANDIDATE_MAPPING",
        unresolved_ambiguities=(
            "Cannot link a participant's reactive-pool count to specific pools or FSPs.",
        ),
    )

    def __init__(self, raw_dir: str | Path):
        self.raw_dir = Path(raw_dir)
        self.review_issues = []

    def extract(self, manifest: SourceManifest) -> dict[str, object]:
        del manifest
        paths = {path.name: path for path in stage_sources(self.raw_dir)}
        fsp = pd.read_excel(
            paths["41591_2025_4182_MOESM3_ESM.xlsx"],
            sheet_name="Supplementary Table 3",
            header=2,
            usecols="A:B",
            engine="openpyxl",
        ).dropna(how="all")
        patients = pd.read_excel(
            paths["41591_2025_4182_MOESM6_ESM.xlsx"],
            sheet_name="Fig 2c",
            header=0,
            usecols="A:C",
            engine="openpyxl",
        ).dropna(how="all")
        patients.columns = ["patient", "gender", "reactive_pools"]
        return {"fsp": fsp, "patients": patients}

    def normalize(self, extracted: dict[str, object], manifest: SourceManifest) -> EventBCorpus:
        fsp = pd.DataFrame(extracted["fsp"])
        patients = pd.DataFrame(extracted["patients"])
        if list(fsp.columns) != ["Immunogenic FSPS", "Peptide pools"]:
            raise ValueError("Nous-209 FSP table columns changed")
        normalized_pools = fsp["Peptide pools"].astype(str).str.replace(" ", "", regex=False)
        expected_deconvolved_pools = {f"Pool{index}" for index in range(1, 17)} - {"Pool15"}
        if len(fsp) != 115 or set(normalized_pools) != expected_deconvolved_pools:
            raise ValueError("Nous-209 deconvolved FSP denominator changed")
        if len(patients) != 37 or set(pd.to_numeric(patients.reactive_pools)) - set(range(1, 17)):
            raise ValueError("Nous-209 patient reactive-pool denominator changed")

        rows: dict[str, list[dict]] = {entity: [] for entity in SCHEMAS}

        def add(entity: str, row: dict, provenance: dict) -> None:
            row["provenance_id"] = provenance["provenance_id"]
            rows[entity].append(row)
            rows["provenance"].append(provenance)

        add(
            "studies",
            {
                "study_id": STUDY_ID,
                "title": "Nous-209 neoantigen vaccine for cancer prevention in Lynch syndrome",
                "publication_ids": json.dumps([f"DOI:{DOI}"]),
                "trial_id": NCT,
                "cancer_type": "Lynch syndrome carriers without active invasive cancer",
                "vaccine_platform": "heterologous GAd prime and MVA boost",
                "vaccination_schedule": "GAd prime followed by MVA boost at week 8",
                "source_urls": json.dumps([BASE_URL + name for name in SOURCE_FILES]),
                "source_paths": json.dumps([doc.local_path for doc in manifest.documents]),
                "source_checksums": json.dumps(SOURCE_FILES, sort_keys=True),
                "source_manifest_id": manifest.manifest_id,
            },
            _prov(
                "studies",
                STUDY_ID,
                document=f"DOI:{DOI}",
                table="main article",
                row="trial identity",
                column="cohort",
                fragment="NCT05078866; 45 enrolled; 37 immunogenicity evaluable",
                method="manual_primary_source_curation",
                origin=ValueOrigin.MANUALLY_CURATED.value,
            ),
        )

        pool_ids = {}
        for pool in [f"Pool{index}" for index in range(1, 17)]:
            pool_id = f"antigen:nous209:{str(pool).lower()}"
            pool_ids[str(pool)] = pool_id
            add(
                "antigens",
                {
                    "antigen_id": pool_id,
                    "study_id": STUDY_ID,
                    "mutant_sequence": pd.NA,
                    "component_type": "MULTI_FSP_ASSAY_POOL",
                    "hla_alleles": json.dumps([]),
                    "hla_evidence_type": "NOT_ASSESSED",
                },
                _prov(
                    "antigens",
                    pool_id,
                    document="41591_2025_4182_MOESM3_ESM.xlsx",
                    table="Supplementary Table 3",
                    row=pool,
                    column="Peptide pools",
                    fragment=str(pool),
                ),
            )

        for source_index, source in fsp.reset_index(drop=True).iterrows():
            fsp_name = str(source["Immunogenic FSPS"]).strip()
            pool = str(source["Peptide pools"]).replace(" ", "").strip()
            antigen_id = f"antigen:nous209:fsp:{fsp_name.lower()}"
            gene = fsp_name.removeprefix("FSP_").rsplit("_", 1)[0]
            add(
                "antigens",
                {
                    "antigen_id": antigen_id,
                    "study_id": STUDY_ID,
                    "gene": gene,
                    "protein_change": fsp_name,
                    "mutant_sequence": pd.NA,
                    "component_type": "DECONVOLVED_IMMUNOGENIC_FSP_IDENTIFIER",
                    "hla_alleles": json.dumps([]),
                    "hla_evidence_type": "NOT_ASSESSED",
                },
                _prov(
                    "antigens",
                    antigen_id,
                    document="41591_2025_4182_MOESM3_ESM.xlsx",
                    table="Supplementary Table 3",
                    row=source_index + 4,
                    column="Immunogenic FSPS",
                    fragment=f"{fsp_name}; {pool}",
                ),
            )
            relationship_id = _id("rel", antigen_id, "CONTAINED_WITHIN", pool_ids[pool])
            add(
                "entity_relationships",
                {
                    "relationship_id": relationship_id,
                    "study_id": STUDY_ID,
                    "source_entity_type": "antigens",
                    "source_entity_id": antigen_id,
                    "target_entity_type": "antigens",
                    "target_entity_id": pool_ids[pool],
                    "relationship_type": "CONTAINED_WITHIN",
                },
                _prov(
                    "entity_relationships",
                    relationship_id,
                    document="41591_2025_4182_MOESM3_ESM.xlsx",
                    table="Supplementary Table 3",
                    row=source_index + 4,
                    column="Peptide pools",
                    fragment=f"{fsp_name} in {pool}",
                ),
            )

        for source_index, source in patients.reset_index(drop=True).iterrows():
            source_patient = _patient_number(source.patient)
            patient_id = f"{STUDY_ID}:patient_{source_patient}"
            reactive_pools = int(source.reactive_pools)
            add(
                "patients",
                {
                    "patient_id": patient_id,
                    "source_patient_id": source_patient,
                    "study_id": STUDY_ID,
                    "cancer_type": "Lynch syndrome carrier; no active invasive cancer",
                    "disease_stage": "cancer interception",
                    "treatment_context": "Nous-209 monotherapy",
                    "hla_alleles": json.dumps([]),
                    "tumor_context": "MSI/FSP risk; no patient-specific tumor mutation mapping",
                },
                _prov(
                    "patients",
                    patient_id,
                    document="41591_2025_4182_MOESM6_ESM.xlsx",
                    table="Fig 2c",
                    row=source_index + 2,
                    column="patient",
                    fragment=f"Pt {source_patient}; {reactive_pools} reactive pools",
                ),
            )
            vaccine_id = f"vaccine:{STUDY_ID}:patient_{source_patient}"
            add(
                "vaccines",
                {
                    "vaccine_id": vaccine_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "vaccine_platform": "GAd-209-FSP prime and MVA-209-FSP boost",
                    "formulation": "viral vectors encoding 209 shared FSPs",
                    "vaccination_dates": json.dumps([]),
                    "relative_schedule": "week 0 prime; week 8 boost",
                    "candidate_count": 209,
                    "mhc_class_intent": MHCClass.BOTH.value,
                    "concurrent_therapy": "none reported for interception cohort",
                },
                _prov(
                    "vaccines",
                    vaccine_id,
                    document=f"DOI:{DOI}",
                    table="main article",
                    row=source_patient,
                    column="vaccine",
                    fragment="Nous-209 shared 209-FSP construct",
                    method="manual_primary_source_curation",
                    origin=ValueOrigin.MANUALLY_CURATED.value,
                ),
            )
            assay_id = f"assay:{STUDY_ID}:patient_{source_patient}:week9_pool_count"
            add(
                "assays",
                {
                    "assay_id": assay_id,
                    "patient_id": patient_id,
                    "study_id": STUDY_ID,
                    "vaccine_id": vaccine_id,
                    "assay_type": AssayType.ELISPOT.value,
                    "sample_type": "PBMC",
                    "timepoint": "WEEK_9_PEAK",
                    "relative_to_vaccine": "POST_BOOST",
                    "stimulation_protocol": "16 pools covering 209 FSPs",
                    "positivity_threshold": "predefined protocol author call",
                    "quantitative_result": reactive_pools,
                    "result_units": "reactive pools out of 16",
                    "qualitative_result": ResponseLabel.POSITIVE.value,
                    "source_interpretation": "patient-level response; pool identities unavailable",
                    "event_type": BiologicalEvent.EVENT_B_VACCINE_INDUCED_RESPONSE.value,
                    "response_label": ResponseLabel.POSITIVE.value,
                    "explicit_assay_inclusion": True,
                    "review_status": ReviewStatus.ACCEPTED.value,
                },
                _prov(
                    "assays",
                    assay_id,
                    document="41591_2025_4182_MOESM6_ESM.xlsx",
                    table="Fig 2c",
                    row=source_index + 2,
                    column="# of reactive pools",
                    fragment=f"Pt {source_patient}: {reactive_pools}/16 reactive pools",
                ),
            )
        return EventBCorpus(**{entity: _frame(entity, rows[entity]) for entity in SCHEMAS})


def reconcile(corpus: EventBCorpus, review_queue=()) -> dict:
    return {
        "source_reported": {
            "enrolled_patients": 45,
            "immunogenicity_evaluable_patients": 37,
            "positive_patients_at_peak": 37,
            "encoded_fsps": 209,
            "assay_pools": 16,
            "deconvolved_immunogenic_fsp_identifiers": 115,
        },
        "extracted": {
            "patients": int(corpus.patients.patient_id.nunique()),
            "patient_level_event_b_observations": int(len(corpus.assays)),
            "primary_candidate_labels": int(corpus.assays.candidate_id.notna().sum()),
            "pool_entities": int(corpus.antigens.component_type.eq("MULTI_FSP_ASSAY_POOL").sum()),
            "fsp_identifiers": int(
                corpus.antigens.component_type.eq("DECONVOLVED_IMMUNOGENIC_FSP_IDENTIFIER").sum()
            ),
            "review_queue": len(tuple(review_queue)),
        },
        "reconciles": bool(
            corpus.patients.patient_id.nunique() == 37
            and len(corpus.assays) == 37
            and corpus.assays.response_label.eq(ResponseLabel.POSITIVE.value).all()
            and corpus.assays.candidate_id.isna().all()
        ),
        "material_discrepancies": [
            "Patient-by-pool identities are not public; reactive-pool counts remain patient-level "
            "and produce zero candidate labels.",
            "The 115-row deconvolution table lists no FSP for Pool15 and formats Pool7 once as "
            "'Pool 7'; all 16 assay pools are preserved but Pool15 has no FSP membership edge."
        ],
    }
