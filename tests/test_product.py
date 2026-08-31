import json
from pathlib import Path

import pandas as pd
import pytest

from epicurus_neo.contracts import validate_candidate_contract, validate_ranked_contract
from epicurus_neo.gates import apply_deterministic_gate, summarize_gate
from epicurus_neo.product import (
    InferenceConfig,
    merge_rna_evidence,
    normalize_product_candidates,
    patient_summaries,
    run_product_inference,
    score_product_candidates,
)


def test_normalize_derives_standard_prime_mixmhcpred_and_tumor_vaf_evidence():
    frame = pd.DataFrame(
        {
            "patient_id": ["p1"],
            "mutation_id": ["m1"],
            "mutant_peptide": ["ACDEFGHIK"],
            "hla_allele": ["HLA-A*02:01"],
            "tumor_vaf": [0.25],
            "mixmhcpred_rank": [10.0],
            "prime_rank": [20.0],
        }
    )
    out = normalize_product_candidates(frame)
    assert out.loc[0, "dna_vaf"] == pytest.approx(0.25)
    assert out.loc[0, "presentation_score"] == pytest.approx(0.90)
    assert out.loc[0, "recognition_score"] == pytest.approx(0.80)


def test_product_default_selects_at_most_one_route_per_mutation():
    assert InferenceConfig().max_per_mutation == 1

    candidates = normalize_product_candidates(_pvac_fixture(), patient_id="demo")
    selected = score_product_candidates(candidates)
    selected = selected[selected["selected"]]

    assert selected.groupby("mutation_id").size().max() == 1


def _pvac_fixture() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "Gene Name": ["KRAS", "KRAS", "TP53", "BRAF"],
            "Transcript": ["ENST1", "ENST1", "ENST2", "ENST3"],
            "Mutation": ["G12D", "G12D", "R175H", "V600E"],
            "HLA Allele": ["HLA-A*02:01", "HLA-B*07:02", "HLA-A*02:01", "HLA-A*03:01"],
            "MT Epitope Seq": ["KLVVVGADGV", "VVVGADGVGK", "HMTEVVRHC", "KIGDFGLATEK"],
            "WT Epitope Seq": ["KLVVVGAGGV", "VVVGAGGVGK", "HMTEVVRRC", "KIGDFGLATVK"],
            "Best MT Score": [30.0, 120.0, 800.0, 60.0],
            "Corresponding WT Score": [900.0, 140.0, 500.0, 400.0],
            "Tumor DNA VAF": [0.42, 0.42, 0.08, 0.3],
            "Tumor RNA Depth": [80, 80, 50, 100],
            "Tumor RNA VAF": [0.38, 0.38, 0.0, 0.27],
            "Tumor RNA Alt Read Support": [30, 30, 0, 28],
            "Gene Expression": [40.0, 40.0, 0.0, 25.0],
        }
    )


def test_normalize_product_candidates_has_no_research_label_requirement():
    normalized = normalize_product_candidates(_pvac_fixture(), patient_id="demo")
    report = validate_candidate_contract(normalized)
    assert report.ok
    assert "label" not in normalized
    assert normalized["candidate_id"].is_unique
    assert normalized.loc[0, "hla_allele"] == "HLA-A*02:01"


def test_invalid_peptide_is_rejected():
    source = _pvac_fixture().iloc[:1].copy()
    source.loc[0, "MT Epitope Seq"] = "NOT-A-PEPTIDE"
    with pytest.raises(ValueError, match="product contract"):
        normalize_product_candidates(source, patient_id="demo")


def test_rna_evidence_merges_at_mutation_key():
    candidates = normalize_product_candidates(
        _pvac_fixture().drop(columns=["Gene Expression", "Tumor RNA VAF"]),
        patient_id="demo",
    )
    evidence = pd.DataFrame(
        {
            "patient_id": ["demo", "demo", "demo"],
            "mutation_id": ["G12D", "R175H", "V600E"],
            "expression_tpm": [50.0, 0.0, 30.0],
            "rna_vaf": [0.4, 0.0, 0.3],
        }
    )
    merged = merge_rna_evidence(candidates, evidence)
    assert merged.loc[merged.mutation_id == "G12D", "expression_tpm"].eq(50.0).all()


def test_scoring_excludes_no_expression_and_selects_deterministically():
    candidates = normalize_product_candidates(_pvac_fixture(), patient_id="demo")
    config = InferenceConfig(k=2, max_per_mutation=1)
    first = score_product_candidates(candidates, config)
    second = score_product_candidates(candidates.sample(frac=1, random_state=4), config)

    assert validate_ranked_contract(first).ok
    assert first.loc[first.mutation_id == "R175H", "exclusion_reason"].iloc[0] == "NO_RNA_EXPRESSION"
    assert first["selected"].sum() == 2
    assert set(first.loc[first.selected, "candidate_id"]) == set(second.loc[second.selected, "candidate_id"])
    assert first.loc[first.selected, "mutation_id"].nunique() == 2


def test_missing_recognition_is_visible_and_increases_uncertainty():
    candidates = normalize_product_candidates(
        _pvac_fixture().drop(columns=["Corresponding WT Score"]), patient_id="demo"
    )
    scored = score_product_candidates(candidates)
    row = scored.iloc[0]
    assert not bool(row["recognized_evidence_available"])
    assert "MISSING=RECOGNIZED" in row["selection_reason"]
    assert row["epicurus_neo_lower_evidence_score"] < row["epicurus_neo_evidence_score"]


def test_run_product_inference_writes_three_reports(tmp_path: Path):
    source = tmp_path / "pvac.tsv"
    _pvac_fixture().to_csv(source, sep="\t", index=False)
    outputs = run_product_inference(source, tmp_path / "report", patient_id="demo")

    assert set(outputs) == {"csv", "json", "markdown"}
    assert all(Path(path).exists() for path in outputs.values())
    payload = json.loads(Path(outputs["json"]).read_text())
    assert payload["patients"][0]["patient_id"] == "demo"
    assert "not_a_validated_response_probability" in payload["policy"]
    markdown = Path(outputs["markdown"]).read_text()
    assert "not validated probabilities" in markdown


def test_patient_abstains_when_no_core_candidates():
    source = _pvac_fixture().iloc[:1].copy()
    source["Gene Expression"] = 0.01
    source["Best MT Score"] = 10000.0
    candidates = normalize_product_candidates(source, patient_id="demo")
    scored = score_product_candidates(candidates)
    summary = patient_summaries(scored)[0]
    assert summary.abstained
    assert summary.abstention_reason == "NO_CANDIDATE_CLEARS_CORE_EVIDENCE_POLICY"


def test_class_i_12mer_is_not_dropped_but_impossible_lengths_are():
    # A validated immunogenic class-I 12mer exists in the Müller NCI cohort (APARLERRHSAL/B0702);
    # class I presents 8-14mers, so the deterministic validity gate must NOT delete a 12mer, while
    # still removing genuinely impossible lengths (7mer, 15mer+).
    frame = pd.DataFrame({
        "mutant_peptide": ["APARLERRHSAL", "SLYNTVATL", "AAAAAAA", "A" * 16],
        "hla_allele": ["B0702", "A0201", "A0201", "A0201"],
        "mhc_class": ["I", "I", "I", "I"],
    })
    gated = apply_deterministic_gate(frame)
    passed = dict(zip(gated["mutant_peptide"], gated["deterministic_gate_pass"]))
    assert passed["APARLERRHSAL"] is True or passed["APARLERRHSAL"]  # 12mer preserved
    assert passed["SLYNTVATL"]                                        # 9mer preserved
    assert not passed["AAAAAAA"]                                      # 7mer removed
    assert not passed["A" * 16]                                       # 16mer removed


def test_validated_vendor_calls_gate_lost_hla_and_unexpressed_routes():
    source = pd.DataFrame(
        {
            "Gene Symbol": ["ERBB2", "SMARCA4", "MAGEA6", "KRAS"],
            "Variant": ["v1", "v2", "v3", "v4"],
            "Protein Variant": ["p.E717Q", "p.R100Q", "p.A10V", "p.G12D"],
            "Source Variant Type": ["SNV", "SNV", "SNV", "SNV"],
            "HLA": ["A0201", "C0401", "A0201", "A0201"],
            "HLA LOH": ["N", "Y", "N", "N"],
            "Expressed": ["Y", "Y", "N", "Y"],
            "Peptide": ["LLQETELVE", "AAQAAAAAA", "AAVAAAAAA", "VVVGADGVG"],
            "Gene Level Expression TPM": [5420.0, 338.0, 0.0, 20.0],
            "SHERPA Presentation Rank": [0.13, 0.02, 0.03, 1.0],
        }
    )
    candidates = normalize_product_candidates(source, patient_id="sij")
    gated = apply_deterministic_gate(candidates)
    summary = summarize_gate(gated)

    assert summary.input_count == 4
    assert summary.survivor_count == 2
    assert summary.reason_counts == {
        "HLA_LOH_LOST_ALLELE": 1,
        "GENE_NOT_EXPRESSED": 1,
    }
    assert gated.loc[gated.gene_symbol == "ERBB2", "deterministic_gate_pass"].iloc[0]
    assert not gated.loc[gated.gene_symbol == "SMARCA4", "deterministic_gate_pass"].iloc[0]
    assert not gated.loc[gated.gene_symbol == "MAGEA6", "deterministic_gate_pass"].iloc[0]


def test_main_pipeline_runs_validity_gate_before_ranking(tmp_path: Path):
    source = pd.DataFrame(
        {
            "Gene Symbol": ["ERBB2", "SMARCA4"],
            "Variant": ["v1", "v2"],
            "HLA": ["A0201", "C0401"],
            "HLA LOH": ["N", "Y"],
            "Expressed": ["Y", "Y"],
            "Peptide": ["LLQETELVE", "AAQAAAAAA"],
            "Gene Level Expression TPM": [5420.0, 338.0],
            "SHERPA Presentation Rank": [0.13, 0.01],
        }
    )
    input_path = tmp_path / "sherpa.tsv"
    source.to_csv(input_path, sep="\t", index=False)
    outputs = run_product_inference(input_path, tmp_path / "report", patient_id="sij")

    ranked = pd.read_csv(outputs["csv"])
    removed = ranked[ranked.gene_symbol == "SMARCA4"].iloc[0]
    assert not bool(removed["selected"])
    assert removed["exclusion_reason"] == "HLA_LOH_LOST_ALLELE"
    payload = json.loads(Path(outputs["json"]).read_text())
    assert payload["deterministic_gate"]["removed_count"] == 1
    assert "HLA_LOH_LOST_ALLELE" in Path(outputs["markdown"]).read_text()
