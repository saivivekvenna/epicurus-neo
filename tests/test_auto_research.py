import json
from pathlib import Path

import pandas as pd

from epicurus_neo.auto_research import (
    build_failure_report,
    make_hypothesis_prompt,
    write_research_artifacts,
)


def _scored() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "candidate_id": ["a", "b", "c", "d"],
            "patient_id": ["p1", "p1", "p1", "p1"],
            "label": ["negative", "positive", "negative", "positive"],
            "epicurus_score": [0.95, 0.8, 0.7, 0.1],
            "mutant_peptide": ["AAA", "BBB", "CCC", "DDD"],
            "hla_allele": ["HLA-A*02:01"] * 4,
            "gene_symbol": ["KRAS", "TP53", "EGFR", "BRAF"],
            "mutation_tcr_face_count": [0.0, 1.0, 0.0, 2.0],
        }
    )


def test_build_failure_report_finds_top_k_errors():
    report = build_failure_report(_scored(), k=2)
    assert report.positives_total == 2
    assert report.positives_missed_at_k == 1
    assert report.negatives_in_top_k == 1
    assert report.false_negatives[0].candidate_id == "d"
    assert report.false_positives[0].candidate_id == "a"
    assert "false_positive_top_k" in report.numeric_feature_means


def test_make_hypothesis_prompt_blocks_bad_behavior():
    prompt = make_hypothesis_prompt(build_failure_report(_scored(), k=2))
    assert "Do not propose LLM direct peptide scoring" in prompt
    assert "Return YAML" in prompt
    assert "failure report" in prompt.lower()


def test_write_research_artifacts(tmp_path: Path):
    report = build_failure_report(_scored(), k=2)
    report_path, prompt_path = write_research_artifacts(report, output_dir=tmp_path)
    assert json.loads(report_path.read_text())["positives_total"] == 2
    assert "hypotheses" in prompt_path.read_text()

