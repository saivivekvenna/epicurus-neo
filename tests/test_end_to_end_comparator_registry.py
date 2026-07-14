from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).parents[1]
REGISTRY = ROOT / "configs/frozen/end_to_end_comparator_registry_v1.json"


def test_registry_makes_nextneopi_headline_and_prime_component_only():
    registry = json.loads(REGISTRY.read_text())
    assert registry["labels_allowed"] is False
    assert registry["final_primary"]["comparator_id"] == "nextneopi_track_a"
    assert registry["final_primary"]["required_for_headline_verdict"] is True
    assert registry["component_diagnostic"] == {
        "comparator_id": "prime_plain",
        "scope": "COMPONENT_LEVEL_IDENTICAL_REACHABLE_UNIVERSE",
        "eligible_for_headline_verdict": False,
    }
    assert "TIES_NEXTNEOPI" in registry["final_verdicts"]
    assert "TIES_PRIME" not in registry["final_verdicts"]


def test_registry_requires_exact_variant_identity_and_explicit_attempts():
    registry = json.loads(REGISTRY.read_text())
    identity = registry["mutation_identity"]
    assert identity["canonical_key"] == "chrom:POS:REF:ALT"
    assert identity["allow_fuzzy_gene_or_peptide_match"] is False
    assert identity["allow_mnv_component_rescue"] is False
    assert registry["final_primary"]["missing_or_unverified_attempt"].startswith("PROTOCOL_FAILURE")
