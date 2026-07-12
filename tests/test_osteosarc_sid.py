from pathlib import Path

from event_b import osteosarc_sid as sid


ROOT = Path(__file__).resolve().parents[1]


def _peptide_html(aux: str, rows: str = "") -> str:
    return f"""
    <section><h2>Vaccine peptides &amp; ELISPOT experiments</h2>
      <section class="peptide-block">
        <header class="peptide-header">
          <span class="pep-source">JLF V1</span>
          <code class="peptide-seq"><span>ABCDE</span><span class="ep">FGHI</span></code>
          <div class="peptide-aux">9 aa · {aux}</div>
        </header>
        <table class="exp-table"><tbody>{rows}</tbody></table>
      </section>
    </section>
    """


def _experiment_row(result: str = "Negative") -> str:
    return f"""
      <tr><td>2024-05-15</td><td>UNC assay</td><td>JLF.1</td><td>NA</td>
      <td><span class="res res-neg">{result}</span></td><td><div>note</div></td></tr>
    """


def test_declared_experiment_count_accepts_singular_and_plural():
    singular = sid.parse_variant_page(_peptide_html("1 experiment", _experiment_row()), "V1")
    plural = sid.parse_variant_page(
        _peptide_html("2 experiments", _experiment_row() + _experiment_row("Positive")), "V2")
    assert singular["peptide_blocks"][0]["declared_experiment_count"] == 1
    assert len(singular["peptide_blocks"][0]["experiments"]) == 1
    assert plural["peptide_blocks"][0]["declared_experiment_count"] == 2
    assert len(plural["peptide_blocks"][0]["experiments"]) == 2


def test_zero_experiment_block_stays_untested_not_negative():
    page = sid.parse_variant_page(_peptide_html("0 experiments"), "V0")
    ledger = sid.build_assay_ledger({"V0": page})
    assert len(ledger) == 1
    assert ledger[0]["label_state"] == "UNTESTED"
    assert ledger[0]["resolution_state"] == "UNKNOWN"


def test_summary_keeps_resolved_rows_units_and_contradictions_distinct():
    base = {
        "variant_id": "V1", "peptide_seq": "PEPTIDE", "resolution_state": "MUTATION_LONG_PEPTIDE",
        "label_state": "POSITIVE_WEAK",
    }
    ledger = [base, {**base, "label_state": "NEGATIVE"}, {
        **base, "variant_id": "V2", "peptide_seq": "OTHER", "label_state": "POSITIVE_STRONG",
    }, {
        **base, "variant_id": "V3", "peptide_seq": "POOLED", "resolution_state": "POOL",
        "label_state": "POSITIVE",
    }]
    catalog = [
        {"variant_id": "V1", "n_vaccines": 1, "elispot_positive_flag": 1},
        {"variant_id": "V2", "n_vaccines": 1, "elispot_positive_flag": 1},
        {"variant_id": "V3", "n_vaccines": 1, "elispot_positive_flag": 1},
    ]
    summary = sid.summarize(catalog, [], ledger, [], [], [])
    assert summary["resolved_nonpool_rows"] == 3
    assert summary["resolved_nonpool_positive_rows"] == 2
    assert summary["resolved_nonpool_negative_rows"] == 1
    assert summary["resolved_unique_peptide_units"] == 2
    assert summary["resolved_unique_positive_units"] == 2
    assert summary["resolved_unique_negative_units"] == 1
    assert summary["resolved_unique_contradictory_units"] == 1
    # The pool positive remains visible globally but is excluded from the resolved denominator.
    assert summary["positives_unqualified"] == 1


def test_full_cached_reconstruction_passes_frozen_headline_invariants(tmp_path):
    cache = ROOT / "data/raw/osteosarc/site_cache"
    if not (cache / "variants_index.html").exists():
        return
    result = sid.build(offline=True, cache_dir=cache, out_dir=tmp_path)
    summary = result["summary"]
    assert summary["unique_variants"] == 182
    assert summary["vaccine_targeted_variants"] == 44
    assert summary["site_elispot_positive_variants"] == 14
    assert summary["resolved_nonpool_rows"] == 29
    assert summary["resolved_unique_peptide_units"] == 15
    assert summary["resolved_unique_contradictory_units"] == 8
    funnel = list(__import__("csv").DictReader((tmp_path / "reachability_funnel.csv").open()))
    hudson_map2 = next(row for row in funnel if row["gene"] == "MAP2" and row["hudson_recognized"] == "true")
    site_map2 = next(row for row in funnel if row["gene"] == "MAP2" and row["recognized_by"] == "site_elispot_positive")
    assert hudson_map2["target_id"] == "MAP2-chr2-209694772"
    assert hudson_map2["in_vaccine"] == "false"
    assert hudson_map2["has_site_elispot"] == "false"
    assert "differs from the Leu867fs vaccine neo-frame" in hudson_map2["adjudication"]
    assert site_map2["target_id"] == "MAP2-chr2-209694768"
