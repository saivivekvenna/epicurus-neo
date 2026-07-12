"""Reconstruct the public osteosarc.com (Sid) recognition record into an evidence-graded ledger.

Runs the frozen preregistration
(`docs/superpowers/specs/2026-07-12-osteosarc-sid-reconstruction-preregistration.md`): fetches the
182 variant pages + VAF TSVs (cached, hashed, network-free on rerun), parses them structurally with a
stdlib html.parser mini-DOM, joins the local public pVACtools/RSEM/Hudson inputs, and writes the five
canonical CSVs + AUDIT.json + REPORT.md + PROVENANCE.json under
`artifacts/milestone_7_decision/osteosarc_sid_reconstruction/`.

No model is fit, tuned, or compared. The frozen Epicurus config is not touched.

    .venv/bin/python -m scripts.osteosarc_sid_reconstruct            # fetch (cache-first) + build
    .venv/bin/python -m scripts.osteosarc_sid_reconstruct --offline  # rerun from cache only
    .venv/bin/python -m scripts.osteosarc_sid_reconstruct --refresh  # force re-fetch
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from event_b import osteosarc_sid as osid  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--offline", action="store_true", help="use cached pages only; fail if missing")
    ap.add_argument("--refresh", action="store_true", help="force re-fetch of all URLs")
    args = ap.parse_args()

    res = osid.build(offline=args.offline, refresh=args.refresh)
    s = res["summary"]

    print("\n============== osteosarc.com (Sid) — public reconstruction ==============")
    print(f"variants={s['unique_variants']}  vaccine-targeted={s['vaccine_targeted_variants']}  "
          f"site-ELISPOT-positive={s['site_elispot_positive_variants']}  (invariants 182/44/14 PASS)")
    print(f"peptide blocks={s['peptide_blocks']}  ledger rows={s['assay_ledger_rows']} "
          f"(tested={s['site_experiment_rows_tested']} + untested-vaccine={s['untested_vaccine_peptide_blocks']})")
    print(f"resolution: individual={s['individual_peptide_tests']} long-peptide={s['long_peptide_tests']} "
          f"pool={s['pool_tests']}")
    print(f"site positives: strong={s['positives_strong']} weak={s['positives_weak']} "
          f"plain={s['positives_unqualified']}  negatives={s['negatives_total']} "
          f"(defensible={s['defensible_negatives_individual_or_longpeptide']} pool-only={s['pool_negatives_not_perpeptide_defensible']})")
    print(f"Hudson TCR: {s['hudson_mutation_tcr_rows']} clonotype rows / "
          f"{s['hudson_distinct_timepoint_mutation_tests']} (tp,mut) tests; recognized={s['hudson_recognized_genes']}")
    print(f"contradictions={s['contradictions']}  site∩hudson overlap={s['site_vs_hudson_overlap_variants']}")
    print("\nrecognized-target reachability (automated funnel first-failure):")
    for gene, stage in sorted(s["recognized_first_failure_stage"].items()):
        print(f"   {gene:10s} -> {stage}")
    print(f"\nartifacts -> {res['out_dir']}   provenance URLs={res['n_provenance']}")
    print(json.dumps(s, indent=2, default=str)[:0] or "", end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
