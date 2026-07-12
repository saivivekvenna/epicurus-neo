# Lossless peptide recovery — final status & next validation blocker

**Milestone-7 exploratory diagnostic. Post-hoc, not preregistered/blind/independent.**

## What was done

Implemented and ran the frozen input-only lossless peptide generator
(`lossless-peptide-generation-1.0.0`, composed with the frozen router policy
`epicurus-evidence-router-1.0.0`) on osteosarc.com / Sid:

`raw GRCh38 allele → Ensembl VEP MANE/canonical → Ensembl CDS/protein → all standard-AA class-I 8–14
windows spanning the mutant (missense) or novel-frame (frameshift) residue → patient HLA panel read
from pVAC → genuine GfellerLab PRIME 2.1 → union with the original pVAC candidate set (PRIME from the
frozen `_cache_prime.tsv`) on a stable genomic candidate identity → frozen evidence-router route-aware
top-20, selected by `genuine_prime = −PRIME %rank`.`

The generator and runner read **only** the raw allele, the Ensembl reference, the pVAC HLA panel, the
RSEM quant, and genuine PRIME. They import/read no functional-assay, therapeutic, or measured-outcome
table (enforced by `tests/test_lossless_peptide_generation.py`). The three IFNγ/TCR-expanded target
labels are joined **only after** ranking, isolated as evaluation.

## Results (reproduced the disclosed feasibility exactly)

| Variant | Role | Windows | Best genuine-PRIME %rank |
|---|---|---:|---:|
| ASPM p.Gly2179Arg (`ENST00000367409` / `NM_018136.5`) | recovered | **77** | **0.088** |
| MAP2 p.Gly868AlafsTer38 (`ENST00000682079` / `NM_001375505.1`, `c.2603_2630del`) | recovered | **259** | **0.010** |
| DYNC1H1 p.Val314Ile (`ENST00000360184` / `NM_001376.5`) | input-only positive control | 77 | 0.002 |

The DYNC control reproduces the known pVAC MT epitopes (`KRFHATISF`, `ATISFDTDT`), confirming the
generator regenerates the incumbent pipeline's candidates from inputs alone.

Coverage of the 3 targets, pVAC-only vs augmented (pVAC + lossless recovery):

| Stage | pVAC-only | Augmented |
|---|---|---|
| candidate-generation recall | 1/3 | **3/3** |
| rankable recall (peptide+HLA) | 1/3 | **3/3** |
| pure genuine-PRIME top-20 mutation coverage | 1/3 | **3/3** |
| route-aware top-20 mutation coverage | 1/3 | **3/3** |

In the augmented route-aware top-20, the recovered MAP2 (`ALMIKFEEI` #4, `RVVPFTKAL` #10) and ASPM
(`RRVRVRRTL` #9) candidates carry `candidate_source = lossless_recovery`; DYNC1H1 (`KRFHATISF`) is the
pVAC incumbent. **The score that placed them is genuine PRIME itself.**

## Reproducibility

Online and offline runs produce an **identical** mode-invariant content hash
`74bbd8fd6947275d61c51b0ab3f8e88d5f6b133026c667719871e21ec6e183e1` and byte-identical
`RECOVERED_CANDIDATES.csv`. Offline serves the gitignored Ensembl cache
(`data/raw/osteosarc/ensembl_cache/`) and **fails closed** on a cache miss. Per-response Ensembl
URL+SHA-256 are recorded in `PROVENANCE.json` (e.g. MAP2 VEP `f523839e…`, CDS `5e978372…`, protein
`c33179dd…`).

## What this is NOT

This does **not** show Epicurus beats PRIME. It is a *reachability* fix: better candidate **generation**
lets genuine PRIME score targets it previously never received. The design was informed by Sid's audit
and the feasibility scores were inspected before the protocol freeze, so this is a **post-hoc
diagnostic on the one patient that motivated it**, not blind/independent/prospective evidence. No model
was fit or tuned; the frozen Epicurus v0.1 config and the frozen router policy are read-only.

## Next untouched-validation blocker (open)

A superiority/benefit claim requires running this **frozen** generator on a patient whose recognition
labels were **not** used to design or tune it — a genuinely untouched cohort with:

1. a raw multi-caller callset (raw GRCh38 alleles), 2. WES/RNA + class-I HLA typing, 3. lossless
input-only generation (this generator, unchanged), and 4. measured recognition outcomes joined only
after ranking.

None is in hand. Highest-value routes to obtain one (from the acquisition brief and the Sid ledger
entry):

- **More Hudson-lab / RTTP patients** with the IFNγ peptide-expansion assay **plus the stimulation-pool
  composition** (the true recognition denominator) → converts n=3 descriptive into a deployable-patient
  benchmark and gives this generator a prospective test.
- **Miller IPV `PRJNA980652`** (label table behind the STM paywall → needs user download) — dense
  denominator, well-powered.
- **CheckMate 153 dbGaP WES/RNA** to complete inputs for the already-run external cohort.

Until such a run clears, this remains an exploratory reachability diagnostic. See
`../../../NORTH_STAR_HISTORY.md` (M7 row "Lossless input-only peptide recovery") for the ledger entry.
