# External-validation lane — cohort audit & acquisition ranking

Goal: rank cohorts by ability to test a **frozen Epicurus formula vs MHCflurry / NetMHCpan / genuine PRIME
at within-patient top-20**, NOT by row count. Genuine PRIME 2.1 is installed, so any cohort with (peptide,
restricting-HLA) is PRIME-scoreable now. Labels are three-state (POSITIVE / genuine TESTED_NEGATIVE /
UNTESTED); UNTESTED is never coerced to negative.

## END-TO-END patient cohorts (per-patient candidate denominator → patient top-20 selectable)

| Cohort | loader | patients | POS / TESTED_NEG / UNTESTED | HLA | expr | VAF | tested-neg? | PRIME | status |
|---|---|--:|---|:--:|:--:|:--:|:--:|:--:|---|
| **Gartner NCI Testing** | `gartner_nci_corpus.load_gartner_nci` | 26 | 46 / 3,722 / 5,009 | ✅(acquired) | decile | decile | ✅ | tool | **USED as frozen transfer eval** |
| **IMPROVE (SRHgroup)** | `data/raw/improve/data.zip` | 70 | 467 / 17,053 / 0 | ✅ | ✅ | ✅ `VarAlFreq` +`CelPrev` | ✅ | precomp+tool | **USED as untouched external validation** |
| **CD8 multimer 2025** | `cd8_multimer_corpus.load_cd8_multimer` | 26 | 34 / 8,069 / 0 | ✅ | TPM | ❌ | ✅ | tool | **USED as frozen training cohort** |
| Müller/neoranking NCI | `neoranking_corpus.load_neoranking_nci` | 56 | 82 / **0** / 292,413 | ✅ | ✅ | ❌ | ❌ (PU-only) | overlap | same lineage as Gartner; PU top-k only |
| Sijbrandij n=1 (SHERPA) | `data/raw/sijbrandij/…SNV_Indel.tsv` | 1 | full universe, **no T-cell labels** | ✅ | TPM | ✅ | n/a | tool | product/demo input, NOT validation |

The three best independent local cohorts (Gartner, IMPROVE, multimer) are now allocated: multimer trains the
frozen residual; Gartner + IMPROVE are its two untouched external tests. **No further end-to-end cohort is
executable locally without an acquisition** (below).

## RECOGNITION-SCALE only (no per-patient denominator → cannot do a real top-20 today)
Zhao DC (352 pts, vaccinated subset only), CEDAR (assay rows, cross-study transfer already falsified 0.478),
TESLA (locked AUROC check), BigMHC (recognition rows), PRIME TableS4 (leakage registry), and all vaccine
Event-B cohorts (Braun, Hu, Ott, mKRAS, PDAC-NeoVax, NOUS-209). Useful for recognition training/leakage, not
for within-patient decision validation.

## Ranked acquisition targets (to add a NEW untouched end-to-end cohort)
Ranked by decision-benchmark value (full denominator + genuine tested negatives + HLA + expression/VAF +
independence), not row count:

1. **Hu + Ott NeoVax melanoma** → dbGaP **phs001451.v3.p1** (tumor/normal WES + tumor RNA-seq, patients
   1–6, 11–12). Highest value: true end-to-end vaccine patients. Needs regeneration of the candidate
   universe + HLA genotype + VAF + TPM (pVACtools/NetMHCpan). **Controlled access, no redistribution.**
   Destination: `data/raw/hu_melanoma_2021/regenerated/`, `data/raw/ott_melanoma_2017/regenerated/`.
   Schema: per-candidate {patient, mutant_peptide(25mer), hla_allele, expr_TPM, vaf, screening_label}.
2. **Müller full negatives** → figshare `Neopep_data_org.txt.zip` + `Mutation_data_org.txt.zip`
   (bassanilab/NeoRanking, links in `docs/download_checklist.md`). Would give richer per-candidate features;
   but same NCI lineage as Gartner (not independent). Destination `data/raw/neoranking/data/`.
3. **Fukuoka DC 2021** → 17-patient patient×peptide table (Table1-equivalent). Publisher endpoints 404 →
   **blocked**; patient-overlap risk with Zhao (same Morisaki program).
4. **Sijbrandij n=1** → add a matched T-cell assay label file to convert its full universe into validation.
5. **Osteosarc** → source files never supplied (`adapters/osteosarc.py` raises); no path today.

## Note on independence
After this milestone, Gartner/IMPROVE/multimer are "consumed" (one as training, two as frozen eval). A
*new* validation must come from acquisition target #1 (Hu/Ott) to keep the north-star claim honest —
"Epicurus beats incumbents" requires a win on patients none of the pipeline has touched.
