# Miller generalization — reconstruction toolchain (pinned, label-blind)

**Scope:** unblock HLA typing and mutanome enumeration for arbitrary Miller IPV patients
(Hu_315 first) using the *exact* pinned providers that reconstructed Hu_287 — no ad-hoc
substitutes. LOCKED_TEST isolation holds: no recognition label is read for download, HLA,
expression, calling, generation, thresholds, or ranking.

## The problem

`benchmark.miller_download.reconstruction_stages()` resolved tools by `PATH` (+ a couple of
gitignored local binaries). But the Hu_287 reconstruction never used PATH tools for HLA or
enumeration:

- **HLA** was typed by **OptiType inside the pinned micromamba `hla` env**
  (`data/raw/tools/micromamba/envs/hla`: optitype 1.3.5 + razers3 + GLPK), driven by
  `scripts/miller_hu287_hla.sh` via `micromamba run -n hla`.
- **Mutanome enumeration** used the **Ensembl VEP REST endpoint** (not a local VEP cache):
  `event_b.lossless_peptide_generation` classifies each variant's consequence and enumerates
  the reference-protein peptide universe, caching every Ensembl response per-patient under
  `ensembl_cache/`; `bcftools norm` left-aligns/splits indels against GRCh38 first.

Because the stage map looked for `OptiTypePipeline.py`/`razers3`/`vep` on `PATH` and for a
local VEP cache + GENCODE GTF, both stages read a **false `NOT_EVALUABLE`** on any machine
whose pinned toolchain is exactly the one that produced Hu_287.

## The fix (scoped)

`src/benchmark/miller_download.py`:

1. **`pinned_env_tool(name)`** — resolves a tool by the existence of
   `data/raw/tools/micromamba/envs/<env>/bin/<tool>` for the pinned envs
   (`OptiTypePipeline.py`,`razers3`→`hla`; `vep`→`vep`). Fail-closed: returns `None` unless
   both the pinned micromamba launcher **and** the env's `bin/<tool>` exist.
2. **`resolve_tool`** now falls through `PATH → local install → pinned micromamba env`.
3. **HLA `OptiType` method** drops the phantom `hla_reference_dna.fasta` sentinel — OptiType
   bundles its own HLA reference, so the pinned `hla` env is self-contained.
4. **Mutanome** gains the real provider **`VEP-REST+lossless`** (tools: `bcftools`; refs:
   GRCh38 fasta + `.fai`) ahead of the local-cache VEP/pvacseq methods, documented as the
   exact provider that produced the frozen Hu_287 universe. Needs network; no local VEP
   cache/GTF.

No change to the universal score, portfolio policy, or any frozen universe/product manifest
(`miller_download.py` is not a hashed dependency of those freezes).

## Result (real machine, 2026-07-13)

All reconstruction stages `RUNNABLE` except scoring, which stays correctly upstream-blocked
until the re-enumerated universe + class-I HLA exist:

| stage | status | provider |
|---|---|---|
| sra_to_fastq | RUNNABLE | fasterq-dump |
| hla_typing_classI | RUNNABLE | OptiType (pinned `hla` env) |
| wes_alignment | RUNNABLE | bwa-mem+samtools |
| somatic_calling | RUNNABLE | Mutect2 (gatk 4.5.0.0) |
| rna_quant | RUNNABLE | salmon |
| mutanome_enumeration | RUNNABLE | VEP-REST+lossless (Ensembl REST) |
| scoring_prime_epicurus | NOT_EVALUABLE | upstream-blocked (needs universe + HLA) |

## Runners

- **HLA:** `scripts/miller_patient_reconstruct.py <PATIENT_ID> hla` → parameterized
  `scripts/miller_hu287_hla.sh` (`PATIENT_ID` env), waits for the normal MD BAM, extracts the
  GRCh38 MHC region, runs OptiType in the `hla` env, writes `HLA_PROVENANCE.json` (per-file
  sha256 + env identity).
- **Enumeration + freeze:** `scripts/miller_patient_universe.py <PATIENT_ID>` →
  `benchmark.miller_universe_core.freeze` (REST enumeration + label-blind universe freeze;
  refuses Hu_287 to protect its legacy frozen provenance).

## Tests

`tests/test_miller_download.py`: `pinned_env_tool` resolution + fail-closed; `resolve_tool`
fall-through; HLA + mutanome RUNNABLE via the pinned toolchain; OptiType requires no external
FASTA. Bare-machine manifest tests neutralize the pinned launcher so injected `which` alone
governs availability.
