# osteosarc.com (Sid) public reconstruction — data-asset preregistration

**Design source of truth for the osteosarc.com/Sid recognition-cohort RECONSTRUCTION.** This is a *pre-fit
data-asset design*, not a modeling experiment. It freezes the canonical schemas, the label/resolution state
machines, the deterministic dedup keys, and the fail-fast integrity checks **before** any table is built and
**before** any model is fit, tuned, or compared. Paired frozen copy:
`artifacts/milestone_7_decision/osteosarc_sid_reconstruction/PREREGISTERED_PROTOCOL.md`.

> **No fitting, tuning, ranking evaluation, or model comparison is licensed until the reconstructed ledger is
> committed and reviewed.** The frozen Epicurus model is NOT touched by this milestone. The recognized-target
> AUROC / hits@k numbers from commit `dd3efd1` are **descriptive and superseded** — see §11.

North Star (unchanged): Epicurus must maximize the chance that a genuinely recognized neoantigen survives
`WES/RNA/HLA → candidate generation → eligibility → top-20`. PRIME is a comparator/component, never a recall
mechanism: it cannot recover a candidate filtered before ranking. This reconstruction measures **where recognized
targets are actually lost in that funnel** on the one public patient for whom we hold both the full input and a
measured recognition label.

---

## 0. What this corrects (the `dd3efd1` diagnostic)

Commit `dd3efd1` (`scripts/osteosarc_rank.py`) was a useful preliminary probe but is **not** a valid benchmark:

1. It treated the **21-mutation `pvactools_curated_aggregated.tsv`** as the candidate universe and the other 20
   curated mutations as **assumed tested-negatives**. They are not established negatives — most were never
   assayed at all.
2. It reported a **single-positive AUROC** (`(n-rank)/(n-1)`) on n=1 in-universe positive (DYNC1H1). That is a
   descriptive placement, not a benchmark.
3. It **never touched the public osteosarc.com site**, which carries a far richer, real assay record: **182
   somatic variants, 44 vaccine-targeted variants, 14 site ELISPOT-positive variants**, per-variant peptide
   blocks, and experiment-level ELISPOT tables with explicit Positive/Weak/Strong/Negative results — a separate
   label stream from the Hudson IFNγ/TCR expander tables.

This reconstruction builds the evidence-graded ledger the diagnostic lacked, keeping the two label streams
separate and never inventing negatives.

---

## 1. Sources (exact, provenance-hashed at build time)

**Runtime (network, cached under gitignored `data/raw/osteosarc/site_cache/`, URL+content SHA256 recorded):**

| id | URL | shape | role |
|---|---|---|---|
| `site_index` | `https://osteosarc.com/variants/` | 182 `<tr data-*>` rows | variant catalog + headline invariants |
| `site_variant[slug]` | `https://osteosarc.com/variant/<GENE-chrN-POS>/` (trailing slash; bare path 308-redirects) | 182 pages | peptide blocks, ELISPOT experiments, detection, predicted peptides |
| `vafs_long` | `https://osteosarc.com/variants/variant_vafs_long.tsv` | 7,600 data rows, 26 cols | per-sample VAF long table (coordinate-level reachability join) |
| `vafs_cols` | `https://osteosarc.com/variants/variant_vafs_long.columns.tsv` | data dictionary | column semantics (committed verbatim) |

**Local public inputs (already staged in `data/raw/osteosarc/`, gitignored via `data/raw/*`):**

| file | shape | role |
|---|---|---|
| `pvactools_all_epitopes.tsv` | ~14,780 peptide×HLA | 2025.01 (T2) candidate universe, full ensemble scores |
| `pvactools_curated_aggregated.tsv` | 21 mutations | curated shortlist (best peptide, tier, VAF) |
| `rsem.2025.01.genes.results` | gene TPM | matched T2 expression |
| `May_all_expanders.tsv`, `Aug_all_expanders.tsv` | 3,581 / 1,789 TCR rows | Hudson IFNγ expander mutation-level TCR label |
| `bucket_manifest.txt` | 262k-line listing | public `b2.osteosarc.com` file index (locates `{May,Aug}_bys_all_expanders.tsv`) |

`{May,Aug}_bys_all_expanders.tsv` (bystander-pool controls) live at
`https://b2.osteosarc.com/hudson_lab/peptide_expansion/analysis/` and MAY be staged in Phase B; if fetched they
are cached and hashed identically. Absence of a bystander expander is **never** a tested negative.

Every fetch uses a stable descriptive User-Agent, retry-with-backoff, and writes an entry to `PROVENANCE.json`
(`url`, `http_status`, `bytes`, `sha256`, `fetched_at`). Phase B must be **rerunnable network-free from cache**.

---

## 2. Canonical schemas (frozen)

Five committed CSVs. Column order is frozen. `variant_id` = `<GENE>-<chrom>-<pos>` (GRCh38, chr-prefixed,
1-based), the site's own stable key; it is the sole cross-table join key at the variant level, and
`(chrom,pos,ref,alt)` is the coordinate key for reachability. All string fields preserve the site's raw text;
derived/normalized fields are additive and separately named.

### 2.1 `variant_catalog.csv` — one row per site variant (expect 182)
`variant_id, gene, chrom, pos, chr_pos, protein_change_raw, consequence, variant_type, ref, alt, change,
n_vaccines, n_pipelines, elispot_positive_flag, tumor_vaf, detected_pipelines, notdetected_pipelines,
vaccines_targeting, elispot_summary, source_url, source_sha256`

- `n_vaccines`, `n_pipelines`, `elispot_positive_flag` come from index `data-vaccines`, `data-pipelines`,
  `data-elispot` (0/1) and are **cross-checked** against the variant page's pill rows (fail on mismatch).
- `detected_pipelines`/`notdetected_pipelines` = `pill-on`(`✓`)/`pill-off` labels from the Detection section
  (e.g. `LENS 2022; Mutect2 2025; pVACtools 2025; oncoanalyser; DRAGEN`).
- `elispot_summary` = the `pill-elispot` text on the variant page (`strong`/`weak`/blank) when present.

### 2.2 `peptide_inventory.csv` — one row per peptide block (long vaccine peptide) on a variant page
`variant_id, gene, block_index, sources, peptide_seq, minimal_epitope, peptide_len, declared_experiment_count,
parsed_experiment_count, source_url`

- `sources` = all `span.pep-source` tags joined (`mRNA; JLF V2; JLF V3; CeGaT`), i.e. which vaccines/designs
  include this long peptide.
- `peptide_seq` = concatenation of all `code.peptide-seq` spans (raw residues, no highlight markup).
- `minimal_epitope` = the highlighted `span.ep` residues (the mRNA minimal epitope) or blank.
- `peptide_len` = residue count of `peptide_seq`; **must equal** the `NN aa` in `peptide-aux` (fail on mismatch).
- `declared_experiment_count` = the `MM experiments` integer in `peptide-aux`; `parsed_experiment_count` = rows
  actually parsed from this block's `exp-table`. **Integrity check: these must be equal per block** (§4).

### 2.3 `assay_ledger.csv` — one row per ELISPOT experiment row (site vaccine-peptide assay stream)
`experiment_key, variant_id, gene, block_index, peptide_seq, minimal_epitope, jlf_peptide_id, exp_date,
experiment_name, pool_raw, result_raw, result_class, label_state, resolution_state, notes_raw, source_url`

- `label_state` ∈ {`POSITIVE_STRONG`,`POSITIVE_WEAK`,`POSITIVE`,`NEGATIVE`,`AMBIGUOUS`,`UNTESTED`} (§3).
- `resolution_state` ∈ {`INDIVIDUAL_PEPTIDE`,`MUTATION_LONG_PEPTIDE`,`POOL`,`MUTATION_TCR`,`UNKNOWN`} (§3).
- `result_class` = the raw CSS class (`res-strong`/`res-weak`/`res-neg`); `result_raw` = the raw cell text
  (`Positive (Strong**)`, `Negative`, …). **Both are preserved** — the class and text can disagree and neither
  is discarded.
- `notes_raw` = all `td-notes` div text joined with ` ¦ ` (preserves sample-timeline and cross-references).
- A peptide block with `declared_experiment_count == 0` emits **one** ledger row with `label_state=UNTESTED`,
  `resolution_state=UNKNOWN`, empty experiment fields (vaccine inclusion without any experiment — never NEGATIVE).

### 2.4 `hudson_tcr_labels.csv` — Hudson IFNγ expander mutation-level TCR stream (SEPARATE label modality)
`timepoint, mutation_label, gene, protein_change, is_mutation_specific, trb, log2fc_umi, fold_expansion,
baseline_pct, poststim_pct, pool_kind, label_state, resolution_state, source_file, source_sha256`

- One row per (timepoint, TRB, mutation) observation; `resolution_state=MUTATION_TCR` always.
- `pool_kind` ∈ {`tumor`,`bystander`} from the `all` vs `bys` file. Bystander rows are controls.
- `label_state` for a (timepoint, mutation) is `POSITIVE` iff any `is_mutation_specific==TRUE` clonotype exists
  for it; otherwise `UNTESTED` for that mutation-timepoint (absence of an expander is NOT a tested negative
  without protocol proof). Per-TRB rows keep their own `is_mutation_specific`.
- This table is **never merged** into `assay_ledger.csv`: different assay modality (peptide ELISPOT vs
  IFNγ/TCR expansion), different resolution semantics.

### 2.5 `reachability_funnel.csv` — per recognized/vaccine target, its survival through the funnel
`target_id, gene, protein_change, chrom, pos, ref, alt, recognized_by, in_vafs_long, detected_pipelines,
in_pvactools_2025, in_curated_21, in_vaccine, has_site_elispot, site_elispot_best, hudson_recognized,
gene_tpm, first_failure_stage, adjudication`

- Rows for: every Hudson-recognized mutation, every site ELISPOT-positive variant, and the 21 curated mutations.
- Funnel stages, in order: `vafs_long_detected → pvactools_2025_candidate → curated_21 → vaccine_included →
  site_elispot_tested → recognized`. `first_failure_stage` = the first stage a recognized target fails to reach.
- Joins use `(chrom,pos,ref,alt)` or exact normalized `protein_change`, **never gene alone** (MAP2 and DYNC1H1
  each have two distinct coordinates on the site — a gene-only join is provably wrong).

---

## 3. State machines (frozen mappings)

### 3.1 `label_state` from raw ELISPOT result (site stream)
Mapping is on the **raw result text** (case/whitespace-normalized), with `result_class` preserved separately:

| raw text (normalized) | label_state |
|---|---|
| `positive (strong)`, `positive (strong**)`, any `positive (strong…)` | `POSITIVE_STRONG` |
| `positive (weak)`, any `positive (weak…)` | `POSITIVE_WEAK` |
| `positive` (unqualified) | `POSITIVE` |
| `negative` | `NEGATIVE` |
| anything else / unparseable / conflicting text-vs-class not in the table | `AMBIGUOUS` |
| (no experiment row exists for a vaccine-included peptide) | `UNTESTED` |

Footnote markers (`**`) and any parenthetical qualifier are preserved in `result_raw`. A plain `Positive` with
`result_class=res-strong` maps to `POSITIVE` (text governs the state; the class is retained, not overwritten).

### 3.2 `resolution_state` (what was actually tested) — inferred only from explicit evidence
| evidence | resolution_state |
|---|---|
| `pool_raw` is a pool id (`S1..S4`, `L1..L3`, `Pool N`, or any non-blank non-`NA` pool token) | `POOL` |
| experiment tests the long vaccine peptide (peptide block sequence), pool blank, no individual-peptide proof | `MUTATION_LONG_PEPTIDE` |
| notes/JLF-id/table explicitly establish a single minimal-epitope (e.g. class-I 9-mer) individual test | `INDIVIDUAL_PEPTIDE` |
| Hudson IFNγ/TCR expander observation | `MUTATION_TCR` |
| `pool_raw` ∈ {`—`,`NA`,blank} with no other proof | `UNKNOWN` |

**`—`/`NA` do NOT prove individual-peptide testing.** Default for a blank pool is `MUTATION_LONG_PEPTIDE` when the
row is attached to a long-peptide block and the notes do not prove individual testing, else `UNKNOWN`. **A
pool-positive does not make any member peptide an exact individual positive** — pool positivity is recorded at
`resolution_state=POOL` and never propagated to members.

---

## 4. Deterministic dedup keys & fail-fast integrity checks

**Dedup keys (canonical, order-independent):**
- assay ledger `experiment_key` = SHA1 of `variant_id | block_index | jlf_peptide_id | exp_date |
  experiment_name | pool_raw | result_raw`. Must be **unique** across the ledger.
- peptide inventory key = `(variant_id, block_index)`; peptide blocks are indexed in document order.
- hudson key = `(timepoint, pool_kind, trb, mutation_label)`. Must be unique.
- reachability key = `target_id` = `variant_id` (or `gene-chrom-pos-ref-alt` for Hudson-only targets absent
  from the site).

**Fail-fast checks (build aborts, no partial write, if any fail):**
1. `variant_catalog` has exactly **182** rows; index and per-page gene/coords agree.
2. `sum(n_vaccines>0) == 44`; `sum(elispot_positive_flag) == 14` (site headline invariants).
3. For every peptide block: `parsed_experiment_count == declared_experiment_count` and
   `peptide_len == aa_count_in_aux`.
4. `assay_ledger.experiment_key` and `hudson_tcr_labels` key are duplicate-free.
5. Every `label_state`/`resolution_state` value is in the frozen enum.
6. Reachability joins resolve by coordinate/protein, never gene-only; MAP2 and DYNC1H1 dual coordinates are
   disambiguated (assert each recognized target maps to exactly one coordinate).
7. All three Hudson-recognized mutations (`ASPM p.G2179R`, `DYNC1H1 p.V314I`, `MAP2 p.…868fs`) appear in
   `reachability_funnel` with a non-null `first_failure_stage`.

Contradictions (same peptide/mutation, differing result across timepoint/protocol/pool) are **reported**
(`AUDIT.json.contradictions`, keyed by peptide/mutation/timepoint/protocol) and **never overwritten or
collapsed**. Repeated protocols on the same peptide are all retained as distinct ledger rows.

---

## 5. Targeted correctness tests (Phase B, `tests/test_osteosarc_sid.py`)

Pinned against live-verified fixtures (cached HTML), network-free:
- **Invariants:** 182 / 44 / 14 reproduce exactly.
- **`EXOC4-chr7-133274996`:** peptide block `KKSVIRTLSTIDDVEDRENEKGR` (23 aa, sources `JLF V3`+`JLF V2`,
  4 experiments) parses 4 rows; a `Pool 1`/`res-weak` `Positive (Weak)` row → `POSITIVE_WEAK`/`POOL`;
  a `UNC…` `res-strong` `Positive (Strong)` row → `POSITIVE_STRONG`; blank-pool rows stay `UNKNOWN`/long-peptide.
- **`MAP2` (both coordinates):** a peptide block with `0 experiments` → single `UNTESTED` row; minimal-epitope
  `span.ep` extracted; the two MAP2 coordinates are distinct rows.
- **`MYO15B`, `SLC25A12`:** parse both `Negative` (`res-neg`) and positive rows; repeated protocols retained;
  minimal vs long peptide distinguished.
- **Label mapping:** every string in the §3.1 table maps to the stated state; unknown text → `AMBIGUOUS`.
- **No-fabrication:** absence of an experiment/expander never yields `NEGATIVE`.
- **Determinism:** two runs from cache produce byte-identical CSVs (sorted, fixed float formatting).

---

## 6. Deliverables (Phase B)

Code: `src/event_b/osteosarc_sid.py` (stdlib `html.parser`/`urllib`, **no new dependency**),
`scripts/osteosarc_sid_reconstruct.py`, `tests/test_osteosarc_sid.py`.
Artifacts under `artifacts/milestone_7_decision/osteosarc_sid_reconstruction/`:
`variant_catalog.csv, peptide_inventory.csv, assay_ledger.csv, hudson_tcr_labels.csv, reachability_funnel.csv,
AUDIT.json, REPORT.md, PROVENANCE.json`. Committed as a **separate scoped commit** from this design.

---

## 7. Preregistered report questions (answered in `REPORT.md`, not before the ledger is built)

1. Unique variants; peptide blocks; exact **individual** tests; **pool** tests; **mutation-level** (Hudson)
   tests; positives (by strength); **defensible** negatives (individual/mutation only, never pool/untested);
   untested vaccine candidates; contradictions.
2. Do the 14 site-ELISPOT-positive variants overlap the 3 Hudson positives? Report overlap **without forcing
   equivalence** across assay modalities.
3. The true **evidence-supported denominator** for a within-patient ranking evaluation today.
4. The exact reachability **stage** at which each recognized target is lost.
5. Which Epicurus changes the evidence justifies (longitudinal/multi-caller union, no hard TPM drop, evidence
   tiers, abstention, diversity/uncertainty top-20) — proposed, **not fit** to this patient.

---

## 8. Out of scope / guardrails

No model fit, tune, retrain, threshold search, or ACCEPT/REJECT ranking gate in this milestone. The frozen
Epicurus config (`configs/frozen/epicurus_v0_1.json`) is not read or modified. Any improvement idea surfaced in
§7.5 must be validated later on the independent labeled cohorts (multimer/Gartner/IMPROVE/CheckMate), never on
this single patient. All emitted tables are derived public data and are committable.

## 9. Registered deviations
_(none yet — appended chronologically as Phase B proceeds)_
