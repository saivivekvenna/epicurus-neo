# osteosarc.com (Sid Sijbrandij) — recognition diagnostic

> First deployment-grade patient for whom we hold **both** the full North-Star input **and** a *measured*
> T-cell recognition label. Public data (osteosarc.com / Research to the People; `b2://osteosarc-data`,
> no DUA). This is a **method-shortcoming diagnostic, not a tuned model and not an ACCEPT-gate result** —
> n=3 measured positives (1 in the candidate universe) cannot fit anything. Frozen Epicurus v0.1 is applied
> unchanged. See `report.json` for machine-readable numbers, `per_mutation_scores.csv` for the full table.

**Not the same patient as RTTP `SR24-58221`.** That is a different individual (disjoint HLA: A\*11:01/A\*24:02…
via a DUA-gated Personalis report, **no** recognition label). The on-disk `artifacts/osteosarc_audit/` and
`artifacts/milestone_7_decision/osteosarc_product/` were built on the RTTP file mislabeled as "osteosarc" and
should be disregarded — this directory (`osteosarc_sid/`) is the real osteosarc.com patient.

## Inputs (all public)
- **Candidates:** pVACtools on the 2025.01 (T2 recurrence) tumor WGS — 21 curated somatic mutations,
  14,780 peptide×HLA rows, with a full ensemble already scored (NetMHCpan, NetMHCpanEL, MHCflurry EL,
  BigMHC EL/IM, DeepImmuno, …). HLA A\*01:01, B\*08:01, B\*27:05, C\*01:02, C\*07:01.
- **Expression:** matched RSEM gene TPM (same tumor/timepoint).
- **Measured label:** Hudson-Lab IFNγ peptide-expansion assay — PBMCs (May & Aug 2025) stimulated with
  tumor-mutation peptides → sorted IFNγ⁺/⁻ → TCR-seq (MiXCR); mutation-specific clonal expansion.

**Measured recognized neoantigens:** `ASPM p.G2179R` (May+Aug), `DYNC1H1 p.V314I` (May+Aug), `MAP2 …868fs` (May).

## Finding 1 — the wall here is RECALL, not ranking (recall = 1/3)

| Recognized | Timepoints | TPM (matched) | In pVACtools shortlist? | Why not |
|---|---|---:|---|---|
| **DYNC1H1 V314I** | May+Aug | **357.1** | ✅ yes (Tier=Pass) | — high-expression strong binder; curation kept it |
| **MAP2 …868fs** | May | **5.2** | ❌ no | frameshift *is* called at this timepoint but **low expression → dropped by the expression/tier filter** |
| **ASPM G2179R** | May+Aug | 16.5 | ❌ no | the specific variant is **not in the single 2025.01 callset** (only other ASPM variants are) — a **timepoint/heterogeneity miss** |

Two of three experimentally-recognized neoantigens never entered the 21-candidate shortlist. **No ranking
model can recover them** — they are filtered upstream. The two misses are *distinct, real* failure modes of
presentation-first, single-timepoint candidate selection: (a) expression thresholding discards a low-TPM but
genuinely immunogenic neoantigen (MAP2), and (b) a single tumor timepoint misses a neoantigen present
elsewhere in the tumor's history (ASPM).

## Finding 2 — within the shortlist, ranking is already saturated; "recognition" models hurt

Placement of the one in-universe recognized mutation (DYNC1H1) among the 21 curated candidates:

| Method | Rank of DYNC1H1 | 1-pos AUROC | Note |
|---|---:|---:|---|
| pVACtools presentation (Best MT %ile) | **1** | 1.00 | presentation ranked it top |
| NetMHCpan-EL %ile | **1** | 1.00 | presentation |
| **Epicurus v0.1 (ours)** | **1** | 1.00 | matches presentation (it is presentation-dominated) |
| genuine PRIME 2.1 | 2 | 0.95 | |
| MHCflurry EL %ile | 2 | 0.95 | presentation |
| **BigMHC_IM** (immunogenicity) | 3 | 0.90 | dedicated recognition model — worse |
| **DeepImmuno** (immunogenicity) | **18** | 0.15 | dedicated recognition model — **near-worst** |

DYNC1H1 is both the strongest presenter *and* the recognized one, so presentation puts it first and Epicurus
(co-)ties for best. Crucially, the two models that explicitly try to add a **recognition/immunogenicity**
axis — BigMHC_IM and DeepImmuno — do **worse**, DeepImmuno dramatically so. This corroborates the project's
repeated result: bolting an orthogonal "recognition" score onto presentation does not help and can hurt.
(Caveat: n=1 in-universe positive — descriptive only. Several non-recognized candidates, e.g. PIP5K1A, are
near-tied top presenters, so presentation's #1 here is partly luck, not clean separation.)

## Shortcomings, per method
- **pVACtools / presentation-first curation:** recall — drops low-expression (MAP2) and off-timepoint (ASPM)
  neoantigens before ranking. This is the dominant, deployment-relevant failure.
- **BigMHC_IM, DeepImmuno (immunogenicity predictors):** mis-rank the recognized peptide relative to plain
  presentation; DeepImmuno actively anti-ranks it.
- **genuine PRIME / MHCflurry EL / NetMHCpan-EL:** fine on this case (rank 1–2) but cannot separate the one
  recognized peptide from equally-strong-presenting non-recognized ones.
- **Epicurus v0.1 (ours):** ranks the recognized one #1, but only because it inherits presentation; it adds
  no recognition insight here and, like every ranker, is blind to the 2/3 filtered upstream.

## What this says about making ours better (principled, to be validated — NOT fit to these 3 points)
1. **Move the lever from the ranker to the funnel.** Within-shortlist ranking is saturated (presentation ≈
   ceiling; immunogenicity models hurt). The recoverable loss is **recall**: a deployable tool should *not*
   hard-filter candidates on expression (MAP2 at 5.2 TPM was recognized) — carry low-expression candidates
   with an evidence flag instead of dropping them.
2. **Aggregate variants across timepoints/regions.** ASPM was missed by a single-timepoint callset; this
   patient has 4 timepoints. Union the mutanome across timepoints to avoid heterogeneity misses.
3. **High-recall + honest abstention, not false-precision reranking.** Since we cannot yet measure the
   recognition axis, the honest product surfaces well-presented candidates, is explicit that recognition is
   *unobserved*, and abstains rather than inventing a recognition ranking. This patient validates that design
   philosophy (the earlier "deterministic evidence policy" prototype), and refutes adding a learned
   recognition score on top of presentation.

Each idea must be validated on the independent labeled cohorts (multimer / Gartner / IMPROVE / CheckMate),
never on this n=3.

## Caveats
- n=3 measured positives (1 in-universe). Single-positive AUROC / hits@k are descriptive for one patient.
- Negatives assumed = the other curated mutations (tested pool ≈ curated set). The exact stimulation-pool
  composition (true denominator) and any additional patients/labels should be **requested from Hudson Lab /
  RTTP** — that is the path from n=3 descriptive to a real benchmark.
- The MAP2 "expression-filter drop" and ASPM "off-timepoint" attributions are consistent with the callsets
  and TPM but not confirmed against the vendor's exact curation thresholds / pool-design document.
