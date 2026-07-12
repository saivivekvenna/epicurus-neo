# CheckMate 153 — external validation of frozen Epicurus v0.1 vs genuine PRIME

**Verdict: TIE (CONSISTENT_WITH_NO_EFFECT). Sixth straight parity result. The recognition wall is confirmed
and sharpened on a genuinely PRIME-untouched cohort — among peptides that all present well, *every* method
(including genuine PRIME) is near-chance at recognition.**

Run: `python -m scripts.checkmate153_dev` · Artifact: `CM153_EXTERNAL.json`

## Cohort

CheckMate 153 (Alban et al., *Nat Med* 2024, `s41591-024-03240-y`) — combinatorial-tetramer neoantigen
screen in NSCLC (nivolumab). Downloaded openly from the Nature static host (supplementary tables MOESM3/4;
the PMC copy is behind a proof-of-work gate, EuropePMC has no supplement for the NIHMS manuscript).

- **14 patients**, **1,197 candidate 9-mers** (predicted binders chosen for tetramer screening).
- Labels: **162 tetramer-POSITIVE** vs **1,035 tetramer-TESTED_NEGATIVE**.
- Class-I, **HLA-resolved** — every candidate carries its single restricting allele (8 common alleles).
- Screen classes: NonBinder 786 · TETRAMER− 249 · TETRAMER+ 162.

**Why it qualifies as a clean, independent test:**
- **PRIME-untouched.** Tetramer labels published Oct-2024, long after PRIME 2.0's 2023 training set
  (Lausanne/NeoDisc). The frozen Epicurus v0.1 residual was trained only on the CD8-multimer cohort;
  CheckMate 153 was never used to fit it. Only **14/1,197 (1.2%)** peptides near/exact-match PRIME's
  training set (report-only; not dropped from a test cohort).
- **Own model discarded.** The paper's trained `scores` column is ignored. Every feature is recomputed:
  genuine GfellerLab **PRIME 2.1 %rank** (incumbent, 1,197/1,197 scored), **MHCflurry presentation
  percentile** (the `el` feature — NetMHCpan-EL is not installed; MHCflurry is percentile-ranked
  within-patient, so the substitution is scale-free), gene-level **RNA expression** from the study's own
  RNA-seq count table (MOESM3), and DNA **VAF** (carried for future models).

## Results — within-patient mean hits@20 · recall@20 · tested-AUROC

### A. All negatives (predicted-non-binders + tetramer-negatives) — the full candidate universe

| Arm | hits@20 | recall@20 | tested AUROC |
|---|--:|--:|--:|
| **MixMHCpred** (presentation) | **5.357** | 75/162 (0.463) | **0.735** |
| Epicurus v0.1 (frozen) | 5.071 | 71/162 | 0.709 |
| genuine PRIME 2.1 | 5.000 | 70/162 | 0.714 |
| MHCflurry-EL (presentation) | 4.786 | 67/162 | 0.720 |

**GATE — frozen Epicurus v0.1 vs genuine PRIME: Δhits@20 = +0.071, CI[−0.357, +0.571] → TIE.**
Pure **MixMHCpred presentation is the single best ranker** — better than genuine PRIME and better than
Epicurus. PRIME's recognition layer buys nothing over its own presentation backbone here.

### B. Binders-only (TETRAMER+ vs TETRAMER−) — the hard recognition test among predicted binders

Dropping the 786 predicted-non-binder easy negatives (411 candidates: 162 pos / 249 tested-neg):

| Arm | hits@20 | recall@20 | tested AUROC |
|---|--:|--:|--:|
| MixMHCpred | 7.429 | 104/162 | 0.594 |
| **Epicurus v0.1 (frozen)** | 7.286 | 102/162 | 0.586 |
| genuine PRIME 2.1 | 7.214 | 101/162 | **0.597** |
| MHCflurry-EL | 6.714 | 94/162 | 0.532 |

**GATE — frozen Epicurus v0.1 vs genuine PRIME: Δhits@20 = +0.071, CI[−0.214, +0.429] → TIE.**

**The load-bearing finding:** once you restrict to peptides that all present well, **AUROC collapses to
~0.53–0.60 for every method, including genuine PRIME (0.597).** Recognition among presented peptides is
near-chance for the entire field. This is the recognition wall, quantified on a fresh, independent,
PRIME-untouched cohort — presentation carries the whole signal (AUROC 0.73 with easy negatives), and no
method has a real recognition edge over presentation.

## Interpretation

1. **Epicurus ≈ genuine PRIME (tie), sixth straight.** v0.2→v0.5 (development) and now CheckMate 153
   (external) all land at PRIME parity. Δ is nominally +0.07/patient in Epicurus's favour but CI spans 0
   in both analyses; 14 patients is underpowered (wide CIs). No ACCEPT — v0.1 stays the frozen model.
2. **PRIME's recognition component adds ≈0 over MixMHCpred presentation** on this cohort (0.597 vs 0.594
   binders-only; PRIME even trails MixMHCpred with easy negatives). The recognition layer is not paying off.
3. **The wall is at recognition-among-presented, and it is near-chance for everyone.** This is the
   sharpest external confirmation yet of the [[epicurus-thesis]] presentation-solved / recognition-hard split.

## Honest limitations

- **`el` = MHCflurry, not NetMHCpan-EL** (not installed). The arm labelled `netmhcpan_el` in the JSON is
  MHCflurry. Both are EL-type presentation %ranks, percentile-ranked within patient.
- **Expression covers 228/1,197** candidates (only ~10 patients have usable RNA-seq gene mapping in MOESM3);
  where absent, `expr` is neutral, so frozen v0.1 is effectively prime+el for most rows. Weakens the model
  slightly; still a tie. TODO: widen gene→count mapping / recover full expression.
- **Range-restricted denominator** — candidates are predicted binders selected for tetramer, not the raw
  somatic mutanome (like IMPROVE's screened set). A legitimate within-patient test, but not Gartner-dense.
- **Underpowered** — 14 patients. Directional, not decisive.
