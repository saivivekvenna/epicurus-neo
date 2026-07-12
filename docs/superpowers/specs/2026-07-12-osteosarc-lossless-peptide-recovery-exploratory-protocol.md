# Osteosarc/Sid lossless peptide recovery — exploratory protocol

**Status: post-hoc engineering diagnostic, not preregistered validation.** Sid's recognition labels
motivated the router, and Codex inspected an exploratory genuine-PRIME run before this document was
committed. The observed feasibility numbers (77 ASPM peptides, 259 MAP2 peptides; best PRIME ranks
approximately 0.088 and 0.010) are disclosed up front. This protocol freezes a reproducible generator
for future patients; it does not make Sid independent or license a superiority claim.

Paired copy: `artifacts/milestone_7_decision/peptide_recovery/EXPLORATORY_PROTOCOL.md`.

## Question

Can variants already recovered by the multi-caller union but absent from pVACtools be converted into
rankable peptide–HLA candidates using only variant/reference annotations, without reading vaccine or
recognition-assay peptides?

## Frozen inputs and provenance

1. GRCh38 chromosome, position, reference, and alternate allele from the public raw VAF table.
2. Ensembl REST VEP consequence with HGVS, MANE, canonical, and protein fields.
3. Ensembl REST canonical/MANE CDS or protein sequence, cached with URL and SHA-256.
4. Patient HLA-A/B/C panel read from the pre-existing pVAC candidate file; exact expected set is
   asserted, not silently replaced.
5. Genuine GfellerLab PRIME 2.1 + MixMHCpred 3.0 at the already-pinned commits.

The generator and ranking runner must not import or read Hudson labels, the assay ledger, vaccine
peptide blocks, site mRNA minimal peptides, or site ELISPOT outcomes. Recognition labels are joined
only after the selected set is frozen for evaluation.

## Frozen transcript expectations (fail closed)

- ASPM chr1:197102716 C>T: canonical/MANE `ENST00000367409` / `NM_018136.5`,
  `c.6535G>A`, `p.Gly2179Arg`, protein position 2179.
- MAP2 chr2:209694772 GGCTACTGTGTGTTCAATAAGTACACAGT>G: canonical/MANE
  `ENST00000682079` / `NM_001375505.1`, `c.2603_2630del`,
  `p.Gly868AlafsTer38`, protein position 868.
- DYNC1H1 chr14:101980529 G>A is an input-only positive control: its MANE transcript must match the
  existing pVAC transcript `ENST00000360184`; generated windows must include `KRFHATISF`.

Any transcript, HGVS, position, reference-residue, or sequence mismatch aborts instead of falling
back to a convenient isoform.

## Frozen generation algorithm

- Standard amino acids only; class-I lengths 8–14 inclusive.
- Missense: mutate the VEP protein residue and enumerate every 8–14mer containing that residue.
- Frameshift: apply the exact VEP one-based inclusive CDS deletion, translate from the CDS start to
  the first stop, and enumerate every 8–14mer containing at least one residue at or after
  `protein_start`. Stop symbols never enter a peptide.
- Cross every unique peptide with the HLA panel from the pVAC file.
- Extend candidate identity by patient + genomic allele + peptide + HLA; never collapse distinct
  peptide–HLA routes.
- Score genuine PRIME; orient the incumbent as `-prime_rank` (higher is better).
- Union recovered candidates with the original pVAC set, route with the frozen evidence router, and
  select the frozen route-aware k=20 portfolio with its existing mutation/gene/HLA caps.

## Primary-source correction frozen before implementation

MAP2 p.Leu867fs and p.Gly868fs remain distinct genomic alleles. They do **not**, however, imply
disjoint mutant antigen sequences. Ensembl's Gly868 consequence translates to
`...DSQLEDLAHCHHLFKTVRIYQGRVVPFTKALMIKFEEIWPQTFH*`, sharing the downstream mutant ORF and
`RVVPFTKAL` with the Leu867-associated construct. VEP identifies `GYCVFNKYTV` as the **reference**
amino-acid context. It is therefore reasonable—but still an inference—to read the Hudson label
`MAP2.p.GYCVFNKYTV868fs` as reference-context naming rather than a mutant epitope sequence.

## Frozen outputs

- unique generated peptide counts and peptide–HLA counts per rescued variant;
- exact Ensembl URLs/hashes, selected transcript/HGVS fields, and short junction context;
- genuine-PRIME scored counts and best ranks;
- static pVAC-only versus pVAC+recovery pure-PRIME and route-aware top-20 candidate sets;
- recognized-mutation candidate-generation recall, rankable recall, and top-20 coverage;
- online and offline result/hash equivalence.

## Interpretation guardrails

This result may establish mechanical recovery on this known case. It cannot show that Epicurus beats
PRIME: the likely mechanism is that better candidate generation lets genuine PRIME score targets it
previously never received. Any performance claim requires a future, untouched patient cohort with a
raw multi-caller callset, WES/RNA/HLA, lossless generation, and measured recognition outcomes.

