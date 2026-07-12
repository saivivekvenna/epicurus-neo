# EXPLORATORY PROTOCOL — lossless variant-to-peptide recovery

Paired frozen copy of
`docs/superpowers/specs/2026-07-12-osteosarc-lossless-peptide-recovery-exploratory-protocol.md`.

**Post-hoc disclosure:** Sid informed the design and the feasibility scores were inspected before this
protocol was committed (77 ASPM and 259 MAP2 sequences; best genuine-PRIME ranks approximately 0.088
and 0.010). This is an engineering replay, not blind or independent validation.

Frozen path:

`raw GRCh38 allele → Ensembl VEP MANE/canonical consequence → Ensembl CDS/protein → all standard-AA
8–14 mutation-covering windows → patient HLA panel from pVAC → genuine PRIME → existing evidence
router → frozen route-aware top-20`.

Generation reads no Hudson labels, assay ledger, vaccine peptide table, mRNA construct, or ELISPOT
result. Labels enter only after ranking. The full expectations, failure conditions, identity rules,
outputs, MAP2 shared-ORF correction, and interpretation guardrails are frozen in the source document.

