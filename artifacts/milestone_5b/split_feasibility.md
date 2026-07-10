# Leakage-safe split feasibility

| Split | Viable | Evaluation positives | Evaluation tested negatives | Note |
| --- | --- | ---: | ---: | --- |
| Patient holdout | Yes | 59 | 119 | 8 positive evaluation patients |
| Full-study holdout | Yes | 61 | 68 | Braun selected by deterministic seed 0 |
| HLA holdout | Yes | 14 | 60 | Restricted to source-resolved HLA groups |
| Peptide-cluster holdout | Yes | 41 | 125 | Patient links remain grouped |
| Cancer-type holdout | Yes | 61 | 68 | Both classes retained |
| Temporal holdout | No | — | — | Fewer than two source-resolved candidate dates |
| Shared-antigen-group holdout | Yes | 11 | 1 | mKRAS G12A 21-mer held out as one group |

Every viable split has at least one positive patient and explicit tested negatives on both sides.
These are feasibility checks only; no model was fitted.
