# Sid end-to-end benchmark — label-blind generation

`python -m scripts.sid_benchmark_generate --offline` · mode offline · PRIME `7b18d4e110`

_Post-hoc n=1 patient / 3 recognized positives — descriptive only. Generation is over the COMPLETE label-blind eligible universe (guard-enforced); the 3 exact labels are joined only after generation+scoring are frozen._


**Universe:** 200 public → 147 eligible → 137 generator-supported (missense/frameshift), 10 unsupported.

**Generation:** 130 ok, 7 failed, 10 unsupported → 59755 peptide×HLA candidates.

**Coverage guard:** COVERAGE_BELOW_THRESHOLD (coverage 0.8844).


## Mutation-level recognized hits@20 (labels joined post-freeze)

| arm | hits@20 / 3 | recognized ranks (variant → rank) |
|---|--:|---|
| presentation_only_mixmhcpred | 2/3 | {'MAP2-chr2-209694772': 25, 'DYNC1H1-chr14-101980529': 2, 'ASPM-chr1-197102716': 20} |
| genuine_prime | 2/3 | {'MAP2-chr2-209694772': 10, 'DYNC1H1-chr14-101980529': 3, 'ASPM-chr1-197102716': 39} |
| frozen_epicurus_v0_1 | 1/3 | {'MAP2-chr2-209694772': 66, 'DYNC1H1-chr14-101980529': 1, 'ASPM-chr1-197102716': 39} |

## Stage of first loss (per recognized positive)

- `ASPM-chr1-197102716` — **reached scoring (rankable)**
- `DYNC1H1-chr14-101980529` — **reached scoring (rankable)**
- `MAP2-chr2-209694772` — **reached scoring (rankable)**

> The old target-conditioned 'lossless 3/3' is withdrawn (see BENCHMARK_PROTOCOL.md §0). This is the honest label-blind end-to-end result; every missed positive has an explicit stage of first loss. Competitor arms (pVACtools boundary-mismatch, etc.) are reported separately per the protocol.
