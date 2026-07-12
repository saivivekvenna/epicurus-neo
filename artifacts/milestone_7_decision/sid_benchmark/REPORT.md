# Sid end-to-end benchmark — label-blind generation

`python -m scripts.sid_benchmark_generate` · mode online · PRIME `7b18d4e110`

_Post-hoc n=1 patient / 3 recognized positives — descriptive only. Generation is over the complete label-blind eligible INPUT universe; output coverage is guard-measured and incomplete. The 3 exact labels are joined only after generation+scoring are frozen._


**Universe:** 200 public → 147 eligible → 139 generator-supported (missense/frameshift), 8 unsupported.

**Generation:** 137 ok, 2 failed, 8 unsupported → 62540 peptide×HLA candidates.

**Coverage guard:** COVERAGE_BELOW_THRESHOLD (coverage 0.9252).


## Verdict

- **Full all-consequence product claim: NOT_EVALUABLE.** Only 130/147 eligible mutations produced candidates (88.4%); 10 consequence classes are unsupported and 7 supported mutations failed transcript generation.
- **Supported-scope, label-blind diagnostic:** all 3 recognized mutations reached scoring. Genuine PRIME recovered 2/3 in its mutation-level top 20; frozen Epicurus v0.1 recovered 1/3. Therefore this run does **not** prove Epicurus reranking beats PRIME.
- **Separate-boundary reference only:** the pre-existing 2025.01 pVAC+PRIME arm recovered 1/3, but it did not consume the same longitudinal 147-variant input and is not a matched head-to-head competitor.

## Mutation-level recognized hits@20 (labels joined post-freeze)

| arm | hits@20 / 3 | recognized ranks (variant → rank) |
|---|--:|---|
| presentation_only_mixmhcpred | 1/3 | {'MAP2-chr2-209694772': 26, 'ASPM-chr1-197102716': 21, 'DYNC1H1-chr14-101980529': 2} |
| genuine_prime | 2/3 | {'MAP2-chr2-209694772': 10, 'ASPM-chr1-197102716': 41, 'DYNC1H1-chr14-101980529': 3} |
| frozen_epicurus_v0_1 | 1/3 | {'MAP2-chr2-209694772': 72, 'ASPM-chr1-197102716': 44, 'DYNC1H1-chr14-101980529': 1} |

## Stage of first loss (per recognized positive)

- `ASPM-chr1-197102716` — **reached scoring (rankable)**
- `DYNC1H1-chr14-101980529` — **reached scoring (rankable)**
- `MAP2-chr2-209694772` — **reached scoring (rankable)**

> The old target-conditioned 'lossless 3/3' is withdrawn (see BENCHMARK_PROTOCOL.md §0). This is the honest label-blind end-to-end result; every missed positive has an explicit stage of first loss. Competitor arms (pVACtools boundary-mismatch, etc.) are reported separately per the protocol.
