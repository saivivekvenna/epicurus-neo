# Actual Epicurus deliverable verdict

## Result

The same shipped Epicurus product logic was run on every locally reconstructable
patient with a mutation-resolved WES/RNA/HLA-derived universe:

| Patient | Generated | Valid | Eligible | Epicurus top 20 | PRIME plain | PRIME cap-2 |
|---|---:|---:|---:|---:|---:|---:|
| Hu_287 | 3/3 | 3/3 | 3/3 | **3/3** | 0/3 | 2/3 |
| Sid | 3/3 | 3/3 | 3/3 | **1/3** | 2/3 | 2/3 |

Both public CLI runs validated their inputs, completed inference, and produced
exactly the same ordered portfolios as the frozen benchmark. The full repository
suite passes: **717 tests**.

## Why Sid failed

Sid was not a generation or gate failure in this canonical run. All three known
positives were generated, deterministic-valid, and had at least one eligible
route. The final product score ranked them at mutation level as follows:

| Mutation | Product mutation rank | Main limiting evidence |
|---|---:|---|
| DYNC1H1 | **1** | none; strong across all components |
| MAP2 | **38** | translation score 0.114; 0 mutant RNA reads at only 4× depth, TPM 5.2 |
| ASPM | **41** | coverage score 0.089; tumor DNA VAF about 4.5% |

MAP2 and ASPM both had excellent presentation and PRIME-derived recognition
evidence, but the production geometric factorization penalized abundance and
coverage enough to push them outside the 20-slot boundary. This is a ranking
failure, not a hard-filter loss.

## Why an earlier Sid result looked better

The earlier Sid 3/3 was **lossless generation + genuine PRIME on a smaller
recovered universe**, not the same complete shipped-product path. An older
target-conditioned recovery result was already withdrawn as invalid. On the
complete 62,540-route universe:

- genuine PRIME recovers 2/3;
- actual Epicurus recovers 1/3.

Therefore the earlier component result must not be presented as Epicurus's
current end-to-end product performance.

## What is genuinely established

1. Epicurus has a runnable candidate-universe-to-portfolio production path with
   exact CLI/benchmark parity.
2. Hu_287 traverses the complete locally reconstructed WES/RNA/HLA pipeline and
   reaches 3/3 in the actual product output.
3. Sid accounts for all 147 eligible variants (137 generated, 10 documented
   non-enumerable), but the actual product reaches only 1/3.
4. Across these two known patients, Epicurus and cap-2 PRIME each recover four
   of six recognized mutations, distributed differently. There is no aggregate
   superiority result.

## Next full-pipeline test

No peptide-table-only benchmark can resolve this. The next decisive step is to
reconstruct additional Miller patients from their public tumor/normal WES and
tumor RNA, freeze the unchanged product path, and measure the same funnel and
top-20 endpoint. Any new ranking policy should be developed without those
patients and frozen before their labels are joined.
