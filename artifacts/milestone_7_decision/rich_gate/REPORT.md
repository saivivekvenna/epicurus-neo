# Rich-feature dynamic gate — IMPROVE nested Partition eval

`python -m scripts.rich_gate_experiment`

_DEVELOPMENT evidence only (single study IMPROVE); NOT external proof._

Leakage controls: recurrent-peptide quarantine, near-peptide guard >=0.8, within-fold preprocessing, patient-disjoint partitions, excluded label/identity/TME/pipeline-score columns.


## Primary + pairwise (all safe families)

| model | ungated | gated | **Δ hits@20** [CI] | p_better | imp/tie/harm | random-matched Δ | beats random |
|---|--:|--:|--:|--:|--:|--:|:--:|
| primary_histgbt | 1.0714 | 0.8714 | **-0.200** [-0.4143, 0.0143] | 0.029 | 10/42/18 | -0.085 | NO |
| pairwise | 1.0714 | 0.8143 | **-0.257** [-0.4857, -0.0286] | 0.01 | 11/38/21 | -0.117 | NO |

## Single-family ablations (HistGBT)

| family | Δ hits@20 [CI] | imp/tie/harm | random-matched Δ | beats random |
|---|--:|--:|--:|:--:|
| only_wt_rank | -0.043 [-0.2, 0.1143] | 7/53/10 | -0.035 | NO |
| only_dai | -0.129 [-0.2857, 0.0286] | 9/46/15 | -0.175 | yes |
| only_rna | -0.100 [-0.3, 0.0714] | 11/47/12 | -0.088 | NO |
| only_dna_vaf | +0.014 [-0.1571, 0.1714] | 13/47/10 | -0.095 | yes |
| only_stability | -0.143 [-0.3, 0.0] | 7/49/14 | -0.111 | NO |
| only_foreign | -0.171 [-0.3714, 0.0286] | 9/44/17 | -0.079 | NO |
| only_physchem | -0.029 [-0.1857, 0.1286] | 10/49/11 | -0.097 | yes |
| only_mutclass | -0.029 [-0.1714, 0.1286] | 8/51/11 | -0.016 | NO |

## Verdict
NULL/negative on IMPROVE nested folds: Δ -0.200 (CI [-0.4143, 0.0143]); random-matched -0.085. Rich features do not convert to net hits@20 here.
