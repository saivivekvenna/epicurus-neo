# Data Workflow

This workflow keeps the benchmark honest while allowing large public datasets to
be downloaded or mounted locally as needed.

## 1. List Known Data Sources

```bash
epicurus list-datasets
epicurus download-plan --dataset neoranking
epicurus download-plan --dataset gartner_nci
```

Some Figshare pages may require manual browser download because they return a
WAF challenge to command-line clients. Place manually downloaded files under
`data/raw/<dataset>/`.

## 2. Normalize Source Tables

NeoRanking neo-peptide matrix:

```bash
epicurus normalize \
  --kind neoranking-neopep \
  --input data/raw/neoranking/Neopep_data_org.txt \
  --output data/processed/neoranking_neopep.normalized.csv
```

Gartner/NCI-style tables:

```bash
epicurus normalize \
  --kind gartner \
  --input data/raw/gartner_nci/NmersTestingSet.txt \
  --output data/processed/gartner_nmers_test.normalized.csv
```

Then validate:

```bash
epicurus validate-schema data/processed/neoranking_neopep.normalized.csv
```

## 3. Run Pre-Locked-Test Cross Validation

Before touching TESLA or other locked tests, run grouped CV:

```bash
epicurus group-cv \
  data/processed/neoranking_neopep.normalized.csv \
  --group-col study_id \
  --metric-group-col patient_id \
  -k 20
```

If folds are reported as `leakage_blocked`, do not override them for headline
claims. Change the grouping or deduplicate the data instead.

## 4. Run a Frozen Train/Test Evaluation

Only after the feature recipe is frozen:

```bash
epicurus train-eval \
  --train data/processed/train.normalized.csv \
  --test data/processed/tesla.normalized.csv \
  --group-col patient_id \
  -k 20 \
  --write-scored outputs/tesla.scored.csv
```

## 5. Select the Submitted Set

```bash
epicurus select-portfolio \
  outputs/patient_case.scored.csv \
  --score-col epicurus_score \
  -k 20 \
  --max-per-hla 8 \
  --max-per-gene 3 \
  --output outputs/patient_case.top20.csv
```

The diversity limits are optional and must be tuned only on validation data, not
on locked tests.

