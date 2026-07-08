# Download Checklist

Place downloaded files under `data/raw/<dataset>/`. Raw data is ignored by git.

## Required First

### NeoRanking

Repository and instructions:

- https://github.com/bassanilab/NeoRanking

Download:

- `Neopep_data_org.txt` or `Neopep_data_org.txt.zip`  
  https://figshare.com/s/a000b0990465ab3e9d33
- `Mutation_data_org.txt` or `Mutation_data_org.txt.zip`  
  https://figshare.com/s/3c27fa3b705a74bdfa10
- `HLA_allotypes.txt`  
  If present in the same Figshare bundle, put it at `data/raw/neoranking/hla/HLA_allotypes.txt`.

Suggested paths:

```text
data/raw/neoranking/Neopep_data_org.txt.zip
data/raw/neoranking/Mutation_data_org.txt.zip
data/raw/neoranking/hla/HLA_allotypes.txt
```

### Gartner/NCI

Collection:

- https://nih.figshare.com/collections/Datasets_for_Development_of_a_model_for_ranking_candidate_HLA_class_I_neoantigens_based_upon_datasets_of_known_neoepitopes_/4792338

Download in this priority order:

- Training Nmers randomly subsampled to reduce expression bias, 2.02 MB  
  https://nih.figshare.com/articles/dataset/Training_Set_of_Long_Peptides_Screened_by_NCI_Surgery_Branch_for_Reactivity_Against_TIL/11400972
- Test Set of Long Peptides, 4.11 MB  
  https://nih.figshare.com/articles/dataset/Test_Set_of_Long_Peptides_Screened_by_NCI_Surgery_Branch_for_Reactivity_Against_TIL/11400984
- All Long Peptides, 12.55 MB  
  https://nih.figshare.com/articles/dataset/All_Long_Peptides_Screened_by_the_NCI_Surgery_Branch_for_Reactivity_Against_TIL/11400966

Optional/heavy:

- Training Mmps, 1.12 GB  
  https://nih.figshare.com/articles/dataset/All_Mutated_Minimal_Peptides_Screened_by_the_NCI_Surgery_Branch_for_Reactivity_Against_TIL/11400975
- Test Set of Mutated Minimal Peptides, 2.2 GB  
  https://nih.figshare.com/articles/dataset/Test_Set_of_Mutated_Minimal_Peptides_screened_by_NCI_Surgery_Branch_for_Reactivity_Against_TIL/11400987
- All Mutated Minimal Peptides, 7.19 GB  
  https://nih.figshare.com/articles/dataset/Training_Set_of_Mutated_Minimal_Peptides_Screened_by_NCI_Surgery_Branch_for_Reactivity_Against_TIL/11400969

### TESLA

Keep locked until the model recipe is frozen.

- Article / PMC page: https://pmc.ncbi.nlm.nih.gov/articles/PMC7652061/
- Supplemental `mmc5.xlsx` from the Cell article, or Mendeley mirror:
  https://data.mendeley.com/datasets/6x87nx8jtc

Suggested path:

```text
data/raw/tesla/mmc5.xlsx
```

## Useful Next

### CEDAR

- Database: https://cedar.iedb.org/
- Export page: https://cedar.iedb.org/database_export_v3.php

### NEPdb

- https://nep.whu.edu.cn/

### BigMHC

- Data: https://data.mendeley.com/datasets/dvmz6pkzvb/4
- Code: https://github.com/KarchinLab/bigmhc

Direct downloads:

- `datasets.zip`, 144.27 MB  
  https://data.mendeley.com/public-files/datasets/dvmz6pkzvb/files/1da0314f-692d-4c39-b81f-fa2a7dba86bf/file_downloaded
- `manafest.csv`, 34.59 KB  
  https://data.mendeley.com/public-files/datasets/dvmz6pkzvb/files/6de56487-7b7e-459e-8aee-c80475face9d/file_downloaded

Suggested paths:

```text
data/raw/bigmhc/datasets.zip
data/raw/bigmhc/manafest.csv
```

Normalize the immunogenicity splits:

```bash
epicurus normalize --kind bigmhc --input data/raw/bigmhc/datasets.zip --zip-member im_train.csv --output data/processed/bigmhc_im_train.normalized.csv
epicurus normalize --kind bigmhc --input data/raw/bigmhc/datasets.zip --zip-member im_val.csv --output data/processed/bigmhc_im_val.normalized.csv
epicurus normalize --kind bigmhc --input data/raw/bigmhc/datasets.zip --zip-member im_test.csv --output data/processed/bigmhc_im_test.normalized.csv
epicurus normalize --kind bigmhc --input data/raw/bigmhc/manafest.csv --output data/processed/bigmhc_manafest.normalized.csv
```

Add presentation and retrieval features for the BigMHC hard-part benchmark:

```bash
epicurus add-mhcflurry-features --input data/processed/bigmhc_im_train_val.normalized.csv --output data/processed/bigmhc_im_train_val.mhcflurry.csv
epicurus add-mhcflurry-features --input data/processed/bigmhc_im_test.normalized.csv --output data/processed/bigmhc_im_test.mhcflurry.csv
epicurus add-retrieval-features --input data/processed/bigmhc_im_test.mhcflurry.csv --reference data/processed/bigmhc_im_train_val.mhcflurry.csv --output data/processed/bigmhc_im_test.mhcflurry_retrieval.csv
epicurus apply-score-selector --validation data/processed/bigmhc_im_val.mhcflurry_retrieval.csv --target data/processed/bigmhc_im_test.mhcflurry_retrieval.csv --output outputs/benchmarks/bigmhc_im_test_selected_score.scored.csv --selection-output outputs/benchmarks/bigmhc_im_test_selected_score.selection.json --group-col hla_allele --score-col retrieval_max_positive_similarity --score-col retrieval_positive_minus_negative_similarity --score-col retrieval_topk_positive_similarity_mean --score-col retrieval_topk_positive_fraction --score-col retrieval_biochemical_topk_positive_similarity_mean --score-col retrieval_biochemical_max_positive_similarity --score-col mhcflurry_presentation_score --score-col mhcflurry_processing_score --score-col epicurus_transfer_score
```
