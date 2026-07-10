# Manual-source contract: NeoVax melanoma (Hu 2021 + Ott 2017)

Hu 2021 (PMC8273876) and Ott 2017 (PMC5577644) are PMC *author manuscripts*, not in the
open-access bulk subset; their supplement files sit behind a JavaScript proof-of-work
download gate, so automated fetch (curl/urllib) cannot retrieve them. A human browser can.

Download URLs (verified to resolve; open in a browser so the proof-of-work gate clears).
NOTE the path segment is `/articles/instance/<id>/bin/`, NOT `/articles/PMC<id>/bin/`:

  Hu Suppl_DataSet.xlsx  [ESSENTIAL]
    https://pmc.ncbi.nlm.nih.gov/articles/instance/8273876/bin/NIHMS1707651-supplement-Suppl_DataSet.xlsx
  Hu Supplementary_Table_1.docx  [optional]
    https://pmc.ncbi.nlm.nih.gov/articles/instance/8273876/bin/NIHMS1707651-supplement-Supplementary_Table_1.docx
  Ott Supp_5.pdf  [ESSENTIAL]
    https://pmc.ncbi.nlm.nih.gov/articles/instance/5577644/bin/NIHMS892660-supplement-Supp_5.pdf
  Ott Supp_6.pdf  [ESSENTIAL]
    https://pmc.ncbi.nlm.nih.gov/articles/instance/5577644/bin/NIHMS892660-supplement-Supp_6.pdf

Fallback if a direct link hangs: open the article page and grab the same named files from its
"Supplementary Materials" / "Associated Data" section:
  Hu:  https://pmc.ncbi.nlm.nih.gov/articles/PMC8273876/
  Ott: https://pmc.ncbi.nlm.nih.gov/articles/PMC5577644/

Place the downloaded files EXACTLY here (names must match):

  data/raw/hu_melanoma_2021/manual/NIHMS1707651-supplement-Suppl_DataSet.xlsx        [ESSENTIAL]
  data/raw/hu_melanoma_2021/manual/NIHMS1707651-supplement-Supplementary_Table_1.docx [optional]
  data/raw/ott_melanoma_2017/manual/NIHMS892660-supplement-Supp_5.pdf                 [ESSENTIAL]
  data/raw/ott_melanoma_2017/manual/NIHMS892660-supplement-Supp_6.pdf                 [ESSENTIAL]

The ingestion will checksum whatever is placed here, record it in the source manifest +
fetch_record for reproducibility, and refuse to proceed (no fabrication) if a file is missing.

## Fukuoka dendritic-cell vaccine cohort

The 17-patient 2021 source is blocked pending lawful manual access and patient-overlap resolution.
See `configs/source_manifests/fukuoka_dc.yml`. Place verified files only under
`data/raw/fukuoka_dc_2021/manual/`; record checksums before any adapter work. Do not substitute later
case reports or infer peptide negatives from the abstract.
