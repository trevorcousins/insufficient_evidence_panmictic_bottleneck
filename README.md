# insufficient_evidence_panmictic_bottleneck

A recent paper in Science proposed that humans went through a very strong population bottleneck around 1 million years ago. We do not believe the data is in support of this. This repo describes some analysis we did, demonstrating that simpler panmictic models without the bottleneck better fit the data. 

## Real data processing

Download 1kGP VCFs from https://www.internationalgenome.org/category/vcf/<br>
Write allele counts at each positions with write_allele_counts_241011.sh<br>
(we did not polarise the alleles using a chimp or gorilla sequence, which is clearly the wrong thing to do, however mushi can model the rate of mididentification very well)<br>
Write the SFS with write_SFS_submission_nopolarisation_241011.sh which calls write_SFS_nopolar_241011.py (this includes a strict mappabiltiy mask, and a B-map from Murphy et al., 2023)<br>

## Running mushi

We inferred a model with mushi using searching_mushi_model_241011.sh, which calls infer_mushi_model_241011.py<br>

## Running FitCoal

We inferred a model with FitCoal using infer_FitCoal_241011.sh<br>

## Analysis 

Various other bits of analysis are done in plots_for_paper_241025_upload.ipynb<br>

## Preprint

Our preprint is available to read at https://www.biorxiv.org/content/10.1101/2024.10.21.619456v1
