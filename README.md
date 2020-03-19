(Informal readme of adhoc steps, will make more formal later)

Step 0: Download data and preprocess
Use https://github.com/lrgr/sigma to get Nik-Zainal data
Download signature estimation repo

Process `ICGC-BRCA-EU.RELEASE_22.SBS.renamed.tsv` into mutation count matrix using `save_sbs_tsv()` in `preprocess_windows.R`
Estimate signatures using signature estimation tool

Use `sig_est_with_window_tsv()` in `preprocess_windows.R` to obtain a data frame with signature estimates and windowed mutation counts for each sample.
This file is saved in `data/nz_sig_est_and_window_counts.tsv`


Step 1: Feature selection using mutual information



