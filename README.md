(Informal README, will make more formal later)

## Step 0: Initialization

### 0.0: Initialize Conda Environment
(TODO: figure out what packages are necessary in the conda environment. If possible, create a conda environment that can run both this package and SigEstimator)


### 0.1: Download Data and Preprocess
(TODO: revisit this -- is the info below still current? What additional data needs to be downloaded (**e.g. COSMIC sigs**)? Are there additional preprocessing steps that need to happen before the main body of code is run?)

Use https://github.com/lrgr/sigma to get Nik-Zainal data
Download signature estimation repo

Process `ICGC-BRCA-EU.RELEASE_22.SBS.renamed.tsv` into mutation count matrix using `save_sbs_tsv()` in `preprocess_windows.R`
Estimate signatures using signature estimation tool
- `conda activate signature-estimation-py-env`
- `python signature-estimation-py/signature_estimation.py -mf MUTATION_COUNTS_FILE -sf SIGNATURES_FILE -of OUTPUT_FILE`

Use `sig_est_with_window_tsv()` in `preprocess_windows.R` to obtain a data frame with signature estimates and windowed mutation counts for each sample.
This file is saved in `data/nz_sig_est_and_window_counts.tsv`


...

**NOTE 9/30/2020: I am sketching the workflow here but this may not be exactly correct -- I do not have access to my lab notebook since I am currently on Northwestern campus due to internet outage, but later this week I plan to update w/ the exact workflow.**

## Evaluation Protocol Step 1: Initialize Test & Train Sets 
In the *outer directory* run:
``` 
Rscript run_10k_panel_script.R -t <FILETAG> -g <SIGNATURE GROUP (see below)> -n <NUMBER OF TRIALS IN EXPERIMENT>
```
Some explanation on the input fields:
- blah

## Evaluation Protocol Step 2: Computation of Window Scoring Function w/ SLURM
Currently the workflow is to:
- log in to `cbcbsub00`
- do `sbatch SLURM_10k_panel_script.sh 
