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

**NOTE 9/30/2020: IMPORTANTLY, there is a necessary directory structure for the output folder and accompanying CONFIG file (not included in repo since it's different for my laptop and my UMD work station) without which, this workflow will not work properly. TODO: need to figure out best practices for generating the output directory and config file generically.**

## Evaluation Protocol Step 1: Initialize Test & Train Sets 
In the *outer directory* run:
``` 
Rscript run_10k_panel_script.R -t <FILETAG> -g <SIGNATURE GROUP (see below)> -n <NUMBER OF TRIALS IN EXPERIMENT>
```
Some explanation on the input fields:
- the FILETAG is a string that identifies files associated with this run. Later down the line, file names will be parsed with some regexes that are not completely foolproof, so it is possible to mess up the workflow with a badly chosen FILETAG. Try to avoid putting `sig` plus a number, `obj` plus a number, `it` plus a number, or anything that looks like a date. A string of all caps letters e.g. `FIRSTBATCH` should be fine.
- the SIGNATURE GROUP is sort of an artifact of how I batched experiments earlier, and the fact that we use random stratified sampling to obtain a test/train split (i.e. for each signature we guarantee that the number of active samples in the test set is proportional to the number of active samples in the cohort). It's no longer important to explain why the groups are what they are, but just know that **-g needs to be either 3 or 4**. If you input 3, it will generate test/train sets for signatures 2, 3, and 13; if you give it 4 it will generate test/train sets for signatures 1, 5, 8, 18, and 30.
- the NUMBER OF TRIALS is the number of test/train sets it will generate per signature in the group chosen by the -g parameter. **If you set this to 15, the SLURM script will not terminate with the resources alotted to it. Instead, break the experiment into 3 batches (each with its own initialization) and set -n to 5 for each batch.**

Notes:
So here is an example of what I would run to initialize an experiment with 15 trials of all 8 signatures -- I apologize for its jankiness:
```
Rscript run_10k_panel_script.R -t FIRST_BATCH_A -g 3 -n 5
Rscript run_10k_panel_script.R -t SEC_BATCH_A -g 3 -n 5
Rscript run_10k_panel_script.R -t THIRD_BATCH_A -g 3 -n 5

Rscript run_10k_panel_script.R -t FIRST_BATCH_B -g 4 -n 5
Rscript run_10k_panel_script.R -t SEC_BATCH_B -g 4 -n 5
Rscript run_10k_panel_script.R -t THIRD_BATCH_B -g 4 -n 5
```

## Evaluation Protocol Step 2: Computation of Window Scoring Function w/ SLURM
Currently the workflow is to:
- log in to `cbcbsub00`
- do `sbatch SLURM_10k_panel_script.sh <FILETAG> <SIGNATURE>`

Notes: 
The FILETAG should be the same one given to the initialization function. The SIGNATURE should be a *single* signature. So here are the inputs to continue the example from above, whose jankiness I also apologize for:
```
sbatch SLURM_10k_panel_script.sh FIRST_BATCH_A 2
sbatch SLURM_10k_panel_script.sh FIRST_BATCH_A 3
sbatch SLURM_10k_panel_script.sh FIRST_BATCH_A 13


sbatch SLURM_10k_panel_script.sh SEC_BATCH_A 2
sbatch SLURM_10k_panel_script.sh SEC_BATCH_A 3
sbatch SLURM_10k_panel_script.sh SEC_BATCH_A 13


sbatch SLURM_10k_panel_script.sh THIRD_BATCH_A 2
sbatch SLURM_10k_panel_script.sh THIRD_BATCH_A 3
sbatch SLURM_10k_panel_script.sh THIRD_BATCH_A 13


sbatch SLURM_10k_panel_script.sh FIRST_BATCH_B 1
sbatch SLURM_10k_panel_script.sh FIRST_BATCH_B 5
sbatch SLURM_10k_panel_script.sh FIRST_BATCH_B 8
sbatch SLURM_10k_panel_script.sh FIRST_BATCH_B 18
sbatch SLURM_10k_panel_script.sh FIRST_BATCH_B 30


sbatch SLURM_10k_panel_script.sh SEC_BATCH_B 1
sbatch SLURM_10k_panel_script.sh SEC_BATCH_B 5
sbatch SLURM_10k_panel_script.sh SEC_BATCH_B 8
sbatch SLURM_10k_panel_script.sh SEC_BATCH_B 18
sbatch SLURM_10k_panel_script.sh SEC_BATCH_B 30


sbatch SLURM_10k_panel_script.sh THIRD_BATCH_B 1
sbatch SLURM_10k_panel_script.sh THIRD_BATCH_B 5
sbatch SLURM_10k_panel_script.sh THIRD_BATCH_B 8
sbatch SLURM_10k_panel_script.sh THIRD_BATCH_B 18
sbatch SLURM_10k_panel_script.sh THIRD_BATCH_B 30
```

## Evaluation Protocol Step 3: Obtaining 96-Category Mutation Count Matrices from Panel Regions
blah
