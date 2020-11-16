# NOTE TO REVIEWERS (November 15, 2020):
This repository is not in its final form. Below we have supplied steps to replicate the primary experiments performed in our manuscript. At the moment we supply a .zip file containing preprocessed data for ease of use, since there are a few time-consuming initial steps required to run the workflow from scratch. We plan to supply additional documentation to facilitate replicating the remaining experiments, as well as running ScalpelSig on new data, in the coming weeks. Please pardon our dust while we continue to make usability improvements.

# example ScalpelSig workflow
The following instructions will perform panel discovery optimized for assessing the activity of a given signature examined in our paper (i.e. Signatures 2, 3, 8, 13, 18, 30), and evaluate the panel using held out data. The train and test data are taken from a publicly available cohort of 560 breast cancer genomes.

## step 0: preliminary setup

### 0.1: download the Signature Estimation repository
``` 
git clone https://github.com/lrgr/signature-estimation-py.git
```

### 0.2: modify the config file
First, open the config file
``` 
vim GLOBAL_CONFIG.R 
```

Then, replace the line that says 
``` 
GLOBAL_PROJECT_DIR = "CHANGE THIS TEXT TO YOUR WORKING DIRECTORY" 
```
so that it accurately reflects the path of the repository, i.e.
```
GLOBAL_PROJECT_DIR = "path/to/scalpelsig_repo/"
```

### 0.3: download and unzip data
UPDATE THIS TEXT

## step 1: initialize test and train sets
Now we are ready to begin running the experiments. 
```
Rscript initialize_10k_panel_script.R -t EXAMPLE_EXPERIMENT -g 3 -n 3
```

## step 2: compute window scoring function on training data
```
Rscript run_10k_panel_script.R -s 2 -t EXAMPLE_EXPERIMENT_iter1 -o 2 -w 250
Rscript run_10k_panel_script.R -s 2 -t EXAMPLE_EXPERIMENT_iter2 -o 2 -w 250
Rscript run_10k_panel_script.R -s 2 -t EXAMPLE_EXPERIMENT_iter3 -o 2 -w 250
```

## step 3: find mutations in panel windows
```Rscript scripts/sbs_mtxs_from_10k_panel_windows.R -t EXAMPLE_EXPERIMENT```

## step 4: estimate panel signatures
To estimate panel signatures, we use the Signature Estimation package. This requires that we activate the Signature Estimation conda library:

```conda activate signature-estimation-py-env```

Then we call the script to continue the experiment

```Rscript scripts/estimate_10k_panel_signatures.R -t EXAMPLE_EXPERIMENT```

Afterwards, we must deactivate the conda library:

```conda deactivate```

## step 5: evaluate performance

```
Rscript scripts/evaluate_10k_panel_results.R -t EXAMPLE_EXPERIMENT -o EXAMPLE_OUTPUT 
```

## step 6: summarize results across trials
This step is done in an interactive R session. Begin by writing `R` in the command line. When the R session initializes, run the following commands:

```
>source("summarize_results.R")
>ls <- list_results_files()
>ls
```
This will list the files in the results directory. Select the file with the desired file tag. In the case of this example, we want `panel_results_df_EXAMPLE_OUTPUT_<timestamp>.tsv` where `<timestamp>` is replaced with a timestamp given by the previous script. If this file is at position 1 in the list, we would run the following to generate the results summary table:

```
>save_summary_df( paste0(GLOBAL_SCRIPT_OUT, ls[[1]]) )
```

