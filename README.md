# NOTE TO REVIEWERS (November 15, 2020):
This repository is not in its final form. Below we have supplied steps to replicate the primary experiments performed in our manuscript. At the moment we supply a .zip file containing preprocessed data for ease of use, since there are a few time-consuming initial steps required to run the workflow from scratch. We plan to supply additional documentation to facilitate replicating the remaining experiments, as well as running ScalpelSig on new data, in the coming weeks. Please pardon our dust while we continue to make usability improvements.

# example ScalpelSig workflow
The following instructions will perform panel discovery optimized for assessing the activity of a given signature examined in our paper (i.e. Signatures 2, 3, 8, 13, 18, 30), and evaluate the panel using held out data. The train and test data are taken from a publicly available cohort of 560 breast cancer genomes.

## step 0: preliminary setup

### 0.1: download the Signature Estimation repository
Run the following command to download the repository:
``` 
git clone https://github.com/lrgr/signature-estimation-py.git
```
The Signature Estimation package uses Anaconda to manage dependencies. Run the following commands to enter the Signature Estimation project directory, build the conda environment, and then return to the ScalpelSig project directory:
```
cd signature-estimation-py
conda env create -f environment.yml
cd ..
```
We will activate the environment later in the workflow.

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
Rscript initialize_10k_panel_script.R -t EXAMPLE_EXPERIMENT -n 3
```
Here, the `-t` argument gives a file tag, which is used to keep track of experimental trials in the downstream pieces of the workflow. The `-n` argument gives the number of iterations, i.e. the number of trials for each signature. For this example, we have set the number of iterations to 3, though in the paper we perform 15 iterations of each experiment.

## step 2: compute window scoring function on training data
This is the most computationally expensive step in the workflow. When we ran the experiments for the paper, we used a cluster to distribute the computation of this step. We use a SLURM script for this purpose, but the script is configured to our local computing setup at the University of Maryland. We recommend that others take similar measures if they would like to fully replicate the results of the paper. The necessary commands (for one out of six signatures, and one out of two alpha settings of ScalpelSig) on our test example are given below.

```
Rscript run_10k_panel_script.R -s 2 -t EXAMPLE_EXPERIMENT_iter1 -o 2 -w 250
Rscript run_10k_panel_script.R -s 2 -t EXAMPLE_EXPERIMENT_iter2 -o 2 -w 250
Rscript run_10k_panel_script.R -s 2 -t EXAMPLE_EXPERIMENT_iter3 -o 2 -w 250
```
Here, `-t` gives the file tag, which is the same file tag as in the previous step but appended with `_iter<i>` where `<i>` ranges from 1 to the `-n` argument given in the previous step. The `-o` argument selects one of the two parametrizations of alpha shown in the paper -- `-o 1` gives alpha=1, `-o 2` gives alpha=0.5 (the latter is the recommended setting). The `-w` argument gives the number of windows in the panel. In the paper, we use 250 windows in our primary experiments.

## step 3: find mutations in panel windows
This step reads the panel windows discovered in the previous step and records the mutation category counts contained inside.
```Rscript scripts/sbs_mtxs_from_10k_panel_windows.R -t EXAMPLE_EXPERIMENT```
Here, and in all future steps, the `-t` tag requires the file tag given in step 1, there is no longer a need to append `_iter<i>` to it.

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
Note that in this step `-t` can be given multiple file tags, delimited by `,` in order to aggregate the results of multiple experiments into a single output. In this example we only deal with a single file tag though. The `-o` argument gives the output tag for the results file. This performs essentially the same purpose as the file tag, but if you are aggregating the results of multiple experiments it is handy, as it allows you to demarkate groups of trials.

## step 6: summarize results across trials
This step is done in an interactive R session. Begin by writing `R` in the command line. When the R session initializes, run the following commands:

```R
>source("summarize_results.R")
>ls <- list_results_files()
>ls
```
This will list the files in the results directory. Select the file with the desired file tag. In the case of this example, we want `panel_results_df_EXAMPLE_OUTPUT_<timestamp>.tsv` where `<timestamp>` is replaced with a timestamp given by the previous script. If this file is at position 1 in the list, we would run the following to generate the results summary table:

```R
>save_summary_df( paste0(GLOBAL_SCRIPT_OUT, ls[[1]]) )
```

This concludes the experiment.

# Outputs of ScalpelSig

The results of evaluation for the experiment above will be stored at `scalpelsig/out/SCRIPT_OUTS/SUMMARY_panel_results_df_EXAMPLE_OUTPUT_<timestamp>.tsv`. The columns of this summary table can be read as follows:

- Signature - the signature for which the ScalpelSig panels are optimized
- Obj1.R.Spearman - mean Spearman's rank correlation between panel exposures and whole-genome exposures across trials (ScalpelSig with alpha=1)
- Obj2.R.Spearman - mean Spearman's rank correlation between panel exposures and whole-genome exposures across trials (ScalpelSig with alpha=0.5)
- MSK.R.Spearman - mean Spearman's rank correlation between panel exposures and whole-genome exposures across trials (MSK-IMPACT panel)
- WES.R.Spearman - mean Spearman's rank correlation between whole exome exposures and whole-genome exposures across trials (whole exome sequencing)
- Obj1.AUPR - mean AUPR for the binary classification task of distinguishing active from inactive samples given panel exposures across trials (ScalpelSig with alpha=1)
- Obj2.AUPR - mean AUPR for the binary classification task of distinguishing active from inactive samples given panel exposures across trials (ScalpelSig with alpha=0.5)
- MSK.AUPR - mean AUPR for the binary classification task of distinguishing active from inactive samples given panel exposures across trials (MSK-IMPACT panel)
- WES.AUPR - mean AUPR for the binary classification task of distinguishing active from inactive samples given whole exome exposures across trials (whole exome sequencing)
- Percent.Active - the percentage of samples that are in the 'active' class for this signature

The genome windows in the ScalpelSig panel for a given run of the experiment can be found at `scalpelsig/out/SCRIPT_OUTS/PANEL_WINDOWS/`. The windows are given as strings delimited by an underscore, denoting the chromosome, the start coordinate, and the end coordinate of each window in the panel.
