# NOTE TO REVIEWERS (November 15, 2020):
This repository is not in its final form. Below we have supplied steps to replicate the primary experiments performed in our manuscript. At the moment we supply a .zip file containing preprocessed data for ease of use, since there are a few time-consuming initial steps required to run the workflow from scratch. We plan to supply additional documentation to facilitate replicating the remaining experiments, as well as running ScalpelSig on new data, in the coming weeks. Please pardon our dust while we continue to make usability improvements.

# example ScalpelSig workflow
The following instructions will perform panel discovery optimized for assessing the activity of a given signature examined in our paper (i.e. Signatures 2, 3, 8, 13, 18, 30), and evaluate the panel using held out data. The train and test data are taken from a publicly available cohort of 560 breast cancer genomes.

## step 0: preliminary setup

### 0.1: download the Signature Estimation repository
``` git clone https://github.com/lrgr/signature-estimation-py.git```

### 0.2: modify the config file
First, open the config file
``` vim GLOBAL_CONFIG.R ```

Then, replace the line that says 
``` GLOBAL_PROJECT_DIR = "CHANGE THIS TEXT TO YOUR WORKING DIRECTORY" ```
so that it accurately reflects the path of the repository, i.e.
``` GLOBAL_PROJECT_DIR = "path/to/scalpelsig_repo/"```

### 0.3: download and unzip data
UPDATE THIS TEXT

## step 1: initialize test and train sets
Now we are ready to begin running the experiments. 
```Rscript blah```

