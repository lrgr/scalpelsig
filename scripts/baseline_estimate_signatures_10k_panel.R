# make sure you activate the signature estimation library's conda environment before running this, or it won't work
# the command is:
# conda activate signature-estimation-py-env

setwd("..")
source("projection_score.R")

DEBUG_FLAG = TRUE

file_tag = "bp_test"

files = list.files(GLOBAL_SCRIPT_BASELINE_SBS, pattern= paste0(".*", file_tag, ".*") )
print(paste0(Sys.time(), "    found ", length(files), " files containing the tag: ", file_tag))

signatures_file = paste0(GLOBAL_DATA_DIR, "cosmic-signatures.tsv") #location of .tsv containing mutational signature frequencies (e.g. COSMIC signatures)


i = 1
for (f in files) {
	print( paste0(Sys.time(), "    ", i, "/",  length(files)) )

	sbs_infile = paste0(GLOBAL_SCRIPT_BASELINE_SBS, f)

	

	s = sub("baseline_sbs_df_", "", f)
	s = paste0("baseline_sig_est_", s)

	sig_est_outfile = paste0(GLOBAL_SCRIPT_BASELINE_SIG_EST, s)
        system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_infile,
                      " -sf ", signatures_file,
                      " -of ", sig_est_outfile)
        )
	i = i + 1
}
