# make sure you activate the signature estimation library's conda environment before running this, or it won't work
# the command is:
# conda activate signature-estimation-py-env

setwd("../..")
source("projection_score.R")
library(optparse)

option_list = list(
	make_option(c("-t", "--tag"), type="character", default=NULL,
		help="tag to refer downstream scripts to the correct training sets", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

file_tag = opt$tag

if (is.null(file_tag)) {
	stop("No file_tag recieved (command line -t ). Please supply file tag.")
}

DEBUG_FLAG = TRUE

#file_tag = "test_run"

files = list.files(GLOBAL_SIZE_TEST_SBS_DIR, pattern= paste0(".*", file_tag, ".*") )
print(paste0(Sys.time(), "    found ", length(files), " files containing the tag: ", file_tag))

signatures_file = paste0(GLOBAL_DATA_DIR, "cosmic-signatures.tsv") #location of .tsv containing mutational signature frequencies (e.g. COSMIC signatures)


i = 1
for (f in files) {
	print( paste0(Sys.time(), "    ", i, "/",  length(files)) )

	sbs_infile = paste0(GLOBAL_SIZE_TEST_SBS_DIR, f)

	

	s = sub("trunc_sbs_df_", "", f)
	s = paste0("trunc_sig_est_", s)

	sig_est_outfile = paste0(GLOBAL_SIZE_TEST_SIG_EST_DIR, s)
        system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_infile,
                      " -sf ", signatures_file,
                      " -of ", sig_est_outfile)
        )
	i = i + 1
}
