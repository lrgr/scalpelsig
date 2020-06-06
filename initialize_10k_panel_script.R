# initialization & organization for run_10k_panel_script.R

# this will essentially just save a bunch of test/train split data and tag it in a standardized way so that the panel script can look it up

library(optparse)
source("projection_score.R")

option_list = list(
	make_option(c("-t", "--tag"), type="character", default=NULL, 
			help="tag to refer downstream scripts to the correct training sets", metavar="character"),
	make_option(c("-g", "--siggroup"), type="numeric", default=NULL,
			help="specifies which group of signatures to generate test/train sets for. 1 or 2 (see documentation)", metavar="character"),
	make_option(c("-n", "--numiters"), type="numeric", default=NULL,
			help="specifies the number of test/train sets to generate for each signature.", metavar="character")
); 


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);


# control variables

file_tag = opt$tag
sig_group = opt$siggroup
num_iters = opt$numiters


if (is.null(file_tag) | is.null(sig_group) | is.null(num_iters)) {
	print("One of the control variables was not supplied: ")
	print(paste0("file_tag: ", file_tag))
	print(paste0("sig_group: ", sig_group))
	print(paste0("num_iters: ", num_iters))
	stop("Please supply all the control variables to run the script.")
}

##############################
##############################

if (sig_group==1) { sig_vec = c(2, 3, 5, 9, 13, 16) } #class-balanced signatures
if (sig_group==2) { sig_vec = c(1, 8, 18, 30) } # class imbalanced signatures

sig_vec = sort(rep(sig_vec, num_iters))

global_sig_df = load_nz_sig_estimates(norm=TRUE)
samp_names = as.character(global_sig_df$Patient)

iter = 1
for (sig in sig_vec) {
	file_name = paste0(GLOBAL_SCRIPT_TEST_TRAIN_DIR, "test_train_", file_tag, "_it", iter, "_sig", sig, ".rds")	
	test_train = tt_stratified_split(sig, samp_names, .10, global_sig_df)
	saveRDS(test_train, file_name)
	if (iter == num_iters) {
		iter = 1
	} else {
		iter = iter + 1
	}
}

