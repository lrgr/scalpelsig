# initialization & organization for run_10k_panel_script.R

# this will essentially just save a bunch of test/train split data and tag it in a standardized way so that the panel script can look it up

library(optparse)
source("projection_score.R")

option_list = list(
	make_option(c("-t", "--tag"), type="character", default=NULL, 
			help="tag to refer downstream scripts to the correct training sets", metavar="character"),
#	make_option(c("-g", "--siggroup"), type="numeric", default=NULL,
#			help="specifies which group of signatures to generate test/train sets for. 3 or 4 (see documentation)", metavar="integer"),
	make_option(c("-n", "--numiters"), type="numeric", default=NULL,
			help="specifies the number of test/train sets to generate for each signature.", metavar="integer")
	#make_option(c("-s", "--sizepanel"), type="numeric", default=NULL,
			#help="specifies the number of 10 kb windows in the panel.", metavar="integer")
); 


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);


# control variables

file_tag = opt$tag
#sig_group = opt$siggroup
num_iters = opt$numiters
#num_windows = opt$sizepanel


if (is.null(file_tag) | is.null(num_iters) ) {
	print("One of the control variables was not supplied: ")
	print(paste0("file_tag: ", file_tag))
	#print(paste0("sig_group: ", sig_group))
	print(paste0("num_iters: ", num_iters))
	#print(paste0("num_windows: ", num_windows)
	stop("Please supply all the control variables to run the script.")
}

##############################
##############################

#if (sig_group==1) { sig_vec = c(2, 3, 5, 9, 13, 16) } #class-balanced signatures (computed with COSMIC v2 1-30)
#if (sig_group==2) { sig_vec = c(1, 8, 18, 30) } # class imbalanced signatures (computed with COSMIC v2 1-30)
#if (sig_group==3) { sig_vec = c(2, 3, 13) } # NEW class-balanced signatures (computed with COSMIC v2 breast cancer sigs)
#if (sig_group==4) { sig_vec = c(8, 18, 30)} # NEW class-inbalanced signatures (COSMIC v2 breast cancer sigs)

sig_vec = c(2, 3, 8, 13, 18, 30)

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

