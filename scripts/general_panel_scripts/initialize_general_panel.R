setwd("../..")

# initialization of test / train sets for the general panel experiments

# this will essentially just save a bunch of test/train split data and tag it in a standardized way so that the panel script can look it up

library(optparse)
source("projection_score.R")

option_list = list(
        make_option(c("-t", "--tag"), type="character", default=NULL,
                        help="tag to refer downstream scripts to the correct training sets", metavar="character"),
        #make_option(c("-g", "--siggroup"), type="numeric", default=NULL,
        #                help="specifies which group of signatures to generate test/train sets for. 1 or 2 (see documentation)", metavar="integer"),
        make_option(c("-n", "--numiters"), type="numeric", default=NULL,
                        help="specifies the number of test/train sets to generate.", metavar="integer")
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

global_sig_df = load_nz_sig_estimates(norm=TRUE)
samp_names = as.character(global_sig_df$Patient)

iter = 1
for (i in 1:num_iters) {
        file_name = paste0(GLOBAL_GP_TEST_TRAIN_DIR, "gp_test_train_", file_tag, "_it", i, ".rds")
        test_train = test_train_random_split(samp_names, .20)
        saveRDS(test_train, file_name)
}


