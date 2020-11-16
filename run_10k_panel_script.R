# due to the heavy memory requirements of constructing the 10kb window panels,
# this script will find and save a panel for a single signature, and a single
# objective function. The results can then be processed down stream.

library(optparse)
source("projection_score.R")


# PARSER

option_list = list(
	make_option(c("-s", "--signum"), type="numeric", default=NULL,
			help="the signature that the panel should optimize for", metavar="numeric"),
	make_option(c("-t", "--tag"), type="character", default=NULL,
			help="tag to identify output files and input training set", metavar="character"),
	make_option(c("-o", "--objectivefn"), type="numeric", default=NULL, 
			help="controls which objective function the panel will use, can be 1 or 2. 1 corresponds to the alpha=1 setting in the paper, 2 corresponds to alpha=0.5", metavar="numeric"),
	make_option(c("-w", "--windowsinpanel"), type="numeric", default=NULL,
			help="number of windows to be included in the panel", metavar="numeric"),
	make_option(c("-a", "--activitythreshold"), type="numeric", default=0.05,
			help="The percentage of mutations in a sample which determine whether the sample is defined as active for that signature. Default is 0.05 (i.e. if Sig X contributes at least 5% of the total mutations in a sample, that sample is active for Sig X).", metavar="numeric")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);



# control variables for the script

sig_num = opt$signum
obj_fn_num = opt$objectivefn 
windows_in_panel = opt$windowsinpanel
file_tag = opt$tag
act_thresh = opt$activitythreshold

if (is.null(sig_num) | is.null(obj_fn_num) | is.null(windows_in_panel) | is.null(file_tag)) {
	print("One of the control variables was not given: ")
	print(paste0("sig_num: ", sig_num))
	print(paste0("obj_fn_num: ", obj_fn_num))
	print(paste0("windows_in_panel: ", windows_in_panel))
	print(paste0("file_tag: ", file_tag))
	stop("Please supply all the control variables to run this script.")
}

if (is.null(act_thresh)) {
	stop("act_thresh was NULL, which means something is wrong.")
}


#####################################################
#####################################################


timestamp_tag = format(Sys.time(), "%d-%b-%Y_%H-%M")

print(paste0(Sys.time(), "    Running run_10k_panel_script.R - recieved the following inputs:"))
print(paste0("sig_num: ", sig_num))
print(paste0("obj_fn_num: ", obj_fn_num))
print(paste0("windows_in_panel: ", windows_in_panel))
print(paste0("file_tag: ", file_tag))
print(paste0("timestamp_tag: ", timestamp_tag))

print("proceeding with script.")
print("******************************************")

# decide objective function from command line input
if (obj_fn_num == 1) { obj_fn = simple_obj_fn }
if (obj_fn_num == 2) { obj_fn = obj_fn_sqrt }
#if (obj_fn_num == 3) { obj_fn = obj_fn_class_balance }


# load test_train set from previous initialization

test_train_file = paste0(GLOBAL_SCRIPT_TEST_TRAIN_DIR, "test_train_", file_tag, "_sig", sig_num, ".rds")

print(paste0(Sys.time(), "    loading test/train set from ", test_train_file))
test_train = readRDS(test_train_file)
train_set = test_train[[2]]

print(paste0(Sys.time(), "    loading score_mtx_ls"))

score_mtx_ls = load_ps_score_mtx_ls(sig_num, mode="10k")

print(paste0(Sys.time(), "    loading global_sig_df"))

global_sig_df = load_nz_sig_estimates(norm=TRUE)

#samp_names = as.character(global_sig_df$Patient)


print(paste0(Sys.time(), "    computing objective score vector"))

obj_vec = compute_obj_score_ps(sig_num, train_set, score_mtx_ls, obj_fn, global_sig_df, activation_thresh = act_thresh)

panel_windows = names( top_n(obj_vec, windows_in_panel) )
window_file_name = paste0(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, "panel_windows_",  file_tag, "_sig", sig_num, "_obj", obj_fn_num, "_nwin", windows_in_panel, "_", timestamp_tag, ".txt")

print(paste0(Sys.time(), "    saving top ", windows_in_panel, " windows at ", window_file_name))
# CHECKPOINT 1: SAVING PANEL WINDOWS
write(panel_windows, file=window_file_name)

