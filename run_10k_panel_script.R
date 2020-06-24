# due to the heavy memory requirements of constructing the 10kb window panels,
# this script will find and save a panel for a single signature, and a single
# objective function. The results can then be processed down stream.

library(optparse)
source("projection_score.R")


# PARSER

option_list = list(
	make_option(c("-s", "--signum"), type="numeric", default=NULL,
			help="the signature that the panel should optimize for", metavar="character"),
	make_option(c("-t", "--tag"), type="character", default=NULL,
			help="tag to identify output files and input training set", metavar="character"),
	make_option(c("-o", "--objectivefn"), type="numeric", default=NULL, 
			help="controls which objective function the panel will use, can be 1, 2, or 3", metavar="character"),
	make_option(c("-w", "--windowsinpanel"), type="numeric", default=NULL,
			help="number of windows to be included in the panel", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser);



# control variables for the script

sig_num = opt$signum #TODO: pass this from parser
obj_fn_num = opt$objectivefn #TODO: pass this from parser
windows_in_panel = opt$windowsinpanel #TODO: pass this from parser
file_tag = opt$tag

if (is.null(sig_num) | is.null(obj_fn_num) | is.null(windows_in_panel) | is.null(file_tag)) {
	print("One of the control variables was not given: ")
	print(paste0("sig_num: ", sig_num))
	print(paste0("obj_fn_num: ", obj_fn_num))
	print(paste0("windows_in_panel: ", windows_in_panel))
	print(paste0("file_tag: ", file_tag))
	stop("Please supply all the control variables to run this script.")
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
if (obj_fn_num == 3) { obj_fn = obj_fn_class_balance }


# load test_train set from previous initialization
#TODO: make an error message for if the given file_tag does not correspond to a saved test_train set

test_train_file = paste0(GLOBAL_SCRIPT_TEST_TRAIN_DIR, "test_train_", file_tag, "_sig", sig_num, ".rds")

print(paste0(Sys.time(), "    loading test/train set from ", test_train_file))
test_train = readRDS(test_train_file)
train_set = test_train[[2]]

print(paste0(Sys.time(), "    loading score_mtx_ls"))

score_mtx_ls = load_ps_score_mtx_ls(sig_num, mode="10k")

print(paste0(Sys.time(), "    loading global_sig_df"))

global_sig_df = load_nz_sig_estimates(norm=TRUE)

#samp_names = as.character(global_sig_df$Patient)

#test_train = tt_stratified_split(sig_num, samp_names, .10, global_sig_df)
#train_set = test_train[[2]] #TODO: load this from file for standardization

print(paste0(Sys.time(), "    computing objective score vector"))

obj_vec = compute_obj_score_ps(sig_num, train_set, score_mtx_ls, obj_fn, global_sig_df)

panel_windows = names( top_n(obj_vec, windows_in_panel) )
window_file_name = paste0(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, "panel_windows_",  file_tag, "_sig", sig_num, "_obj", obj_fn_num, "_", timestamp_tag, ".txt")

print(paste0(Sys.time(), "    saving top ", panel_windows, " windows at ", window_file_name))
# CHECKPOINT 1: SAVING PANEL WINDOWS
write(panel_windows, file=window_file_name)

####################################################
####################################################

#THE PART BELOW IS NOW PERFORMED BY scripts/sbs_mtxs_from_10k_panel_windows.R

#score_mtx_ls = NULL
#obj_vec=NULL

#gc()

#print(paste0(Sys.time(), "    load_10k_sbs_arr_ls()"))
#sbs_arr_ls = load_10k_sbs_arr_ls()
#print(paste0(Sys.time(), "    finished loading."))

#panel_df = get_panel_sbs_df(panel_windows, sbs_arr_ls)

#sbs_outfile = paste0(GLOBAL_SCRIPT_PANEL_SBS_DIR, "panel_sbs_df_",  file_tag, "_sig", sig_num, "_obj", obj_fn_num, "_", timestamp_tag, ".tsv")
#print(paste0(Sys.time(), "    saving panel sbs matrix to ", sbs_outfile))

#save_panel_sbs_tsv(panel_df, sbs_outfile)

#print(paste0(Sys.time(), "    finished task: run_10k_panel_script.R -s ", sig_num, " -t ", file_tag, " -o ", obj_fn_num, " -w ", panel_windows))
