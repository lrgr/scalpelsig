# compute how well a panel distinguishes signature presence

library(pROC)
library(PRROC)

source("preprocess_windows.R")
source("win_sigprob_analysis.R")

# if laptop: "~/projects/hotspot_signature_panel/data/"
# if workstation: "/fs/cbcb-lab/mdml/users/franzese/projects/signature-panel/signature-panel/data/"

GLOBAL_DATA_DIR = "~/projects/hotspot_signature_panel/data/"
GLOBAL_CHR_MTX_DIR = paste0(GLOBAL_DATA_DIR, "individual_chromosome_matrices/")

GLOBAL_BASELINE_SIG_EST_DIR = paste0(GLOBAL_DATA_DIR, "BASELINE_PANELS/SIGNATURE_ESTIMATES/")

# if laptop: 3
# if workstation: 35

GLOBAL_NCORES = 3


source("GLOBAL_CONFIG.R")


######################
# load data          #
######################

# TODO: documentation - what script creates these files?
load_sig_estimate_df <- function(filename, norm=FALSE, replace_na=FALSE) {
        df = read.csv(file=filename, sep="\t", header=TRUE)
	colnames(df)[1] = "Patient"

        if (norm==TRUE) {
		rs = rowSums(df[,-1])
	        df[,-1] = df[,-1] / rs
		}

	# replace NA values with 0
	if (replace_na==TRUE) {
		df[is.na(df)] = 0
	}
	
	return(df)
}



load_logged_test_sets <- function(sig_num) {
	sig_keyword = paste0("_sig", sig_num, "_")
	f_ls = list.files(GLOBAL_LOGFILE_DIR, pattern=sig_keyword)
	
	ret = list()
	i = 1
	for (f in f_ls) {
		test_train = readRDS(paste0(GLOBAL_LOGFILE_DIR, f))
		test_set = test_train[[1]]
		ret[[i]] = test_set
		i = i + 1
	}
	
	return(ret)
}



################################################



# get scores for panel signature activity in a given test set
get_panel_score_df <- function(sig_num, panel_sig_df, test_set) {
	sig_name = paste0("Signature.", sig_num)

	test_df = data.frame(Patient=test_set, Score=numeric(length(test_set)))
	rownames(test_df) = test_set

	ps_in_panel = test_set[test_set %in% panel_sig_df$Patient] # patients in both test_set and have mutations in panel

	for (p in ps_in_panel) {
		test_df[p, "Score"] = panel_sig_df[panel_sig_df$Patient==p, sig_name]
	}

	return(test_df)
}


# add signature activity labels to panel score df
add_labels_to_panel_df <- function(panel_score_df, sig_activity_labels) {
	panel_score_df$Active = logical(nrow(panel_score_df))
	test_set = as.character(panel_score_df$Patient)

	for (p in test_set) {
		panel_score_df[p, "Active"] = sig_activity_labels[p]
	}
	return(panel_score_df)
}



# get labels for signature activity in all patients
# if a sample has at least 5% signature contribution, we say that it is active (modify thresh parameter to change)
get_sig_activity_labels <- function(sig_num, global_sig_df, thresh=0.05) {
	sig_name = paste0("Signature.", sig_num)

	ret = global_sig_df[ , sig_name] >= thresh
	names(ret) = global_sig_df$Patient
	return(ret)
}


# sig_num : integer for which COSMIC signature to assess
# test_set : vector of sample names
# infile : file name (with path) for signature estimate .tsv corresponding to the panel
compute_panel_auroc <- function(sig_num, test_set, infile, global_sig_df=NULL, debug=FALSE) {
	if (is.null(global_sig_df)) {
		global_sig_df = load_nz_sig_estimates(norm=TRUE)
	}

	sig_activity_labels = get_sig_activity_labels(sig_num, global_sig_df, 0.05)

	panel_sig_df = load_sig_estimate_df(infile, replace_na=TRUE)

	panel_score_df = get_panel_score_df(sig_num, panel_sig_df, test_set)
	panel_score_df = add_labels_to_panel_df(panel_score_df, sig_activity_labels)

	return(as.numeric( auc(roc(panel_score_df$Active, panel_score_df$Score, quiet=TRUE)) ) )
}


compute_panel_aupr <- function(sig_num, test_set, sig_est_infile, global_sig_df=NULL, debug=FALSE) {
	if (is.null(global_sig_df)) {
		global_sig_df = load_nz_sig_estimates(norm=TRUE)
	}

	sig_activity_labels = get_sig_activity_labels(sig_num, global_sig_df, 0.05)
	
	panel_sig_df = load_sig_estimate_df(sig_est_infile, replace_na=TRUE)

	panel_score_df = get_panel_score_df(sig_num, panel_sig_df, test_set)
	panel_score_df = add_labels_to_panel_df(panel_score_df, sig_activity_labels)

	fg = panel_score_df[ panel_score_df$Active==TRUE, "Score"]
	bg = panel_score_df[ panel_score_df$Active==FALSE, "Score"] 

	pr <- pr.curve(scores.class0=fg, scores.class1=bg) # from PRROC package, entries [[2]] and [[3]] give aupr computed by 2 different methods
	return(pr[[3]])
}



# this is used to evaulate the performance of the MSK IMPACT and WES panels against the proper test sets of each
# panel found by our framework.
panel_auroc_logged_test_sets <- function(sig_num, sig_est_infile, global_sig_df=NULL, debug=TRUE) {
	test_set_ls = load_logged_test_sets(sig_num)
	if (debug) { print(paste0("Loaded ", length(test_set_ls), " logged test sets for signature ", sig_num)) }

	auroc_vec = numeric(length(test_set_ls))
	for (i in 1:length(test_set_ls)) {
		auroc_vec[i] = compute_panel_auroc(sig_num, test_set_ls[[i]], sig_est_infile, global_sig_df)
	}

	if (debug) {
		print(paste0("median auroc: ", median(auroc_vec)))
		print(paste0("min : ", min(auroc_vec)))
		print(paste0("max : ", max(auroc_vec)))
	}
	return(auroc_vec)
}

# get auroc for each baseline panel
compute_baseline_auroc <- function(sig_num, test_set, global_sig_df=NULL, eval_mode="auroc", debug=FALSE) {
	if (is.null(global_sig_df)) {
		global_sig_df = load_nz_sig_estimates(norm=TRUE)
	}

	if (eval_mode != "auroc" & eval_mode != "aupr") {
		stop(paste0("Error in compute_baseline_auroc(): eval_mode most be either \'auroc\' or \'aupr\'. Recieved ", eval_mode, " instead."))
	}
	
	baseline_sig_est_files = list.files(GLOBAL_SCRIPT_BASELINE_SIG_EST)
	if (debug) { print(paste0("found ", length(baseline_sig_est_files), " baseline panels.")) }

	auc_vec = numeric(length(baseline_sig_est_files))

	sig_activity_labels = get_sig_activity_labels(sig_num, global_sig_df, 0.05)

	for (i in 1:length(auc_vec)) {
		if (debug) { print(paste0(i, "/", length(auc_vec))) }
		panel_sig_infile = paste0(GLOBAL_SCRIPT_BASELINE_SIG_EST, baseline_sig_est_files[i])
		panel_sig_df = load_sig_estimate_df(panel_sig_infile, replace_na=TRUE)
		
		panel_score_df = get_panel_score_df(sig_num, panel_sig_df, test_set)
		panel_score_df = add_labels_to_panel_df(panel_score_df, sig_activity_labels)

		if (eval_mode=="auroc") {
			curr_score = as.numeric( auc(roc(panel_score_df$Active, panel_score_df$Score, quiet=TRUE)) )
		} else if (eval_mode=="aupr") {
			fg = panel_score_df[ panel_score_df$Active==TRUE, "Score"]
			bg = panel_score_df[ panel_score_df$Active==FALSE, "Score"] 
			pr <- pr.curve(scores.class0=fg, scores.class1=bg) # from PRROC package, entries [[2]] and [[3]] give aupr computed by 2 different methods
			
			curr_score = pr[[3]]
		}

		auc_vec[i] = curr_score
		names(auc_vec)[i] = panel_sig_infile
	}

	return(auc_vec)
}




############### main ##################

main <- function() {
	# load inputs
	global_sig_df = load_nz_sig_estimates(norm=TRUE)
	#panel_sig_df = load_sig_estimate_df("~/projects/hotspot_signature_panel/data/BASELINE_PANELS/SIGNATURE_ESTIMATES/test_sigest.tsv")

	test_set = as.character(global_sig_df$Patient)[1:56] # 10% chunk of patients

	baseline_aucs = compute_baseline_auroc(2, test_set, global_sig_df)

	return(baseline_aucs)	

	# extract panel scores for test set
	#panel_score_df = get_panel_score_df(2, panel_sig_df, test_set)

	# compute signature activity labels
	#sig_activity_labels = get_sig_activity_labels(2, global_sig_df, 0.05)

	# add labels to test set
	#panel_score_df = add_labels_to_panel_df(panel_score_df, sig_activity_labels)

	# compute ROC
	#roc_obj = roc(panel_score_df$Active, panel_score_df$Score)
}



