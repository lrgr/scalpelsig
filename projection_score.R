# this will implement the projection score idea

source("preprocess_windows.R")
source("assess_panel.R")


GLOBAL_OUT_DIR = "~/projects/hotspot_signature_panel/out/"
GLOBAL_PANEL_DIR = paste0(GLOBAL_OUT_DIR, "SIGNATURE_PANELS/")
GLOBAL_PANEL_SBS_DIR = paste0(GLOBAL_PANEL_DIR, "SBS_MATRICES/")
GLOBAL_PANEL_SIG_EST_DIR = paste0(GLOBAL_PANEL_DIR, "SIGNATURE_ESTIMATES/")
GLOBAL_LOGFILE_DIR = paste0(GLOBAL_OUT_DIR, "LOG_FILES/")


source("GLOBAL_CONFIG.R")

registerDoParallel(cores=GLOBAL_NCORES)


# given a training set, find windows of the genome which maximize this objective function:

# SUM_(over active samples) alpha_r    -    SUM_(over inactive samples) alpha_r

# to do this what needs to happen?
# 1) compute alpha_r for all signatures, regions, and samples
#	a) make an SBS count vector for each region & sample,
#	b) compute the scalar projection of this vector with each signature.
# 	c) save these somehow

# 2) mark samples with their signature activity
# 	assess_panel.R function get_sig_activity_labels() does this

# 3) given train and test sets, find the regions which maximize the objective function on the training set
# 	a) sum the alpha_r's across positive patients, sum alpha_r's across negative patients, take difference

# 4) evaluate AUROC on test set


###################
# loading data    #
###################


# load outputs of 1a)
load_sbs_array <- function(chrom, win_size=10^5) {
	return(readRDS(paste0(GLOBAL_CHR_MTX_DIR, "samp_window_sbs_array_winsize", win_size, "_chr", chrom, ".rds")))
}


# load outputs of 1b & c)
load_ps_score_mtx_ls <- function(sig_num) {
	return(readRDS(paste0(GLOBAL_CHR_MTX_DIR,  "projection_score_mtx_ls_sig", sig_num, ".rds")) )
}



#############################################
# quicker way to convert panel to sbs tsv   #
#############################################

# just find the relevant windows in the sbs array and sum the vectors together, then save it in the same format
# that sbs_df_from_mut_df() or whatever uses





#####################################
# 1a) compute 96SBS count array     #
#####################################


# this is for non-overlapping windows
# (i.e. [n_sample x n_window x n_mutation_categores])
# create a  array containing the 96-sbs mutation count for each sample and each region
get_sbs_counts_chrom_windows <- function(mut_df, chrom, win_size=10^6, debug=TRUE) {
	name_ls = list()
	
	samp_names = as.character(unique(mut_df$Patient))

	n_samp = length(samp_names)
	name_ls[[1]] = samp_names

	windows = windows_in_chromosome(mut_df, chrom, win_size) # function in preprocess_windows.R
	n_window = length(windows)
	win_names = window_names(windows, chrom)
	name_ls[[2]] = win_names

	n_categories = 96
	name_ls[[3]] = mut_types() # function in preprocess_windows.R

	chrom_df = mut_df[mut_df$Chromosome==chrom, c("Patient", "Chromosome", "Start.Position", "SBS_96")]

	arr = array(dim=c(n_samp, n_window, n_categories), dimnames=name_ls)


	if (debug) {print(paste0(Sys.time(), "    Starting main loop"))}

	for (s in 1:n_samp) {
		if (debug & (s==1 | s==2 | s==3 | s %% 50 == 0)) {
			print(paste0(s, "/", n_samp))
		}
		
		curr_samp = samp_names[s]
		samp_df = chrom_df[chrom_df$Patient == curr_samp,]
		
		for (w in 1:n_window) {
			curr_window = windows[[w]]
			curr_window_name = win_names[w]
			arr[curr_samp, curr_window_name, ] = get_96sbs_vec(samp_df, curr_window)
		}
	}


	if (debug) {print(paste0(Sys.time(), "    finished"))}	
	return(arr)
}



parallel_get_sbs_counts_chrom_windows <- function(mut_df, chrom, win_size, debug=TRUE) {
	name_ls = list()
	
	samp_names = as.character(unique(mut_df$Patient))

	n_samp = length(samp_names)
	name_ls[[1]] = samp_names

	windows = windows_in_chromosome(mut_df, chrom, win_size) # function in preprocess_windows.R
	n_window = length(windows)
	win_names = window_names(windows, chrom)
	name_ls[[2]] = win_names

	n_categories = 96
	name_ls[[3]] = mut_types() # function in preprocess_windows.R

	chrom_df = mut_df[mut_df$Chromosome==chrom, c("Patient", "Chromosome", "Start.Position", "SBS_96")]

	#arr = array(dim=c(n_samp, n_window, n_categories), dimnames=name_ls)
	
	acomb <- function(...) abind(..., along=3)
	if (debug) {print(paste0(Sys.time(), "    Starting main loop"))}

	arr <- foreach (s=1:n_samp, .combine='acomb', .multicombine=TRUE) %dopar% {
		curr_samp = samp_names[s]
		if (debug & (s==1 | s==2 | s==3 | s %% 50 == 0)) {
			print(paste0(s, "/", n_samp))
		}
		curr_samp = samp_names[s]

		#samp_name_ls = list()
		#samp_name_ls[[1]] = curr_samp
		#samp_name_ls[[2]] = win_names
		#samp_name_ls[[3]] = name_ls[[3]]
		#samp_arr = array(dim=c(1, n_window, n_categories), dimnames=samp_name_ls)
		samp_df = chrom_df[chrom_df$Patient == curr_samp,]
		
		samp_window_sbs_mtx = matrix(nrow=n_window, ncol=n_categories)
		#samp_window_sbs_mtx <- foreach(w=1:n_window, .combine=rbind) %do% {
		#	curr_window = windows[[w]]
		#	curr_window_name = win_names[w]
		#	get_96sbs_vec(samp_df, curr_window)
		#	#arr[curr_samp, curr_window_name, ] = get_96sbs_vec(samp_df, curr_window)
		#}
		for (w in 1:n_window) {
			curr_window = windows[[w]]
			curr_window_name = win_names[w]
			samp_window_sbs_mtx[w, ] = get_96sbs_vec(samp_df, curr_window)
		}
		rownames(samp_window_sbs_mtx) = win_names
		colnames(samp_window_sbs_mtx) = mut_types()
		samp_window_sbs_mtx
	}
	arr = aperm(arr, c(3,1,2)) # rearrange array so that samples are in first dimension
	dimnames(arr)[[1]] <- samp_names

	if (debug) {print(paste0(Sys.time(), "    finished")) }
	return(arr)
}





get_96sbs_vec <- function(samp_df, window, mut_categories=NULL) {
	if (is.null(mut_categories)) {
		mut_categories=mut_types()
	}

	start = window[1]
	end = window[2]
	muts = samp_df[samp_df$Start.Position >= start & samp_df$Start.Position < end,]$SBS_96

	vec = integer(96)
	names(vec) = mut_categories
	for (i in 1:length(mut_categories)) {
		cat = mut_categories[i]
		vec[i] = sum(muts==cat)
	}

	return(vec)
}



save_96sbs_arrays <- function(win_size=10^5, mut_df=NULL, debug=TRUE) {
	if (is.null(mut_df)) {
		mut_df = load_nz_mutation_df()
	}
	
	chrom_vec = c(1:22, "X", "Y")

	for (chrom in chrom_vec) {
		if (debug) { print(paste0("Starting get_sbs_counts_chrom_windows() for chromosome ", chrom)) }
		sbs_arr = get_sbs_counts_chrom_windows(mut_df, chrom, win_size=win_size) 
		filename = paste0(GLOBAL_CHR_MTX_DIR, "samp_window_sbs_array_winsize", win_size, "_chr", chrom, ".rds") 
		if (debug) { print(paste0("Saving sbs_array at: ", filename)) }
		saveRDS(sbs_arr, filename)
	}
}


parallel_save_96sbs_arrays <- function(win_size=10^5, mut_df=NULL, debug=TRUE) {
	if (is.null(mut_df)) {
		mut_df = load_nz_mutation_df()
	}
	
	chrom_vec = c(1:22, "X", "Y")

	for (chrom in chrom_vec) {
		if (debug) { print(paste0("Starting get_sbs_counts_chrom_windows() for chromosome ", chrom)) }
		sbs_arr = parallel_get_sbs_counts_chrom_windows(mut_df, chrom, win_size=win_size) 
		filename = paste0(GLOBAL_CHR_MTX_DIR, "samp_window_sbs_array_winsize", win_size, "_chr", chrom, ".rds") 
		if (debug) { print(paste0("Saving sbs_array at: ", filename)) }
		saveRDS(sbs_arr, filename)
	}
}




######################################
# compute scalar projection          #
######################################

# scalar projection of a in the direction of b
scalar_projection <- function(a, b) {
	unit_b = b / magnitude(b)
	ret = dot_prod(a, unit_b)
	return( ret )
}

magnitude <- function(v) {
	sum_sq = sum(v^2)
	return( sum_sq^(1/2) )
}

dot_prod <- function(v1, v2) {
	return( sum(v1 %*% v2) )
}

do_f <- function(f, ls, s) {
	ret = list()
	for (i in 1:length(ls)) {
		ret[[i]] = f(ls[[i]], s)
	}
	return(ret)
}



###############################################
# compute score for each window and sample    #
###############################################


score_sbs_array <- function(score_fn, sig_vec, sbs_array) {
	samp_names = dimnames(sbs_array)[[1]]

	window_names = dimnames(sbs_array)[[2]]

	score_mtx = matrix(nrow=length(samp_names), ncol=length(window_names))
	rownames(score_mtx) = samp_names
	colnames(score_mtx) = window_names

	for (i in 1:length(samp_names)) {
		for (j in 1:length(window_names)) {
			sbs_vec = sbs_array[i, j, ]
			score = score_fn(sbs_vec, sig_vec)
			score_mtx[i,j] = score
		}
	}

	return(score_mtx)
}

ps_sbs_array <- function(sig_vec, sbs_array) {
	return( score_sbs_array(scalar_projection, sig_vec, sbs_array) )
}




##################################################################
# 1b & c) compute projection score for each window and patient   #
##################################################################

projection_scores_by_sig <- function(sig_vec, win_size=10^5, debug=TRUE) {
	
	chrom_vec = c(1:22, "X", "Y") 

	score_mtx_ls = list()
	for (i in 1:length(chrom_vec)) {
		chrom = chrom_vec[i]

		if (debug) { print(paste0("projection_scores_by_sig() starting chromosome ", chrom)) }

		chr_sbs_arr = load_sbs_array(chrom, win_size)
		chr_score_mtx = ps_sbs_array(sig_vec, chr_sbs_arr)
		score_mtx_ls[[i]] = chr_score_mtx
	}

	return(score_mtx_ls)
}

save_projection_scores <- function(debug=TRUE) {
	sig_df = load_COSMIC_signatures()

	for (s in 1:nrow(sig_df)) {
		sig_vec = as.numeric(sig_df[s, c(-1)])
		if (debug) { print(paste0("running projection_scores_by_sig() for signature ", s)) }
		score_mtx_ls = projection_scores_by_sig(sig_vec)

		outfile = paste0(GLOBAL_CHR_MTX_DIR, "projection_score_mtx_ls_sig", s, ".rds")
		saveRDS(score_mtx_ls, outfile)
	}
}


############################
# test/train framework     #
############################


test_train_random_split <- function(sample_names, test_set_proportion) {
	num_in_test = as.integer(test_set_proportion * length(sample_names))

	test_inds = sample.int(length(sample_names), size=num_in_test)

	test_set = sample_names[test_inds]
	train_set = sample_names[-test_inds]

	ls = list()
	ls[[1]] = test_set
	ls[[2]] = train_set
	return(ls)
}


tt_stratified_split <- function(sig_num, sample_names, test_set_proportion, global_sig_df=NULL, activity_thresh=0.05, debug=FALSE) {
	if (debug) {print(paste0("Stratified split of test and train for Signature ", sig_num))}
	label_vec = get_sig_activity_labels(sig_num, global_sig_df, activity_thresh)
	
	n_samp = length(label_vec)
	n_pos_samp = sum(label_vec) # gives number of samples labeled TRUE

	pos_prop = n_pos_samp / n_samp
	neg_prop = 1 - pos_prop

	if (debug) { print(paste0("num samples: ", n_samp, "    num active samples: ", n_pos_samp, 
				  "    proportion of active samples: ", pos_prop)) }
	
	num_in_test = as.integer(test_set_proportion * length(sample_names))
	num_pos_in_test = as.integer(pos_prop * num_in_test)
	num_neg_in_test = num_in_test - num_pos_in_test

	pos_samp_names = names(which(label_vec==TRUE))
	neg_samp_names = names(which(label_vec==FALSE))

	pos_test = sample(pos_samp_names, num_pos_in_test)
	neg_test = sample(neg_samp_names, num_neg_in_test)

	if (debug) {print(paste0("test set size: ", length(pos_test) + length(neg_test), 
				 "     num active: ", length(pos_test),
				 "     num inactive: ", length(neg_test)))}

	pos_train = setdiff(pos_samp_names, pos_test)
	neg_train = setdiff(neg_samp_names, neg_test)

	test_set = c(pos_test, neg_test)
	train_set = c(pos_train, neg_train)

	ls = list()
	ls[[1]] = test_set
	ls[[2]] = train_set
	return(ls)
}


# get labels for signature activity in all patients
# if a sample has at least 5% signature contribution, we say that it is active (modify thresh parameter to change)
get_sig_activity_labels <- function(sig_num, global_sig_df=NULL, thresh=0.05) {
	if (is.null(global_sig_df)) {
		global_sig_df = load_nz_sig_estimates(norm=TRUE)
	}
	
	sig_name = paste0("Signature.", sig_num)

        ret = global_sig_df[ , sig_name] >= thresh
	names(ret) = global_sig_df$Patient
	return(ret)
}

# returns a vector which contains the objective function score for each region of the genome
compute_obj_score_ps <- function(sig_num, train_samps, score_mtx_ls=NULL, obj_fn=NULL, global_sig_df=NULL, debug=TRUE) {
	
	if (is.null(score_mtx_ls)) {
		if (debug) {print("compute_obj_score_ps() recieved no score_mtx_ls, loading...")}
		score_mtx_ls = load_ps_score_mtx_ls(sig_num)
	}

	if (is.null(obj_fn)) {
		if (debug) {print("compute_obj_score_ps() recieved no obj_fn, using default (simple_obj_fn).")}
		obj_fn = simple_obj_fn
	}

	if (is.null(global_sig_df)) {
		if (debug) {print("compute_obj_score_ps() recieved no global_sig_df, loading...")}
		global_sig_df = load_nz_sig_estimates(norm=TRUE)
	}

	activity_vec = get_sig_activity_labels(sig_num, global_sig_df, 0.05)

	train_act_vec = activity_vec[train_samps]
	pos_samps = names(train_act_vec[train_act_vec==TRUE])
	neg_samps = names(train_act_vec[train_act_vec==FALSE])

	n_windows = 0
	window_names = character(0)
	for (i in 1:length(score_mtx_ls)) {
		n_windows = n_windows + ncol(score_mtx_ls[[i]])
		window_names = c(window_names, colnames(score_mtx_ls[[i]]) )
	}

	obj_vec = numeric(n_windows)
	names(obj_vec) = window_names

	# main loop, compute objective function on each window
	for (i in 1:length(score_mtx_ls)) {
		
		score_mtx = score_mtx_ls[[i]]
		for (j in 1:ncol(score_mtx)) {
			window = colnames(score_mtx)[j]
			obj_vec[window] = obj_fn(window, pos_samps, neg_samps, score_mtx)
		}
	}
	
	return(obj_vec)
}

# where score_mtx is [samples x regions]
simple_obj_fn <- function(region, pos_samps, neg_samps, score_mtx) {
	t1 = sum(score_mtx[pos_samps, region])
	t2 = sum(score_mtx[neg_samps, region])

	return(t1 - t2)
}

obj_fn_sqrt <- function(region, pos_samps, neg_samps, score_mtx) {
	t1 = sum(score_mtx[pos_samps, region]^(1/2))
	t2 = sum(score_mtx[neg_samps, region]^(1/2))

	return(t1 - t2)
}

obj_fn_class_balance <- function(region, pos_samps, neg_samps, score_mtx) {
	lambda = length(pos_samps) / length(neg_samps)

	t1 = sum(score_mtx[pos_samps, region]^(1/2))
	t2 = sum(score_mtx[neg_samps, region]^(1/2))

	return(t1 - (lambda * t2))
}



# estimate pval by comparing to baseline panel AUROC
est_pval <- function(obj_auroc, sig_num, test_set, global_sig_df=NULL, debug=FALSE) {
	if (debug) {print("running compute_baseline_auroc()")}
	baseline_aurocs = compute_baseline_auroc(sig_num, test_set, global_sig_df)

	if (debug) {print(paste0("Estimating p-val from ", length(baseline_aurocs), " random panels."))}
	if (debug) {print(paste0("Summary statistics for random panel AUROCs - ",
				"max: ", round(max(baseline_aurocs), 4), "    ",
			 	"med: ", round(median(baseline_aurocs), 4)) ) }
	n_better = sum(baseline_aurocs >= obj_auroc)

	pval = n_better / length(baseline_aurocs)
	if (debug) {print(paste0("estimated pval: ", pval))}
	return(pval)
}



###########################
# random helper functions #
###########################


# save panel sbs df as a tsv
save_panel_sbs_tsv <- function(sbs_df, outfile) {
	        write.table(sbs_df, file=outfile, sep="\t", quote=FALSE)
}







############# main ################

sig_specific_main <- function(sig_num, mut_df=NULL, global_sig_df=NULL, debug=TRUE) {
	#ret = matrix(nrow = 30, ncol=3)
	#rownames(ret) = 1:30
	#colnames(ret) = c("simple_obj_fn", "obj_fn_sqrt", "obj_fn_class_balance")
	
	if (is.null(mut_df)) {
		if (debug) { print(paste0(Sys.time(),"    loading mut_df")) }
		mut_df = load_nz_mut_df_with_sigprob()
	}

	if (is.null(global_sig_df)) {
		if (debug) { print(paste0(Sys.time(), "    loading global_sig_df")) }
		global_sig_df = load_nz_sig_estimates(norm=TRUE)
	}

	samp_names = as.character(global_sig_df$Patient)

	if (debug) { print(paste0(Sys.time(), "    splitting test and train")) }
	#test_train = test_train_random_split(samp_names, .10)
	test_train = tt_stratified_split(sig_num, samp_names, .10, global_sig_df, debug=debug)
	test_set = test_train[[1]]
	train_set = test_train[[2]]
	timestamp_tag = format(Sys.time(), "%d-%b-%Y_%H-%M")

	test_train_outfile = paste0(GLOBAL_LOGFILE_DIR, "test_train_sig", sig_num, "_", timestamp_tag, ".rds")
	if (debug) {print(paste0("saving test & train set to : ", test_train_outfile))}
	saveRDS(test_train, test_train_outfile)



	if (debug) { print(paste0(Sys.time(), "    starting main loop")) }
	#for (sig_num in 1:30) {
		# if test set contains all positive or all negative examples, go to next signature
		activity_vec = get_sig_activity_labels(sig_num, global_sig_df, 0.05)
		test_activity = activity_vec[test_set]
		check = sum(test_activity)
		if (check==0 | check==length(test_set)) {
			if (debug) {print("bad test set")}
			next
		}
		
		if (debug) { print(paste0(Sys.time(), "    sig_num: ", sig_num, " loading score_mtx_ls")) }
		score_mtx_ls = load_ps_score_mtx_ls(sig_num)

		if (debug) { print(paste0(Sys.time(), "    compute_obj_score_ps() with simple_obj_fn")) }
		obj_v1 = compute_obj_score_ps(sig_num, train_set, score_mtx_ls, obj_fn=simple_obj_fn, global_sig_df=global_sig_df)
		if (debug) { print(paste0(Sys.time(), "    compute_obj_score_ps() with obj_fn_sqrt")) }
		obj_v2 = compute_obj_score_ps(sig_num, train_set, score_mtx_ls, obj_fn=obj_fn_sqrt, global_sig_df=global_sig_df)
		if (debug) { print(paste0(Sys.time(), "    compute_obj_score_ps() with obj_fn_class_balance")) }
		obj_v3 = compute_obj_score_ps(sig_num, train_set, score_mtx_ls, obj_fn=obj_fn_class_balance, global_sig_df=global_sig_df)

		if (debug) { print(paste0(Sys.time(), "    getting panel windows")) }
		panel_1_windows = names( top_n(obj_v1, 27) )
		panel_2_windows = names( top_n(obj_v2, 27) )
		panel_3_windows = names( top_n(obj_v3, 27) )

		if (debug) {print(paste0(Sys.time(), "    generating panel dfs")) }
		panel_1_df = select_window_muts(panel_1_windows, mut_df)
		panel_2_df = select_window_muts(panel_2_windows, mut_df)
		panel_3_df = select_window_muts(panel_3_windows, mut_df)

		sbs_outfile_1 = paste0(GLOBAL_PANEL_SBS_DIR, "ps_obj1_sig", sig_num, "_panel_sbs_", timestamp_tag, ".tsv")
		sbs_outfile_2 = paste0(GLOBAL_PANEL_SBS_DIR, "ps_obj2_sig", sig_num, "_panel_sbs_", timestamp_tag, ".tsv")
		sbs_outfile_3 = paste0(GLOBAL_PANEL_SBS_DIR, "ps_obj3_sig", sig_num, "_panel_sbs_", timestamp_tag, ".tsv")

		if (debug) {print(paste0(Sys.time(), "    saving to:")); print(sbs_outfile_1); print(sbs_outfile_2); print(sbs_outfile_3)}
		save_panel_sbs_tsv(sbs_df_from_mut_df(panel_1_df), sbs_outfile_1)
		save_panel_sbs_tsv(sbs_df_from_mut_df(panel_2_df), sbs_outfile_2)
		save_panel_sbs_tsv(sbs_df_from_mut_df(panel_3_df), sbs_outfile_3)

		sig_est_outfile_1 = paste0(GLOBAL_PANEL_SIG_EST_DIR, "ps_obj1_sig", sig_num, "_est_", timestamp_tag, ".tsv")
		sig_est_outfile_2 = paste0(GLOBAL_PANEL_SIG_EST_DIR, "ps_obj2_sig", sig_num, "_est_", timestamp_tag, ".tsv")
		sig_est_outfile_3 = paste0(GLOBAL_PANEL_SIG_EST_DIR, "ps_obj3_sig", sig_num, "_est_", timestamp_tag, ".tsv")

		if (debug) {print("running signature estimator for panel 1")}
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_outfile_1, 
			      " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
			      " -of ", sig_est_outfile_1)
		)

		if (debug) {print("running signature estimator for panel 2")}
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_outfile_2, 
			      " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
			      " -of ", sig_est_outfile_2)
		)

		if (debug) {print("running signature estimator for panel 3")}
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_outfile_3, 
			      " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
			      " -of ", sig_est_outfile_3)
		)

		if (debug) {print(paste0(Sys.time(), "    running compute_panel_auroc()"))}
		auroc_1 = compute_panel_auroc(sig_num, test_set, sig_est_outfile_1, global_sig_df=global_sig_df)
		auroc_2 = compute_panel_auroc(sig_num, test_set, sig_est_outfile_2, global_sig_df=global_sig_df)
		auroc_3 = compute_panel_auroc(sig_num, test_set, sig_est_outfile_3, global_sig_df=global_sig_df)
		best_auroc = max(c(auroc_1, auroc_2, auroc_3))

		if (debug) {print(paste0(Sys.time(), "    results for signature ", sig_num, ": ")); print(auroc_1); print(auroc_2); print(auroc_3)}

		if (debug) {print(paste0(Sys.time(), "    running est_pval()")) }
		pval = est_pval(best_auroc, sig_num, test_set, global_sig_df, debug=debug)
		

		if (debug) { print("saving results file")}
		safe_ret = matrix(c(auroc_1, auroc_2, auroc_3, best_auroc, pval))
		rownames(safe_ret) = c("simple_obj_fn", "obj_fn_sqrt", "obj_fn_class_balance", "best", "est pval")
		colnames(safe_ret) = c("test_set_auroc")
		safe_ret_outfile = paste0(GLOBAL_PANEL_DIR, "ps_panel_results_sig", sig_num, "_", timestamp_tag, ".txt")
		write.table(safe_ret, file=safe_ret_outfile)

		#ret[sig_num, ] = c(auroc_1, auroc_2, auroc_3)	
	#}

	#auc_mtx_outfile = paste0(GLOBAL_OUT_DIR, "ps_auc_mtx_", timestamp_tag, ".txt")
	#write.table(ret, file=auc_mtx_outfile)
	
}




main <- function(mut_df=NULL, global_sig_df=NULL, debug=TRUE) {
	ret = matrix(nrow = 30, ncol=3)
	rownames(ret) = 1:30
	colnames(ret) = c("simple_obj_fn", "obj_fn_sqrt", "obj_fn_class_balance")
	
	if (is.null(mut_df)) {
		if (debug) { print(paste0(Sys.time(),"    loading mut_df")) }
		mut_df = load_nz_mut_df_with_sigprob()
	}

	if (is.null(global_sig_df)) {
		if (debug) { print(paste0(Sys.time(), "    loading global_sig_df")) }
		global_sig_df = load_nz_sig_estimates(norm=TRUE)
	}

	samp_names = as.character(global_sig_df$Patient)

	if (debug) { print(paste0(Sys.time(), "    splitting test and train")) }
	test_train = test_train_random_split(samp_names, .10)
	test_set = test_train[[1]]
	train_set = test_train[[2]]
	timestamp_tag = format(Sys.time(), "%d-%b-%Y_%H-%M")




	if (debug) { print(paste0(Sys.time(), "    starting main loop")) }
	for (sig_num in 1:30) {
		# if test set contains all positive or all negative examples, go to next signature
		activity_vec = get_sig_activity_labels(sig_num, global_sig_df, 0.05)
		test_activity = activity_vec[test_set]
		check = sum(test_activity)
		if (check==0 | check==length(test_set)) {
			if (debug) {print("bad test set")}
			next
		}
		
		if (debug) { print(paste0(Sys.time(), "    sig_num: ", sig_num, " loading score_mtx_ls")) }
		score_mtx_ls = load_ps_score_mtx_ls(sig_num)

		if (debug) { print(paste0(Sys.time(), "    compute_obj_score_ps() with simple_obj_fn")) }
		obj_v1 = compute_obj_score_ps(sig_num, train_set, score_mtx_ls, obj_fn=simple_obj_fn, global_sig_df=global_sig_df)
		if (debug) { print(paste0(Sys.time(), "    compute_obj_score_ps() with obj_fn_sqrt")) }
		obj_v2 = compute_obj_score_ps(sig_num, train_set, score_mtx_ls, obj_fn=obj_fn_sqrt, global_sig_df=global_sig_df)
		if (debug) { print(paste0(Sys.time(), "    compute_obj_score_ps() with obj_fn_class_balance")) }
		obj_v3 = compute_obj_score_ps(sig_num, train_set, score_mtx_ls, obj_fn=obj_fn_class_balance, global_sig_df=global_sig_df)

		if (debug) { print(paste0(Sys.time(), "    getting panel windows")) }
		panel_1_windows = names( top_n(obj_v1, 27) )
		panel_2_windows = names( top_n(obj_v2, 27) )
		panel_3_windows = names( top_n(obj_v3, 27) )

		if (debug) {print(paste0(Sys.time(), "    generating panel dfs")) }
		panel_1_df = select_window_muts(panel_1_windows, mut_df)
		panel_2_df = select_window_muts(panel_2_windows, mut_df)
		panel_3_df = select_window_muts(panel_3_windows, mut_df)

		sbs_outfile_1 = paste0(GLOBAL_PANEL_SBS_DIR, "ps_obj1_sig", sig_num, "_panel_sbs_", timestamp_tag, ".tsv")
		sbs_outfile_2 = paste0(GLOBAL_PANEL_SBS_DIR, "ps_obj2_sig", sig_num, "_panel_sbs_", timestamp_tag, ".tsv")
		sbs_outfile_3 = paste0(GLOBAL_PANEL_SBS_DIR, "ps_obj3_sig", sig_num, "_panel_sbs_", timestamp_tag, ".tsv")

		if (debug) {print(paste0(Sys.time(), "    saving to:")); print(sbs_outfile_1); print(sbs_outfile_2); print(sbs_outfile_3)}
		save_panel_sbs_tsv(sbs_df_from_mut_df(panel_1_df), sbs_outfile_1)
		save_panel_sbs_tsv(sbs_df_from_mut_df(panel_2_df), sbs_outfile_2)
		save_panel_sbs_tsv(sbs_df_from_mut_df(panel_3_df), sbs_outfile_3)

		sig_est_outfile_1 = paste0(GLOBAL_PANEL_SIG_EST_DIR, "ps_obj1_sig", sig_num, "_est_", timestamp_tag, ".tsv")
		sig_est_outfile_2 = paste0(GLOBAL_PANEL_SIG_EST_DIR, "ps_obj2_sig", sig_num, "_est_", timestamp_tag, ".tsv")
		sig_est_outfile_3 = paste0(GLOBAL_PANEL_SIG_EST_DIR, "ps_obj3_sig", sig_num, "_est_", timestamp_tag, ".tsv")

		if (debug) {print("running signature estimator for panel 1")}
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_outfile_1, 
			      " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
			      " -of ", sig_est_outfile_1)
		)

		if (debug) {print("running signature estimator for panel 2")}
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_outfile_2, 
			      " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
			      " -of ", sig_est_outfile_2)
		)

		if (debug) {print("running signature estimator for panel 3")}
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_outfile_3, 
			      " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
			      " -of ", sig_est_outfile_3)
		)

		if (debug) {print(paste0(Sys.time(), "    running compute_panel_auroc()"))}
		auroc_1 = compute_panel_auroc(sig_num, test_set, sig_est_outfile_1, global_sig_df=global_sig_df)
		auroc_2 = compute_panel_auroc(sig_num, test_set, sig_est_outfile_2, global_sig_df=global_sig_df)
		auroc_3 = compute_panel_auroc(sig_num, test_set, sig_est_outfile_3, global_sig_df=global_sig_df)

		if (debug) {print(paste0(Sys.time(), "    results for signature ", sig_num, ": ")); print(auroc_1); print(auroc_2); print(auroc_3)}

		if (debug) { 
			print("saving incomplete results file")
			safe_ret = matrix(c(auroc_1, auroc_2, auroc_3))
			rownames(safe_ret) = c("simple_obj_fn", "obj_fn_sqrt", "obj_fn_class_balance")
			colnames(safe_ret) = c("test_set_auroc")
			safe_ret_outfile = paste0(GLOBAL_PANEL_DIR, "ps_panel_inprogress_results_sig", sig_num, "_", timestamp_tag, ".txt")
			write.table(safe_ret, file=safe_ret_outfile)
		}

		ret[sig_num, ] = c(auroc_1, auroc_2, auroc_3)	
	}

	auc_mtx_outfile = paste0(GLOBAL_OUT_DIR, "ps_auc_mtx_", timestamp_tag, ".txt")
	write.table(ret, file=auc_mtx_outfile)
}



run_this <- function() {
	mut_df = load_nz_mut_df_with_sigprob()
	global_sig_df = load_nz_sig_estimates(norm=TRUE)
	
	for (s in c(2,2,2,2,2, 3,3,3,3,3, 5,5,5,5,5, 9,9,9,9,9, 13,13,13,13,13, 16,16,16,16,16)) {
		sig_specific_main(s, mut_df, global_sig_df)
	}
}
