# script for statistical analysis of regional signature probabilities

# if laptop: "~/projects/hotspot_signature_panel/data/"
# if workstation: "/fs/cbcb-lab/mdml/users/franzese/projects/signature-panel/signature-panel/data/"

GLOBAL_DATA_DIR = "~/projects/hotspot_signature_panel/data/"
GLOBAL_CHR_MTX_DIR = paste0(GLOBAL_DATA_DIR, "individual_chromosome_matrices/")

# if laptop: 3
# if workstation: 35

GLOBAL_NCORES = 3

#####################
# load data         #
#####################

load_nz_aggregated_window_sigprob <- function() {
	filename = paste0(GLOBAL_DATA_DIR, "nz_aggregated_window_sigprob.tsv")
	df = read.csv(filename, sep = "\t", header=TRUE)
	return(df)
}


load_nz_mutwindow_chr_10k <- function(chr) {
	filename = paste0(GLOBAL_CHR_MTX_DIR, "nz_mutwindow_10000_sigprob_chr_", chr, ".tsv")
	
	df = data.table::fread(file=filename, header=FALSE)
	sigs = df[[1]]
	regions = scan(file=filename, what=character(), nlines=1)

	mtx = as.matrix(df[ , -1]) # cut off column containing signature names
	rownames(mtx) = sigs
	colnames(mtx) = regions
	
	return(mtx)
}

load_nz_mutwindow_chr_100k <- function(chr) {
	filename = paste0(GLOBAL_CHR_MTX_DIR, "nz_mutwindow_1e+05_sigprob_chr_", chr, ".tsv")

	df = data.table::fread(file=filename, header=FALSE)
	sigs = df[[1]]
	regions = scan(file=filename, what=character(), nlines=1)

	mtx = as.matrix(df[ , -1]) # cut off column containing signature names
	rownames(mtx) = sigs
	colnames(mtx) = regions
	
	return(mtx)	
}

load_fc_mutwindow_chr_10k <- function(chr) {
	filename = paste0(GLOBAL_CHR_MTX_DIR, "fc_mtx_mutwindow_10k_chr_", chr, ".tsv")

	df = data.table::fread(file=filename, header=FALSE)
	sigs = df[[1]]
	regions = scan(file=filename, what=character(), nlines=1)

	mtx = as.matrix(df[ , -1]) # cut off column containing signature names
	rownames(mtx) = sigs
	colnames(mtx) = regions
	
	return(mtx)	
}


load_fc_mutwindow_chr_100k <- function(chr) {
	filename = paste0(GLOBAL_CHR_MTX_DIR, "fc_mtx_mutwindow_100k_chr_", chr, ".tsv")

	df = data.table::fread(file=filename, header=FALSE)
	sigs = df[[1]]
	regions = scan(file=filename, what=character(), nlines=1)

	mtx = as.matrix(df[ , -1]) # cut off column containing signature names
	rownames(mtx) = sigs
	colnames(mtx) = regions
	
	return(mtx)	
}



# fread from library(data.table) works for reading but it's over 3 million rows and it crashes when you attempt access of entries
# don't use these, it seems they are too big for read.csv() to handle them without getting stuck.
#load_nz_mutwindow_10k <- function() {
#	df = read.csv("~/projects/hotspot_signature_panel/data/nz_mutwindow_10000_sigprob.tsv",
#		      sep = "\t", header=TRUE)
#	return(df)
#}

#load_nz_mutwindow_1k <- function() {
#	df = read.csv("~/projects/hotspot_signature_panel/data/nz_mutwindow_1000_sigprob.tsv",
#		      sep = "\t", header=TRUE)
#	return(df)
#}

##############

# mutation df
load_nz_mut_df_with_sigprob <- function() {
	filename = paste0(GLOBAL_DATA_DIR, "nz_mut_df_with_sig_prob.tsv")
	df = read.csv(file=filename, sep="\t", header=TRUE)
	return(df)
}

# signature df
load_nz_sig_estimates <- function(norm=FALSE) {
	filename = paste0(GLOBAL_DATA_DIR, "nz_signature_estimation.tsv")
	df = read.csv(file=filename, sep="\t", header=TRUE)
	colnames(df)[1] = "Patient"

	if (norm==TRUE) {
		rs = rowSums(df[,-1])
		df[,-1] = df[,-1] / rs
	}

	return(df)
}


#####################
# bookkeeping       #
#####################

print_loop_progress <- function(curr, end, interval= -1) {
	if (interval== -1) {
		# then use default setting
		interval = round(end / 10)
	}
	if (curr==1 | curr==2 | curr==3 | curr %% interval == 0 ) {
		print(paste0(curr, "/", end))
	}
}


##############################
# general helper functions   #
##############################

top_n <- function(x, n=30) {
	if (length(x) < n) {
		n = length(x)
	}
	return( sort(x, decreasing=TRUE)[1:n] )
}

top_n_no_overlap <- function(x, n, debug=TRUE) {
	s = sort(x, decreasing=TRUE)
	if (length(s)==0) {
		return(numeric(0))
	}	

	result_indices = numeric(0)
	all_regs = names(s)
	regs_so_far = character(0)
	for (i in 1:length(all_regs)) {
		# only want regions that don't overlap what we've seen so far
		# if (debug) { print(all_regs[i]) }
		if(is_overlapping(all_regs[i], regs_so_far) == FALSE) {
			result_indices = c(result_indices, i)
			regs_so_far = names(s[result_indices])
		}

		if (length(result_indices) >= n) {
			break
		}
	}

	return(s[result_indices])
}

# don't use lmao
whichpart <- function(x, n=30) {
	nx = length(x)
	p = nx-n
	xp = sort(x, partial=p)[p]
	return(which(x > xp))
}


#######################################
# compare to background score         #
#######################################

expected_window_score <- function(df, win_num) {
	#print(win_num)
	muts_in_win = sum(df[ , win_num]) # number of mutations in the window
	#print("uhh")
	total_score_vec = rowSums(df) # summed signature scores across all windows
	total_muts = sum(df) # number of mutations total
	bg_prob = total_score_vec / total_muts # normalize

	expected_win_score = muts_in_win * bg_prob # expected score if there was no regional variation
	return(expected_win_score)
}

# matrix with fold change from expected values for each window
fold_change_mtx <- function(df, debug=TRUE) {
	ret = matrix(nrow=nrow(df), ncol=ncol(df))
	rownames(ret) = rownames(df)
	colnames(ret) = colnames(df)

	for (i in 1:ncol(df)) {
		if (debug) { print_loop_progress(i, ncol(df)) }	
		expected = expected_window_score(df, i)
		fold_change_vec = df[ , i] / expected
		ret[ , i] = fold_change_vec
	}
	return(ret)
}

# same as above function but computes sums row-wise, which maybe should optimize
library(doParallel)
registerDoParallel(cores=GLOBAL_NCORES)

par_fc_mtx <- function(mtx, debug=TRUE) {
	#ret = matrix(nrow=nrow(df), ncol=ncol(df))

	ret <- foreach (i = 1:ncol(mtx), .combine=cbind) %dopar% {
		#if (debug) { print_loop_progress(i, ncol(df)) }	
		expected = expected_window_score(mtx, i)
		fold_change_vec = mtx[ , i] / expected
		fold_change_vec
	}
	
	rownames(ret) = rownames(mtx)
	colnames(ret) = colnames(mtx)
	
	return(ret)
}
		



# report windows with highest fold change for each signature
report_top_windows <- function(fold_change_df, n=5) {
	ls = list()
	
	for (i in 1:nrow(fold_change_df)) {
		curr_row = fold_change_df[i, ] # get row of matrix for signature i
		ls[rownames(fold_change_df)[i]] = list(top_n(curr_row, n)) # append top n entries in row to ls
	}
	return(ls)
}

# report top n windows with highest fold change for a given signature
active_sig_windows <- function(fold_change_df, signature, n=50) {
	sig_row = fold_change_df[signature, ] # get row of matrix for given signature
	return(top_n_no_overlap(sig_row, n))
}



##############################################################################
# assessing relationship between regional and global signature abundance     #
##############################################################################


# given list of windows, acquire truncated mutation df that contains only the mutations within the windows

select_window_muts <- function(window_strs, mut_df) {
	#df = muts_in_window(window_strs[1], mut_df)
	
	#for (i in 2:length(window_strs)) {
	#	curr_muts = muts_in_window(window_strs[i], mut_df)
	#	df = rbind(df, curr_muts)
	#}
	#return(df)

	return(muts_in_windows(window_strs, mut_df))
}


### helper functions ###

#library(stringr)
# given the name of a window, return a list containing its chromosome, start base, and end base
parse_window_str <- function(window_str) {
	vec = strsplit(window_str, "_")[[1]]
	chr_str = vec[1]
	new_chr_str = substr(chr_str, 4, nchar(chr_str)) # remove first 3 letters "chr" from chromosome string
	vec[1] = new_chr_str
	return(vec)
}

# given name of window and mutation df, return truncated df containing mutations in that window
muts_in_window <- function(window_str, mut_df) {
	info = parse_window_str(window_str)
	chr = info[1]
	start = info[2]
	end = info[3]

	# df of all mutations that fall in the window
	df = mut_df[mut_df$Chromosome==chr & mut_df$Start.Position >= start & mut_df$Start.Position < end, ] 
	return(df)
}

### smarter way ###

parse_window_strs <- function(strs) {
	ls = list()
	for (i in 1:length(strs)) {
		v <- parse_window_str(s)
		ls[[i]] = v
	}
	return(ls)
}

# get rows of mutations which fall in the given window
window_muts_help <- function(window_str, mut_df) {
	window_str_vec = parse_window_str(window_str)

	chr = window_str_vec[1]
	start = window_str_vec[2]
	end = window_str_vec[3]

	row_indices = which(mut_df$Chromosome==chr & mut_df$Start.Position >= start & mut_df$Start.Position < end)
	return(row_indices)
}

muts_in_windows <- function(window_strs, mut_df) {
	row_indices = integer()
	

	for (s in window_strs) {
		new_rows = window_muts_help(s, mut_df)
		row_indices = union(row_indices, new_rows)
	}

	return( mut_df[row_indices, ] )
}




# given a string representing a target region t
# and a char vector of reference regions r_vec
# determine whether t overlaps any of the reference regions

is_overlapping <- function(t, r_vec, debug=FALSE) {
	t_info = parse_window_str(t)
	t_chr = t_info[[1]]
	t_start = as.integer(t_info[[2]])
	t_end = as.integer(t_info[[3]])

	for (r in r_vec) {
		r_info = parse_window_str(r)
		r_chr = r_info[[1]]
		r_start = as.integer(r_info[[2]])
		r_end = as.integer(r_info[[3]])

		if ((t_chr == r_chr) & overlap_help(t_start, t_end, r_start, r_end) ) {
			if (debug) { 
				print(paste0("Overlap found between ", t, " ; ", r)) 
				return(r)
			}
			return(TRUE)
		}
	}

	return(FALSE)
}

overlap_help <- function(t_start, t_end, r_start, r_end) {
	res = abs(t_start - r_start)
	window_size = max(t_end - t_start, r_end - r_start)
	return(res < window_size)
}

# number of elements in the given vector of regions that overlap with some other region in the vector
count_overlaps <- function(r_vec, debug=TRUE) {
	count = 0
	#debug_ls = list()
	for (i in 1:length(r_vec)) {
		r = r_vec[i]
		if (is_overlapping(r, r_vec[-i]) != FALSE) {
			
			#the_thing = is_overlapping(r, r_vec[-i], TRUE)
			#thing_ind = which(r_vec==the_thing)[1]

			#debug_ls[[count+1]] = c(r, the_thing)

			#check = is_overlapping(the_thing, r)
			#if (check==FALSE) {
			#	print("ok what is going on here")
			#	print(paste0("r: ", r))
			#	print(paste0("the thing: ", the_thing))
			#	return(NULL)
			#}
			
			count = count + 1
		}
	}
	#if (debug) { print(debug_ls); return(debug_ls) }
	return(count)
}


#debug_fn <- function(t, ls) {
#	count = 0
#	for (i in 1:length(ls[-t])) {
#		if (setequal(ls[[t]], ls[-t][[i]])) {
#			count = count + 1
#		}
#	}
#	return(count)
#}

########################################################
# identify signature score for each patient in a panel #
########################################################


get_10k_panel_muts <- function(signature, mut_df, n=50, debug=TRUE) {
	panel_scores = cross_chrom_top_windows_10k(signature, n=50, debug=TRUE)
	panel_windows = names(panel_scores)
	panel_df = select_window_muts(panel_windows, mut_df)
	return(panel_df)
}

# sum the probability of the given signature across mutations in the panel
patient_panel_sig_scores <- function(panel_df, signature, norm=TRUE) {
	patient_vec = unique(panel_df$Patient)
	sig_name = paste0("Signature.", signature)

	score_vec = numeric(length(patient_vec))
	names(score_vec) = patient_vec

	for (i in 1:length(patient_vec)) {
		p = patient_vec[i]
		patient_muts = panel_df[panel_df$Patient==p,]
		n_muts = nrow(patient_muts)
		raw_score = sum(patient_muts[ , sig_name])
		if (norm) {
			score_vec[i] = raw_score / n_muts #normalized
		} else {
			score_vec[i] = raw_score
		}
	}
	return(score_vec)
}

# find the exposures of patients that have mutations in the panel DF
patient_panel_exposures <- function(panel_df, signature, sig_df) {
	patient_vec = unique(panel_df$Patient)
	sig_name = paste0("Signature.", signature)

	exposure_vec = numeric(length(patient_vec))
	names(exposure_vec) = patient_vec

	for (i in 1:length(patient_vec)) {
		p = patient_vec[i]
		exposure_vec[i] = sig_df[sig_df$Patient==p, sig_name]
	}

	return(exposure_vec)
}




#########################################################################
# save fold-change matrices in chromosome chunks for the mutwindows     #
#########################################################################

save_10k_fc_mtxs <- function() {
	chr_vec = c(1:22, "X", "Y")
	for (chr in chr_vec) {
		print(paste0("loading mutwindow_chr_10k for chromosome: ", chr))
		mutwindow_mtx = load_nz_mutwindow_chr_10k(chr)
		
		print(paste0("running par_fc_mtx for chromosome: ", chr))
		ret = par_fc_mtx(mutwindow_mtx)

		filename = paste0(GLOBAL_CHR_MTX_DIR, "fc_mtx_mutwindow_10k_chr_", chr, ".tsv")
		print(paste0("writing table for chromosome: ", chr))
		print(paste0("filename: ", filename))
		write.table(ret, file=filename, sep="\t", quote=FALSE)
	}
}

save_100k_fc_mtxs <- function() {
	chr_vec = c(1:22, "X", "Y")
	for (chr in chr_vec) {
		print(paste0("loading mutwindow_chr_100k for chromosome: ", chr))
		mutwindow_mtx = load_nz_mutwindow_chr_100k(chr)
		
		print(paste0("running par_fc_mtx for chromosome: ", chr))
		ret = par_fc_mtx(mutwindow_mtx)

		filename = paste0(GLOBAL_CHR_MTX_DIR, "fc_mtx_mutwindow_100k_chr_", chr, ".tsv")
		print(paste0("writing table for chromosome: ", chr))
		print(paste0("filename: ", filename))
		write.table(ret, file=filename, sep="\t", quote=FALSE)
	}
}


corr_plot <- function(signature, panel_df, outfile, normalize_sigscore=TRUE, normalize_exposure=TRUE, add_labels=FALSE) {
	sig_df = load_nz_sig_estimates(norm=normalize_exposure)

	x = patient_panel_sig_scores(panel_df, signature, norm=normalize_sigscore)
	y = patient_panel_exposures(panel_df, 2, sig_df)

	inds = which(!is.na(x))
	x = x[inds]
	y = y[inds]


	my_xlab = paste0("Panel Sig", signature, " Score")
	if (normalize_sigscore) {
		my_xlab = paste0("Normalized ", my_xlab)
	}

	my_ylab = paste0("Signature ", signature, " Exposure")
	if (normalize_exposure) {
		my_ylab = paste0("Normalized ", my_ylab)
	}

	png(outfile)
	plot(x, y, main=paste0("Panel Signature Score vs. Global Exposure for Sig ", signature), xlab=my_xlab, ylab = my_ylab)
	if (add_labels) {
		text(x, y+(max(y) / 100), labels=names(x))
	}
	dev.off()
	print(paste0("Saved result at ", outfile))
	print(cor(x,y, method="spearman"))
	results = list()
	results$x = x
	results$y = y
	return(results)
}





#################################################
# new objective function for finding windows    #
#################################################

# use the SCALAR PROJECTION 
# given a 96-vector of mutation counts t
# and a 96-vector mutation signature s
# the SCALAR PROJECTION alpha of t onto s can be seen as a proxy for exposures, it gives us some notion
# of, if we were to explain mutations in this region SOLELY with signature s, how well would we do?
# NOTE: maybe we should use laplace smoothing or something.
# to compute alpha, use the following:

# alpha = t * cos(theta)
# (where theta is the angle between t and s)

# or equivalently 
# alpha = t \dot s
# since mutation signatures are unit vectors

# then we can define an objective function which rewards regions with high exposures from samples where the signature
# is (globally) active, and penalizes high exposures from samples where the signature is (globally) inactive

# NOTATION: alpha_P := alpha when you construct t out of the mutation counts in a panel P

#OBJ : SUM_(over positive samples) alpha_P    -    SUM_(over negative samples) alpha_P

# we would ideally like to find P of a certain size that maximizes the objective, but the search space is not very nice
# (SUM_(over each region r in P) alpha_r != alpha_P)
# So we will either have to take a greedy approach or do some kind of randomized search mumbo jumbo



# possible modifications of OBJ:
# take sqrt of alpha to dissuade overfitting to samples with many mutations
# (or possibly a piecewise function where f(x) = x if x < 1, sqrt(x) o/w

# use a scaling factor to mitigate class imbalance (i.e. multiply the second term by |positive samples| / |negative samples|)









###################################
# run from scratch                #
###################################


#cross_chrom_top_windows_10k <- function(signature, n=50, debug=TRUE) {
#	chrom_vec = c(1:22, "X", "Y")
#
#	candidates = numeric(0)
#	for (chr in chrom_vec) {
#		if (debug) { print(paste0("retrieving windows from chromosome ", chr)) }
#		fc_mtx = load_fc_mutwindow_chr_10k(chr)
#		topn = active_sig_windows(fc_mtx, signature, n)
#		candidates = c(candidates, topn)
#	}
#	return(top_n(candidates, n))
#}

cross_chrom_top_windows_help <- function(signature, n, mode=c("10k", "100k"), debug=TRUE) {
	chrom_vec = c(1:22, "X", "Y")

	candidates = numeric(0)
	for (chr in chrom_vec) {
		if (debug) { print(paste0("retrieving windows from chromosome ", chr)) }
		if (mode=="10k") {
			if (debug) { print("using 10k fn") }
			fc_mtx = load_fc_mutwindow_chr_10k(chr)
		} else {
			if (debug) { print("using 100k fn") }
			fc_mtx = load_fc_mutwindow_chr_100k(chr)
		}
		topn = active_sig_windows(fc_mtx, signature, n)
		candidates = c(candidates, topn)
	}
	return(top_n(candidates, n))
}

get_fc_ls <- function(mode="100k", debug=TRUE) {
	chrom_vec = c(1:22, "X", "Y")
	fc_ls = list()
	for (i in 1:length(chrom_vec)) {
		if (mode=="10k") {
			if (debug) { print("using 10k fn") }
			fc_ls[[i]] = load_fc_mutwindow_chr_10k(chrom_vec[i])
		} else {
			if (debug) { print("using 100k fn") }
			fc_ls[[i]] = load_fc_mutwindow_chr_100k(chrom_vec[i])
		}
	}

	return(fc_ls)
}

cross_chrom_top_windows_all_sigs <- function(n, mode=c("10k", "100k"), fc_ls=NULL, debug=TRUE) {
	chrom_vec = c(1:22, "X", "Y")

	if (is.null(fc_ls)) {
		fc_ls = list()
		for (i in 1:length(chrom_vec)) {
			if (mode=="10k") {
				if (debug) { print("using 10k fn") }
				fc_ls[[i]] = load_fc_mutwindow_chr_10k(chrom_vec[i])
			} else {
				if (debug) { print("using 100k fn") }
				fc_ls[[i]] = load_fc_mutwindow_chr_100k(chrom_vec[i])
			}
		}
	}

	res_ls = list()
	for (sig_num in 1:30) {
		candidates = numeric(0)
		for (i in 1:length(chrom_vec)) {
			#if (debug) { print(paste0("running active_sig_windows() for sig: ", sig_num, " ; chrom: ", chrom_vec[i])) }
			topn = active_sig_windows(fc_ls[[i]], sig_num, n)
			candidates = c(candidates, topn)
		}
		res_ls[[sig_num]] = top_n(candidates, n)
	}

	return(res_ls)
}


cross_chrom_top_windows_100k <- function(signature, n=50, debug=TRUE) {
	return(cross_chrom_top_windows_help(signature, n, mode="100k", debug=debug))
}

cross_chrom_top_windows_10k <- function(signature, n=50, debug=TRUE) {
	return(cross_chrom_top_windows_help(signature, n, mode="10k", debug=debug))
}


get_panel_df_10k <- function(signature, n=50, mut_df=NULL, debug=TRUE) {
	if (is.null(mut_df)) {
		mut_df = load_nz_mut_df_with_sigprob()
	}
	panel_windows = cross_chrom_top_windows_10k(signature, n, debug)
	panel_df = select_window_muts(names(panel_windows), mut_df)
	return(panel_df)
}

get_panel_df_100k <- function(signature, n=50, mut_df=NULL, debug=TRUE) {
	if (debug) { print("Loading mut_df") }
	if (is.null(mut_df)) {
		mut_df = load_nz_mut_df_with_sigprob()
	}
	if (debug) { print("running cross_chrom_top_windows_100k()") }
	panel_windows = cross_chrom_top_windows_100k(signature, n, debug)
	if (debug) { print("running select_window_muts()") }
	panel_df = select_window_muts(names(panel_windows), mut_df)
	return(panel_df)
}

get_all_sig_panels_100k <- function(n=50, mut_df=NULL, fc_ls=NULL, debug=TRUE) {
	if (debug) { print("Loading mut_df") }
	if (is.null(mut_df)) {
		mut_df = load_nz_mut_df_with_sigprob()
	}
	p_window_ls = cross_chrom_top_windows_all_sigs(n, mode="100k", fc_ls=fc_ls, debug=debug)

	p_df_ls = list()
	for (i in 1:length(p_window_ls)) {
		if (debug) { print(paste0("Running select_window_muts for Sig ", i)) }
		p_df_ls[[i]] = select_window_muts(names(p_window_ls[[i]]), mut_df)
	}

	return(p_df_ls)
}


save_corr_plots_100k <- function(n=27, mut_df=NULL, fc_ls=NULL, debug=TRUE) {
	if (is.null(mut_df)) {
		mut_df = load_nz_mut_df_with_sigprob()
	}
	p_window_ls = cross_chrom_top_windows_all_sigs(n, mode="100k", fc_ls=fc_ls, debug=debug)

	for (i in 1:length(p_window_ls)) {
		if (debug) { print(paste0("Running select_window_muts for Sig ", i)) }
		panel_df = select_window_muts(names(p_window_ls[[i]]), mut_df)
		
		my_outfile = paste0("out/sig", i, "_n", n, "_100kpanel_scatter_nosignorm.png")
		
		corr_plot(i, panel_df, my_outfile, FALSE, TRUE, FALSE)
	}
}


# get mutations in top 10 fold change windows for Sig 2
not_sure_of_name_yet <- function() {
	mut_df = load_nz_mut_df_with_sigprob()
	fc_chr1 = load_fc_mutwindow_chr_100k(1)
	top10_sig2 = report_top_windows(fc_chr1, 20)$Signature.2
	df = select_window_muts(names(top10_sig2), mut_df)
	return(df)
}
