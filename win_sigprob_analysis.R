# script for statistical analysis of regional signature probabilities


GLOBAL_DATA_DIR = "/fs/cbcb-lab/mdml/users/franzese/projects/signature-panel/signature-panel/data/"
GLOBAL_CHR_MTX_DIR = paste0(GLOBAL_DATA_DIR, "individual_chromosome_matrices/")
GLOBAL_NCORES = 35

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


# fread from library(data.table) works for reading but it's over 3 million rows and it crashes when you attempt access of entries
# don't use these, it seems they are too big for read.csv() to handle them without getting stuck.
load_nz_mutwindow_10k <- function() {
	df = read.csv("~/projects/hotspot_signature_panel/data/nz_mutwindow_10000_sigprob.tsv",
		      sep = "\t", header=TRUE)
	return(df)
}

load_nz_mutwindow_1k <- function() {
	df = read.csv("~/projects/hotspot_signature_panel/data/nz_mutwindow_1000_sigprob.tsv",
		      sep = "\t", header=TRUE)
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
	sort(x, decreasing=TRUE)[1:n]
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
