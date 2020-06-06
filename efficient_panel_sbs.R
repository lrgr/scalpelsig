# blah

library("foreach")


load_10k_sbs_arr_ls <- function(debug=FALSE) {
	chrom_vec = c(1:22, "X", "Y")
	
	file_prefix = paste0(GLOBAL_CHR_MTX_DIR, "samp_window_sbs_arrays/", "samp_window_sbs_array_winsize10000_chr")
	
	arr_list = list()
	for (i in 1:length(chrom_vec)) {
		if (debug) { print(paste0(Sys.time(), "    load_10k_sbs_arr_ls() chrom ", i, "/", length(chrom_vec))) }
		sbs_arr_file = paste0(file_prefix, chrom_vec[i], ".rds")
		arr_list[[ chrom_vec[i] ]] = readRDS(sbs_arr_file)
	}
	
	return(arr_list)
}





# given a set of windows and precomputed sbs arrays, returns a matrix with sbs counts for each patient in the panel
# the returned matrix has dimensions [<num_samples>, 96]
# panel_window_names : a character vector containing the name of each window in the panel
# sbs_array_or_arr_ls : either a combined sbs array with shape [<num_samples>, <num_windows>, 96] OR
#                       a list of sbs arrays (with same shape) indexed by chromosome
get_panel_sbs_mtx <- function(panel_window_names, sbs_array_or_arr_ls) {
	if (class(sbs_array_or_arr_ls) == "array") {
		return(get_panel_sbs_mtx_from_single_arr(panel_window_names, sbs_array_or_arr_ls))
	} else if (class(sbs_array_or_arr_ls) == "list") {
		return(get_panel_sbs_mtx_from_arr_ls(panel_window_names, sbs_array_or_arr_ls))
	} else {
		stop("get_panel_sbs_mtx() was given something that was neither an sbs_array nor an sbs_array_ls.")
	}
}

get_panel_sbs_df <- function(panel_window_names, sbs_array_or_arr_ls) {
	return(as.data.frame(get_panel_sbs_mtx(panel_window_names, sbs_array_or_arr_ls)))
}


# helper function for get_panel_sbs_mtx()
# handles the case where a single combined array is provided
#get_panel_sbs_mtx_from_single_arr <- function(panel_window_names, sbs_array) {
#	window_sbs_ls <- foreach(i = 1:length(panel_window_names)) %do% {
#		curr_window = panel_window_names[i]
#		win_sbs_mtx = sbs_array[ , curr_window, ] # [<num_samples>, 96] sbs matrix for the current window
#		win_sbs_mtx
#	}
#
#	mtx_add <- function(x) Reduce("+", x)
#	panel_sbs_mtx = mtx_add(window_sbs_ls)
#
#	return(panel_sbs_mtx)
#}

get_panel_sbs_mtx_from_single_arr <- function(panel_window_names, sbs_array) {
	n_samp = dim(sbs_array)[1]
	samp_names = dimnames(sbs_array)[[1]]

	n_cat = dim(sbs_array)[3]
	cat_names = dimnames(sbs_array)[[3]]

	ret_sbs_mtx = matrix(0, nrow = n_samp, ncol = n_cat)
	rownames(ret_sbs_mtx) = samp_names
	colnames(ret_sbs_mtx) = cat_names

	for (i in 1:length(panel_window_names)) {
		curr_window = panel_window_names[i]
		print(i)
		win_sbs_mtx = sbs_array[ , curr_window, ] # [<num_samples>, 96] sbs matrix for the current window
		print("done")
		ret_sbs_mtx = ret_sbs_mtx + win_sbs_mtx
	}
	return(ret_sbs_mtx)
}


# helper function for get_panel_sbs_mtx()
# handles the case where a list of sbs_arrays is provided
get_panel_sbs_mtx_from_arr_ls <- function(panel_window_names, sbs_array_ls, debug=FALSE) {
	n_samp = dim(sbs_array_ls[[1]])[1]
	samp_names = dimnames(sbs_array_ls[[1]])[[1]]
	
	n_cat = dim(sbs_array_ls[[1]])[3]
	cat_names = dimnames(sbs_array_ls[[1]])[[3]]
	#if (debug) { saveRDS(panel_window_names, file="DEBUG_LOG_panel_win_names.rds") }
	#window_sbs_ls <- foreach(i = 1:length(panel_window_names)) %do% {
	ret_sbs_mtx = matrix(0, nrow = n_samp, ncol = n_cat)
	rownames(ret_sbs_mtx) = samp_names
	colnames(ret_sbs_mtx) = cat_names

	for(i in 1:length(panel_window_names)) {
		if (debug) { print(paste0( "panel window ", i , "/", length(panel_window_names))) }
		curr_window = panel_window_names[i]
		if (debug) { print(paste0("curr_window: ", curr_window)) }
		chrom = chrom_from_win_str(curr_window)
		if (debug) { print(paste0("chrom: ", chrom)) }
		chrom_sbs_arr = sbs_array_ls[[chrom]]

		win_sbs_mtx = chrom_sbs_arr[ , curr_window, ] # [<num_samples>, 96] sbs matrix for the current window
		#win_sbs_mtx
		ret_sbs_mtx = ret_sbs_mtx + win_sbs_mtx
	}
	if (debug) { print("finished loop") }

	#mtx_add <- function(x) Reduce("+", x) # a function for element-wise addition of multiple matrices
	#panel_sbs_mtx = mtx_add(window_sbs_ls) # element-wise addition of the matrices for each window

	return(ret_sbs_mtx)
}


# helper function for get_10k_panel_sbs_mtx()
# given a window str, returns the chromosome that this window belongs to
chrom_from_win_str <- function(win_str, ret_as_num=TRUE) {
	s = gsub("chr", "", win_str) # remove 'chr' from the string
	s = gsub("_.*", "", s) # remove the underscore after the chromosome name and everything trailing it
	if (ret_as_num==TRUE) {
		if (s == "X") { n = 23 }
		else if (s == "Y") { n = 24 }
		else { n = as.numeric(s) }
		return(n)
	} else {
		return(s)
	}
}
