# blah

library("foreach")

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
get_panel_sbs_mtx_from_single_arr <- function(panel_window_names, sbs_array) {
	window_sbs_ls <- foreach(i = 1:length(panel_window_names)) %do% {
		curr_window = panel_window_names[i]
		win_sbs_mtx = sbs_array[ , curr_window, ] # [<num_samples>, 96] sbs matrix for the current window
		win_sbs_mtx
	}

	mtx_add <- function(x) Reduce("+", x)
	panel_sbs_mtx = mtx_add(window_sbs_ls)

	return(panel_sbs_mtx)
}


# helper function for get_panel_sbs_mtx()
# handles the case where a list of sbs_arrays is provided
get_panel_sbs_mtx_from_arr_ls <- function(panel_window_names, sbs_array_ls) {
	window_sbs_ls <- foreach(i = 1:length(panel_window_names)) %do% {
		curr_window = panel_window_names[i]
		chrom = chrom_from_win_str(curr_window)
		chrom_sbs_arr = sbs_array_ls[[chrom]]

		win_sbs_mtx = chrom_sbs_arr[ , curr_window, ] # [<num_samples>, 96] sbs matrix for the current window
		win_sbs_mtx
	}

	mtx_add <- function(x) Reduce("+", x) # a function for element-wise addition of multiple matrices
	panel_sbs_mtx = mtx_add(window_sbs_ls) # element-wise addition of the matrices for each window

	return(panel_sbs_mtx)
}


# helper function for get_10k_panel_sbs_mtx()
# given a window str, returns the chromosome that this window belongs to
chrom_from_win_str <- function(win_str) {
	s = gsub("chr", "", win_str) # remove 'chr' from the string
	s = gsub("_.*", "", s) # remove the underscore after the chromosome name and everything trailing it
	return(s)
}
