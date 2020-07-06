source("GLOBAL_CONFIG.R")



##### parsing panel file names ######

parse_sig_num <- function(file_name) {
	sig_chunk = regmatches(file_name, regexpr("_sig[0-9]+_", file_name))
	sig_chunk = sub("_sig", "", sig_chunk)
	sig_chunk = sub("_", "", sig_chunk)
	sig_num = as.numeric(sig_chunk)
	return(sig_num)
}

parse_obj_fn <- function(file_name) {
	obj_chunk = regmatches(file_name, regexpr("_obj[0-9]+_", file_name))
	obj_chunk = sub("_obj", "", obj_chunk)
	obj_chunk = sub("_", "", obj_chunk)
	obj_num = as.numeric(obj_chunk)
	return(obj_num)
}


parse_iter <- function(file_name) {
	it_chunk = regmatches(file_name, regexpr("_it[0-9]+_", file_name))
	it_chunk = sub("_it", "", it_chunk)
	it_chunk = sub("_", "", it_chunk)
	it_num = as.numeric(it_chunk)

	if (is.na(it_num)) {
		print(file_name)
		stop("The file printed above introduced NA by coercion.")
	}
	return(it_num)
}


# given a list of panel files, detect the set of signatures being evaluated
detect_sig_nums <- function(file_pool) {
	return(sort(unique(parse_sig_num(file_pool))))
}

detect_obj_nums <- function(file_pool) {
	return(sort(unique(parse_obj_fn(file_pool))))
}


###### parsing window names ######

parse_chr_from_win_name <- function(win_name) {
	chr_chunk = regmatches(win_name, regexpr("chr[0-9,X,Y]+", win_name))
	chr = sub("chr", "", chr_chunk)
	if (chr == "") { stop(paste0("In parse_chr_from_win_name(), win_name: ", win_name, " yielded empty string as chromosome.")) }
	return(chr)
}

parse_start_from_win_name <- function(win_name) {
	start_chunk = sub("chr[0-9,X,Y]+_", "", win_name) # remove "chrNN_"
	#start_chunk = regmatches(start_chunk, regexpr("[0-9]+_", start_chunk))
	start_chunk = sub("_.*", "", start_chunk) # remove everything after trailing underscore
	start_bp = as.numeric(start_chunk)
	return(start_bp)
}

parse_end_from_win_name <- function(win_name) {
	end_chunk = sub(".*_.*_", "", win_name) # remove everything before the two underscores
	end_bp = as.numeric(end_chunk)
	if (is.na(end_bp)) {
		print(win_name) 
		stop("The above win_name introduced NA by coercion")
	}
	return(end_bp)
}





############## functions for window df construction and other forms of data organization ####################


# given a pool of panel window files
# creates a data frame with one entry for each window
# the data frame records, for each window, what signature, what obj function, what random trial did the panel come from
construct_window_df <- function(file_pool) {
	i = 1
	panel_ls = list()
	for (f in file_pool) {
		panel_windows = scan(paste0(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, f), what=character(), quiet=TRUE)
		panel_ls[[i]] = panel_windows
		i = i + 1
	}
	all_windows = unlist(panel_ls)
	
	# initialize data frame
	n = length(all_windows)
	
	Window.Name = character(n)
	Chromosome = character(n) # what chromosome is the window on
	Start = numeric(n) # start coordinate in bp
	End = numeric(n) # end coordinate in bp
	Panel.Rank = numeric(n) # windows are sorted by obj fn score, lower rank indicates better scoring window
	Signature = numeric(n) # which signature was this panel designed for
	Obj.Fn = numeric(n) # which objective function did the panel use
	Trial.Num = numeric(n) # same as "iter" / "iteration"
	Panel.File = character(n) # file name containing the windows in this panel

	win_index = 1 # iterates from 1 to n (total number of windows in all the panels)
	for (f in file_pool) {
		curr_panel = scan(paste0(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, f), what=character(), quiet=TRUE)
		curr_sig = parse_sig_num(f)
		curr_obj = parse_obj_fn(f)
		curr_trial = parse_iter(f)

		curr_panel_rank = 1
		for (w in curr_panel) {
			Window.Name[win_index] = w
			Chromosome[win_index] = parse_chr_from_win_name(w)
			Start[win_index] = parse_start_from_win_name(w)
			End[win_index] = parse_end_from_win_name(w)
			Panel.Rank[win_index] = curr_panel_rank
			Signature[win_index] = curr_sig
			Obj.Fn[win_index] = curr_obj
			Trial.Num[win_index] = curr_trial
			Panel.File[win_index] = f

			curr_panel_rank = curr_panel_rank + 1
			win_index = win_index + 1
		}
	}

	results_df = data.frame(Signature, Obj.Fn, Window.Name, Panel.Rank, Trial.Num, Chromosome, Start, End, Panel.File)
	results_df = results_df[order(Signature, Obj.Fn, Panel.File, Panel.Rank), ]
	return(results_df)
}



# given a df of windows (from construct_window_df() ), return a list of dfs, one for each individual panel
# in the win_df
panel_dfs_from_win_df <- function(win_df) {
	ret = list()
	panel_files = unique(win_df$Panel.File)

	i = 1
	for (f in panel_files) {
		ret[[i]] = win_df[win_df$Panel.File == f, ]
		i = i + 1
	}
	return(ret)
}


# convenience function for subsetting a window df (e.g. from construct_window_df() or load_canonical_win_df())
# input: df - the window df to be subsetted
# input: sig - the signature or set of signatures to be included in the subset - if no input is given, allows all sigs
# input: obj - the objective function or set of obj fns to be included - if no input is given, allows all obj fns
# 
# output: a window df subsetted according to the given arguments
subset_win_df <- function(df, sig=NULL, obj=NULL) {
	if (is.null(sig)) {
		sig_df = df
	} else {
		sig_df = df[df$Signature %in% sig, ]
	}

	if (is.null(obj)) {
		obj_df = sig_df
	} else {
		obj_df = sig_df[sig_df$Obj.Fn %in% obj, ]
	}

	return(obj_df)
}



################ functions for analysis of the data ######################


# given a df of windows (from construct_window_df(), or a subset of this output),
# breaks the df into individual panels,
# then finds windows that appear in at least <commonality_thresh> % of the panels
# returns these windows as a character vector
# by default, commonality_thresh is set to 1, so it finds windows shared by all panels in the input df
common_panel_windows <- function(df, commonality_thresh=1) {
	p_dfs = panel_dfs_from_win_df(df) # break df into its constituent panels
	total_panels = length(p_dfs)

	all_windows = unique(as.character(df$Window.Name)) # get all unique window names within the df
	n_windows = length(all_windows)
	
	# initialize vector to count the number of panels that each window appears in
	panel_count = numeric(n_windows)
	names(panel_count) = all_windows

	# for each window, count the number of panels in which it appears
	for (w in all_windows) {
		for (p in p_dfs) {
			if (w %in% as.character(p$Window.Name)) {
				panel_count[w] = panel_count[w] + 1
			}
		}
	}
	
	t = commonality_thresh * total_panels # number of panels a window has to be in to be included

	# return the windows that appear in at least t panels.
	return(names(panel_count[panel_count >= t]) )
}


# given a window name and a df of windows, find the mean rank of that window among panels in the df
window_mean_rank <- function(win_name, df) {
	win_entries = df[df$Window.Name == win_name, ]
	ranks = win_entries$Panel.Rank
	return(mean(ranks))
}






############# DEPRECATED functions ################
# these functions may be used in other parts of the code (e.g. the experiments section),
# but moving forward their use should be halted


# DEPRECATED: use common_panel_windows() instead - default setting does the same thing as this
# given a list of panel dfs (e.g. from panel_dfs_from_win_df() ),
# find the set of windows that is shared by all the panels in the list
panel_df_intersection <- function(panel_df_ls) {
	i = 1
	for (df in panel_df_ls) {
		df_windows = df$Window.Name
		
		if (i == 1) {
			# in first iteration, just include all windows of the first panel
			shared = df_windows
		} else {
			# in subsequent iterations, take the intersection with previous windows and the current panel
			shared = intersect(shared, df_windows)
		}

		i = i + 1
	}
	return(shared)
}



# DEPRECATED - use subset_win_df() instead
get_panel_windows_by_sig <- function(sig_num, file_pool) {
	# sanity check to make sure the requested signature is in the file pool
	sigs_in_pool = detect_sig_nums(file_pool)
	if ( !(sig_num %in% sigs_in_pool) ) {
		print("Requested signature: ", sig_num, " was not found in the provided pool of files.")
		print("Detected signatures in given file pool:")
		print(sigs_in_pool)
		stop("Please provide a signature from within the file pool.")
	}

	# get panels from the given pool that examine the given signature
	sig_files = regmatches(file_pool, regexpr(pattern=paste0(".*_sig", sig_num, "_.*"), file_pool))
	
	# go through panel files and extract window names, add them to a list
	sig_panels = list()
	i = 1
	for (f in sig_files) {
		panel_windows = scan(paste0(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, f), what=character(), quiet=TRUE)
		sig_panels[[i]] = panel_windows
		i = i + 1
	}
	return(sig_panels)
}





################### wrapper functions for experiments & cursory investigations ######################


### sets up the 'canonical' panel window df ###
load_canonical_win_df <- function() {	
	tags = c("B1_BC_SIGS_BAL", "B2_BC_SIGS_BAL", "B3_BC_SIGS_BAL", "B1_BC_SIGS_IMBAL", "B2_BC_SIGS_IMBAL", "B3_BC_SIGS_IMBAL")

	file_pool = character(0)
	for (t in tags) {
		new_files = list.files(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, pattern=paste0(".*", t, ".*"))
		file_pool = c(file_pool, new_files)
	}

	win_df = construct_window_df(file_pool)
	
	# cut the panels to 2.5 MB (rather than 2.7 MB)
	win_df = win_df[win_df$Panel.Rank <= 250, ]

	return(win_df)
}





shared_windows_sig1_resp_obj <- function(win_df=NULL) {
	if (is.null(win_df)) {
		win_df = load_canonical_win_df()
	}
	
	s1_df = win_df[win_df$Signature == 1, ]

	shared_ls = list()
	for (curr_obj in c(1,2,3)) {
		o_df = s1_df[s1_df$Obj.Fn == curr_obj, ]
		p_dfs = panel_dfs_from_win_df(o_df)
		shared_ls[[curr_obj]] = panel_df_intersection(p_dfs)
	}

	return(shared_ls)
}


count_shared_windows_resp_sig_and_obj <- function(win_df=NULL) {
	if (is.null(win_df)) {
		win_df = load_canonical_win_df()
	}

	sig_vec = unique(win_df$Signature)
	obj_vec = c(1,2,3)
	n = length(sig_vec) * length(obj_vec)

	Signature = numeric(n)
	Obj.Fn = numeric(n)
	t100.Shared = numeric(n)
	t90.Shared = numeric(n)
	t80.Shared = numeric(n)
	t50.Shared = numeric(n)

	i = 1
	for (s in sig_vec) {
		for (o in obj_vec) {
			sub_df = subset_win_df(win_df, sig=s, obj=o)
			t100.Shared[i] = length(common_panel_windows(sub_df, 1))
			t90.Shared[i] = length(common_panel_windows(sub_df, .9))
			t80.Shared[i] = length(common_panel_windows(sub_df, .8))
			t50.Shared[i] = length(common_panel_windows(sub_df, .5))

			Signature[i] = s
			Obj.Fn[i] = o
			i = i + 1
		}
	}

	results_df = data.frame(Signature, Obj.Fn, t50.Shared, t80.Shared, t90.Shared, t100.Shared)
	results_df = results_df[order(Signature, Obj.Fn), ]

	return(results_df)
}

# uses obj 2 panels by default, can use different objective function panels by changing the my_obj argument
count_shared_windows_sig_pairs <- function(win_df=NULL, my_obj=2) {
	if (is.null(win_df)) {
		win_df = load_canonical_win_df()
	}

	t_vec = c(.5, .75, .9, 1)

	sig_vec = unique(win_df$Signature)
	sig_pairs = combn(sig_vec, 2, simplify=FALSE)

	n = length(sig_pairs)

	Sig.A = numeric(n)
	Sig.B = numeric(n)
	t50.Shared = numeric(n)
	t75.Shared = numeric(n)
	t90.Shared = numeric(n)
	t100.Shared = numeric(n)

	results_df = data.frame(Sig.A, Sig.B, t50.Shared, t75.Shared, t90.Shared, t100.Shared)	

	i = 1
	for (p in sig_pairs) {
		sub_df1 = subset_win_df(win_df, sig=p[1], obj=my_obj)
		sub_df2 = subset_win_df(win_df, sig=p[2], obj=my_obj)

		t_ind = 1
		res_vec = numeric(length(t_vec))
		names(res_vec) = t_vec
		for (t in t_vec) {
			w1 = common_panel_windows(sub_df1, t)
			w2 = common_panel_windows(sub_df2, t)

			t_shared = intersect(w1, w2)
			res_vec[t_ind] = length(t_shared)
			
			t_ind = t_ind + 1
		}
		
		results_df[i, "Sig.A"] = p[1]
		results_df[i, "Sig.B"] = p[2]
		results_df[i, 3:6] = res_vec
		i = i + 1
	}

	return(results_df)
}


main <- function() {


	tags = c("B1_BC_SIGS_BAL", "B2_BC_SIGS_BAL", "B3_BC_SIGS_BAL", "B1_BC_SIGS_IMBAL", "B2_BC_SIGS_IMBAL", "B3_BC_SIGS_IMBAL")

	file_pool = character(0)
	for (t in tags) {
		new_files = list.files(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, pattern=paste0(".*", t, ".*"))
		file_pool = c(file_pool, new_files)
	}

	win_df = construct_window_df(file_pool)

	return(win_df)
	
}

