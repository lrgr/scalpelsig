source("preprocess_windows.R")
source("win_sigprob_analysis.R")

# if laptop: "~/projects/hotspot_signature_panel/data/"
# if workstation: "/fs/cbcb-lab/mdml/users/franzese/projects/signature-panel/signature-panel/data/"

GLOBAL_DATA_DIR = "~/projects/hotspot_signature_panel/data/"
GLOBAL_CHR_MTX_DIR = paste0(GLOBAL_DATA_DIR, "individual_chromosome_matrices/")

# if laptop: 3
# if workstation: 35

GLOBAL_NCORES = 3

source("GLOBAL_CONFIG.R")




################################################################
##### BASELINE 1: SAMPLE WINDOWS STARTING AT EACH MUTATION #####
################################################################

# this was the original randomized baseline, but it actually turned out to be a reasonable randomized algorithm.
# how it works: 
# we sample n_windows mutations from all the mutations in the data set
# we take windows (e.g. 100kb long) starting at the genome coordinates of the sampled mutations (reroll if overlap)
# the panel consists of the mutations that fall within these windows
# 
# when I first did this with ~1000 panels, some of them were actually quite good. 
# at first this is kind of weird, but when you think about it we are more likely to sample regions with high
# mutation density, since we sample uniformly from mutation sites.
# A cool idea would be to turn this into a more refined randomized algorithm, but I don't know if I'll have time.
# But anyway, we decided this was too conservative of a baseline, so below this block of code I wrote a dumber
# random baseline framework.





# get all region strings
get_all_mutwindows <- function(mode) {
	chrom_vec = c(1:22, "X", "Y")
	ret = character(0)
	
	for (chr in chrom_vec) {
		if (mode=="100k") {
			filename = paste0(GLOBAL_CHR_MTX_DIR, "fc_mtx_mutwindow_100k_chr_", chr, ".tsv")
		} else {
			filename = paste0(GLOBAL_CHR_MTX_DIR, "fc_mtx_mutwindow_10k_chr_", chr, ".tsv")
		}
		regions = scan(file=filename, what=character(), nlines=1)
		ret = c(ret, regions)
	}

	return(ret)
}

# select n random non-overlapping windows
get_random_windows <- function(n, window_vec, check_overlap=TRUE) {
	ret = character(0)
	perm = sample(window_vec, length(window_vec)) #random permutation of the windows


	if (check_overlap) {
	for (i in 1:length(window_vec)) {
		new = perm[i]
		if (!is_overlapping(new, ret, debug=FALSE)) {
			ret = c(ret, new)
			if (length(ret) >= n) {
				break
			}
		}
	}

	} else {
		ret = perm[1:n]
	}
	
	return(ret)
}


# save panel sbs df as a tsv
save_panel_sbs_tsv <- function(sbs_df, outfile) {
	write.table(sbs_df, file=outfile, sep="\t", quote=FALSE)
}


############## main ################

# randomly sample WINDOWS_IN_PANEL mutwindows from the genome
# put all the mutations 
main <- function(WINDOWS_IN_PANEL, NUM_PANELS, OUTPREFIX, WIN_MODE="100k", START_INDEX=1, mut_df=NULL, verbose=1) {
	registerDoParallel(cores=GLOBAL_NCORES)

	if ((WIN_MODE != "100k") & (WIN_MODE != "10k")) { stop(paste0("WIN_MODE argument of main() was set to ", WIN_MODE, " but it can only be '100k' or '10k' ")) }
	SBS_OUTPREFIX = paste0(OUTPREFIX, "SBS_MATRICES/random_panel_ALLSAMPLE_sbs_m", WIN_MODE, "_n", WINDOWS_IN_PANEL, "_iter_")

	if (is.null(mut_df)) {
		if (verbose) { print("loading mut_df") }
		mut_df = load_nz_mut_df_with_sigprob()
	}

	if (verbose) { print("loading all window names") }
	possible_windows = get_all_mutwindows(mode=WIN_MODE)

	ret_ls <- foreach (i = START_INDEX:(NUM_PANELS+START_INDEX)) %dopar% {
		if (verbose>=1) { print(paste0("Starting random panel ", i, "/", NUM_PANELS)) }

		if (verbose>=2) { print(paste0("randomly selecting ", WINDOWS_IN_PANEL, " windows")) }
		selected_windows = get_random_windows(WINDOWS_IN_PANEL, possible_windows)

		if (verbose>=2) { print("getting mutations in panel") }
		panel_df = select_window_muts(selected_windows, mut_df) # this is a big runtime bottleneck
		
		sbs_outfile = paste0(SBS_OUTPREFIX, i, ".tsv")
		if (verbose>=2) { print(paste0("saving SBS matrix to ", sbs_outfile)) }
		panel_sbs = sbs_df_from_mut_df(panel_df)
		panel_sbs
	}

	for (i in 1:length(ret_ls) ) {
		panel_sbs_df = ret_ls[[i]]
		sbs_outfile = paste0(SBS_OUTPREFIX, i+START_INDEX, ".tsv")
		save_panel_sbs_tsv(panel_sbs_df, sbs_outfile)
	}	
}

#GLOBAL_BASELINE_PANEL_DIR = paste0(GLOBAL_DATA_DIR, "BASELINE_PANELS/")

#main(27, 1000, GLOBAL_BASELINE_PANEL_DIR, verbose=1)



############################################################
##### BASELINE 2: UNIFORM SAMPLING OF NONEMPTY WINDOWS #####
############################################################


load_combined_100k_arr <- function() {
	filename = paste0(GLOBAL_DATA_DIR, "COMBINED_sbs_array_winsize1e+05.rds")
	arr = readRDS(filename)
	return(arr)
}



# okay so here we will just chop the genome into blocks of 10kb or 100kb, figure out the ones that have at least
# one mutation in them, and then sample uniformly from those


# find windows with at least one mutation 

find_active_window_ls <- function(sbs_array) {
	n_windows = dim(sbs_array)[2]
	win_names = dimnames(sbs_array)[[2]]

	active_windows = list()
	n_active = 1
	foreach(i=1:n_windows) %do% {
		#win_sums[i] = sum(sbs_array[ , i, ]) 
		if ( sum(sbs_array[ , i, ]) > 0 ) {
			active_windows[[n_active]] = win_names[i]
			n_active = n_active+1
		}
	}
	return(active_windows)
}

find_active_window_vec_from_arr_ls <- function(sbs_arr_ls) {
	active_windows = character(0)
	for (chr_arr in sbs_arr_ls) {
		curr = as.character(find_active_window_ls(chr_arr))
		active_windows = c(active_windows, curr)
	}
	return(active_windows)
}


active_100k_windows <- function(combined_sbs_array=NULL, debug=TRUE) {
	if (is.null(combined_sbs_array)) {
		if (debug) { print(paste0(Sys.time(), "    active_100k_windows() loading combined_sbs_array")) }
		combined_sbs_array = load_combined_100k_arr()
	}

	active_windows = find_active_window_ls(combined_sbs_array)
	return(as.character(active_windows))
}

active_10k_windows <- function(sbs_arr_ls = NULL, debug=TRUE) {
	if (is.null(sbs_arr_ls)) {
		if (debug) { print(paste0(Sys.time(), "    loading 10kb sbs array.")) }
		sbs_arr_ls = load_10k_sbs_arr_ls()
		if (debug) { print(paste0(Sys.time(), "    done.")) }
	}
	
	active_windows = find_active_window_vec_from_arr_ls(sbs_arr_ls)
	return(active_windows)
}


#TODO write modified MAIN() to use these windows 

window_main <- function(WINDOWS_IN_PANEL, NUM_PANELS, WIN_MODE="100k", START_INDEX=1, OUTDIR=NULL, mut_df=NULL, sbs_array=NULL, verbose=1) {
	registerDoParallel(cores=GLOBAL_NCORES)

	if (is.null(OUTDIR)) {
		OUTDIR = GLOBAL_BASELINE_SBS_DIR
	}	

	if ((WIN_MODE != "100k") & (WIN_MODE != "10k")) { stop(paste0("WIN_MODE argument of main() was set to ", WIN_MODE, " but it can only be '100k' or '10k' ")) }
	SBS_OUTPREFIX = paste0(OUTDIR, "random_panel_ACTIVEWINDOWS_sbs_m", WIN_MODE, "_n", WINDOWS_IN_PANEL, "_iter_")

	if (is.null(mut_df)) {
		if (verbose) { print("loading mut_df") }
		mut_df = load_nz_mut_df_with_sigprob()
	}
	
	if (verbose) { print("loading all window names") }
	#possible_windows = get_all_mutwindows(mode=WIN_MODE)
	if (WIN_MODE == "100k") {
		possible_windows = active_100k_windows(sbs_array)
	}


	ret_ls <- foreach (i = START_INDEX:(NUM_PANELS+START_INDEX)) %dopar% {
		if (verbose>=1) { print(paste0("Starting random panel ", i, "/", NUM_PANELS)) }

		if (verbose>=2) { print(paste0("randomly selecting ", WINDOWS_IN_PANEL, " windows")) }
		selected_windows = get_random_windows(WINDOWS_IN_PANEL, possible_windows, check_overlap=FALSE)	

		if (verbose>=2) { print("getting mutations in panel") }
		panel_df = select_window_muts(selected_windows, mut_df) # this is a big runtime bottleneck
		
		#sbs_outfile = paste0(SBS_OUTPREFIX, i, ".tsv")
		#if (verbose>=2) { print(paste0("saving SBS matrix to ", sbs_outfile)) }
		panel_sbs = sbs_df_from_mut_df(panel_df)

		if (verbose >=3) { 
			print(paste0("Iteration ", i, ": "))
			print(selected_windows)
			print(dim(panel_df))
			print(dim(panel_sbs))
		}	

		panel_sbs
	}

	saveRDS(ret_ls, "ret_ls.rds")

	success_count = 0
	for (i in 1:length(ret_ls) ) {
		panel_sbs_df = ret_ls[[i]]
		
		if (!is.null(panel_sbs_df)) {
			sbs_outfile = paste0(SBS_OUTPREFIX, i+START_INDEX, ".tsv")
			save_panel_sbs_tsv(panel_sbs_df, sbs_outfile)
			success_count = success_count + 1
		}
	}
	print(paste0(success_count, "/", NUM_PANELS, " panels generated successfully."))	
}


#main(27, 1000, GLOBAL_BASELINE_PANEL_DIR, verbose=1)

#window_main(27, 50, GLOBAL_BASELINE_PANEL_DIR, verbose=1)


