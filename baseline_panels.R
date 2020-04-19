source("preprocess_windows.R")
source("win_sigprob_analysis.R")

# if laptop: "~/projects/hotspot_signature_panel/data/"
# if workstation: "/fs/cbcb-lab/mdml/users/franzese/projects/signature-panel/signature-panel/data/"

GLOBAL_DATA_DIR = "~/projects/hotspot_signature_panel/data/"
GLOBAL_CHR_MTX_DIR = paste0(GLOBAL_DATA_DIR, "individual_chromosome_matrices/")

# if laptop: 3
# if workstation: 35

GLOBAL_NCORES = 3

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
get_random_windows <- function(n, window_vec) {
	ret = character(0)
	perm = sample(window_vec, length(window_vec)) #random permutation of the windows

	for (i in 1:length(window_vec)) {
		new = perm[i]
		if (!is_overlapping(new, ret, debug=FALSE)) {
			ret = c(ret, new)
			if (length(ret) >= n) {
				break
			}
		}
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
main <- function(WINDOWS_IN_PANEL, NUM_PANELS, OUTPREFIX, START_INDEX=1, mut_df=NULL, verbose=1) {
	registerDoParallel(cores=GLOBAL_NCORES)
	SBS_OUTPREFIX = paste0(OUTPREFIX, "SBS_MATRICES/random_panel_ALLSAMPLE_sbs_m100k_n", WINDOWS_IN_PANEL, "_iter_")

	if (is.null(mut_df)) {
		if (verbose) { print("loading mut_df") }
		mut_df = load_nz_mut_df_with_sigprob()
	}

	if (verbose) { print("loading all window names") }
	possible_windows = get_all_mutwindows(mode="100k")

	foreach (i = START_INDEX:(NUM_PANELS+START_INDEX)) %dopar% {
		if (verbose>=1) { print(paste0("Starting random panel ", i, "/", NUM_PANELS)) }

		if (verbose>=2) { print(paste0("randomly selecting ", WINDOWS_IN_PANEL, " windows")) }
		selected_windows = get_random_windows(WINDOWS_IN_PANEL, possible_windows)

		if (verbose>=2) { print("getting mutations in panel") }
		panel_df = select_window_muts(selected_windows, mut_df)
		
		sbs_outfile = paste0(SBS_OUTPREFIX, i, ".tsv")
		if (verbose>=2) { print(paste0("saving SBS matrix to ", sbs_outfile)) }
		panel_sbs = sbs_df_from_mut_df(panel_df)	
		save_panel_sbs_tsv(panel_sbs, sbs_outfile)
	}
}

GLOBAL_BASELINE_PANEL_DIR = paste0(GLOBAL_DATA_DIR, "BASELINE_PANELS/")

#main(27, 1000, GLOBAL_BASELINE_PANEL_DIR, verbose=1)
