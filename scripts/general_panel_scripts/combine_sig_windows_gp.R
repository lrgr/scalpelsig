# the previous script in the workflow (find_windows_general_panel.R) computes windows for several panels for individual signatures,
# and saves the names of those windows to text files.
# This script will look at those panels, take the top windows from them, and integrate them into a single general panel.


setwd("../..")
source("projection_score.R")

library(optparse)

option_list = list(
        make_option(c("-t", "--tag"), type="character", default=NULL,
                        help="tag to refer downstream scripts to the correct training sets", metavar="character"),
        make_option(c("-n", "--numiters"), type="numeric", default=NULL,
                        help="number of random interations", metavar="integer"),
        make_option(c("-s", "--sizepanel"), type="numeric", default=NULL,
                        help="number of windows in the panel total", metavar="integer"),
	make_option(c("-v", "--verbose"), type="numeric", default=2,
			help="Level of verbosity, can be set 1-4", metavar="integer")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

file_tag = opt$tag
num_iters = opt$numiters
num_windows = opt$sizepanel
verbose = opt$verbose

if (is.null(file_tag)) {
        stop("No file_tag recieved (command line -t ). Please supply a file_tag.")
}
if (is.null(num_iters)) {
        stop("No num_iters recieved (command line -n ). Please supply num_iters.")
}

if (is.null(num_windows)) {
        stop("No num_windows recieved (command line -s ). Please supply num_windows.")
}

#print(paste0(Sys.time(), "    loading 10kb sbs array"))
#sbs_arr_ls_10k = load_10k_sbs_arr_ls()
#print(paste0(Sys.time(), "    done."))
#print(paste0(Sys.time(), "    saving sbs matrices for panels with file_tag: ", file_tag))

files = list.files(GLOBAL_GP_IND_WIN_DIR, pattern= paste0(".*", file_tag, ".*") )
print(paste0(Sys.time(), "    found ", length(files), " files containing the tag: ", file_tag))


# adds x windows from new_panel to windows_so_far, such that the new windows are not already in windows_so_far 
# input - windows_so_far: a character vector containing the names of windows already in the general panel
# input - new_panel: a character vector containing names of windows in a sorted panel for a new individual signature
# input - x: number of windows to add onto windows_so_far
add_x_unique <- function(windows_so_far, new_panel, x) {
	n = length(windows_so_far)
	goal_len = n + x
	
	new_wins = new_panel[1:x]
	
	p = union(windows_so_far, new_wins)

	for (i in (x+1):length(new_panel)) {
		if (length(p) >= goal_len) {
			break
		}
		
		if (! (new_panel[i] %in% p) ) {
			p = c(p, new_panel[i])
		}
	}
	if (length(p) != goal_len) {
		print(paste0("Warning: in add_x_unique(), the returned panel had ", length(p), " windows - wanted ", goal_len, "."))
	}
	return(p)
}


for (o in 2) {
	# separate the panel files based on objective function
	obj_file_set = regmatches(files, regexpr(paste0(".*_obj", o, "_.*"), files))

	for (i in 1:num_iters) {
        	# get the windows of each individual signature that correspond to this test/train set
        	it_files = regmatches(obj_file_set, regexpr(paste0(".*_it", i, "_.*"), obj_file_set))

		n_sigs = length(it_files)
	
	
		wins_per_segment = floor(num_windows / n_sigs)

		if (verbose >= 2) {
			print(paste0("Obj ", o, ", Iteration ", i, ": ", n_sigs, " files found."))
			print(paste0("==> ", wins_per_segment, " windows per signature in the general panel."))
			if (verbose >= 3) { print(paste0("(", num_windows - (n_sigs * wins_per_segment), " unoccupied)")) }
			if (verbose >= 3) { print(it_files) }
		}

        	gp_windows = character(0)

        	for (f in it_files) {
                	panel_windows = scan(paste0(GLOBAL_GP_IND_WIN_DIR, f), what=character(), quiet=TRUE)

			if (verbose >= 3) {
				print(paste0("Length of current panel: ", length(gp_windows)))
				print(paste0("Adding ", wins_per_segment, " panel windows from ", f))
			}		

			gp_windows = add_x_unique(gp_windows, panel_windows, wins_per_segment)
		
			if (verbose >= 3) {
				print(paste0("Length of panel after add: ", length(gp_windows)))
			}
        	}

		s = sub("_sig[0-9]+", "", f) # trim signature indicator from last file on the list
		s = sub("indsig_panel_windows_", "", s) # trim file prefix
	
		outfile = paste0(GLOBAL_GP_WINDOWS_DIR, "gp_panel_windows_", s)

		if (verbose >= 2) { print(paste0("Final panel length for Obj ", o, ", Iteration ", i, ": ", length(gp_windows))) }

		if (verbose >= 1) { print(paste0(Sys.time(), "    saving panel to file: ", outfile)) } 
		write(gp_windows, file=outfile)
	}
}




