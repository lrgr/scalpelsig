setwd("../..")
library(optparse)
source("GLOBAL_CONFIG.R")

option_list = list(
        make_option(c("-t", "--tag"), type="character", default=NULL,
                        help="tag to refer downstream scripts to the correct training sets", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

file_tag = opt$tag

if (is.null(file_tag)) {
        stop("No file_tag recieved (command line -t ). Please supply a file_tag.")
}

win_num_vec = c(10, 50, 100, 150, 200, 250)


files = list.files(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, pattern= paste0(".*", file_tag, ".*") )
print(paste0(Sys.time(), "    found ", length(files), " files containing the tag: ", file_tag))

i = 1
for (f in files) {
	print(paste0( Sys.time(), "    ", i, "/", length(files)) )
	
	panel_windows = scan(paste0(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, f), what=character(), quiet=TRUE)

	# sanity check, make sure original panel is big enough to truncate
	if (length(panel_windows) < max(win_num_vec)) {
		stop(paste0("There are only ", length(panel_windows), " windows in ", f, ", but I am trying to create a truncated panel with ", max(win_num_vec), " windows."))
	}	

	# truncate the panel window vector, save truncated vectors
	for (n in win_num_vec) {
		trunc_panel = panel_windows[1:n]

		s = sub("panel_windows_", "", f)		

		outfile = paste0(GLOBAL_SIZE_TEST_WINDOWS_DIR, "trunc_panel_windows_size", n, "_", s)
		
		write(trunc_panel, outfile)
	}
	i = i + 1
}
