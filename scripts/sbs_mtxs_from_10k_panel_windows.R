setwd("..")
source("projection_score.R")

library(optparse)

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

#file_tag = "SEC_BATCH_BAL"

print(paste0(Sys.time(), "    loading 10kb sbs array"))
sbs_arr_ls_10k = load_10k_sbs_arr_ls()
print(paste0(Sys.time(), "    done."))
print(paste0(Sys.time(), "    saving sbs matrices for panels with file_tag: ", file_tag))

files = list.files(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, pattern= paste0(".*", file_tag, ".*") )
print(paste0(Sys.time(), "    found ", length(files), " files containing the tag: ", file_tag))

i = 1
for (f in files) {
	print( paste0(Sys.time(), "    ", i, "/",  length(files)) )
	panel_windows = scan(paste0(GLOBAL_SCRIPT_PANEL_WINDOWS_DIR, f), what=character(), quiet=TRUE)
	panel_df = get_panel_sbs_df(panel_windows, sbs_arr_ls_10k)

	s = sub("panel_windows_", "", f) #trim 'panel_windows_' prefix from file name
	s = sub(".txt", "", s) #trim ".txt" from file name
	outfile = paste0(GLOBAL_SCRIPT_PANEL_SBS_DIR, "panel_sbs_df_", s, ".tsv")
	
	save_panel_sbs_tsv(panel_df, outfile)
	gc()
	i = i + 1
}
