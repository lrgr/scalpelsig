setwd("../..")
source("projection_score.R")
library(optparse)

option_list = list(
	make_option(c("-t", "--tags"), type="character", default=NULL,
		help="tags, separated by comma, of panels to be evaluated (e.g. \"first_batch,second_batch,third_batch\")", metavar="character"),
	make_option(c("-e", "--evalmode"), type="character", default=NULL,
		help="evaluation mode for panels. Can be set to \"auroc\" or \"aupr\".", metavar="character")
);

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

file_tags = opt$tags
if (is.null(file_tags)) {
        stop("No file_tags recieved (command line -t ). Please supply file_tags.")
}

EVAL_MODE = opt$evalmode
if (is.null(EVAL_MODE)) {
	stop("No EVAL_MODE recieved (command line -e ). Please supply an EVAL_MODE (either \"auroc\" or \"aupr\")" )
}


DEBUG_FLAG = TRUE


# given a file name as a string fn, returns a list with information about the data that generated the file
parse_trunc_sig_est_file_name <- function(fn) {
	ret = list()
	
	size_chunk = regmatches(fn, regexpr("_size[0-9]+_", fn))
	size_chunk = sub("_size", "", size_chunk)
	size_chunk = sub("_", "", size_chunk)
	size = as.numeric(size_chunk)
	ret[["size"]] = size

	sig_chunk = regmatches(fn, regexpr("_sig[0-9]+_", fn)) # chunk of the string that looks like '_sigNN_'
	sig_chunk = sub("_sig", "", sig_chunk) # remove '_sig' from front of chunk
	sig_chunk = sub("_", "", sig_chunk) # remove trailing underscore
	sig_num = as.numeric(sig_chunk) 
	ret[["sig_num"]] = sig_num

	obj_chunk = regmatches(fn, regexpr("_obj[0-9]+_", fn))
	obj_chunk = sub("_obj", "", obj_chunk)
	obj_chunk = sub("_", "", obj_chunk)
	obj_num = as.numeric(obj_chunk) 
	ret[["obj_num"]] = obj_num

	it_chunk = regmatches(fn, regexpr("_it[0-9]+_", fn))
	it_chunk = sub("_it", "", it_chunk)
	it_chunk = sub("_", "", it_chunk)
	it_num = as.numeric(it_chunk)
	ret[["it_num"]] = it_num

	s = sub(".tsv", "", fn)
	timestamp = sub(".*obj[0-9]+_", "", s)
	ret[["timestamp"]] = timestamp

	s = sub("panel_sig_est_", "", fn)
	tag = sub("_it[0-9]+_.*", "", s)
	ret[["tag"]] = tag

	return(ret)
}


tag_ls = strsplit(file_tags, ",")[[1]]

files = character(0)
for (file_tag in tag_ls) {
	curr_files = as.character(list.files(GLOBAL_SIZE_TEST_SIG_EST_DIR, pattern= paste0(".*", file_tag, ".*") ))
	files = c(files, curr_files)
	print(paste0(Sys.time(), "    found ", length(curr_files), " files containing the tag: ", file_tag))
}
print(paste0(Sys.time(), "    found ", length(files), " files in total."))

global_sig_df = load_nz_sig_estimates(norm=TRUE)

n = length(files)
Signature = numeric(n)
Obj.Fn = numeric(n)
Iteration = numeric(n)
File.Tag = character(n)
Timestamp.Tag = character(n)
Eval.Mode = character(n)
Eval.Result = numeric(n)
File.Name = character(n)

Panel.Size = numeric(n)

MSK.IMPACT.Result = numeric(n)
WES.Result = numeric(n)

Est.Pval = numeric(n)
Baseline.Med = numeric(n)
Baseline.Max = numeric(n)
BP.Max.File = character(n)

i = 1
# loop through each file containing the 
for (f in files) {
	print(paste0(Sys.time(), "    ", i, "/", length(files))) 
	#print(f)
	info = parse_trunc_sig_est_file_name(f)
	
	Signature[i] = info[["sig_num"]]
	Obj.Fn[i] = info[["obj_num"]]
	Iteration[i] = info[["it_num"]]
	File.Tag[i] = info[["tag"]]
	Timestamp.Tag[i] = info[["timestamp"]]
	Eval.Mode[i] = EVAL_MODE
	File.Name[i] = f

	Panel.Size[i] = info[["size"]]

	#print(info)

	sig_est_outfile = paste0(GLOBAL_SIZE_TEST_SIG_EST_DIR, f)	

	s = sub("trunc_sig_est_size[0-9]+_", "", f) # strip "panel_sig_est_" from front of the string
	
	t = sub("_[0-9]+-.*", "", s) # remove the timestamp tag and trailing '.tsv'
	t = sub("_obj[0-9]+", "", t) # remove '_objNN' from filename
	tt_file = paste0(GLOBAL_SCRIPT_TEST_TRAIN_DIR, "test_train_", t, ".rds")
	test_train = readRDS(tt_file)
	test_set = test_train[[1]]
	
	#u = sub(".*it[0-9]*_", "", t) # strip everything before 'sig<number>' in the filename

	sig_num = info[["sig_num"]]

	msk_impact_sig_est = paste0(GLOBAL_PANEL_SIG_EST_DIR, "msk_test_panel_sig_est.tsv")
	wes_sig_est = paste0(GLOBAL_PANEL_SIG_EST_DIR, "gencode_exon_panel_sig_est.tsv")

	# COMPUTE AUROC / AUPR OF PANEL
	if (EVAL_MODE=="auroc") {
		result = compute_panel_auroc(sig_num, test_set, sig_est_outfile, global_sig_df)
		
		# benchmark panel results
		msk_result = compute_panel_auroc(sig_num, test_set, msk_impact_sig_est, global_sig_df)
		wes_result = compute_panel_auroc(sig_num, test_set, wes_sig_est, global_sig_df)
	} else if (EVAL_MODE=="aupr") {
		result = compute_panel_aupr(sig_num, test_set, sig_est_outfile, global_sig_df)		

		# benchmark panel results
		msk_result = compute_panel_aupr(sig_num, test_set, msk_impact_sig_est, global_sig_df)
		wes_result = compute_panel_aupr(sig_num, test_set, wes_sig_est, global_sig_df)
	} else {
		stop("eval_mode was something other than \'auroc\' or \'aupr\'")
	}

	Eval.Result[i] = result
	MSK.IMPACT.Result[i] = msk_result
	WES.Result[i] = wes_result

	# random baseline computation
	#print("computing baseline vec")
	#baseline_vec = compute_baseline_auroc(sig_num, test_set, global_sig_df, eval_mode=EVAL_MODE)
	#print("done.")
	
	#baseline_med = median(baseline_vec)
	#baseline_max = max(baseline_vec)
	
	#max_index = which(baseline_vec == baseline_max)
	#max_bp = names(baseline_vec)[max_index]
	#BP.Max.File[i] = max_bp

	#n_better = sum(baseline_vec >= result)
	#pval = n_better / length(baseline_vec)
	
	#Est.Pval[i] = pval
	#Baseline.Med[i] = baseline_med
	#Baseline.Max[i] = baseline_max

	i = i + 1
}

# results df without random baseline
results_df = data.frame(Panel.Size, Eval.Result, MSK.IMPACT.Result, WES.Result, Signature, Obj.Fn, Iteration, Eval.Mode, File.Tag, Timestamp.Tag, File.Name)

# results df with random baseline
#results_df = data.frame(Eval.Result, MSK.IMPACT.Result, WES.Result, Est.Pval, Baseline.Med, Baseline.Max, Signature, Obj.Fn, Eval.Mode, Iteration, File.Tag, Timestamp.Tag, File.Name, BP.Max.File)

results_df = results_df[order(Signature, Obj.Fn, -Panel.Size, -Eval.Result), ]

results_timestamp = format(Sys.time(), "%d-%b-%Y_%H-%M")

results_df_outfile = paste0(GLOBAL_SCRIPT_OUT, "size_test_results_df_", results_timestamp, ".tsv")

write.table(results_df, file=results_df_outfile, sep="\t", quote=FALSE)
