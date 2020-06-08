setwd("..")
source("projection_score.R")


DEBUG_FLAG = TRUE

EVAL_MODE = "auroc"


# given a file name as a string fn, returns a list with information about the data that generated the file
parse_sig_est_file_name <- function(fn) {
	ret = list()

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

	s = sub("panel_sbs_df_", "", fn)
	tag = sub("_it[0-9]+_.*", "", s)
	ret[["tag"]] = tag

	return(ret)
}


file_tag = "test_run"
files = list.files(GLOBAL_SCRIPT_PANEL_SIG_EST_DIR, pattern= paste0(".*", file_tag, ".*") )
print(paste0(Sys.time(), "    found ", length(files), " files containing the tag: ", file_tag))

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

i = 1
for (f in files) {
	print(paste0(i, "/", length(files))) 
	#print(f)
	info = parse_sig_est_file_name(f)
	
	Signature[i] = info[["sig_num"]]
	Obj.Fn[i] = info[["obj_num"]]
	Iteration[i] = info[["it_num"]]
	File.Tag[i] = info[["tag"]]
	Timestamp.Tag[i] = info[["timestamp"]]
	Eval.Mode[i] = EVAL_MODE
	File.Name[i] = f

	#print(info)

	sig_est_outfile = paste0(GLOBAL_SCRIPT_PANEL_SIG_EST_DIR, f)	

	s = sub("panel_sig_est_", "", f) # strip "panel_sig_est_" from front of the string
	
	t = sub("_[0-9]+-.*", "", s) # remove the timestamp tag and trailing '.tsv'
	t = sub("_obj[0-9]+", "", t) # remove '_objNN' from filename
	tt_file = paste0(GLOBAL_SCRIPT_TEST_TRAIN_DIR, "test_train_", t, ".rds")
	test_train = readRDS(tt_file)
	test_set = test_train[[1]]
	
	#u = sub(".*it[0-9]*_", "", t) # strip everything before 'sig<number>' in the filename

	sig_num = info[["sig_num"]]

	if (EVAL_MODE=="auroc") {
		result = compute_panel_auroc(sig_num, test_set, sig_est_outfile, global_sig_df)
	} else if (EVAL_MODE=="aupr") {
		result = compute_panel_aupr(sig_num, test_set, sig_est_outfile, global_sig_df)
	} else {
		stop("eval_mode was something other than \'auroc\' or \'aupr\'")
	}

	Eval.Result[i] = result

	i = i + 1
}

results_df = data.frame(Eval.Result, Signature, Obj.Fn, Eval.Mode, Iteration, File.Tag, Timestamp.Tag, File.Name)
results_df = results_df[order(Signature, Obj.Fn, Iteration), ]

results_timestamp = format(Sys.time(), "%d-%b-%Y_%H-%M")

results_df_outfile = paste0(GLOBAL_SCRIPT_OUT, "panel_results_df_", results_timestamp, ".tsv")

write.table(results_df, file=results_df_outfile, sep="\t", quote=FALSE)
