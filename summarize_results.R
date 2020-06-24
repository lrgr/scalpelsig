# functions to interpret / summarize the outputs of evaluate_10k_panel_results.R

source("GLOBAL_CONFIG.R")

# example:
# load_panel_results_df(paste0(GLOBAL_SCRIPT_OUT, "panel_results_df_15-Jun-2020_13-57.tsv"))
load_panel_results_df <- function(infile) {
	return(read.csv(infile, sep="\t"))
}

list_results_files <- function() {
	return(list.files(GLOBAL_SCRIPT_OUT, pattern=".tsv"))
}


# given results data frame df,
# get medians of the relevant chunks of the results
# i.e., subset the df by Signature, and within each signature by objective function
# take the median result (auroc/aupr) for the obtained panel, MSK-IMPACT, and WES of each subset, and record them
median_results_df <- function(df, verbose=2) {
	sigs = sort(unique(df$Signature))
	n_sigs = length(sigs)	

	if (verbose >= 2) { 
		print(paste0("Found ", length(sigs), " distinct signatures in given df:")) 
		print(sigs)
	}
	

	Obj1.Score = numeric(n_sigs)
	Obj2.Score = numeric(n_sigs)
	Obj3.Score = numeric(n_sigs)
	
	MSK.IMPACT.Score = numeric(n_sigs)
	WES.Score = numeric(n_sigs)
	Signature = numeric(n_sigs)
	Eval.Mode = character(n_sigs)


	# iterate through signatures
	i = 1
	for (s in sigs) {
		sig_df = df[df$Signature==s, ] # get entries in df with current signature

		# sanity check to make sure that all the observations in a signature are either AUROC or AUPR
		# (if both are within a single signature, then there is a problem)
		curr_eval_mode = as.character(sig_df$Eval.Mode)
		if (length( unique(curr_eval_mode) ) != 1) {
			print("Sanity check failed:")
			print(paste0("Signature ", s, " contained ", length(unique(curr_eval_mode)), " distinct entries for Eval.Mode"))
			print(unique(curr_eval_mode))
			stop()
		}
		em = curr_eval_mode[1]
		Eval.Mode[i] = em

		obj1_df = sig_df[sig_df$Obj.Fn==1, ] # get entries in df for obj fn 1, 2, 3
		obj2_df = sig_df[sig_df$Obj.Fn==2, ]
		obj3_df = sig_df[sig_df$Obj.Fn==3, ]

		obj1_res = obj1_df$Eval.Result
		obj2_res = obj2_df$Eval.Result
		obj3_res = obj3_df$Eval.Result

		# the obj1, 2, and 3 dfs have the same test sets, so the benchmark panels repeat their results
		# so it is sufficient to just take 1 of the obj dfs.
		msk_impact_res = obj1_df$MSK.IMPACT.Result
		wes_res = obj1_df$WES.Result

		# get medians across each subset of results
		obj1_med = median(obj1_res)
		obj2_med = median(obj2_res)
		obj3_med = median(obj3_res)

		mski_med = median(msk_impact_res)
		wes_med = median(wes_res)
		
		# place scores into appropriate vectors

		Obj1.Score[i] = obj1_med
		Obj2.Score[i] = obj2_med
		Obj3.Score[i] = obj3_med

		MSK.IMPACT.Score[i] = mski_med
		WES.Score[i] = wes_med
		
		Signature[i] = s

		i = i + 1	
	}
	
	results_df = data.frame(Signature, Obj1.Score, Obj2.Score, Obj3.Score, MSK.IMPACT.Score, WES.Score, Eval.Mode)
	
	return(results_df)
}


save_summary_df <- function(results_df_infile, outfile=NULL) {
	print(paste0("Loading panel results df from ", results_df_infile))
	df = load_panel_results_df(results_df_infile)

	print("results df dimensions: ")
	print(dim(df))
	
	summary_df = median_results_df(df)
	
	
	if (is.null(outfile)) {
		# use default outfile
		# this assumes that the infile was taken from GLOBAL_SCRIPT_OUT
		file = sub(GLOBAL_SCRIPT_OUT, "", results_df_infile) # remove file path from infile

		outfile = paste0(GLOBAL_SCRIPT_OUT, "SUMMARY_", file)
	}
	
	print(paste0("writing summary df to ", outfile))
	write.table(summary_df, file=outfile, sep="\t", quote=FALSE)
}

print_guide <- function() {
	print("How to use this file:") 
	print("1) list_results_files()")
	print("2) run save_summary_df() on paste0(GLOBAL_SCRIPT_OUT, <file name from step 1>)")
}

print_guide()
