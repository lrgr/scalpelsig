# functions to interpret / summarize the outputs of evaluate_10k_panel_results.R
setwd("../..")
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
# i.e., subset the df by Signature, and within each signature by panel size, and within each size by objective function
# take the median result (auroc/aupr) for the obtained panel, MSK-IMPACT, and WES of each subset, and record them
st_median_results_df <- function(df, verbose=2) {
	if ( !("Panel.Size" %in% colnames(df)) ) {
		stop("The data frame supplied to st_median_results_df() does not have a \'Panel.Size\' column. Please supply a data frame that does.")
	}
	
	sigs = sort(unique(df$Signature))
	sizes = sort(unique(df$Panel.Size))
	n_sigs = length(sigs)	
	n_sizes = length(sizes)

	if (verbose >= 2) { 
		print(paste0("Found ", length(sigs), " distinct signatures in given df:")) 
		print(sigs)
		print(paste0("Found ", length(sizes), " distinct panel sizes in given df:"))
		print(sizes)
	}
	

	Obj1.Score = numeric(n_sigs * n_sizes)
	Obj2.Score = numeric(n_sigs * n_sizes)
	Obj3.Score = numeric(n_sigs * n_sizes)
	
	MSK.IMPACT.Score = numeric(n_sigs * n_sizes)
	WES.Score = numeric(n_sigs * n_sizes)
	Signature = numeric(n_sigs * n_sizes)
	Eval.Mode = character(n_sigs * n_sizes)
	Panel.Size = numeric(n_sigs * n_sizes)

	Obj1.Raw.Sp = numeric(n_sigs * n_sizes)
	Obj2.Raw.Sp = numeric(n_sigs * n_sizes)
	Obj3.Raw.Sp = numeric(n_sigs * n_sizes)

	MSK.Raw.Sp = numeric(n_sigs * n_sizes)
	WES.Raw.Sp = numeric(n_sigs * n_sizes)

	Obj1.Norm.Sp = numeric(n_sigs * n_sizes)
	Obj2.Norm.Sp = numeric(n_sigs * n_sizes)
	Obj3.Norm.Sp = numeric(n_sigs * n_sizes)

	MSK.Norm.Sp = numeric(n_sigs * n_sizes)
	WES.Norm.Sp = numeric(n_sigs * n_sizes)


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

		for ( curr_size in sizes ) {
			# get all observations w/in the signature that have the current Panel.Size
			size_df = sig_df[sig_df$Panel.Size==curr_size, ]

			Eval.Mode[i] = em

			obj1_df = size_df[size_df$Obj.Fn==1, ] # get entries in df for obj fn 1, 2, 3
			obj2_df = size_df[size_df$Obj.Fn==2, ]
			obj3_df = size_df[size_df$Obj.Fn==3, ]

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
			Panel.Size[i] = curr_size

			# correlation measures
			if ("Raw.Spearman" %in% colnames(df)) {
				obj1_rsp = median(obj1_df$Raw.Spearman)
				obj2_rsp = median(obj2_df$Raw.Spearman)
				obj3_rsp = median(obj3_df$Raw.Spearman)

				msk_rsp = median(obj1_df$MSK.R.Spearman)
				wes_rsp = median(obj1_df$WES.R.Spearman)

				Obj1.Raw.Sp[i] = obj1_rsp
				Obj2.Raw.Sp[i] = obj2_rsp
				Obj3.Raw.Sp[i] = obj3_rsp

				MSK.Raw.Sp[i] = msk_rsp
				WES.Raw.Sp[i] = wes_rsp

				obj1_nsp = median(obj1_df$Norm.Spearman)
				obj2_nsp = median(obj2_df$Norm.Spearman)
				obj3_nsp = median(obj3_df$Norm.Spearman)
				
				msk_nsp = median(obj1_df$MSK.N.Spearman)
				wes_nsp = median(obj1_df$WES.N.Spearman)

				Obj1.Norm.Sp[i] = obj1_nsp
				Obj2.Norm.Sp[i] = obj2_nsp
				Obj3.Norm.Sp[i] = obj3_nsp

				MSK.Norm.Sp[i] = msk_nsp
				WES.Norm.Sp[i] = wes_nsp
				
			}


			i = i + 1
		}
	}
	
	if ("Raw.Spearman" %in% colnames(df)) {
		results_df = data.frame(Signature, Panel.Size, Obj1.Raw.Sp, Obj2.Raw.Sp, Obj3.Raw.Sp, MSK.Raw.Sp, WES.Raw.Sp, Obj1.Norm.Sp, Obj2.Norm.Sp, Obj3.Norm.Sp, MSK.Norm.Sp, WES.Norm.Sp, Obj1.Score, Obj2.Score, Obj3.Score, MSK.IMPACT.Score, WES.Score, Eval.Mode)
	} else {
		results_df = data.frame(Signature, Panel.Size, Obj1.Score, Obj2.Score, Obj3.Score, MSK.IMPACT.Score, WES.Score, Eval.Mode)
	}	

	return(results_df)
}


save_st_summary_df <- function(results_df_infile, outfile=NULL) {
	print(paste0("Loading panel results df from ", results_df_infile))
	df = load_panel_results_df(results_df_infile)

	print("results df dimensions: ")
	print(dim(df))
	
	summary_df = st_median_results_df(df)
	
	
	if (is.null(outfile)) {
		# use default outfile
		# this assumes that the infile was taken from GLOBAL_SCRIPT_OUT
		file = sub(GLOBAL_SCRIPT_OUT, "", results_df_infile) # remove file path from infile

		outfile = paste0(GLOBAL_SCRIPT_OUT, "SIZE_TEST_SUMMARY_", file)
	}
	
	print(paste0("writing summary df to ", outfile))
	write.table(summary_df, file=outfile, sep="\t", quote=FALSE)
}

print_guide <- function() {
	print("How to use this file:") 
	print("1) start an R session")
	print("2) call source(\'summarize_panel_size_results.R\')")
	print("3) list_results_files()")
	print("4) run save_st_summary_df() on paste0(GLOBAL_SCRIPT_OUT, <file name from step 1>)")
}

print_guide()
