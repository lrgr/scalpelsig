source("projection_score.R")

DEBUG_FLAG = TRUE


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

	s = sub("panel_sig_est_", "", fn)
	tag = sub("_it[0-9]+_.*", "", s)
	ret[["tag"]] = tag

	return(ret)
}

# 0 -> exposures we extracted using cosmic sigs {1, 2, 3, 5, 6, 8, 10, 13, 17, 18, 20, 26, 30}
# 1 -> exposures we extracted using cosmic sigs {1, 2, 3, 5, 6, 8, 13, 17, 18, 20, 26}
# 2 -> exposures taken directly from the supplementary material of Staaf et al 2019
STAAF_EXPOSURES_CHOICE = 2

# 0 -> exposures extracted using cosmic sigs {1, 2, 3, 5, 6, 8, 10, 13, 17, 18, 20, 26, 30}
# 1 -> exposures extracted using cosmic sigs {1, 2, 3, 5, 6, 8, 13, 17, 18, 20, 26}
PANEL_EXPOSURES_CHOICE = 0



staaf_sig_df = load_staaf_sig_estimates(norm=TRUE, alt_setting=STAAF_EXPOSURES_CHOICE)
no_norm_staaf_sig_df = load_staaf_sig_estimates(norm=FALSE, alt_setting=STAAF_EXPOSURES_CHOICE)


n = 5
Signature = numeric(n)
Panel.AUPR = numeric(n)

MSK.IMPACT.AUPR = numeric(n)
Percent.Active = numeric(n)

Norm.Spearman = numeric(n)
MSK.N.Spearman = numeric(n)

Raw.Spearman = numeric(n)
R.Sp.pval = numeric(n)
MSK.R.Spearman = numeric(n)
MSK.R.Sp.pval = numeric(n)


STAAF.EXP.ARG = numeric(n)
PAN.EXP.ARG = numeric(n)
act_thresh = 0.05

i = 1
# loop through each file containing the 
for (sig_num in c(2,3,8,13,18)) {
	print(paste0(Sys.time(), "    ", "evaluating ft panel for signature ", sig_num)) 

	STAAF.EXP.ARG[i] = STAAF_EXPOSURES_CHOICE
	PAN.EXP.ARG[i] = PANEL_EXPOSURES_CHOICE	

	Signature[i] = sig_num

	if (PANEL_EXPOSURES_CHOICE==0) {
		sig_est_file = paste0(GLOBAL_DATA_DIR, "trained_panels/", "ft_panel_sig", sig_num, "_sig_est.tsv")
		msk_staaf_sig_est = paste0(GLOBAL_PANEL_SIG_EST_DIR, "msk_staaf_panel_sig_est.tsv")
	} else {
		sig_est_file = paste0(GLOBAL_DATA_DIR, "trained_panels/", "ft_panel_sig", sig_num, "_sig_est_WITHOUT10_30.tsv")
		msk_staaf_sig_est = paste0(GLOBAL_PANEL_SIG_EST_DIR, "msk_staaf_panel_sig_est_WITHOUT10_30.tsv")
	}

	test_set = staaf_sig_df$Patient # all patients in Staaf et al dataset

	Percent.Active[i] = get_percent_active(sig_num, staaf_sig_df, act_thresh)

	# COMPUTE AUROC / AUPR OF PANEL
	result = compute_panel_aupr(sig_num, test_set, sig_est_file, staaf_sig_df, activation_thresh=act_thresh)		

	# benchmark panel results
	msk_result = compute_panel_aupr(sig_num, test_set, msk_staaf_sig_est, staaf_sig_df, activation_thresh=act_thresh)
	


	Panel.AUPR[i] = result
	MSK.IMPACT.AUPR[i] = msk_result


	#spearman computation

        panel_sp_norm = compute_panel_spearman(sig_num, test_set, sig_est_file, staaf_sig_df)
        msk_sp_norm = compute_panel_spearman(sig_num, test_set, msk_staaf_sig_est, staaf_sig_df)
        

        Norm.Spearman[i] = panel_sp_norm
        MSK.N.Spearman[i] = msk_sp_norm
        

	panel_sp_nn_htest = compute_panel_spearman_htest(sig_num, test_set, sig_est_file, no_norm_staaf_sig_df)
	msk_sp_nn_htest = compute_panel_spearman_htest(sig_num, test_set, msk_staaf_sig_est, no_norm_staaf_sig_df)

	panel_sp_nonorm = panel_sp_nn_htest$estimate
	panel_pval = panel_sp_nn_htest$p.value

	msk_sp_nonorm = msk_sp_nn_htest$estimate
	msk_pval = msk_sp_nn_htest$p.value

        #panel_sp_nonorm = compute_panel_spearman(sig_num, test_set, sig_est_file, no_norm_staaf_sig_df)
        #msk_sp_nonorm = compute_panel_spearman(sig_num, test_set, msk_staaf_sig_est, no_norm_staaf_sig_df)

        Raw.Spearman[i] = panel_sp_nonorm
        MSK.R.Spearman[i] = msk_sp_nonorm

	R.Sp.pval[i] = panel_pval
	MSK.R.Sp.pval[i] = msk_pval
        i = i + 1
}

results_timestamp = format(Sys.time(), "%d-%b-%Y_%H-%M")
	
# results df without random baseline
results_df = data.frame(Signature, Raw.Spearman, R.Sp.pval, MSK.R.Spearman, MSK.R.Sp.pval, Norm.Spearman, MSK.N.Spearman, Panel.AUPR, MSK.IMPACT.AUPR, Percent.Active, STAAF.EXP.ARG, PAN.EXP.ARG)

if (STAAF_EXPOSURES_CHOICE==0) { results_df_outfile = paste0(GLOBAL_OUT_DIR, "ft_panel_results_df_", results_timestamp, ".tsv") }
if (STAAF_EXPOSURES_CHOICE==1) { results_df_outfile = paste0(GLOBAL_OUT_DIR, "ft_panel_results_df_without10_30_", results_timestamp, ".tsv") }
if (STAAF_EXPOSURES_CHOICE==2) { results_df_outfile = paste0(GLOBAL_OUT_DIR, "ft_panel_results_df_SUPP_EXPOSURES_", results_timestamp, ".tsv") }


results_df = results_df[order(Signature), ]

write.table(results_df, file=results_df_outfile, sep="\t", quote=FALSE)
