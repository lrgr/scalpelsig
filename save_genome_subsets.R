# need to find the subset of mutations which fall under a large set of genomic regions
# (goal is to have a .tsv for whole exome, and for MSK-IMPACT panel)
# this is very inefficient with naive df subsetting, since there are so many query regions
# so this script will split the mut_df into chromosomes and use binary search to compute this subset more efficiently


source("preprocess_windows.R")
source("assess_panel.R")
source("GLOBAL_CONFIG.R")


###########################
# HELPER FUNCTIONS        #
###########################

note <- function(text) {
	s = paste0(Sys.time(), "    ", text)
	print(s)
}

print_iter <- function(current, max, initial=1:5, num_updates=NULL) {
	#print(initial)
	#print(current %in% initial)
	chrom_flag = FALSE
	if (is.null(num_updates)) {
		# want to print each update if iterating over chromosomes 
		if (max == 24) { 
			num_updates = 24
			chrom_flag = TRUE
		} else {
			num_updates = 20
		}
	}

	if (chrom_flag) {
		chrom_vec = c(1:22, "X", "Y")
		s = paste0(Sys.time(), "    ", "chrom ", chrom_vec[current], "/24")
		print(s)
	} else if ((current %in% initial) | (current %% (floor(max / num_updates)) == 0)) {
		s = paste0(Sys.time(), "    ", current, "/", max)
		print(s)
	}
}



save_panel_sbs_tsv <- function(sbs_df, outfile) {
	write.table(sbs_df, file=outfile, sep="\t", quote=FALSE)
}


#############################
# LOAD DATA                 #
#############################

# this loads the saved output of this script for gencode
load_gencode_exon_panel_muts <- function() {
	filename = paste0(GLOBAL_DATA_DIR, "gencode_exon_panel_mut_df.tsv")
	df = read.csv(file=filename, sep="\t", header=TRUE)
	return(df)
}

load_msk_panel_muts <- function() {
	filename = paste0(GLOBAL_DATA_DIR, "msk_region_panel_mut_df.tsv")
	df = read.csv(file=filename, sep="\t", header=TRUE)
	return(df)
}


load_gencode_exons <- function(unprocessed=FALSE) {
	filename = paste0(GLOBAL_DATA_DIR, "gencode_exons.tsv")
	df = read.csv(file=filename, sep="\t", header=TRUE)
	if (unprocessed) {
		return(df)
	} else {
		return(chr_dfs_from_gencode_df(df))
	}
}

load_msk_regions <- function(unprocessed=FALSE) {
	filename = paste0(GLOBAL_DATA_DIR, "msk_impact_regions.tsv")
	df = read.csv(file=filename, sep="\t", header=TRUE)
	if (unprocessed) {
		return(df)
	} else {
		return(chr_dfs_from_gencode_df(df))
	}
}



load_ft_regions <- function(sig_num, unprocessed=FALSE) {
	df = coord_df_from_windows(paste0(GLOBAL_DATA_DIR, "trained_panels/windows/ft_panel_windows_sig", sig_num, ".txt"))
	if (unprocessed) { return(df) }
	else { return(chr_dfs_from_gencode_df(df)) }
}







window_info_from_str <- function(window_str) {
        v = strsplit(window_str, "_")[[1]]
        info = list()
        info[["chr"]] = v[1]
        info[["start"]] = as.numeric(v[2])
        info[["end"]] = as.numeric(v[3])
        return(info)
}


# turn a list of windows into df of genome coordinates

coord_df_from_windows <- function(window_file) {
	panel_windows = scan(window_file, what=character(), quiet=TRUE) # read strings from panel window file
	n = length(panel_windows)
	
	# initialize vectors for dataframe
	Chromosome = character(n)
	Start = numeric(n)
	End = numeric(n)
	CHROM.NUMBER = numeric(n) # this is to make sorting easier, won't get into the final result df

	for (i in 1:n) {
		window_info = window_info_from_str(panel_windows[i])
		Chromosome[i] = window_info[["chr"]]
		Start[i] = window_info[["start"]]
		End[i] = window_info[["end"]]
		
		# get the 
		chr_str = window_info[["chr"]]
		chr_str = sub("chr", "", chr_str)
		if (chr_str=="X") { chr_num = 23 }
		else if (chr_str=="Y") { chr_num = 24 } 
		else { chr_num = as.numeric(chr_str) }
		CHROM.NUMBER[i] = chr_num
	}

	df = data.frame(Chromosome, Start, End, CHROM.NUMBER)
	df = df[order(CHROM.NUMBER, Start, End), ]
	df = df[ , 1:3] # remove CHROM.NUMBER column
	return(df)
}




# formatting gencode dataframe into chromosomes

# break into individual chromosomes
# regions need to be sorted by start coord
# some regions have same start site, ties in the order are broken with end coord
chr_dfs_from_gencode_df <- function(gencode_reg_df) {
	chrom_vec = c(1:22, "X", "Y")
	ret = list()
	for (i in 1:length(chrom_vec)) {
		curr_chrom = paste0("chr", chrom_vec[i])
		chrom_regs = gencode_reg_df[ gencode_reg_df$Chromosome==curr_chrom, ]
		
		chrom_regs = chrom_regs[ with(chrom_regs, order(chrom_regs$Start, chrom_regs$End)), ] # orders the regions by start site, breaking ties with end site
		
		ret[[curr_chrom]] = chrom_regs
	}

	return(ret)
}




#################################
# BINARY SEARCH ON REGIONS      #
#################################

# mut_coord should be given as an integer, the bp location of the mutation on the chromosome
# target_region_starts : sorted int vector, start of each region
# target_region_ends : int vector, target_region_ends[i] is the end of region[i] in target_region_ends
mutIsInRegions <- function(mut_coord, target_region_starts, target_region_ends, debug=FALSE) {
	left = 1
	right = length(target_region_starts)
	
	while (left <= right) {
		mid = floor((left + right) / 2)
		curr_start = target_region_starts[mid]
		curr_end = target_region_ends[mid]
		curr_region = c(curr_start, curr_end)

		if (mut_in_region(mut_coord, curr_region)) {
			if (debug) { 
				print("mutIsInRegions() FOUND overlap.")
				print(paste0("mut_coord: ", mut_coord, "    region start: ", curr_start, "    region_end: ", curr_end)) 
			} 
			return(TRUE)
		} else if (curr_start < mut_coord) {
			left = mid+1
		} else if (curr_end > mut_coord) {
			right = mid-1
		} else {
			stop("Error in mutIsInRegions, this shouldn't happen")
		}
	}
	return(FALSE)
}

# mut_coord is integer
# target_region is pair of integers
mut_in_region <- function(mut_coord, target_region) {
	r_start = target_region[1]
	r_end = target_region[2]
	if ((r_start <= mut_coord) & (r_end >= mut_coord)) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}



# mut_coord should be given as an integer, the bp location of the mutation on the chromosome
# target_regions should be a list of integer pairs, sorted by the first entry in the pair
# (the integer pairs are the start and end bp coordinate of each region)
mutIsInRegionsOLD <- function(mut_coord, target_regions) {
	left = 1
	right = length(target_regions)
	
	while (left <= right) {
		mid = floor((left + right) / 2)
		curr_region = target_regions[[mid]]
		curr_start = curr_region[1]
		curr_end = curr_region[2]

		if (mut_in_region(mut_coord, curr_region)) {
			return(TRUE)
		} else if (curr_start < mut_coord) {
			left = mid+1
		} else if (curr_end > mut_coord) {
			right = mid-1
		} else {
			stop("Error in mutIsInRegions, this shouldn't happen")
		}
	}
	return(FALSE)
}

# mut_coord is integer
# target_region is pair of integers
mut_in_regionOLD <- function(mut_coord, target_region) {
	r_start = target_region[[1]]
	r_end = target_region[[2]]
	if ((r_start <= mut_coord) & (r_end >= mut_coord)) {
		return(TRUE)
	} else {
		return(FALSE)
	}
}



################################
# SUBSET MUT DF                #
################################

# save_prefix should include information about which regions you are subsetting by, e.g. 'gencode_exome' or 'msk_impact'
region_subset_mut_df <- function(chrom_reg_ls, mut_df, save_prefix=NULL, save_df=TRUE, debug=TRUE) {
	chrom_vec = c(1:22, "X", "Y")

	ret_ls = list()
	if (debug) { note("region_subset_mut_df() beginning main loop") }

	for (i in 1:length(chrom_vec)) {
		if (debug) { print_iter(i, length(chrom_vec)) }
		chrom_df = mut_df[mut_df$Chromosome == chrom_vec[i],]
		chrom_target_regions_df = chrom_reg_ls[[i]]

		result = chrom_region_subset_mut_df(chrom_target_regions_df, chrom_df)

		if (save_df==TRUE) {
			if (!is.null(save_prefix)) {
				filename = paste0(save_prefix, "_region_subset_muts_chr", chrom_vec[i], ".tsv")
				if (debug) { note(paste0("region_subset_mut_df() saving to: ", filename)) }
				write.table(result, file=filename, sep="\t", quote=FALSE)
			}
		}
		
		ret_ls[[i]] = result
	}


	# combine remaining dfs into one panel df
	panel_df = foreach(i=1:length(ret_ls), .combine=rbind) %do% {
		ret_ls[[i]]
	}

	return(panel_df)
}

chrom_region_subset_mut_df <- function(chrom_target_regions, chrom_df, debug=TRUE) {
	mut_in_target = logical(nrow(chrom_df))

	foreach(i=1:nrow(chrom_df)) %do% {
		mut=chrom_df[i, "Start.Position"]
		
		if (debug) { print_iter(i, nrow(chrom_df)) }
		target_reg_starts = chrom_target_regions$Start
		target_reg_ends = chrom_target_regions$End

		mut_in_target[i] = mutIsInRegions(mut, target_reg_starts, target_reg_ends)
		FALSE
	}

	return(chrom_df[mut_in_target, ])
}





########################
# UNIT TESTS           #
########################


test_mutIsInRegions <- function() {
	#even_region_list = list()
	#R_sucks = 1
	#for (i in 1:20) {
	#	if (i %% 2 == 0) {
	#		even_region_list[[R_sucks]] = c(i, i)
	#		R_sucks = R_sucks + 1
	#	}
	#}

	even_region_vec = (1:10) * 2

	ret = logical(20)

	for(i in 1:20) {
		ret[i] = mutIsInRegions(i, even_region_vec, even_region_vec)
	}

	names(ret) = 1:20
	print(ret)
	correct = (1:20) %% 2 == 0
	if (sum(ret == correct) == 20) {
		print("TEST PASSED")
	}
}

# this isn't actually a unit test it just runs it on the smallest chromosome and you can manually inspect
test_chrom_region_subset_mut_df <- function(mut_df=NULL, gencode_chr_ls=NULL) {
	if (is.null(mut_df)) {
		mut_df = load_nz_mut_df_with_sigprob()
	}
	if (is.null(gencode_chr_ls)) {
		gencode_chr_ls = load_gencode_exons()
	}
	
	chr22_target_df = gencode_chr_ls[[22]]
	chr22_mut_df = mut_df[mut_df$Chromosome==22,]
	test = chrom_region_subset_mut_df(chr22_target_df, chr22_mut_df)
	return(test)
}



#############################
# POST PROCESSING           #
#############################

# post processing:
# 1. convert panel df into 96 sbs df
# 2. extract signatures from 96 sbs df and save them

# how to do this:
# 1. use sbs_df_from_mut_df() in preprocess_windows.R
# 2. call python script - NOTE: must use the signature-estimation-py conda env

post_process_panel_df <- function(panel_df, outfile_tag, alt_setting=0) {
	sbs_df = sbs_df_from_mut_df(panel_df)
	
	sbs_tsv_outfile = paste0(GLOBAL_PANEL_SBS_DIR, outfile_tag, "_panel_sbs.tsv")
	save_panel_sbs_tsv(sbs_df, sbs_tsv_outfile)

	if (alt_setting==0) {
		sig_est_outfile = paste0(GLOBAL_PANEL_SIG_EST_DIR, outfile_tag, "_panel_sig_est.tsv")
	
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_tsv_outfile,
			                                    " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
							    " -of ", sig_est_outfile)
		      )
	} else if (alt_setting==1) {
		sig_est_outfile = paste0(GLOBAL_PANEL_SIG_EST_DIR, outfile_tag, "_panel_sig_est_WITHOUT10_30.tsv")

		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_tsv_outfile,
							" -sf ", GLOBAL_DATA_DIR, "cosmic-sigs-without10_30.tsv",
							" -of ", sig_est_outfile)
			)
	}
}


post_process_ft_mut_df <- function(ft_mut_df, sig_num, alt_setting=0) {
	sbs_df = sbs_df_from_mut_df(ft_mut_df)
	
	sbs_tsv_outfile = paste0(GLOBAL_DATA_DIR, "trained_panels/", "ft_panel_sig", sig_num, "_sbs.tsv")
	save_panel_sbs_tsv(sbs_df, sbs_tsv_outfile)

	if (alt_setting==0) { 
		sig_est_outfile = paste0(GLOBAL_DATA_DIR, "trained_panels/", "ft_panel_sig", sig_num, "_sig_est.tsv")
	
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_tsv_outfile,
	        	                            " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
						    " -of ", sig_est_outfile)
			)
	} else if (alt_setting==1) {
		sig_est_outfile = paste0(GLOBAL_DATA_DIR, "trained_panels/", "ft_panel_sig", sig_num, "_sig_est_WITHOUT10_30.tsv")

		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_tsv_outfile,
							" -sf ", GLOBAL_DATA_DIR, "cosmic-sigs-without10_30.tsv",
							" -of ", sig_est_outfile)
			)
	}
}


SPECIAL_post_process_panel_df <- function(outfile_tag, alt_setting=0) {
	#sbs_df = sbs_df_from_mut_df(panel_df)
	
	sbs_tsv_outfile = paste0(GLOBAL_PANEL_SBS_DIR, outfile_tag, "_panel_sbs.tsv")
	#save_panel_sbs_tsv(sbs_df, sbs_tsv_outfile)

	if (alt_setting==0) {
		sig_est_outfile = paste0(GLOBAL_PANEL_SIG_EST_DIR, outfile_tag, "_panel_sig_est.tsv")
	
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_tsv_outfile,
			                                    " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
							    " -of ", sig_est_outfile)
		      )
	} else if (alt_setting==1) {
		sig_est_outfile = paste0(GLOBAL_PANEL_SIG_EST_DIR, outfile_tag, "_panel_sig_est_WITHOUT10_30.tsv")

		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_tsv_outfile,
							" -sf ", GLOBAL_DATA_DIR, "cosmic-sigs-without10_30.tsv",
							" -of ", sig_est_outfile)
			)
	}
}





SPECIAL_post_process_ft_mut_df <- function(sig_num, alt_setting=0) {
	#sbs_df = sbs_df_from_mut_df(ft_mut_df)
	
	sbs_tsv_outfile = paste0(GLOBAL_DATA_DIR, "trained_panels/", "ft_panel_sig", sig_num, "_sbs.tsv")
	#save_panel_sbs_tsv(sbs_df, sbs_tsv_outfile)

	if (alt_setting==0) { 
		sig_est_outfile = paste0(GLOBAL_DATA_DIR, "trained_panels/", "ft_panel_sig", sig_num, "_sig_est.tsv")
	
		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_tsv_outfile,
	        	                            " -sf ", GLOBAL_DATA_DIR, "cosmic-signatures.tsv ",
						    " -of ", sig_est_outfile)
			)
	} else if (alt_setting==1) {
		sig_est_outfile = paste0(GLOBAL_DATA_DIR, "trained_panels/", "ft_panel_sig", sig_num, "_sig_est_WITHOUT10_30.tsv")

		system(paste0("python signature-estimation-py/signature_estimation.py -mf ", sbs_tsv_outfile,
							" -sf ", GLOBAL_DATA_DIR, "cosmic-sigs-without10_30.tsv",
							" -of ", sig_est_outfile)
			)
	}
}




# want to see this panel's auroc for test sets used in assessing the signature panels
# so we have to dig up the log files (which contain the test sets)
get_comprable_aurocs <- function(panel_sig_est_file) {
	
}





# don't actually use this -
# this is a quick and dirty way of joining the dfs, I just don't have time to reprocess it right now

# load each chromosome csv and join them into one dataframe

join_gencode_exon_tsvs <- function() {
	chrom_vec = c(1:22, "X", "Y")

	df <- foreach(i=1:length(chrom_vec), .combine=rbind) %do% {
		print_iter(i, length(chrom_vec))
		load_gencode_exon_test(chrom_vec[i])
	}

	outfile = paste0(GLOBAL_DATA_DIR, "gencode_exon_panel_mut_df.tsv")
	write.table(df, file=outfile, sep="\t", quote=FALSE)
}



# load one of the test savefiles

load_gencode_exon_test <- function(chr) {
	filename = paste0("gencode_exon_panel_test_region_subset_muts_chr", chr, ".tsv")
	df = read.csv(file=filename, sep="\t", header=TRUE)
}


join_msk_chrom_tsvs <- function() {
	chrom_vec = c(1:22, "X", "Y")

	df <- foreach(i=1:length(chrom_vec), .combine=rbind) %do% {
		print_iter(i, length(chrom_vec))
		load_msk_chrom_tsv(chrom_vec[i])
	}

	outfile = paste0(GLOBAL_DATA_DIR, "msk_region_panel_mut_df.tsv")
	write.table(df, file=outfile, sep="\t", quote=FALSE)
}

load_msk_chrom_tsv <- function(chr) {
	filename = paste0("msk_impact_panel_test_region_subset_muts_chr", chr, ".tsv")
	df = read.csv(file=filename, sep="\t", header=TRUE)
	return(df)
}


######################################

main <- function(mut_df=NULL, gencode_chr_ls=NULL, msk_chr_ls = NULL) {
	if (is.null(mut_df)) {
		mut_df = load_nz_mut_df_with_sigprob()
	}
	if (is.null(gencode_chr_ls)) {
		gencode_chr_ls = load_gencode_exons()
	}
	if (is.null(msk_chr_ls)) {
		msk_chr_ls = load_msk_regions()
	}

	#registerDoParallel(cores=GLOBAL_NCORES)
	
	#region_subset_mut_df(gencode_chr_ls, mut_df, save_prefix="gencode_exon_panel_test")
	region_subset_mut_df(msk_chr_ls, mut_df, save_prefix="msk_impact_panel_test")
}

ft_main <- function(mut_df=NULL) {
	if (is.null(mut_df)) {
		mut_df = load_staaf_mutation_df()
	}

	#for (sig_num in c(2,3,8,13,18)) {
	#	ftp_chr_ls = load_ft_regions(sig_num)
	#	panel_mut_df = region_subset_mut_df(ftp_chr_ls, mut_df)
	#	post_process_ft_mut_df(panel_mut_df, sig_num)
	#}

	msk_chr_ls = load_msk_regions()
	msk_mut_df = region_subset_mut_df(msk_chr_ls, mut_df)
	post_process_panel_df(msk_mut_df, outfile_tag = "msk_staaf")
}


