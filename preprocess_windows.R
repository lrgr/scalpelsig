# process list of mutations into windowed mutation counts



##################
# loading data   #
##################



load_nz_mutation_df <- function() {
	df = read.csv(file="~/tests/sigma_test/sigma/data/mutations/raw/ICGC-BRCA-EU.RELEASE_22.SBS.renamed.tsv",
		 sep="\t", header=TRUE)
	return(df)
}

load_nz_mut_df_with_sigprob <- function() {
	df = read.csv(file="~/projects/hotspot_signature_panel/data/nz_mut_df_with_sig_prob.tsv",
		      sep="\t", header=TRUE)
	return(df)
}

load_nz_sig_estimates <- function(norm=FALSE) {
	df = read.csv(file="~/projects/hotspot_signature_panel/data/nz_signature_estimation.tsv",
		      sep="\t", header=TRUE)
	colnames(df)[1] = "Patient"
	
	if (norm==TRUE) {
		rs = rowSums(df[,-1])
		df[,-1] = df[,-1] / rs
	}

	return(df)
}

load_COSMIC_signatures <- function() {
	df = read.csv(file="~/projects/hotspot_signature_panel/data/cosmic-signatures.tsv",
		      sep="\t", header=TRUE)
	colnames(df)[1] = "Signature"
	
	v = character(nrow(df))
	for (i in 1:nrow(df)) {
		v[i] = paste0("Signature.", i)
	}
	df[,1] = v

	return(df)
}


load_nz_sbs_mtx <- function(norm=FALSE) {
	df = read.csv(file="~/projects/hotspot_signature_panel/data/nz_sbs.tsv",
		      sep="\t", header=TRUE)
	
	if (norm==TRUE) {
		rs = rowSums(df)
		df = df / rs
	}
	return(df)
}


load_aggregated_window_sigprob <- function() {
	df = read.csv(file="~/projects/hotspot_signature_panel/data/nz_aggregated_window_sigprob.tsv",
		      sep="\t", header=TRUE)
	return(df)
}


# don't use this
# this is named differently than the mutation df
#load_nz_sbs <- function() {
#	df = read.csv(file="~/projects/hotspot_signature_panel/data/test_runs/counts.ICGC-BRCA-EU_BRCA_22.SBS-96.tsv",
#		      sep="\t", header=TRUE)
#	return(df)
#}



###############################
# bookkeeping functions       #
###############################

samples <- function(df) {
	return(unique(df$Patient))
}

signature_names <- function(n=30) {
	ret = character(n)
	for (i in 1:n) {
		ret[i] = paste0("Signature.", i)
	}
	return(ret)
}

mut_types <- function(dots=FALSE) {
	sbs_96 = character(96)
	i = 1
	
	if (dots==F) {
		sbs_ls = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
	} else {
		sbs_ls = c("C.A", "C.G", "C.T", "T.A", "T.C", "T.G")
	}
	base_ls = c("A", "C", "G", "T")
	for (sbs in sbs_ls) {
		for (left_base in base_ls) {
			for (right_base in base_ls) {
				if (dots==F) {
					sbs_96[i] = paste0(left_base, "[", sbs, "]", right_base)
				} else {
					sbs_96[i] = paste0(left_base, ".", sbs, ".", right_base)
				}
				i = i + 1
			}
		}
	}
	return(sbs_96)
}

print_loop_progress <- function(curr, end, interval = -1) {
	if (interval == -1) {
		# then use default setting
		interval = round(end / 10)
	}
	if (curr==1 | curr==2 | curr==3 | curr %% interval == 0) {
		print(paste0(curr, "/", end))
	}
}


########################
# processing data      #
########################

# given a dataframe in mutation list format (e.g. from load_nz_mutation_df() )
# output a dataframe in sbs_96 format (i.e. each row is a sample, each column is a trinucleotide context mutation category)
sbs_df_from_mut_df <- function(df, debug=FALSE) {
	sample_names = samples(df)
	sbs_names = mut_types()

	sbs_mtx = matrix(nrow = length(sample_names), ncol=length(sbs_names))
	rownames(sbs_mtx) = sample_names
	colnames(sbs_mtx) = sbs_names

	i = 1
	for (s in sample_names) {
		sample_muts = df[df$Patient==s,]
		for (m in sbs_names) {
			count = length(which(sample_muts$SBS_96==m))
			sbs_mtx[s, m] = count
		}
	}
	
	return(as.data.frame(sbs_mtx))
}



#### helper functions for window_df_from_mut_df() ####

sample_df <- function(df, sample) {
	return(df[df$Patient==sample, ])
}

windows_in_chromosome <- function(df, chrom, win_size=10^6) {
	c_df = df[df$Chromosome==chrom,]
	last_mut_position = max(c_df$Start.Position)
	v = seq(0, last_mut_position + (win_size -1), win_size)
	pairs = list()
	for (i in 1:(length(v))) {
		pairs[[i]] = c(v[i], v[i] + (win_size))
	}
	return(pairs)
}

window_names <- function(window_ls, chrom) {
	v = character(length(window_ls))
	i = 1
	for (w in window_ls) {
		v[i] = paste0("chr", chrom, "_", w[1], "_", w[2])
		i = i + 1
	}

	return(v)
}

mut_windows_in_chr <- function(df, chrom, win_size=10^4) {
	c_df = df[df$Chromosome==chrom,]
	mut_positions = unique(c_df$Start.Position)
	pairs = list()
	for (i in 1:length(mut_positions)) {
		pairs[[i]] = c(mut_positions[i], mut_positions[i] + (win_size))
	}
	return(pairs)
}

chrom_window_mtx <- function(df, chrom, debug=FALSE) {
	samp_names = samples(df)
	windows = windows_in_chromosome(df, chrom)

	ret = matrix(nrow = length(samp_names), ncol = length(windows))
	rownames(ret) = samp_names
	colnames(ret) = window_names(windows, chrom)

	chrom_df = df[df$Chromosome==chrom,]

	sample_count = 1
	for (s in samp_names) {
		if (debug & (sample_count==1 | sample_count %% 50 == 0 | sample_count == length(samp_names))) { 
			print(paste0("Sample ", sample_count, "/", length(samp_names)) ) 
		}
		sample_count = sample_count + 1
		s_df = sample_df(chrom_df, s)
		i = 1 #counts windows
		for (w in windows) {
			sw_muts = s_df[s_df$Start.Position >= w[1] & s_df$Start.Position < w[2], ]
			ret[s, i] = dim(sw_muts)[1]
			i = i + 1
		}
	}
	return(ret)
}

#### end helper functions for window_df_from_mut_df() ####



# given a dataframe in mutation list format
# output a dataframe in windowed format (each row is a sample, each column is the number of mutations in a given genome window)
window_df_from_mut_df <- function(df, debug=FALSE) {
	samp_names = samples(df)
	chrom_list = c(1:22, "X", "Y")

	ret = matrix(nrow=length(samp_names), ncol=0)
	for (chrom in chrom_list) {
		if (debug) { print(paste0("Computing window matrix for chromosome ", chrom)) }
		chrom_mtx = chrom_window_mtx(df, chrom, debug=debug)
		ret = cbind(ret, chrom_mtx)
	}

	return( as.data.frame(ret) )
}

############################

agg_sigprob_chrom_window_mtx <- function(df, chrom, win_size=10^6, mut_windows=FALSE, debug=FALSE) {
	sig_names = signature_names()

	if (mut_windows==FALSE) {
		windows = windows_in_chromosome(df, chrom, win_size)
	} else {
		windows = mut_windows_in_chr(df, chrom, win_size)
	}
	#print(windows)

	ret = matrix(nrow = length(sig_names), ncol = length(windows))
	rownames(ret) = sig_names
	colnames(ret) = window_names(windows, chrom)

	chrom_df = df[df$Chromosome==chrom,]

	i = 1 # counts windows
	for (w in windows) {
		# signature probabilities for each mutation within the window
		w_sigprob_df = chrom_df[chrom_df$Start.Position >= w[1] & chrom_df$Start.Position < w[2], sig_names]
		agg = colSums(w_sigprob_df, na.rm=TRUE)
		if(debug) { print_loop_progress(i, length(windows), 2000) }
		ret[ ,i] = agg
		i = i + 1
	}
	return(ret)
}


aggregate_window_sigprob_mtx <- function(mut_sigprob_df, win_size=10^6, mut_windows=FALSE, debug=TRUE) {
	signatures = signature_names() 
	chrom_list = c(1:22, "X", "Y")

	ret = matrix(nrow=length(signatures), ncol=0)
	for (chrom in chrom_list) {
		if (debug) { print(paste0("Computing aggregate window matrix for chromosome ", chrom)) }
		chrom_mtx = agg_sigprob_chrom_window_mtx(mut_sigprob_df, chrom, win_size=win_size, mut_windows=mut_windows, debug=debug)
		ret = cbind(ret, chrom_mtx)
	}

	return(ret)
}


library(doParallel)
registerDoParallel(cores=3)

# does same as above function but parallelized
par_aggregate_window_sigprob_mtx <- function(mut_sigprob_df, win_size=10^6, mut_windows=FALSE, debug=TRUE) {
	signatures = signature_names()
	chrom_list = c(1:22, "X", "Y")

	#ret = matrix(nrow=length(signatures), ncol=0)
	ret <- foreach (i = 1:24, .combine=cbind) %dopar% {
		chrom = chrom_list[i]
		chrom_mtx = agg_sigprob_chrom_window_mtx(mut_sigprob_df, chrom, win_size=win_size, mut_windows=mut_windows, debug=debug)
	 }
	return(ret)
}



par_chr_mtx_sigprob_ls <- function(mut_sigprob_df, win_size=10^6, mut_windows=FALSE, debug=TRUE) {
	signatures = signature_names()
	chrom_list = c(1:22, "X", "Y")

	#ret = matrix(nrow=length(signatures), ncol=0)
	ret <- foreach (i = 1:24) %dopar% {
		chrom = chrom_list[i]
		chrom_mtx = agg_sigprob_chrom_window_mtx(mut_sigprob_df, chrom, win_size=win_size, mut_windows=mut_windows, debug=debug)
	}
	return(ret)
}



###########################################################
# estimate signature probabilities for each mutation      #
###########################################################

# run-from-scratch wrapper for the two functions below
est_mutwise_signature_probabilities <- function(normalize=TRUE, debug=TRUE) {
	norm_exp_mtx = load_nz_sig_estimates(norm=TRUE)
	norm_sbs_mtx = load_nz_sbs_mtx(norm=TRUE)
	cosmic_signatures = load_COSMIC_signatures()
	if (debug) {
		print(paste0("norm_exp_mtx dimensions: ", dim(norm_exp_mtx)))
		print(paste0("norm_sbs_mtx dimensions: ", dim(norm_sbs_mtx)))
		print(paste0("cosmic_signatures dimensions: ", dim(cosmic_signatures)))
	}

	ret = get_sample_sig_prob_matrices(norm_exp_mtx, norm_sbs_mtx, cosmic_signatures, norm=normalize, debug=debug)
	return(ret)
}

# compute a sample-specific [n signatures x n categories] matrix X
# where x_ij = p(Signature=i | Mutation category = j)
# which is the probability that a mutation of a given category was emitted by the given signature

# input: norm_exposure_vec - a <n signatures> vector with the NORMALIZED exposure for each signature of the sample
# input: norm_mutation_vec - a <n categories> vector with the NORMALIZED number of mutations of each category in a samp
# input: cosmic_signatures - a [n signatures x n categories] data frame which contains the cosmic signatures
sample_prob_sig_matrix <- function(norm_exposure_vec, norm_mutation_vec, cosmic_signatures, norm) {
	# compute p(s|c) via bayes rule
	# p(s|c) = p(c|s) * p(s) / p(c)
	# p(c|s) = the probability that cosmic signature s emits mutation category c
	# p(s) = exposure for sig s in sample / # mutations in sample
	# p(c) = # category c muts in sample / # mutations in sample

	sig_names = cosmic_signatures[,"Signature"]
	cat_names = mut_types(dots=TRUE)

	psc_mtx = matrix(nrow=length(sig_names), ncol=length(cat_names))
	rownames(psc_mtx) = sig_names
	colnames(psc_mtx) = mut_types(dots=FALSE)

	for (s in 1:length(sig_names)) {
		for (c in 1:length(cat_names)) {
			p_c_s = cosmic_signatures[,-1][[s, c]]
			p_s = norm_exposure_vec[[s]]
			p_c = norm_mutation_vec[[c]]
			p_s_c = (p_c_s * p_s) / p_c

			psc_mtx[s, c] = p_s_c
		}
	}
	if (norm==TRUE) {
		cs = colSums(psc_mtx)
		psc_mtx = t( t(psc_mtx) / cs) # divide column-wise
	}

	return(psc_mtx)
}


# get a matrix for each sample
get_sample_sig_prob_matrices <- function(norm_exposure_mtx, norm_sbs_mtx, cosmic_signatures, norm, debug=TRUE) {
	n_sig = nrow(cosmic_signatures)
	n_cat = ncol(norm_sbs_mtx)
	n_samp = nrow(norm_sbs_mtx)

	sig_names = cosmic_signatures[, "Signature"]
	cat_names = mut_types(dots=FALSE)
	samp_names = rownames(norm_sbs_mtx)
	name_ls = list()
	name_ls[[1]] = sig_names
	name_ls[[2]] = cat_names
	name_ls[[3]] = samp_names

	arr = array(dim=c(n_sig, n_cat, n_samp), dimnames=name_ls)

	for (s in 1:n_samp) {
		if (debug) {
			if (s==1 | s %% 50 == 0) { print(paste0(s, "/", n_samp)) }
		}
		exp_vec = norm_exposure_mtx[s,]
		sbs_vec = norm_sbs_mtx[s,]
		
		if (debug) {
			check_samp_name = exp_vec[1]
			check2_samp_name = rownames(sbs_vec)
			if (check_samp_name != samp_names[s] | check2_samp_name != samp_names[s]) { 
				print("Sample name mismatch in exp_vec or sbs_vec.")
				print(paste0("s: ", s))
				print(paste0("samp_names[s]: ", samp_names[s]))
				print(paste0("exp_vec name: ", exp_vec[1]))
				print(paste0("sbs_vec name: ", rownames(sbs_vec)))
				stop()
			}
		}
		exp_vec = exp_vec[-1] # remove element with sample name

		mtx = sample_prob_sig_matrix(exp_vec, sbs_vec, cosmic_signatures, norm)
		arr[ , , s] = mtx
	}

	return(arr)
}


################ append signature probabilities to mutation dataframe #################
# (mutation dataframe is the output of load_nz_mutation_df() )

append_signature_probabilities <- function(mut_df, sig_prob_arr, debug=TRUE) {
	mtx = matrix(nrow = nrow(mut_df), ncol=nrow(sig_prob_arr[ , , 1]))
	colnames(mtx) = rownames(sig_prob_arr[ , , 1])

	for (i in 1:nrow(mut_df)) {
		if (debug) {
			if (i==1 | i == 10000 | i == 50000 | i %% 250000) {
				print(paste0(i, '/', nrow(mut_df)))
			}	
		}
		p = mut_df[i, "Patient"]
		c = mut_df[i, "SBS_96"]
		sig_prob_vec = sig_prob_arr[ , c, p]
		mtx[i, ] = sig_prob_vec
	}

	ret = cbind(mut_df, mtx)
	return(ret)
}




#############################
# save processed files      #
#############################

save_sbs_tsv <- function(outfile="nz_sbs.tsv") {
	mut_df = load_nz_mutation_df()
	sbs_df = sbs_df_from_mut_df(mut_df)
	write.table(sbs_df, file=outfile, sep="\t", quote=FALSE)
}

save_window_tsv <- function(outfile="nz_window_counts.tsv") {
	mut_df = load_nz_mutation_df()
	window_df = window_df_from_mut_df(mut_df)
	write.table(window_df, file=outfile, sep="\t", quote=FALSE)
}

# combine the signature estimates with the window mutation counts into a single saved .tsv
sig_est_with_window_tsv <- function(outfile="data/nz_sig_est_and_window_counts.tsv") {
	mut_df = load_nz_mutation_df()
	sig_df = load_nz_sig_estimates()

	window_df = window_df_from_mut_df(mut_df)
	
	ret = cbind(sig_df, window_df)
	write.table(ret, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
}

save_mut_df_with_sig_probs <- function(outfile="data/nz_mut_df_with_sig_prob.tsv", as_tsv=TRUE) {
	mut_df = load_nz_mutation_df()
	sig_prob_arr = est_mutwise_sig_probabilities(normalize=TRUE)

	ret = append_signature_probabilities(mut_df, sig_prob_arr)
	
	if (as_tsv==FALSE & outfile=="data/nz_mut_df_with_sig_prob.tsv") {
		outfile="data/nz_mut_df_with_sig_prob.rds"
	}
	if (as_tsv==FALSE) {
		saveRDS(ret, file=outfile)
	} else {
		write.table(ret, file=outfile, sep="\t", quote=FALSE, row.names=FALSE)
	}
}

save_aggregated_window_sigprob <- function(outfile="data/nz_aggregated_window_sigprob.tsv") {
	mut_sigprob_df = load_nz_mut_df_with_sigprob()
	ret = aggregate_window_sigprob_mtx(mut_sigprob_df)

	write.table(ret, file=outfile, sep="\t", quote=FALSE)
}

par_save_aggregated_window_sigprob <- function(outfile="", window_size=10^6, mut_windows=FALSE) {
	if (nchar(outfile)==0) {
		if (mut_windows==FALSE) {
			outfile=paste0("data/nz_aggregated_window_", window_size, "_sigprob.tsv")
		} else {
			outfile=paste0("data/nz_mutwindow_", window_size, "_sigprob.tsv")
		}
		print(paste0("running with default outfile: ", outfile))
	}

	mut_sigprob_df = load_nz_mut_df_with_sigprob()
	print("starting par_aggregate_window_sigprob_mtx")
	ret = par_aggregate_window_sigprob_mtx(mut_sigprob_df, win_size = window_size, mut_windows = mut_windows)

	print("done. writing table.")
	write.table(ret, file=outfile, sep="\t", quote=FALSE)
}

par_save_individual_chr_sigprob <- function(outprefix="", window_size=10000, mut_windows=TRUE) {
	if (nchar(outprefix)==0) {
		if (mut_windows==FALSE) {
			outprefix=paste0("data/individual_chromosome_matrices/nz_window_", window_size, "_sigprob_chr_")
		} else {
			outprefix=paste0("data/individual_chromosome_matrices/nz_mutwindow_", window_size,  "_sigprob_chr_")
		}
		print(paste0("running with default outprefix: ", outprefix))
	}

	chrom_vec = c(1:22, "X", "Y")	
	
	mut_sigprob_df = load_nz_mut_df_with_sigprob()
	print("starting par_chr_mtx_sigprob_ls")
	ret = par_chr_mtx_sigprob_ls(mut_sigprob_df, win_size=window_size, mut_windows=mut_windows)
	
	#ret = list()

	#for (i in 1:length(chrom_vec)) {
	#	mtx = matrix(nrow=10, ncol=5)
#		mtx[,] = chrom_vec[i]
#		ret[[i]] = mtx
#	}
	
	print("done. writing tables.")


	for (i in 1:length(ret)) {
		outfile = paste0(outprefix, chrom_vec[i], ".tsv")
		print(paste0("writing to outfile: ", outfile))
		write.table(ret[[i]], file=outfile, sep="\t", quote=FALSE)
	}
	print("done.")
}
