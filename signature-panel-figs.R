library(ggplot2)
library(dplyr)
# script to generate figures for signature panel project

PATH_TO_RESULTS_FILES = "C:/Users/nickf/Desktop/ubuntu_stuff/"
SIZE_TEST_SUMMARY_FILE = "C:/Users/nickf/Desktop/ubuntu_stuff/size_test_results_medians_26-Jun-2020.tsv"

# load results df
load_panel_size_result_summary <- function() {
  return(read.csv(SIZE_TEST_SUMMARY_FILE, sep="\t"))
}

process_summary_df_for_ggplot <- function(df) {
  o1_df = df[ , c("Signature", "Panel.Size", "Obj1.Score")]
  o2_df = df[ , c("Signature", "Panel.Size", "Obj2.Score")]
  o3_df = df[ , c("Signature", "Panel.Size", "Obj3.Score")]
  
  colnames(o1_df)[3] = "AUPR"
  colnames(o2_df)[3] = "AUPR"
  colnames(o3_df)[3] = "AUPR"
  
  o1_df$Panel.Size = o1_df$Panel.Size * 10
  o2_df$Panel.Size = o2_df$Panel.Size * 10
  o3_df$Panel.Size = o3_df$Panel.Size * 10
  
  Objective.Function = rep(1, nrow(o1_df))
  o1_df = cbind(o1_df, Objective.Function)
  
  Objective.Function = rep(2, nrow(o2_df))
  o2_df = cbind(o2_df, Objective.Function)

  Objective.Function = rep(3, nrow(o3_df))
  o3_df = cbind(o3_df, Objective.Function)
  
  processed_df = rbind(o1_df, o2_df, o3_df)
  processed_df$Objective.Function = as.factor(processed_df$Objective.Function) 
  return(processed_df)
}

load_ps_summary_df_for_ggplot <- function() {
  return(process_summary_df_for_ggplot(load_panel_size_result_summary()))
}


get_msk_benchmark <- function(sig_num) {
  sum_df = load_panel_size_result_summary()
  msk = sum_df[sum_df$Signature==sig_num, "MSK.IMPACT.Score"][1]
  return(msk)
}

get_wes_benchmark <- function(sig_num) {
  sum_df = load_panel_size_result_summary()
  wes = sum_df[sum_df$Signature==sig_num, "WES.Score"][1]
  return(wes)
}

make_ps_line_plot <- function(sig_num) {
  gg_df = load_ps_summary_df_for_ggplot()
  sig_df = gg_df[gg_df$Signature==sig_num, ]
  
  msk_bench = get_msk_benchmark(sig_num)
  wes_bench = get_wes_benchmark(sig_num)
  
  p = sig_df %>% ggplot( aes(x=Panel.Size, y=AUPR, group=Objective.Function, color=Objective.Function)) + geom_line() +
                  ggtitle(paste0("Signature ", sig_num, " Panel Performance")) +
                  xlab("Panel Size (Kb)") +
                  scale_x_continuous(breaks=c(100, 500, 1000, 1500, 2000, 2500)) +
                  geom_hline( aes(yintercept=wes_bench, linetype="Whole Exome (~30,000 Kb)"), color="magenta") +
                  geom_hline( aes(yintercept=msk_bench, linetype="MSK-IMPACT (~2,500 Kb)"), color="cyan") +
                  scale_linetype_manual(name="Benchmarks", values=c(2,2), guide=guide_legend(override.aes = list(color= c("cyan", "magenta"))))
  
  return(p)
}

save_all_ps_line_plots <- function(file_prefix="panel_size_line_plot", out_dir = NULL, file_format="pdf") {
  df = load_panel_size_result_summary()
  sigs = unique(df$Signature)
  
  for (s in sigs) {
    p = make_ps_line_plot(s)
    
    fn = paste0(file_prefix, "_sig", s, ".", file_format)
    
    ggsave(filename=fn, plot=p, path=out_dir, width = 7.5, height=5)
  }
}