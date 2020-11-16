setwd('..')
source('baseline_panels.R')
source('projection_score.R')

file_tag = "bp_test"

NUM_PANELS = 1000
NUM_WIN_IN_PANEL = 250

print(paste0(Sys.time(), "    finding active 10kb windows."))
possible_windows = active_10k_windows()
print(paste0(Sys.time(), "    done."))

outprefix = paste0(GLOBAL_SCRIPT_BASELINE_WINDOWS, "baseline_panel_windows_", file_tag)

print(paste0(Sys.time(), "    sampling ", NUM_PANELS, " random panels with ", NUM_WIN_IN_PANEL, " windows each."))
for (i in 1:NUM_PANELS) {
	curr_outfile = paste0(outprefix, "_it", i, "_nwin", NUM_WIN_IN_PANEL, ".txt")
	selected_windows = get_random_windows(NUM_WIN_IN_PANEL, possible_windows, check_overlap=FALSE)
	write(selected_windows, curr_outfile)
}
print(paste0(Sys.time(), "    done."))
