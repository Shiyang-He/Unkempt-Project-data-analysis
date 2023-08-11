argv=commandArgs(TRUE)
library(tailfindr)
#library(rbokeh)

df <- find_tails(fast5_dir = argv[2], 
	save_dir = argv[1],
	csv_filename = 'rna_tails.csv',
	num_cores = 8)
#	save_plots = TRUE,
#	plot_debug_traces = TRUE,
#	plotting_library = 'rbokeh')
