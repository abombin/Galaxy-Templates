library(ggplot2)
library(jsonlite)


extract_plots = function(json_file) {
  params_json = read_json(json_file)
  object_list = readRDS(params_json$input_rds)
  plots_list = object_list[["plots"]]
  pdf("plots.pdf")   # Open a new PDF device
  for (p in plots_list) {
    print(p)            # Print each plot (forces new page automatically)
  }
  dev.off()  
}

# run tool
if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  extract_plots(args[1])
}