library(ggplot2)
library(Seurat)

convert_to_a3 = function(object_list, output_rds) {

  object_list = readRDS(object_list)
  so_object_list = object_list[["object"]]
  plots_list = object_list[["plots"]]
  
  # convert to A3
  edit_so_list = list()
  for ( i in 1:length(so_object_list)) {
    cur_name = names(so_object_list)[i]
    cur_so = so_object_list[[i]]
    cur_so[["RNA"]] = as(object = cur_so[["RNA"]], Class = "Assay")
    edit_so_list[[cur_name]] = cur_so
  }
  
  saveRDS(edit_so_list, file = output_rds)
  
  #save plots
  
  pdf("plots.pdf")   # Open a new PDF device
  for (p in plots_list) {
    print(p)            # Print each plot (forces new page automatically)
  }
  dev.off()  
  
}

convert_to_a3("processRawData_result.rds", "processRawData_result.rds")