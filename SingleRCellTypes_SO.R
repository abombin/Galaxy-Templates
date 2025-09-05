#!/usr/bin/env Rscript

SingleRCellTypes_SO <- function(CombinedSO,RDS, Metadata, species="Human", useClusters=FALSE, 
                                 ClusterColumn=NULL, Legend_Dot_Size=2, 
                                 doFineTuning=FALSE, Number_of_cells_per_partition=400) {

## --------- ##
## Libraries ##
## --------- ##
    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    library(grid)
    library(gridExtra)
    library(gridBase)
    library(cowplot)
    

## -------------------------------- ##
## User-Defined Template Parameters ##
## -------------------------------- ##



# Visualization Parameters:
Reduction <- "umap"


## -------------------------------- ##
## Errors                           ##
## -------------------------------- ##

## -------------------------------- ##
## Functions                        ##
## -------------------------------- ##

 
## --------------- ##
## Main Code Block ##
## --------------- ##

# Loading Seurat object
    cat("1. Reading Seurat Object from dataset: seurat_object.rds\n\n")
# auto removed:     fs <- CombinedSO$fileSystem()
# auto removed:     path <- fs$get_path("seurat_object.rds", 'r')
    so <- CombinedSO

# Loading CellDex annotation
#cell.dex <- list(HPCA, BP, mousernaseq, immgen)
    cat("2. Reading pre-saved CellDex Data from dataset: CellDexDatabase.rds\n\n")
# auto removed:     fs2 <- RDS$fileSystem()
# auto removed:     path2 <- fs2$get_path("CellDexDatabase.rds", 'r')
    CellDexData <- RDS

Reduction = "umap"
if (useClusters)
  {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData,
                            use.clusters = ClusterColumn)
  } else {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData)

  }
    gc()

## Print figures
    print(anno$p1)
    print(anno$p2)

Reduction = "tsne"
if (useClusters)
  {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData,
                            use.clusters = ClusterColumn)
  } else {
      anno <- annotateCellTypes(object = so,
                            species = species,
                            reduction.type = Reduction,
                            legend.dot.size = Legend_Dot_Size,
                            do.finetuning = doFineTuning,
                            local.celldex = CellDexData)

  }
    gc()

## Print figures
    print(anno$p1)
    print(anno$p2)

## Save the annotated Seurat object
# auto removed:     output <- new.output()
# auto removed:     output_fs <- output$fileSystem()
return(anno$object)

# auto removed:     return(NULL)
}

# Command line interface
# Check if script is being run from command line
if (!interactive()) {
  # Set up logging
  log_file <- "tool_stdout.log"
  sink(log_file, append = TRUE, split = TRUE)
  
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 1) {
    stop("Usage: Rscript SingleRCellTypes_SO.R <config.json>")
  }
  
  json_file <- args[1]
  
  # Check if JSON file exists
  if (!file.exists(json_file)) {
    stop(paste("JSON file not found:", json_file))
  }
  
  # Load required libraries
  library(jsonlite)
  
  # Read parameters from JSON file
  cat("Reading parameters from:", json_file, "\n")
  params <- fromJSON(json_file)
  
  # Validate required parameters
  required_params <- c("CombinedSO", "RDS", "Metadata")
  missing_params <- setdiff(required_params, names(params))
  if (length(missing_params) > 0) {
    stop(paste("Missing required parameters:", paste(missing_params, collapse = ", ")))
  }
  
  # Set default values for optional parameters
  if (is.null(params$species)) params$species <- "Human"
  if (is.null(params$useClusters)) params$useClusters <- FALSE
  if (is.null(params$ClusterColumn)) params$ClusterColumn <- NULL
  if (is.null(params$Legend_Dot_Size)) params$Legend_Dot_Size <- 2
  if (is.null(params$doFineTuning)) params$doFineTuning <- FALSE
  if (is.null(params$Number_of_cells_per_partition)) params$Number_of_cells_per_partition <- 400
  
  # Load the data files specified in the JSON
  cat("Loading CombinedSO from:", params$CombinedSO, "\n")
  cat("Available memory before loading:", round(as.numeric(system("free -m | grep '^Mem:' | awk '{print $7}'", intern=TRUE)), 2), "MB\n")
  CombinedSO <- readRDS(params$CombinedSO)
  gc() # Force garbage collection
  cat("Available memory after loading CombinedSO:", round(as.numeric(system("free -m | grep '^Mem:' | awk '{print $7}'", intern=TRUE)), 2), "MB\n")
  
  cat("Loading RDS from:", params$RDS, "\n")
  RDS <- readRDS(params$RDS)
  gc() # Force garbage collection
  cat("Available memory after loading RDS:", round(as.numeric(system("free -m | grep '^Mem:' | awk '{print $7}'", intern=TRUE)), 2), "MB\n")
  
  cat("Loading Metadata from:", params$Metadata, "\n")
  Metadata <- read.csv(params$Metadata)
  
  # Call the main function with parameters from JSON
  result <- SingleRCellTypes_SO(
    CombinedSO = CombinedSO,
    RDS = RDS,
    Metadata = Metadata,
    species = params$species,
    useClusters = params$useClusters,
    ClusterColumn = params$ClusterColumn,
    Legend_Dot_Size = params$Legend_Dot_Size,
    doFineTuning = params$doFineTuning,
    Number_of_cells_per_partition = params$Number_of_cells_per_partition
  )
  
  # Save output if specified in JSON
  if (!is.null(params$output_file)) {
    cat("Saving result to:", params$output_file, "\n")
    saveRDS(result, params$output_file)
  }
  
  cat("Analysis completed successfully!\n")
  
  # Close logging
  sink()
}

