# Color by Metadata [CCBR] [scRNA-seq] (0fac593c-01e9-4c02-96d0-cc4876c4dfef): v163
ColorByMetadata_CellTypes <- function(seurat_object,
                                     samples_to_include = "E1,E2,E3,H1,H2,H3,T1,T2,T3",
                                     metadata_to_plot = "HPCA_main,BP_encode_main",
                                     number_of_columns_for_final_image = 0,
                                     show_labels = FALSE,
                                     dot_size = 0.01,
                                     legend_text_size = 1,
                                     legend_position = "right",
                                     columns_to_summarize = c(),
                                     summarization_cut_off = 5,
                                     save_the_entire_dataset = FALSE,
                                     use_cite_seq = FALSE) {
    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    library(grid)
    library(gridExtra)

    ## -------------------------------- ##
    ## Process Parameters               ##
    ## -------------------------------- ##
    
    # Convert comma-separated strings to vectors
    samples_to_include = strsplit(samples_to_include, ",")[[1]]
    metadata_to_plot = strsplit(metadata_to_plot, ",")[[1]]
    print(samples_to_include)
    print(metadata_to_plot)

    samples_to_include = paste0('c("', paste(samples_to_include, collapse = '","'), '")')
    metadata_to_plot = paste0('c("', paste(metadata_to_plot, collapse = '","'), '")')

    if (columns_to_summarize == "None") {
        columns_to_summarize = NULL
    } else {
        columns_to_summarize = strsplit(columns_to_summarize, ",")[[1]]
        columns_to_summarize = paste0('c("', paste(columns_to_summarize, collapse = '","'), '")')
    }

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    # Note: All parameters are now function arguments
    # seurat_object is passed as parameter



    ## --------------- ##
    ## Error Messages ##
    ## -------------- ##


    ## --------- ##
    ## Functions ##
    ## --------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    ## Load SO
    cat("Reading Seurat Object\n\n")
    SO <- seurat_object
    print(SO)

    All_Reduction_Types <- c("umap", "tsne", "pca")
    for (Reduction_Type in All_Reduction_Types) {
        results <- plotMetadata(
            object = SO,
            samples.to.include = samples_to_include,
            metadata.to.plot = metadata_to_plot,
            columns.to.summarize = columns_to_summarize,
            summarization.cut.off = summarization_cut_off,
            reduction.type = Reduction_Type,
            use.cite.seq = use_cite_seq,
            show.labels = show_labels,
            legend.text.size = legend_text_size,
            legend.position = legend_position,
            dot.size = dot_size
        )

        ## Print Graphic output

        # Calculate number of metadata columns for layout
        n_metadata_cols <- length(metadata_to_plot)
        if (number_of_columns_for_final_image == 0) {
            n <- ceiling(n_metadata_cols^0.5)
        } else {
            n <- number_of_columns_for_final_image
        }
        do.call("grid.arrange", c(results$plot, ncol = n))
    }

    ## Save dataset if requested

    if (save_the_entire_dataset) {
        # auto removed:     output <- new.output()
        # auto removed:     output_fs <- output$fileSystem()
        return(results$object)
        # auto removed:     return(NULL)
    } else {
        return(results$object@meta.data)
    }
}

# Command Line Interface
# Check if script is being run from command line
if (!interactive()) {
  # Load required library
  if (!require("jsonlite", quietly = TRUE)) {
    install.packages("jsonlite")
    library("jsonlite")
  }
  
  # Get command line arguments
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) != 1) {
    cat("Usage: Rscript ColorByMetadata_CellTypes.R <config.json>\n")
    cat("Example: Rscript ColorByMetadata_CellTypes.R color_by_metadata_config.json\n")
    quit(status = 1)
  }
  
  config_file <- args[1]
  
  # Check if config file exists
  if (!file.exists(config_file)) {
    cat("Error: Configuration file", config_file, "not found\n")
    quit(status = 1)
  }
  
  # Load configuration
  tryCatch({
    params <- fromJSON(config_file)
    cat("Loading configuration from:", config_file, "\n")
  }, error = function(e) {
    cat("Error reading JSON file:", e$message, "\n")
    quit(status = 1)
  })
  
  # Set default values for missing parameters
  if (is.null(params$samples_to_include)) params$samples_to_include <- "E1,E2,E3,H1,H2,H3,T1,T2,T3"
  if (is.null(params$metadata_to_plot)) params$metadata_to_plot <- "HPCA_main,BP_encode_main"
  if (is.null(params$number_of_columns_for_final_image)) params$number_of_columns_for_final_image <- 0
  if (is.null(params$show_labels)) params$show_labels <- FALSE
  if (is.null(params$dot_size)) params$dot_size <- 0.01
  if (is.null(params$legend_text_size)) params$legend_text_size <- 1
  if (is.null(params$legend_position)) params$legend_position <- "right"
  if (is.null(params$columns_to_summarize)) params$columns_to_summarize <- c()
  if (is.null(params$summarization_cut_off)) params$summarization_cut_off <- 5
  if (is.null(params$save_the_entire_dataset)) params$save_the_entire_dataset <- FALSE
  if (is.null(params$use_cite_seq)) params$use_cite_seq <- FALSE
  
  # Load input data
  tryCatch({
    cat("Loading Seurat Object from:", params$seurat_object, "\n")
    cat("Available memory before loading:", round(as.numeric(system("free -m | grep '^Mem:' | awk '{print $7}'", intern=TRUE)), 2), "MB\n")
    seurat_obj <- readRDS(params$seurat_object)
    cat("Available memory after loading:", round(as.numeric(system("free -m | grep '^Mem:' | awk '{print $7}'", intern=TRUE)), 2), "MB\n")
    cat("Seurat object loaded successfully\n")
    print(seurat_obj)
  }, error = function(e) {
    cat("Error loading Seurat object:", e$message, "\n")
    quit(status = 1)
  })
  
  # Run the main function
  cat("Running ColorByMetadata_CellTypes analysis...\n")
  result <- ColorByMetadata_CellTypes(
    seurat_object = seurat_obj,
    samples_to_include = params$samples_to_include,
    metadata_to_plot = params$metadata_to_plot,
    number_of_columns_for_final_image = params$number_of_columns_for_final_image,
    show_labels = params$show_labels,
    dot_size = params$dot_size,
    legend_text_size = params$legend_text_size,
    legend_position = params$legend_position,
    columns_to_summarize = params$columns_to_summarize,
    summarization_cut_off = params$summarization_cut_off,
    save_the_entire_dataset = params$save_the_entire_dataset,
    use_cite_seq = params$use_cite_seq
  )
  
  # Save output
  output_file <- if (!is.null(params$output_file)) params$output_file else "color_by_metadata_result.rds"
  saveRDS(result, output_file)
  cat("Results saved to:", output_file, "\n")
  
  # Memory cleanup
  cat("Available memory after processing:", round(as.numeric(system("free -m | grep '^Mem:' | awk '{print $7}'", intern=TRUE)), 2), "MB\n")
}
