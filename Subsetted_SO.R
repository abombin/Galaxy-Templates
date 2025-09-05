# Subset Seurat Object [CCBR] [scRNA-seq] (ec3f23f9-bcba-4f3a-8a08-8ba611fbb6c7): v149
Subsetted_SO <- function(seurat_object,
                        category_to_filter = "HPCA_main",
                        seed = 10,
                        samples_to_include = "E1,E2,E3,H1,H2,H3,T1,T2,T3",
                        sample_name = "orig_ident", 
                        keep_or_remove = TRUE,
                        values_to_filter = c("T_cells"),
                        greater_less_than = "greater than",
                        cut_off = 0.5,
                        plot_as_interactive_plot = FALSE,
                        use_cite_seq_data = FALSE,
                        legend_position = "right",
                        legend_symbol_size = 2,
                        number_of_legend_columns = 1,
                        dot_size = 0.1,
                        dot_size_highlighted_cells = 0.5,
                        colors = c("aquamarine3", "salmon1", "lightskyblue3", "plum3", "darkolivegreen3", "goldenrod1", "burlywood2", "gray70", "firebrick2", "steelblue", "palegreen4", "orchid4", "darkorange1", "yellow", "sienna", "palevioletred1", "gray60", "cyan4", "darkorange3", "mediumpurple3", "violetred2", "olivedrab", "darkgoldenrod2", "darkgoldenrod", "gray40", "palegreen3", "thistle3", "khaki1", "deeppink2", "chocolate3", "paleturquoise3", "wheat1", "lightsteelblue", "salmon", "sandybrown", "darkolivegreen2", "thistle2", "gray85", "orchid3", "darkseagreen1", "lightgoldenrod1", "lightskyblue2", "dodgerblue3", "darkseagreen3", "forestgreen", "lightpink2", "mediumpurple4", "lightpink1", "thistle", "navajowhite", "lemonchiffon", "bisque2", "mistyrose", "gray95", "lightcyan3", "peachpuff2", "lightsteelblue2", "lightyellow2", "moccasin", "gray80", "antiquewhite2", "lightgrey")) {
    ## --------- ##
    ## Libraries ##
    ## --------- ##

    library(nidapFunctions)
    nidapLoadPackages("SCWorkflow")
    library(cowplot)
    library(grid)
    library(gridExtra)
    library(scales)

    ### GlobalCodeFix
    library("Seurat")
    library("ggplot2")
    library("gridExtra")
    library("grid")
    library("gridBase")
    library("cowplot")
    library("RColorBrewer")
    library("colorspace")
    library("tidyverse")

    library(plotly)

    ## -------------------------------- ##
    ## User-Defined Template Parameters ##
    ## -------------------------------- ##

    # Note: seurat_object is passed as parameter
    # All other parameters are now function arguments

    ## --------------- ##
    ## Error Messages ##
    ## -------------- ##


    ## --------- ##
    ## Functions ##
    ## --------- ##

    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    ## Load SO.
    cat("Reading Seurat Object\n\n")
    so <- seurat_object
    print(so)

   samples = strsplit(samples_to_include, ",")[[1]]
   print(samples)

   if (length(colors) == 1) {
       colors = strsplit(colors, ",")[[1]]
   }

    if ("orig.ident" %in% colnames(so@meta.data)) {
        Idents(so) <- so$orig.ident
    } else if ("orig_ident" %in% colnames(so@meta.data)) {
        Idents(so) <- so$orig_ident
    } else {
        print("orig.ident or orig_ident not found")
    }

    so <- subset(so, idents = samples)
    ## HOTFIX ENDS

    print("subsetted so")

    reduction <- "umap"
    results_umap <- filterSeuratObjectByMetadata(
        object = so,
        samples.to.include = samples,
        sample.name = sample_name,
        category.to.filter = category_to_filter,
        values.to.filter = values_to_filter,
        keep.or.remove = keep_or_remove,
        greater.less.than = greater_less_than,
        seed = seed,
        cut.off = cut_off,
        legend.position = legend_position,
        reduction = reduction,
        plot.as.interactive.plot = plot_as_interactive_plot,
        legend.symbol.size = legend_symbol_size,
        colors = colors,
        dot.size = dot_size,
        number.of.legend.columns = number_of_legend_columns,
        dot.size.highlighted.cells = dot_size_highlighted_cells,
        use.cite.seq.data = use_cite_seq_data
    )

    reduction <- "tsne"
    results_tsne <- filterSeuratObjectByMetadata(
        object = so,
        samples.to.include = samples,
        sample.name = sample_name,
        category.to.filter = category_to_filter,
        values.to.filter = values_to_filter,
        keep.or.remove = keep_or_remove,
        greater.less.than = greater_less_than,
        seed = seed,
        cut.off = cut_off,
        legend.position = legend_position,
        reduction = reduction,
        plot.as.interactive.plot = plot_as_interactive_plot,
        legend.symbol.size = legend_symbol_size,
        colors = colors,
        dot.size = dot_size,
        number.of.legend.columns = number_of_legend_columns,
        dot.size.highlighted.cells = dot_size_highlighted_cells,
        use.cite.seq.data = use_cite_seq_data
    )
    ## If interactive plot requested, then ...
    if (plot_as_interactive_plot == TRUE) {
        gp1 <- ggplotly(results$plot1)
        gp2 <- ggplotly(results$plot2)
        p <- subplot(gp1, gp2, nrows = 2)
        print(p)
    } else {
        ## Else, print non-interactive plot.
        print(plot_grid(results_umap$plot1, results_umap$plot2, nrow = 1))
        print(plot_grid(results_tsne$plot1, results_tsne$plot2, nrow = 1))
    }

    ## Return the subsetted Seurat object.
    # auto removed:     output <- new.output()
    # auto removed:     output_fs <- output$fileSystem()
    return(results_umap$object)

    # auto removed:     return(NULL)
}

filterSeuratObjectByMetadata <- function(object,
                                         samples.to.include,
                                         sample.name,
                                         category.to.filter,
                                         values.to.filter,
                                         keep.or.remove = TRUE,
                                         greater.less.than = "greater than",
                                         seed = 10,
                                         cut.off = 0.5,
                                         legend.position = "top",
                                         reduction = "umap",
                                         plot.as.interactive.plot = FALSE,
                                         legend.symbol.size = 2,
                                         colors = c(
                                             "aquamarine3",
                                             "salmon1",
                                             "lightskyblue3",
                                             "plum3",
                                             "darkolivegreen3",
                                             "goldenrod1",
                                             "burlywood2",
                                             "gray70",
                                             "firebrick2",
                                             "steelblue",
                                             "palegreen4",
                                             "orchid4",
                                             "darkorange1",
                                             "yellow",
                                             "sienna",
                                             "palevioletred1",
                                             "gray60",
                                             "cyan4",
                                             "darkorange3",
                                             "mediumpurple3",
                                             "violetred2",
                                             "olivedrab",
                                             "darkgoldenrod2",
                                             "darkgoldenrod",
                                             "gray40",
                                             "palegreen3",
                                             "thistle3",
                                             "khaki1",
                                             "deeppink2",
                                             "chocolate3",
                                             "paleturquoise3",
                                             "wheat1",
                                             "lightsteelblue",
                                             "salmon",
                                             "sandybrown",
                                             "darkolivegreen2",
                                             "thistle2",
                                             "gray85",
                                             "orchid3",
                                             "darkseagreen1",
                                             "lightgoldenrod1",
                                             "lightskyblue2",
                                             "dodgerblue3",
                                             "darkseagreen3",
                                             "forestgreen",
                                             "lightpink2",
                                             "mediumpurple4",
                                             "lightpink1",
                                             "thistle",
                                             "navajowhite",
                                             "lemonchiffon",
                                             "bisque2",
                                             "mistyrose",
                                             "gray95",
                                             "lightcyan3",
                                             "peachpuff2",
                                             "lightsteelblue2",
                                             "lightyellow2",
                                             "moccasin",
                                             "gray80",
                                             "antiquewhite2",
                                             "lightgrey"
                                         ),
                                         dot.size = 0.1,
                                         number.of.legend.columns = 1,
                                         dot.size.highlighted.cells = 0.5,
                                         use.cite.seq.data = FALSE) {
    ## --------------- ##
    ## Parameters      ##
    ## --------------- ##





    ## --------------- ##
    ## Functions       ##
    ## --------------- ##

    # Drawing TSNE/UMAP/PCA plot
    .drawtsne <- function(SO, reduction, scale.col, col.grad) {
        SO.clus <- SO@meta.data[[category.to.filter]]

        plot1 <- DimPlot(SO, reduction = reduction, group.by = "ident")
        class(plot1$data$ident) <- "numeric"

        if (reduction == "tsne") {
            clusmat <- data.frame(
                umap1 = plot1$data$tSNE_1,
                umap2 = plot1$data$tSNE_2,
                clusid = as.numeric(SO@meta.data[[category.to.filter]])
            )
        } else if (reduction == "umap") {
            clusmat <- data.frame(
                umap1 = plot1$data$UMAP_1,
                umap2 = plot1$data$UMAP_2,
                clusid = as.numeric(SO@meta.data[[category.to.filter]])
            )
        } else if (reduction == "pca") {
            clusmat <- data.frame(
                umap1 = plot1$data$PC_1,
                umap2 = plot1$data$PC_2,
                clusid = as.numeric(SO@meta.data[[category.to.filter]])
            )
        } else if (reduction == "protein_tsne") {
            clusmat <- data.frame(
                umap1 = plot1$data$protein_tsne_1,
                umap2 = plot1$data$protein_tsne_2,
                clusid = as.numeric(SO@meta.data[[category.to.filter]])
            )
        } else if (reduction == "protein_umap") {
            clusmat <- data.frame(
                umap1 = plot1$data$protein_umap_1,
                umap2 = plot1$data$protein_umap_2,
                clusid = as.numeric(SO@meta.data[[category.to.filter]])
            )
        } else {
            clusmat <- data.frame(
                umap1 = plot1$data$protein_pca_1,
                umap2 = plot1$data$protein_pca_2,
                clusid = as.numeric(SO@meta.data[[category.to.filter]])
            )
        }

        # Preparing
        clusmat %>%
            group_by(clusid) %>%
            summarise(
                umap1.mean = mean(umap1),
                umap2.mean = mean(umap2)
            ) -> umap.pos
        title <- as.character(category.to.filter)
        clusmat %>% dplyr::arrange(clusid) -> clusmat

        plot2 <- ggplot(clusmat, aes(x = umap1, y = umap2)) +
            theme_bw() +
            theme(legend.title = element_blank()) +
            geom_point(
                aes(colour = clusid),
                alpha = 0.5,
                shape = 20,
                size = dot.size
            ) +
            scale_color_gradientn(
                colors = brewer.pal(n = 5, name = col.grad),
                values = scale.col
            ) +
            guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()
            ) +
            ggtitle(title) +
            xlab("umap-1") +
            ylab("umap-2")
        return(plot2)
    }

    .distinctColorPalette <- function(k = 1, seed) {
        current.color.space <- our.color.space@coords
        # Set iter.max to 20 to avoid convergence warnings.
        set.seed(seed)
        km <- kmeans(current.color.space, k, iter.max = 20)
        colors <- unname(hex(LAB(km$centers)))
        return(colors)
    }
    ## --------------- ##
    ## Main Code Block ##
    ## --------------- ##

    # Checking if samples are selected
   #samples = strsplit(samples_to_include, ",")[[1]]

    # if (length(samples) == 0) {
    #     samples <- unique(object@meta.data[[sample.name[1]]])
    # }

    ## Replace dots in metadata column names with underscores.
    colnames(object@meta.data) <- gsub("\\.", "_", colnames(object@meta.data))
    new.sample.name <- gsub("\\.", "_", sample.name[1])

    ## If you have protien data, then ...
    if (use.cite.seq.data) {
        reduction <- paste("protein_", reduction, sep = "")
    }


    ## Original color-picking code.
    rand.col <- 2e3
    set.seed(seed)
    our.color.space <- colorspace::RGB(runif(rand.col), runif(rand.col), runif(rand.col))
    our.color.space <- as(our.color.space, "LAB")


    ## User-selected metadata column is used to set idents.
    Filter.orig <- object@meta.data[[category.to.filter[1]]]
    colname <- category.to.filter[1]

    ident.of.interest <- as.factor(object@meta.data[[colname]])
    names(ident.of.interest) <- names(object@active.ident)
    object@active.ident <- as.factor(vector())
    object@active.ident <- ident.of.interest

    ## Get colors from user parameter and add more if the default list is too short.
    if (class(object@meta.data[[category.to.filter[1]]]) != "numeric") {
        col.length <- length(levels(as.factor(Filter.orig)))
        if (length(colors) < col.length) {
            more.cols <- .distinctColorPalette(col.length - length(colors), 10)
            colors <- c(colors, more.cols)
        }
        names(colors) <- levels(as.factor(Filter.orig))

        ## Keep or remove cells based on user input values.
        if (keep.or.remove) {
            subset.value <- values.to.filter
            meta.col <- unique(object@meta.data[[category.to.filter[1]]])
            print("Missing values:")
            print(setdiff(subset.value, meta.col))
            subset.value <- intersect(meta.col, subset.value)
        } else {
            meta.col <- unique(object@meta.data[[colname]])
            vals.to.remove <- values.to.filter
            subset.value <- setdiff(meta.col, vals.to.remove)
        }

        ## Subset Seurat object.
        SO.sub <- subset(object, idents = subset.value)

        ## Log output of tables of cell types by samples before and after filtes.
        print("Breakdown of filtered data:")
        print(table(
            object@meta.data[[category.to.filter[1]]],
            object@meta.data[[new.sample.name]]
        ))

        cat("\n")
        print("After Filtering:")
        print(table(
            SO.sub@meta.data[[category.to.filter[1]]],
            SO.sub@meta.data[[new.sample.name]]
        ))

        ## Set filter for the subsetted SO.
        SO.sub@meta.data[[colname]] <-
            as.factor(as.character(SO.sub@meta.data[[colname]])) # Relevel Factors

        filter.sub <- SO.sub@meta.data[[colname]]

        # Set colors for unfiltered and filtered data by sample name:
        filt.length <- length(levels(as.factor(filter.sub)))
        idx <- vector("list", filt.length)
        names(idx) <- levels(as.factor(filter.sub))
        for (i in 1:filt.length) {
            id <- Filter.orig %in% levels(as.factor(filter.sub))[i]
            idx[[i]] <- rownames(object@meta.data)[id]
        }
        cols2 <- colors[levels(as.factor(filter.sub))]

        ## Make before and after plots.
        title <-
            paste0(
                "filtered by ",
                category.to.filter[1] # ,
                ##             " and split by ",
                ##             category.to.filter[2]
            )
        plot1 <- DimPlot(
            object,
            reduction = reduction,
            group.by = colname,
            pt.size = dot.size,
            raster = FALSE
        ) +
            theme_classic() +
            scale_color_manual(values = colors) +
            theme(legend.position = legend.position) +
            guides(colour = guide_legend(
                ncol = number.of.legend.columns,
                override.aes = list(size = legend.symbol.size)
            )) +
            ggtitle(colname)
        plot2 <- DimPlot(
            object,
            reduction = reduction,
            cells.highlight = idx,
            cols.highlight = rev(cols2[1:filt.length]),
            sizes.highlight = dot.size.highlighted.cells,
            raster = FALSE
        ) +
            theme_classic() +
            theme(legend.position = legend.position) +
            guides(colour = guide_legend(
                ncol = number.of.legend.columns,
                reverse = TRUE,
                override.aes = list(size = legend.symbol.size)
            )) +
            ggtitle(title)

        ## Else, filter on numeric data with a user defined threshold and direction.
    } else {
        filter.direction <- greater.less.than
        meta.col <- unique(object@meta.data[[category.to.filter]])
        value <- cut.off
        if (filter.direction == "greater than") {
            SO.sub <- subset(object, subset = category.to.filter > cut.off)
        } else {
            SO.sub <- subset(object, subset = category.to.filter < cut.off)
        }


        clusid <- object@meta.data[[category.to.filter]]
        maxclus <- max(clusid)
        clusmid <- 0.01 / maxclus
        min <- min(clusid)
        midpt.1 <- 0.99 * value
        midpt.0 <- value
        midpt.2 <- 1.01 * value
        max <- max(clusid)
        col.points <- c(min, midpt.1, midpt.0, midpt.2, max)
        col.points <- scales::rescale(col.points, c(0, 1))



        plot1 <- .drawtsne(object, reduction, col.points, "RdBu")

        clusid <- scales::rescale(SO.sub@meta.data[[category.to.filter]], to = c(0, 1))
        clus.quant <- quantile(clusid[clusid > 0], probs = c(0, .25, .5, .75, 1))
        min <- clus.quant[1]
        midpt.1 <- clus.quant[3]
        midpt.3 <- clus.quant[2]
        midpt.4 <- clus.quant[4]
        max <- clus.quant[5]
        col.points.2 <- c(min, midpt.3, midpt.1, midpt.4, max)

        plot2 <- .drawtsne(SO.sub, reduction, col.points.2, "Blues")
    }



    result.list <- list(
        "object" = SO.sub,
        "plot1" = plot1,
        "plot2" = plot2
    )
    return(result.list)
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
    cat("Usage: Rscript Subsetted_SO.R <config.json>\n")
    cat("Example: Rscript Subsetted_SO.R subsetted_so_config.json\n")
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
    params <- read_json(config_file)
    cat("Loading configuration from:", config_file, "\n")
  }, error = function(e) {
    cat("Error reading JSON file:", e$message, "\n")
    quit(status = 1)
  })
  
  # Set default values for missing parameters
  if (is.null(params$category_to_filter)) params$category_to_filter <- "HPCA_main"
  if (is.null(params$seed)) params$seed <- 10
  if (is.null(params$samples_to_include)) params$samples_to_include <- "E1,E2,E3,H1,H2,H3,T1,T2,T3"
  if (is.null(params$sample_name)) params$sample_name <- "orig_ident"
  if (is.null(params$keep_or_remove)) params$keep_or_remove <- TRUE
  if (is.null(params$values_to_filter)) params$values_to_filter <- c("T_cells")
  if (is.null(params$greater_less_than)) params$greater_less_than <- "greater than"
  if (is.null(params$cut_off)) params$cut_off <- 0.5
  if (is.null(params$plot_as_interactive_plot)) params$plot_as_interactive_plot <- FALSE
  if (is.null(params$use_cite_seq_data)) params$use_cite_seq_data <- FALSE
  if (is.null(params$legend_position)) params$legend_position <- "right"
  if (is.null(params$legend_symbol_size)) params$legend_symbol_size <- 2
  if (is.null(params$number_of_legend_columns)) params$number_of_legend_columns <- 1
  if (is.null(params$dot_size)) params$dot_size <- 0.1
  if (is.null(params$dot_size_highlighted_cells)) params$dot_size_highlighted_cells <- 0.5
  if (is.null(params$colors)) params$colors <- c("aquamarine3", "salmon1", "lightskyblue3", "plum3", "darkolivegreen3", "goldenrod1", "burlywood2", "gray70", "firebrick2", "steelblue", "palegreen4", "orchid4", "darkorange1", "yellow", "sienna", "palevioletred1", "gray60", "cyan4", "darkorange3", "mediumpurple3", "violetred2", "olivedrab", "darkgoldenrod2", "darkgoldenrod", "gray40", "palegreen3", "thistle3", "khaki1", "deeppink2", "chocolate3", "paleturquoise3", "wheat1", "lightsteelblue", "salmon", "sandybrown", "darkolivegreen2", "thistle2", "gray85", "orchid3", "darkseagreen1", "lightgoldenrod1", "lightskyblue2", "dodgerblue3", "darkseagreen3", "forestgreen", "lightpink2", "mediumpurple4", "lightpink1", "thistle", "navajowhite", "lemonchiffon", "bisque2", "mistyrose", "gray95", "lightcyan3", "peachpuff2", "lightsteelblue2", "lightyellow2", "moccasin", "gray80", "antiquewhite2", "lightgrey")
  
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
  cat("Running Subsetted_SO analysis...\n")
  result <- Subsetted_SO(
    seurat_object = seurat_obj,
    category_to_filter = params$category_to_filter,
    seed = params$seed,
    samples_to_include = params$samples_to_include,
    sample_name = params$sample_name,
    keep_or_remove = params$keep_or_remove,
    values_to_filter = params$values_to_filter,
    greater_less_than = params$greater_less_than,
    cut_off = params$cut_off,
    plot_as_interactive_plot = params$plot_as_interactive_plot,
    use_cite_seq_data = params$use_cite_seq_data,
    legend_position = params$legend_position,
    legend_symbol_size = params$legend_symbol_size,
    number_of_legend_columns = params$number_of_legend_columns,
    dot_size = params$dot_size,
    dot_size_highlighted_cells = params$dot_size_highlighted_cells,
    colors = params$colors
  )
  
  # Save output
  output_file <- if (!is.null(params$output_file)) params$output_file else "subsetted_so_result.rds"
  saveRDS(result, output_file)
  cat("Results saved to:", output_file, "\n")
  
  # Memory cleanup
  cat("Available memory after processing:", round(as.numeric(system("free -m | grep '^Mem:' | awk '{print $7}'", intern=TRUE)), 2), "MB\n")
}
