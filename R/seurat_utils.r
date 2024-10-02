# Author: Edoardo Filippi
# mail: efilippi@uni-mainz.de

#### MARKER ANALYSIS ####

#' Generate Multiple Plots for a Seurat Object
#'
#' This function generates various plots (e.g., dot plot, heatmap) for a Seurat object, based on the specified parameters.
#'
#' @param seurat_object A Seurat object containing the single-cell RNA-seq data.
#' @param which A character vector specifying the types of plots to generate. 
#'        Options are \code{"dotplot"} and \code{"heatmap"}. Default is \code{c("dotplot", "heatmap")}.
#' @param clusters A vector specifying the clusters to include in the plots. 
#'        If \code{NA}, all clusters will be included. Default is \code{NA}.
#' @param assay A character string specifying the assay to use from the Seurat object. 
#'        Default is \code{"RNA"}.
#' @param cluster_column A character string specifying the column name in the Seurat object's metadata 
#'        that contains the cluster assignments. Default is \code{"harmony_clusters"}.
#' @param name A character string specifying the prefix to use for the saved plot files. 
#'        Default is an empty string \code{""}.
#' @param markers A logical value indicating whether to include marker genes in the plots. 
#'        Default is \code{FALSE}.
#' @param extension_plot A character string specifying the file extension for the saved plots. 
#'        Default is \code{".png"}.
#' @param maxn_genes An integer specifying the maximum number of genes to include in the plots. 
#'        Default is \code{100}.
#' @param n_genes An integer specifying the number of top genes to include in each plot. 
#'        Default is \code{25}.
#' @param maxn_genes_per_plot An integer specifying the maximum number of genes to include per plot. 
#'        Default is \code{100}.
#'
#' @return This function saves the generated plots to files and does not return a value.
#'
#' @examples
#' # Generate dot plots and heatmaps for a Seurat object
#' plot_heatmap(seurat_object = my_seurat, which = c("dotplot", "heatmap"), clusters = c(1, 2, 3))
#'
#' @export
plot_heatmap <- function(seurat_object, which = c("dotplot", "heatmap"), clusters = FALSE, assay = "RNA",
                   cluster_column = "harmony_clusters", name = "", markers = FALSE,
                   extension_plot = ".png", maxn_genes = 100, n_genes = 25, maxn_genes_per_plot = 100, sorting_method = "abs") {
  
  # Check arguments
  if (!inherits(seurat_object, "Seurat")) stop("seurat_object must be a Seurat object")
  if (!is.character(which) 
      || !all(which %in% c("dotplot", "heatmap"))) stop("which argument must be a character vector containing 'dotplot' and/or 'heatmap'")
  if (!isFALSE(clusters) 
      || (!is.integer(clusters) && !is.vector(clusters))) stop("clusters argument must be an integer")
  if (!is.character(assay)) stop("assay argument must be a string")
  if (!is.character(cluster_column)) stop("cluster_column argument must be a string (or int)")
  if (!is.character(name)) stop("name argument must be a character")
  if (!isFALSE(markers) 
      && ((!is.character(markers) && !is.vector(markers))
      || (!is.integer(markers) && !is.vector(markers)))) stop("markers argument must be a character vector or False")
  if (!is.character(extension_plot)) stop("extension_plot argument must be a string")
  if (!is.numeric(maxn_genes) || maxn_genes <= 0) stop("maxn_genes argument must be a positive numeric value")
  if (!is.numeric(n_genes) || n_genes <= 0) stop("n_genes argument must be a positive numeric value")
  if (!is.numeric(maxn_genes_per_plot) || maxn_genes_per_plot <= 0) stop("maxn_genes_per_plot argument must be a positive numeric value")


  # Messages
  message(paste0("Parameters: which: ", paste(which, collapse = ", "),
                 " - clusters: ", clusters,
                 " - assay: ", assay,
                 " - cluster_column: ", cluster_column,
                 " - name: ", name,
                 " - markers: ", paste(markers, collapse = ", "),
                 " - extension_plot: ", extension_plot,
                 " - maxn_genes: ", maxn_genes,
                 " - n_genes: ", n_genes,
                 " - maxn_genes_per_plot: ", maxn_genes_per_plot))
  

  # Set up output dir, make it second level if the name is an absolute path 
  if (grepl("^[A-Za-z]:/", name))
    output_dir <- set_up_output(paste0(name, "heatmap/")) 
  else   
    output_dir <- set_up_output(paste0(output_folder, "plot_heatmap_", name, "/"))
  
  # Arguments check
  if (!isFALSE(clusters) && is.vector(clusters)) {
    message("Heatmap: clusters provided")
    clusters <- as.character(clusters)
  } # Nota -> Na is considered vector of length 1
  else clusters <-  as.character(unique(seurat_object@meta.data[[cluster_column]]))
  if (!is.vector(which) || !is.character(which)) stop(paste0("Error: the variable which must be a vector with the desired plot names"))
  if (!isFALSE(markers)) compute_markers <- FALSE
  else compute_markers <- TRUE
  
  # Definition of function used to return a list of markers correctly formatted
  find_markers_local <- function() {
    
    # Initialize emty list where the lists of markers will be kept
    cluster_markers <- list()
    
    # Iterate through all the clusters and compute the markers
    for (cluster_number in clusters) {
      cluster_markers[[cluster_number]] <- FindMarkers(new_seurat_object, ident.1 = cluster_number, 
                                                       ident.2 = clusters[!(clusters == cluster_number)])
      message(paste0("markers computed for cluster: ", cluster_number))
    }
    
    # create a dataframe that contains the n_genes top markers per cluster, oly keeping the ones present in variable features
    markers_df <- as.data.frame(lapply(cluster_markers, function(markers_df_for_cluster) {
      if (sorting_method == "abs") markers_df_for_cluster <- markers_df_for_cluster[order(abs(markers_df_for_cluster$avg_log2FC), decreasing = TRUE), ]
      else  if (sorting_method == "min") markers_df_for_cluster <- markers_df_for_cluster[order(markers_df_for_cluster$avg_log2FC), ]
      else  if (sorting_method == "max") markers_df_for_cluster <- markers_df_for_cluster[order(abs(markers_df_for_cluster$avg_log2FC), decreasing = TRUE), ]
      else  if (sorting_method == "pvalue") markers_df_for_cluster <- markers_df_for_cluster[order(abs(markers_df_for_cluster$p_val_adj)), ]
      else  stop("Heatmap plot, deg: invalid parameter: sorting_method")

      # Select the top n_genes markers, after having checked that they are present in the variable features
      top_markers <- row.names(markers_df_for_cluster) %>%
        .[. %in% rownames(seurat_object@assays[[assay]]$scale.data)] %>%
        .[1:n_genes]
      # return(top_markers)

    })
    )
    
    # Linearise the dataframe
    markers_vector <- unlist(c(markers_df), use.names = FALSE)
    
    return(markers_vector)
  }
  
  # Loops thrugh the which argument and plots all the ones that are needed using the selected method 
  for (which_element in which){
    
    if (which_element == "dotplot") {
      n_genes <- 5
      
      # Create new seurat object
      new_seurat_object <- create_object_from_cluster_id(seurat_object, clusters, assay = assay,
                                                         clusters_column = cluster_column, save = FALSE, new_idents = cluster_column)
      # Old code
      if (FALSE) {
        new_seurat_object <- CreateSeuratObject(seurat_object[[assay]],
                                                assay = "RNA",
                                                meta.data = seurat_object@meta.data)
        Idents(new_seurat_object) <- cluster_column
      }
      
      # Compute the markers to plot if needed
      if (compute_markers) markers <- find_markers_local()
      
      # Create a DotPlot
      save_plot(DotPlot(
        new_seurat_object,
        features = unique(markers),
        group.by = cluster_column,
        idents = clusters,
        cols = c("white", "blue"),
        dot.scale = 3 
      ), paste0(output_dir, "dotplot", extension_plot), x = 20, y = 5)
    } 
    if (which_element == "heatmap") {
      
      # Create new seurat object
      new_seurat_object <- create_object_from_cluster_id(seurat_object, clusters, assay = assay,
                                                         clusters_column = cluster_column, save = FALSE)
      
      # Compute the markers to plot if needed, otherwise filter to keep the ones that are given (take the top n)
      if (compute_markers) markers_filtered <- find_markers_local()
      else { 
        if (!any(unique(markers) %in% rownames(seurat_object@assays[[assay]]$scale.data))) {
          warning("non of the genes:", paste(markers, collapse = ", "), "is present in the seurat object given")
          return()
        }
        markers_filtered <- unique(markers) %>%
        .[. %in% rownames(seurat_object@assays[[assay]]$scale.data)]  %>%
        .[1:maxn_genes]
      }

      # To avoid that te maxgenes per plot parameter is smaller than th number of genes (throws an error)
      if (maxn_genes_per_plot > length(markers_filtered)) 
        maxn_genes_per_plot <- length(markers_filtered)
      
      # Loop to create the heatmaps
      for (i in seq(1, length(markers_filtered), by = maxn_genes_per_plot)) {
        
        subset_markers_filtered <- markers_filtered[i:min(i + maxn_genes_per_plot, length(markers_filtered))]
        
        # Create a heatmap
        save_plot(DoHeatmap(
          new_seurat_object,
          features = subset_markers_filtered,
          group.by = cluster_column,
          slot = "scale.data", # scale data contains only the variable features 
          disp.min = 0 # to modify
        ) + scale_fill_gradientn(colors = c("powderblue", "white", "red")), 
        paste0(output_dir, "heatmap_", i, "_", min(i + maxn_genes_per_plot - 1, length(markers_filtered)), extension_plot), x = 15, y = 20)
        
      }
      # Create a heatmap
      # save_plot(DoHeatmap(
      #  new_seurat_object,
      #  features = markers_filtered,
      #  group.by = cluster_column,
      #  slot = "scale.data", # scale data contains only the variable features 
      #  disp.min = 0 # to modify
      # ) + scale_fill_gradientn(colors = c("powderblue","white", "red")), paste0(output_dir, "heatmap", extension_plot), x = 15, y = 20)
      # r colors: https://www.datanovia.com/en/blog/awesome-list-of-657-r-color-names/
      
      # Save plotted features in a dataframe
      openxlsx::write.xlsx(as.data.frame(markers_filtered), file = paste0(output_dir, "markers.xlsx"), sheetName = "marker_genes", append = TRUE)
      message(paste0("results saved in: ", output_dir, "markers.xlsx"))
      
    }
    if (which_element == "barplot") {
      n_genes <- 50
      
      # Create new seurat object
      new_seurat_object <- create_object_from_cluster_id(seurat_object, clusters, assay = assay,
                                                         clusters_column = cluster_column, save = FALSE, new_idents = cluster_column)
      
      # Compute the markers to plot if needed
      if (compute_markers) markers <- find_markers_local(n_genes)
      
      barplot(t(as.matrix(df$Value)), beside = TRUE, col = rainbow(length(unique(df$Cluster))),
              names.arg = unique(df$Cluster), xlab = "Cluster", ylab = "Expression Value",
              main = "Grouped Bar Plot by Cluster")
      
      # Create a DotPlot
      save_plot(DoHeatmap(
        new_seurat_object,
        features = unique(markers),
        group.by = cluster_column 
      ), paste0(output_dir, "_barplot", extension_plot), x = 12, y = 12)
    }
  }
  
}

# Function to find differentially expressed markers and plot feature plot
find_and_plot_markers <- function(seurat_object, cluster_id = "all", reduction_name = "pca", 
                                  save_data = TRUE, cluster_column = "",
                                  name = "", assay = "RNA", method = "default",
                                  subset_id = "control", o2 = NA, condition = "PD",
                                  nothreshold = FALSE,
                                  control = "non_PD", condition_column = "subject_pathology",
                                  extension_plot = ".png", feature_plots_top9_deg = FALSE, ...) {
  
  # Set up output dir
  output_dir <- set_up_output(paste0(output_folder, "markers_", name, "/"))
  message("Running differentially expressed markers analysis")

  # Setting parameters
  if (cluster_column != "" && cluster_column %in% names(seurat_object@meta.data)) Idents(seurat_object) <- cluster_column
  else if (cluster_column != "" && !(cluster_column %in% names(seurat_object@meta.data))) 
    stop("invalid idents column selected")
  
  # Find differentially expressed markers with the specified method
  if (method == "microglia_only") {
    microglia <- c("0", "1", "2")
    message("running for microglia")
    markers <- FindMarkers(seurat_object, ident.1 = cluster_id, ident.2 = microglia[!(microglia == cluster_id)])
  } else if (method == "condition") {
    message("running for condition")
    if (nothreshold) 
      markers <- FindMarkers(seurat_object, ident.1 = control, ident.2 = condition, group.by = condition_column,
                             logfc.threshold = 0, min.pct = 0) #,  test.use="DESeq2")
    else 
      markers <- FindMarkers(seurat_object, ident.1 = control, ident.2 = condition, group.by = condition_column) #,  test.use="DESeq2")
  } else if (method == "condition_vf") {
    message("running for condition")
    if (nothreshold) 
      markers <- FindMarkers(seurat_object, ident.1 = control, ident.2 = condition, group.by = condition_column,
                             logfc.threshold = 0, min.pct = 0, 
                             slot = "scale.data", features = row.names(seurat_object@assays[[assay]]$scale.data)) #,  test.use="DESeq2")
    else 
      markers <- FindMarkers(seurat_object, ident.1 = control, ident.2 = condition, group.by = condition_column, 
                             slot = "scale.data", features = row.names(seurat_object@assays[[assay]]$scale.data)) #,  test.use="DESeq2")
    
  } else if (method == "condition_and_clusters" && !is.na(subset_id)) {
    message("running for condition and clusters")
    markers <- FindMarkers(seurat_object, ident.1 = control, group.by = condition_column,
                           subset.ident = cluster_id)#, test.use="DESeq2")
  } else if (method == "condition_and_clusters_vf" && !is.na(subset_id)) {
    message("running for condition and clusters")
    markers <- FindMarkers(seurat_object, ident.1 = control, group.by = condition_column,
                           subset.ident = cluster_id, 
                           slot = "scale.data", features = row.names(seurat_object@assays[[assay]]$scale.data))#, test.use="DESeq2")
  } else if (method == "default") {
    markers <- FindMarkers(seurat_object, ident.1 = cluster_id)
  } else if (method == "default_vf") {
    markers <- FindMarkers(seurat_object, ident.1 = cluster_id, 
                           slot = "scale.data", features = row.names(seurat_object@assays[[assay]]$scale.data)) 
  } else if (method == "default_LR") {
    markers <- FindMarkers(seurat_object, ident.1 = cluster_id, test.use = "LR")
  } else if (method == "condition_DESeq2") {
    message("running for condition with DESeq2")
    
    # Delete the features that have an expression of 0 everywhere, (hwo is it possile????) rowsumns on the condition > 0, if rowsums is 0 then there are no cells that satisfy it
    seurat_object <- subset(seurat_object, features = rownames(seurat_object@assays$RNA$data)[rowSums(seurat_object@assays$RNA$data > 0) > 0])
    
    if (nothreshold) 
      markers <- FindMarkers(seurat_object, ident.1 = control, ident.2 = condition, group.by = condition_column,
                             logfc.threshold = 0, min.pct = 0, test.use = "DESeq2")
    else 
      markers <- FindMarkers(seurat_object, ident.1 = control, ident.2 = condition, group.by = condition_column, test.use = "DESeq2")
    
  } else if (method == "condition_and_clusters_DESeq2" && !is.na(subset_id)) {
    message("running for condition and clusters with DESeq2")
    seurat_object <- subset(seurat_object, features = rownames(seurat_object@assays$RNA$data)[rowSums(seurat_object@assays$RNA$data > 0) > 0])
    
    markers <- FindMarkers(seurat_object, ident.1 = control, group.by = condition_column,
                           subset.ident = cluster_id, test.use = "DESeq2")
  } else {
    stop("no condition met")
  }
  

  
  if (save_data) {
    if (!dir.exists(output_dir)) dir.create(output_dir)
    
    markers_excel <- markers# [1:1000, ]
    markers_excel$gene <- rownames(markers_excel)
    
    if (method == "default" || method == "microglia_only") {
      
      # Check if the group is numeric (for some reason if it is numeric it adds a g at the end)
      if (is.na(as.numeric(cluster_id))) expr <-
          AggregateExpression(object = seurat_object)[[assay]][, cluster_id]
      else       expr <- 
          AggregateExpression(object = 
                                seurat_object)[[assay]][, paste0("g", as.character(cluster_id))]
      
      # Add the average expression column
      markers_excel <- cbind(markers_excel, "average expression" = 
                               expr[names(expr) %in% row.names(markers_excel)])
    }
    
    openxlsx::write.xlsx(markers_excel, file = paste0(output_dir, "expressed_markers_", cluster_id, "_", name, ".xlsx"), sheetName = "marker_genes", append = TRUE)
    message(paste0("results saved in: ", output_dir, "expressed_markers_", cluster_id, ".xlsx"))
  }
  if (feature_plots_top9_deg) save_plot(FeaturePlot(seurat_object, features = rownames(markers[1:9, ]), cols = c("lightgrey", "blue"),
                             reduction = reduction_name), #  min.cutoff="q15"),
                 paste0(output_dir, "feature_plot_", cluster_id, extension_plot), x = 10, y = 10)
  return()
  # Create a feature plot for the top markers
  if (!is.na(o2)) save_plot(FeaturePlot(o2, features = marker_genes, cols = c("lightgrey", "blue"),
                                        reduction = reduction_name), #  min.cutoff="q15"),
                            paste0(output_dir, "feature_plot", cluster_id, extension_plot), x = 10, y = 10)#
  
}

# Volcano plot from gene table
volcano_plot <- function(markers_dir, count_threshold = 0, extension_plot = ".png") {
  
  # Check packages
  library(EnhancedVolcano)

  # Set up output dir
  output_dir <- set_up_output(paste0(output_folder, "markers_", markers_dir, "/"), message)

  # List excel files
  excel_files <- list.files(output_dir, pattern = "\\.xlsx$", full.names = TRUE) # TODO: add correct pattern (take only marksers tables)
  
  # Plot volcano for all
  for (excel_file in excel_files) {
    
    # Load source
    source <- as.data.frame(readxl::read_excel(excel_file))
    
    # Filter dataframe
    source <- source[source$pct.1 > count_threshold & source$pct.2 > count_threshold, ]
    
    # Create dataframe to plot
    source_df <- source[, c("avg_log2FC", "p_val_adj")]
    rownames(source_df) <- unlist(source["gene"])

    # Save plot
    save_plot(EnhancedVolcano(source_df,
                              lab = rownames(source_df),
                              x = "avg_log2FC",
                              y = "p_val_adj") + 
                labs(subtitle = tools::file_path_sans_ext(basename(excel_file))),
              paste0(output_dir, "volcano_plot_", tools::file_path_sans_ext(basename(excel_file)), as.character(count_threshold), "", extension_plot), x = 10, y = 10)
  }
}

# Save a gene expression table for selected gene list
gene_table <- function(seurat_object, gene_list, message = "resutls of gene analysis", name = "gabriel",
                       assay = "RNA", method = "subjects", name2 = "", cluster_column = "microglia_clusters") {
  
  # Set up output dir
  output_dir <- set_up_output(paste0(output_folder, "specific_gene_analysis_", name2, "/"), message)
  
  # Initialize df
  Idents(seurat_object) <- cluster_column
  row_names <- gene_list
  col_names <- c(as.character(unique(Idents(seurat_object))))
  df <- as.data.frame(matrix(NA, nrow = length(row_names), ncol = length(col_names),
                             dimnames = list(row_names, col_names)))
  features <- intersect(row.names(seurat_object@assays[[assay]]$data), gene_list)
  
  # Aggregate expression
  expr <- AggregateExpression(object = seurat_object, features = features)[[assay]]
  df[features, ] <- expr[features, ]
  
  # Initialize df
  Idents(seurat_object) <- "subject_pathology"
  col_names <- c(as.character(unique(Idents(seurat_object))))
  df[, col_names] <- NA
  
  # Diff expression and adding to a table
  for (group in col_names) {
    if (sum(seurat_object@meta.data$subject_pathology == group) < 3) next
    markers <- FindMarkers(seurat_object, ident.1 = group, features = features,
                           logfc.threshold = 0, min.pct = 0)
    log_2fc_column <- as.vector(markers$avg_log2FC)
    names(log_2fc_column) <- row.names(markers)
    
    df[features, group] <- log_2fc_column[features] # can i add it also in another way??
  }
  
  # Save df
  df <- cbind(row.names(df), df)
  colnames(df)[1] <- "gene"
  openxlsx::write.xlsx(df, file = paste0(output_dir, method, name, "_analysis.xlsx"))
}

violin_plot <- function(seurat_object, 
                        gene_list = FALSE, 
                        name = "", 
                        n_min = 3,
                        markers_analysis = FALSE, 
                        extension_plot = ".png", 
                        cluster = FALSE, 
                        ngenes_to_plot = 10) {

  if (isFALSE(markers_analysis) && isFALSE(gene_list)) stop("you must provide either a list of genes or a marker analysis as source")
  # Here markers are only used to load info on the expression data 
  # Set up output dir
  if (grepl("^[A-Za-z]:/", name))
    output_dir <- set_up_output(paste0(name, "violin_plots/")) 
  else   
    output_dir <- set_up_output(paste0(output_folder, "violin_plots_", name, "/"))


  # Load DEGs table
  if (isFALSE(marker_analysis)) {
    if(!isFALSE(gene_list)) message("gene_list AND marker source given, gene_list overridden by source")
    tryCatch({
        if (isFALSE(cluster)) {
            message(paste0("loading... ", output_folder, "markers_", markers_analysis,
                "/expressed_markers_all_", markers_analysis, ".xlxs"))
            markers_table <- openxlsx::read.xlsx(paste0(output_folder, "markers_", markers_analysis,
                "/expressed_markers_all_", markers_analysis, ".xlsx"))
        } else {
            message(paste0("loading... ", output_folder, "markers_", markers_analysis,
                "/expressed_markers_", cluster, "_", markers_analysis, ".xlxs"))
            markers_table <- openxlsx::read.xlsx(paste0(output_folder, "markers_", markers_analysis,
                "/expressed_markers_", cluster, "_", markers_analysis, ".xlsx"))
        }
    }, error = function(e) {
      stop("Could not load the deg results to compute histogram in wgcna, probably the table does not exist?: \n", e)
    })
    # sort the top n genes by absolute log2fc
    gene_df <- tibble::column_to_rownames(as.data.frame(markers_table), var = "gene")
    gene_list <- rownames(gene_df[order(abs(gene_df$avg_log2FC), decreasing = TRUE), ])[1:ngenes_to_plot]
  }
  else message("no DEGs table given")

  # Subset count matrix to include only needed genes
  count_matrix <- Matrix::Matrix(GetAssay(seurat_object, assay = "RNA")$data, sparse = TRUE)
  count_matrix <- t(count_matrix[row.names(count_matrix) %in% gene_list, ])
  # TODO: convert in sparse matrix?
  sorting_dataframe <- seurat_object@meta.data$subject_pathology
  # Iterate through the genes and create violin-plot
  purrr::walk2(data.frame(count_matrix), colnames(count_matrix), function(expr_list, gene) {

    # Matching with condition 
    expr_df <- data.frame(cbind(expr_list, sorting_dataframe))
    c1 <- nrow(expr_df)
    expr_df <- expr_df[expr_df$expr_list > 0, ]
    c2 <- nrow(expr_df)
    expr_df$expr_list <- as.numeric(expr_df$expr_list)
    expr_df$sorting_dataframe <- factor(expr_df$sorting_dataframe) #, levels = c("PD", "genetic_PD", "non_PD"))
    
    # Prepare text
    if (!markers_analysis == "") {
      
      # Read values
      qvalue <- format(as.numeric(markers_table_pd[markers_table_pd$gene == gene, "p_val_adj"], scientific = TRUE, digits = 4))
      # qvalue_gpd <- format(as.numeric(markers_table_gpd[markers_table_gpd$gene == gene, "p_val_adj"], scientific = TRUE, digits = 4))
      log2fc <- format(as.numeric(markers_table_pd[markers_table_pd$gene == gene, "avg_log2FC"], scientific = TRUE, digits = 4))
      # log2fc_gpd <- format(as.numeric(markers_table_gpd[markers_table_gpd$gene == gene, "avg_log2FC"], scientific = TRUE, digits = 4))
      
      # Prepare text
      text <- paste0("cells which show expression: ", c2, "/", c1, "\n",
                     # "genetic PD adj_pvalue: ", qvalue_gpd, " - log2FC: ", log2fc_gpd, "\n",
                     "adj_pvalue: ", qvalue, " - log2FC: ", log2fc)
      # message(text)
    }
    else 
      text <- paste0("cells which show expression: ", c2, "/", c1)
    
    # Create violin plot
    if (c2 > n_min) 
      save_plot(ggplot(expr_df, aes(x = sorting_dataframe, y = expr_list, colour = sorting_dataframe)) + #, fill=sorting_dataframe)) +
                  geom_violin(alpha = 0.2) + geom_jitter(alpha = 0.8, size = 0.2) + # scale_fill_manual(values = c("red", "blue", "green")) +
                  labs(x = NULL, y = "Expression Value") +
                  ggtitle(paste0("Violin Plot for ", gene)) + 
                  # theme(plot.margin = margin(b=1.5, unit="cm")) +
                  # labs(caption = str_wrap(text, width = 60)) +
                  labs(caption = text) +
                  theme(legend.position = "none"),
                paste0(output_dir, "volin_plot_", gene, extension_plot))
    else message(gene, " expressed in too little cells (<", n_min, ")")

  })
}

# Saves a table with the average expression vales of the data
avg_expression <- function(seurat_object, gene_list, name = "") {
  
  # Destination Folder
  output_dir <- set_up_output(paste0(output_folder, "avg_expression_", name, "/"), message)
  
  # Extract count matrix
  count_matrix <- Matrix(GetAssay(seurat_object, assay = "RNA")$data, sparse = TRUE)
  count_matrix <- data.frame(count_matrix[row.names(count_matrix) %in% gene_list, ])
  
  # saves the gene list as dataframes and assigns rownames
  means <- as.data.frame(gene_list)
  row.names(means) <- gene_list
  
  # Iterate though the different conditions and compute means
  # TODO: use apply?
  for (pathology in unique(seurat_object@meta.data$subject_pathology)) {
    cells <- row.names(seurat_object@meta.data[seurat_object@meta.data$subject_pathology == pathology, ])

    means <- cbind(means, as.data.frame(rowMeans(count_matrix[, colnames(count_matrix) %in% cells])))
    names(means)[length(means)] <- pathology
  }

  # Save excel
  openxlsx::write.xlsx(means, file = paste0(output_dir, "means_", name, ".xlsx"))
  message(paste0("results saved in: ", output_dir, "means_", name, ".xlsx"))
}
  
#### CLUSTERING ####

#' Performs clustering, selecting the number of dimensons to consider based on a threshold on the standard deviation, and the visualization, the results of the visualization are passed as parameters
#'
#'
#' This function performs clustering on a Seurat object using a specified dimensional reduction method 
#' (e.g., PCA) and saves the clustering results in the object's metadata. The function also optionally saves 
#' the Seurat object to a file.
#'
#' @param seurat_object A Seurat object containing the single-cell RNA-seq data.
#' @param reduction A character string specifying the dimensional reduction method to use for clustering. 
#'        Default is \code{"pca"}.
#' @param dimensions An integer specifying the number of dimensions to use for clustering. 
#'        Default is \code{15}.
#' @param desired_resolution A numeric value specifying the resolution parameter for clustering. 
#'        Higher values lead to more clusters. Default is \code{0.6}.
#' @param save A logical value indicating whether to save the Seurat object after clustering. 
#'        Default is \code{FALSE}.
#' @param column_name A character string specifying the name of the metadata column to store the cluster 
#'        assignments. Default is \code{"seurat_clusters"}.
#'
#' @return The Seurat object with updated clustering information in the metadata.
#'
#' @details This function first builds a nearest neighbor graph using the specified reduction method and 
#' dimensions. It then performs clustering on the graph and saves the resulting cluster assignments 
#' in the specified metadata column of the Seurat object. If the \code{column_name} is not 
#' \code{"seurat_clusters"}, the default \code{"seurat_clusters"} column is removed from the metadata.
#'
#' @examples
#' # Perform clustering using the first 20 principal components
#' seurat_object <- clustering(seurat_object = my_seurat, reduction = "pca", 
#'                             dimensions = 20, desired_resolution = 0.8)
#'
#' @export
clustering <- function(seurat_object, reduction = "pca", 
                       dimensions = 15, desired_resolution = 0.6, save = FALSE, column_name = "seurat_clusters") {
  
  # Message
  cat("Clustering.... \n")
  message(paste0("Parameters: reduction: ", reduction,
                 " - dimensions: ", dimensions,
                 " - desired_resolution: ", desired_resolution,
                 " - save: ", save,
                 " - column_name: ", column_name))
  
  # build the graph
  seurat_object <- FindNeighbors(seurat_object, reduction = reduction, dims = 1:dimensions, verbose = FALSE)
  # graph.name assiged according to: https://github.com/satijalab/seurat/issues/2995
  
  # Message
  message("Graph computed")
  message(paste0("running FindClusters, graph used:", seurat_object[[reduction]]@assay.used, "_snn",
                 " - saved in metadata column: ", column_name))
  
  # Compute the clusters
  seurat_object <- FindClusters(seurat_object, resolution = desired_resolution, 
                                cluster.name = column_name, 
                                graph.name = paste0(seurat_object[[reduction]]@assay.used, "_snn"), verbose = FALSE)
  
  # Delete seurat_clusters column if created by mistake
  if (column_name != "seurat_clusters")  seurat_object$seurat_clusters <- NULL 
  
  # Saving the object if needed
  if (save) saveRDS(seurat_object, file = paste0(output_folder, "after_clustering.rds"))
  message("Clustering done")
  
  return(seurat_object)
}

# Visualization, given a reduction to use and the de1sired stdeviation to select the desired number of features 
#' UMAP Visualization and Plotting for a Seurat Object
#'
#' This function generates UMAP visualizations for a Seurat object, creating and saving multiple plots 
#' based on the specified parameters. It can also optionally run UMAP dimensionality reduction if required.
#'
#' @param seurat_object A Seurat object containing the single-cell RNA-seq data.
#' @param reduction_name A character string specifying the name to save the UMAP reduction under. 
#'        Default is \code{"umap_"}.
#' @param reduction A character string specifying the dimensional reduction method to use as input for UMAP. 
#'        Default is \code{"pca"}.
#' @param dimensions An integer specifying the number of dimensions to use for UMAP. 
#'        Default is \code{16}.
#' @param cluster_column A character string specifying the column name in the Seurat object's metadata 
#'        that contains the cluster assignments. Default is \code{"seurat_clusters"}.
#' @param plots A character vector specifying the names of the plot files to save. 
#'        If \code{NA}, default names will be generated based on the \code{name} and \code{extension_plot} parameters.
#' @param save A logical value indicating whether to save the Seurat object after generating the plots. 
#'        Default is \code{FALSE}.
#' @param run_umap A logical value indicating whether to run UMAP before generating the plots. 
#'        Default is \code{TRUE}.
#' @param name A character string specifying the prefix for the plot names. 
#'        Default is \code{""}.
#' @param message A character string specifying a custom message to append to an information log. 
#'        Default is \code{"results"}.
#' @param extension_plot A character string specifying the file extension for the saved plots. 
#'        Default is \code{".png"}.
#'
#' @return The Seurat object with the UMAP reduction added to it (if \code{run_umap} is \code{TRUE}).
#'
#' @details This function handles UMAP dimensionality reduction and generates several plots based on the specified 
#' reduction. The plots can be customized and saved with different naming conventions. The \code{daniela} parameter 
#' controls a specific set of plots, while the default generates a different set. If the UMAP has already been run, 
#' \code{run_umap} can be set to \code{FALSE} to skip rerunning it.
#'
#' @examples
#' # Generate UMAP plots and save them to files
#' seurat_object <- visualization_UMAP(seurat_object = my_seurat, dimensions = 20, 
#'                                     save = TRUE, name = "sample_umap")
#'
#' @export
visualization_UMAP <- function(seurat_object, reduction_name = "umap_",
                                          reduction = "pca", dimensions = 16, 
                                          cluster_column = "seurat_clusters",
                                          plots = FALSE, save = FALSE, 
                                          run_umap = TRUE, name = "", extension_plot = ".png",
                                          daniela = FALSE) {

  # Check that cluster column is actually in the metadada
  if (!(cluster_column %in% colnames(seurat_object@meta.data))) stop(paste("Error: cluster_column", cluster_column, "does not exist in the Seurat object's metadata."))
  if (isFALSE(run_umap) && !reduction_name %in% names(seurat_object@reductions)) stop("compute_umap (or run_umap) is set to false and the selected reduction does not exist")
  # Setting plot names
  if (nchar(name) == 0) name <- reduction_name
  if (isFALSE(plots)) {
    plots <- c(
      paste0(name, "_by_pathology", extension_plot),
      paste0(name, "_clusters", extension_plot),
      paste0(name, "_split_pathology", extension_plot),
      paste0(name, "_by_subject", extension_plot)
    )
  }
  
  # Setting up output directory + info
  output_dir <- paste0(output_folder, "plots_dim_red_", name, "/")

  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Function info
  cat(paste0("Visualisation and UMAP... \n"))
  
  message(paste0("Parameters: reduction_name: ", reduction_name,
                 " - reduction: ", reduction,
                 " - dimensions: ", dimensions,
                 " - cluster_column: ", cluster_column,
                 # " - plots: ", plots,
                 " - save: ", save,
                 " - run_umap: ", run_umap,
                 " - name: ", name,
                 " - extension_plot: ", extension_plot))
  
  # UMAP
  if (run_umap) seurat_object <- RunUMAP(seurat_object, dims = 1:dimensions, reduction = reduction, 
                                         reduction.name = reduction_name, verbose = FALSE)
  
  # Visualisation
  # TODO: delete daniela or/ and add some other way of controlling plots
  if (daniela) {
    save_plot(DimPlot(seurat_object, reduction = reduction_name, group.by = "Condition"), paste0(output_dir, plots[1]))
    save_plot(DimPlot(seurat_object, reduction = reduction_name, group.by = cluster_column, label = TRUE),  paste0(output_dir, plots[2]))
    save_plot(DimPlot(seurat_object, reduction = reduction_name, split.by = "Condition"), paste0(output_dir, plots[3]), x = 14, y = 7)
    save_plot(DimPlot(seurat_object, reduction = reduction_name, group.by = "Sample", label = TRUE), paste0(output_dir, plots[4]))
  } else {
    save_plot(DimPlot(seurat_object, reduction = reduction_name, group.by = "subject_pathology"), paste0(output_dir, plots[1]))
    save_plot(DimPlot(seurat_object, reduction = reduction_name, group.by = cluster_column, label = TRUE) + theme(legend.position = "none"),  paste0(output_dir, plots[2]))
    save_plot(DimPlot(seurat_object, reduction = reduction_name, split.by = "subject_pathology"), paste0(output_dir, plots[3]), x = 14, y = 7)
    save_plot(DimPlot(seurat_object, reduction = reduction_name, group.by = "subject", label = TRUE) + theme(legend.position = "none"), paste0(output_dir, plots[4]))
  }
  
  # Save the seurat object if requested
  if (save) {
    saveRDS(seurat_object, file = paste0(output_folder, "after_cluster_visualization.rds"))
    message(paste0("seurat object saved in ´", paste0(output_folder, "after_cluster_visualization.rds")))}
  
  return(seurat_object)
}
#### SEURAT OBJECTS MANIPULATION #####

# Find the files that respect a pattern in a root directory and returns the info
preparation_for_data_loading <- function(root, pattern, path_to_patient_info, source = "celescope") {
  
  if (source == "celescope") {
    # location of files and creation of objects 
    count_matrix_files <- find_matching_matrices_paths(root, pattern)
    subjects_info <- load_conditions(path_to_patient_info)
    message(paste0(capture.output(subjects_info), collapse = "\n"))

    count_matrix_files <- match_names(subjects_info, count_matrix_files)
    message("all file names correctly found? ", nrow(subjects_info) == length(count_matrix_files))

    return(list(count_matrix_files = count_matrix_files, subjects_info = subjects_info))

  } else if (source == "textfile") {

    count_matrix_files <- find_matching_matrices_paths(root, pattern, source = "textfile")
    subjects_info <- load_conditions(path_to_patient_info)
    message(paste0(capture.output(subjects_info), collapse = "\n"))

  return(list(count_matrix_files = count_matrix_files, subjects_info = subjects_info))
  } 
  
}

#' Prepare Data for Loading by Finding Matching Files
#'
#' This function searches for files in a root directory that match a specified pattern, 
#' loads patient information, and returns the paths to the matching files along with the patient information.
#'
#' @param root A character string specifying the root directory where the search for files should be performed.
#' @param pattern A character string specifying the pattern to match file names against.
#' @param path_to_patient_info A character string specifying the path to the file containing patient information.
#' @param source A character string specifying the source type. Options are \code{"celescope"} or \code{"textfile"}.
#'        This determines how the matching files are identified. Default is \code{"celescope"}.
#'
#' @return A list containing two elements:
#' \item{count_matrix_files}{A vector of file paths that match the specified pattern.}
#' \item{subjects_info}{A data frame containing the loaded patient information.}
#'
#' @details The function first finds all files in the root directory that match the given pattern. It then loads 
#' the patient information from the provided path. If the \code{source} is \code{"celescope"}, it verifies that 
#' all file names are correctly matched with the patient information. The function returns a list with the file 
#' paths and the patient information.
#'
#' @examples
#' # Prepare data for loading from celescope source
#' data_info <- preparation_for_data_loading(root = "/data/sequencing_results", 
#'                                           pattern = "matrix.mtx.gz", 
#'                                           path_to_patient_info = "/data/patient_info.csv")
#'
#' @export
seurat_objects_and_quality_control <- function(count_matrix_files, subjects_info, save = FALSE,
                                               normalization = FALSE, source = "celescope", ...) {


  # Setting up output dir 
  output_dir <- set_up_output(paste0(output_folder, "quality_control/"), message)
  
  # Define function and needed objects
  seurat_objects <- list()
  plot_doublets <- function(plot, points) {
    # https://github.com/satijalab/seurat/blob/master/R/visualization.R

    # Funzione che serve per recuperare i punti 
    GetXYAesthetics <- function(plot, geom = "GeomPoint", plot.first = TRUE) {
      geoms <- sapply(
        X = plot$layers,
        FUN = function(layer) {
          return(class(x = layer$geom)[1])
        }
      )
      # handle case where raster is set to True
      if (geom == "GeomPoint" && "GeomScattermore" %in% geoms) {
        geom <- "GeomScattermore"
      }
      geoms <- which(x = geoms == geom)
      if (length(x = geoms) == 0) {
        stop("Cannot find a geom of class ", geom)
      }
      geoms <- min(geoms)
      if (plot.first) {
        # x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
        x <- as_label(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)
        # y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
        y <- as_label(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)
      } else {
        x <- as_label(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)
        y <- as_label(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)
      }
      return(list("x" = x, "y" = y))
    }
    xynames <- GetXYAesthetics(plot = plot)
    
    
    label_data <- plot$data[points, ]
    # label_data$labels <- labels
    # Per fre gli scatter usano una funzione anche da loro definita che si chiama SingleCorPlot
    # secondo me dovrebe essere possibile ottenere lo stesso risultato semplicemente con geom scatter

    plot <- plot + geom_point( #https://ggplot2.tidyverse.org/reference/geom_point.html
      mapping = aes_string(x = xynames$x, y = xynames$y),
      data = label_data,  
      colour = "blue")
    
    # da LabelPoints
    # plot <- plot + geom_text_repel(
    #  mapping = aes_string(x = xynames$x, y = xynames$y, label = 'labels'),
    #  data = label_data)
    
    return(plot)
  }
    
  # Main iteration
  for (file_path in count_matrix_files) {
    
    # Subject info
    object_name <- remove_parts(file_path, ...)

    if (object_name %in% subjects_info$subject) {
      pathology <- subjects_info[subjects_info$subject == object_name, "condition"]
    } else {
      message("skipping subject: ", object_name, "check that the parts_to_remove parameter and the subject name in the subject description file are correctly set")
      next
    }

    message(paste0("loading ", object_name, " pathology: ", pathology))
    
    # functions definition (here to defne them in the right environment)
    doublet_finder <- function(seurat_object) {
      
      # remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')
      library("DoubletFinder")
      
      # Source -> https://github.com/chris-mcginnis-ucsf/DoubletFinder
      ## Pre-process Seurat object (standard) --------------------------------------------------------------------------------------
      seurat_object <- NormalizeData(seurat_object)
      seurat_object <- FindVariableFeatures(seurat_object)
      seurat_object <- ScaleData(seurat_object)
      seurat_object <- RunPCA(seurat_object)
      seurat_object <- RunUMAP(seurat_object, dims = 1:10)
      seurat_object <- FindNeighbors(seurat_object, dims = 1:10)
      seurat_object <- FindClusters(seurat_object, resolution = 0.5)
      
      ## Fit linear model for expected number of doublets--------------------------------------------------------------------------------------
      num_cells <- c(1000, 2000, 3000, 4000, 5000)
      doublet_rate <- c(0.008, 0.016, 0.023, 0.031, 0.039)
      linear_model <- lm(doublet_rate ~ num_cells)
      perc_doublet <- predict(linear_model, newdata = data.frame(num_cells = nrow(seurat_object@meta.data)))
      
      ## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
      # sweep.res.list_ <- paramSweep(seurat_object, PCs = 1:10, sct = FALSE)
      # sweep.stats_ <- summarizeSweep(sweep.res.list_, GT = FALSE)
      # bcmvn_ <- find.pK(sweep.stats)
      
      ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
      # homotypic.prop <- modelHomotypic(seurat_object@meta.data$seurat_clusters)           ## ex: annotations <- seurat_object@meta.data$ClusteringResults
      n_exp_poi <- round(perc_doublet * nrow(seurat_object@meta.data))                    ## Assuming 7.5% doublet formation rate - tailor for your dataset
      # nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))                              ## adjuste pvalue for hih confidenc of identifying doublets
      
      ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
      seurat_object <- doubletFinder(seurat_object, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = n_exp_poi, reuse.pANN = FALSE, sct = FALSE)
      
      ## Save Plots ---------------------------------------------------------------------------------------------------------------
      if (!dir.exists(paste0(output_dir, "/exdimred/"))) dir.create(paste0(output_dir, "/exdimred/"))
      save_plot(DimPlot(seurat_object, group.by = "seurat_clusters"),
                paste0(output_dir, "/exdimred/", object_name, "", extension_plot))
      
      if (!dir.exists(paste0(output_dir, "/exdimred_doublets/"))) dir.create(paste0(output_dir, "/exdimred_doublets/"))
      save_plot(DimPlot(seurat_object, group.by = names(seurat_object@meta.data)[ncol(seurat_object@meta.data)]),
                paste0(output_dir, "exdimred_doublets/", object_name, "", extension_plot))
      
      # Scatter plot
      if (!dir.exists(paste0(output_dir, "/scatter_doublets/"))) dir.create(paste0(output_dir, "/scatter_doublets/"))
      
      return(seurat_object@meta.data)
    }
    
    empty_drops <- function(seurat_object) {
      
      e_out <- emptyDrops(counts(count_data), lower = 100, test.ambient = TRUE)
      
    }
    
    # Call functions for quality control
    cell_removal <- function(seurat_object, method = "DoubletFinder") {
      
      if (method == "fixed_percentage") {
        seurat_object <- subset(seurat_object, subset = nFeature_RNA < n_features)
      }
      if (method == "fixed_value") {
        seurat_object <- subset(seurat_object, subset = nFeature_RNA < n_features)
      }
      if (method == "model_percentage") {
        num_cells <- c(1000, 2000, 3000, 4000, 5000)
        doublet_rate <- c(0.008, 0.016, 0.023, 0.031, 0.039)
        linear_model <- lm(doublet_rate ~ num_cells)
        perc_doublet <- predict(linear_model, newdata = data.frame(num_cells = nrow(seurat_object@meta.data)))
        n_features <- quantile(df$column_name, probs = perc_doublet, na.rm = TRUE)
        seurat_object <- subset(seurat_object, subset = nFeature_RNA < n_features)
      }
      if (method == "DoubletFinder") {
        
        
        metadata_df <- doublet_finder(seurat_object)
        
        cells_to_keep <- rownames(
          metadata_df[
            metadata_df[
              , ncol(metadata_df)
            ] == "Singlet", ])
        
        cells_to_discard <- rownames(
          metadata_df[
            metadata_df[
              , ncol(metadata_df) # Takes last column
            ] == "Doublet", ]) # selects rownames when the value in the last column equals to doublet
        
        save_plot(
          plot_doublets(
            plot = FeatureScatter(
              seurat_object, 
              feature1 = "nCount_RNA", 
              feature2 = "percent.mt"), 
            points = cells_to_discard),           
          plotname = paste0(output_dir, "scatter_doublets/", 
                            object_name, 
                            "_features_scatter_mt_plot", extension_plot),
          x = 10, y = 7)
        
        save_plot(
          plot_doublets(
            plot = FeatureScatter(
              seurat_object, 
              feature1 = "nCount_RNA", 
              feature2 = "nFeature_RNA"), 
            points = cells_to_discard),           
          plotname = paste0(output_dir, "scatter_doublets/", 
                            object_name, 
                            "_features_scatter_plot", extension_plot),
          x = 10, y = 7)
        
        seurat_object <- subset(seurat_object, cells = cells_to_keep)
      }
      
      return(seurat_object)
    }

    # Load the data
    if (source == "celescope") count_data <- Read10X(data.dir = file_path)
    else if (source == "textfile") count_data <-  read.table(file_path, header = TRUE)

    seurat_object <- CreateSeuratObject(counts = count_data)
    
    # add metadata
    seurat_object$subject_pathology <- pathology
    
    # Calculate percentage mithocondrial
    seurat_object$percent.mt <- PercentageFeatureSet(seurat_object, 
                                                     pattern = "^MT-")
    # Scatter plot
    save_plot(FeatureScatter(seurat_object, feature1 = "nCount_RNA", 
                             feature2 = "percent.mt") +
                FeatureScatter(seurat_object, feature1 = "nCount_RNA", 
                               feature2 = "nFeature_RNA"),
              paste0(output_dir, object_name, "_features_scatter_plot", extension_plot), x = 10, y = 7)
    
    # Violin plot
    save_plot(VlnPlot(seurat_object, 
                      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                      ncol = 3), 
              paste0(output_dir, object_name, "_features_violin_plot", extension_plot))
    
    # Originale: 5% mt, >200, <2500 nFeatures
    # Subset seurat object according to quality filters

    seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & percent.mt < 10, # avevo usato 5 o 10 in % mt?
                            features = rownames(seurat_object@assays$RNA$counts)[rowSums(seurat_object@assays$RNA$counts > 0) > 3]) 

    # With doublet finder or other method
    serurat_object <- cell_removal(seurat_object)
                            
    # Normalization 
    if (normalization) seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
    
    # Adding to seurat_objects list
    cat(paste0(object_name, " loaded, normalized and quality control performed! \n"))
    seurat_objects <- c(seurat_objects, seurat_object)
    
  }
  
  # Saving seurat object if necessary
  if (save) {
    message("saving seurat oject...")
    saveRDS(seurat_objects, file = paste0(output_folder, "after_loading.rds"))
    message("saved in: ", paste0(output_folder, "after_loading.rds"))
  }
  return(seurat_objects)
}

#' Merge a List of Seurat Objects
#'
#' This function merges a list of Seurat objects into a single Seurat object without performing any additional processing.
#'
#' @param seurat_objects A list of Seurat objects to be merged.
#' @param subjects_info A character vector containing the identifiers for each Seurat object, used to label the cells in the merged object.
#' @param save A logical value indicating whether to save the merged Seurat object to a file. 
#'        Default is \code{FALSE}.
#'
#' @return A single merged Seurat object containing the combined data from the input list of Seurat objects.
#'
#' @details The function merges multiple Seurat objects provided in the \code{seurat_objects} list into a single Seurat object. 
#' The \code{subjects_info} vector is used to assign unique cell identifiers corresponding to each Seurat object. 
#' If \code{save} is \code{TRUE}, the merged Seurat object is saved to a file named \code{"after_merging.rds"} in the output folder.
#'
#' @examples
#' # Merge a list of Seurat objects without saving the result
#' merged_object <- merge_seurat_objects(seurat_objects = list(seurat_obj1, seurat_obj2), 
#'                                       subjects_info = c("Subject1", "Subject2"))
#'
#' @export
merge_seurat_objects <- function(seurat_objects, subjects_info, save = FALSE) {
  
  # merge seurat objects
  cat("Merging seurat objects... \n")
  merged_seurat <- merge(x = seurat_objects[[1]], y = seurat_objects[-1], 
                         add.cell.ids = subjects_info)
  
  # message
  cat(paste0(length(seurat_objects), " seurat objects merged! \n"))
  
  # Save if necessary
  if (save) saveRDS(merged_seurat, file = paste0(output_folder, "after_merging.rds"))
  return(merged_seurat)
}

# To create metadata according to the object information
#' Create and Update Metadata for a Seurat Object
#'
#' This function creates new metadata columns for a Seurat object by extracting information from 
#' the sample identifiers. It can also save the updated Seurat object to a file.
#'
#' @param seurat_object A Seurat object containing the single-cell RNA-seq data.
#' @param option An integer specifying the option for metadata creation. 
#'        Currently, only \code{option = 1} is implemented. Default is \code{1}.
#' @param save A logical value indicating whether to save the updated Seurat object to a file. 
#'        Default is \code{TRUE}.
#'
#' @return The Seurat object with updated metadata, including new columns for \code{subject} and \code{barcode}.
#'
#' @details The function extracts the \code{subject} and \code{barcode} information from the row names of the 
#' Seurat object's metadata, which are assumed to follow a specific naming convention with values separated by 
#' underscores. It then adds these as new columns in the metadata. If \code{save} is \code{TRUE}, 
#' the updated Seurat object is saved to a file named \code{"after_loading.rds"} in the output folder.
#'
#' @examples
#' # Create metadata for a Seurat object and save the updated object
#' seurat_obj <- create_metadata(seurat_object = my_seurat)
#'
#' @export
create_metadata <- function(seurat_object, option = 1, save = TRUE) {
  # TODO: make the column have factors
  
  # Message
  cat("Creating seurat object metadata... \n")
  
  # Retrieves the sample values
  sample_values <- rownames(seurat_object@meta.data)
  
  # Splitting the values based on the underscore character
  split_values <- strsplit(sample_values, "_")
  
  # Extracting subject and barcode separately
  subject <- sapply(split_values, function(x) paste0(x[1], "_", x[2]))
  barcode <- sapply(split_values, function(x) x[length(x)])
  
  # Creating a data frame with the separated values
  seurat_object@meta.data$subject <- subject
  seurat_object@meta.data$barcode <- barcode
  
  # Saves the object if needed
  if (save) saveRDS(seurat_object, file = paste0(output_folder, "after_loading.rds"))
  
  message("Done")
  return(seurat_object)
}

#' Subset a Seurat Object Based on Cluster IDs
#'
#' This function subsets a Seurat object based on specified cluster IDs, and optionally saves the subsetted object to a file.
#'
#' @param seurat_object A Seurat object containing the single-cell RNA-seq data.
#' @param clusters A vector specifying the cluster IDs to subset from the Seurat object.
#' @param assay A character string specifying the assay to use. Default is \code{"RNA"}.
#' @param clusters_column A character string specifying the column in the Seurat object's metadata that contains 
#'        the cluster assignments. Default is \code{"seurat_clusters"}.
#' @param save A logical value indicating whether to save the subsetted Seurat object to a file. 
#'        Default is \code{FALSE}.
#' @param filename A character string specifying the name of the file to save the subsetted Seurat object. 
#'        Default is \code{"subset_of_object.rds"}.
#' @param new_idents A character string specifying the new identity class for the cells. 
#'        If \code{NA}, the function will retain the original \code{"seurat_clusters"} identities. Default is \code{NA}.
#' @param variable_features A logical value indicating whether to find and select variable features in the subsetted object. 
#'        Default is \code{FALSE}.
#'
#' @return The subsetted Seurat object.
#'
#' @details The function subsets the input Seurat object based on the specified cluster IDs. 
#' If \code{variable_features} is \code{TRUE}, it identifies the top 2000 variable features in the subset. 
#' The subsetted object can be saved to a file if \code{save} is \code{TRUE}. The function also ensures that 
#' necessary packages are installed before proceeding.
#'
#' @examples
#' # Subset a Seurat object by clusters 1 and 2, and save the result
#' subset_obj <- create_object_from_cluster_id(seurat_object = my_seurat, clusters = c(1, 2), save = TRUE)
#'
#' @export
create_object_from_cluster_id <- function(seurat_object, clusters, assay = "RNA", clusters_column = "seurat_clusters", 
                                          save = FALSE, filename = "subset_of_object.rds", new_idents = FALSE, 
                                          variable_features = FALSE) {
  
  # Check arguments
  if (!is.vector(clusters)) stop("clusters must be a vector")
  if (isFALSE(new_idents)) new_idents <- clusters_column

  # Set corrrect identity column
  Idents(seurat_object) <- clusters_column
  subset_seurat_object <- subset(seurat_object, idents = clusters)
  # https://stackoverflow.com/questions/71542822/conditional-subsetting-of-seurat-object/72888084#72888084?newreg=32dfabd1d23c470d95a7cb56ff7933ae
  
  # Variable feature if needed
  if  (variable_features) seurat_object <- FindVariableFeatures(seurat_object, selection.method  = "vst", 
                                                                nfeatures = 2000, verbose = FALSE)
  # Save and exit
  if (save) saveRDS(subset_seurat_object, file = paste0(output_folder, filename))
  message("subset created")
  return(subset_seurat_object)
}

# Plot a selection of markers from a df (saved in different columns)
#' Plot Markers from a DataFrame
#'
#' This function generates various types of plots (e.g., feature plots or heatmaps) based on markers provided in a DataFrame.
#' The markers are matched with the features present in the Seurat object, and plots are generated for visualizing the expression of these markers across the clusters or cells.
#'
#' @param seurat_object A Seurat object containing the data to be plotted.
#' @param markers_location A file path to an Excel file containing the markers to be plotted, where each column contains a set of marker genes.
#' @param reduction_name A character string specifying the reduction to use for plotting (default is "umap_microglia_harmony_reduction").
#' @param name A character string that will be used to name the output directory and the plots.
#' @param many_plot A logical value indicating whether to generate many plots, including a heatmap (default is TRUE).
#' @param feature_plot A logical value indicating whether to generate feature plots (default is FALSE).
#' @param extension_plot A character string specifying the file extension for the saved plots (default is ".png").
#' @param heatmap_by_column A logical value indicating whether to generate a separate heatmap for each column in the marker DataFrame (default is FALSE).
#' @param cluster_column A character string specifying the column in the Seurat object metadata that contains the cluster identities (default is "microglia_clusters").
#' @param subplot_n An integer specifying the number of features to display in each subplot (default is 9).
#' @param max_feature_plots An integer specifying the maximum number of feature plots to generate (default is 11).
#' @param max_genes_plot_heatmap An integer specifying the maximum number of genes to display in the heatmap (default is 100).
#' @param column_list A logical value or character vector. If it is a character vector, the function will only consider columns in the marker DataFrame that match the vector elements (default is FALSE).
#'
#' @details
#' This function takes a Seurat object and a DataFrame of markers (in an Excel file) and plots the markers in various ways, including feature plots and heatmaps. The markers must be present in the Seurat object, and plots are saved in the specified output directory.
#'
#' The function supports two main types of plots:
#' - **Feature plots**: For each column in the marker DataFrame, the function generates feature plots, displaying the expression of marker genes across the cells, colored according to their expression levels.
#' - **Heatmaps**: If `many_plot` is TRUE, the function generates heatmaps to display the expression of marker genes across clusters or cells.
#'
#' @return The function saves the generated plots to the specified output directory but does not return a value.
#'
#' @examples
#' \dontrun{
#' plot_markers_from_df(seurat_object = my_seurat_obj,
#'                      markers_location = "path/to/markers.xlsx",
#'                      reduction_name = "umap",
#'                      name = "marker_plots",
#'                      feature_plot = TRUE,
#'                      heatmap_by_column = TRUE,
#'                      cluster_column = "cell_clusters")
#' }
#' 
#' @export
plot_markers_from_df <- function(seurat_object, 
    markers_location, 
    reduction_name = "umap_microglia_harmony_reduction", 
    name  = "",
    many_plot = TRUE, 
    feature_plot = FALSE, 
    extension_plot = ".png",
    heatmap_by_column = FALSE, 
    cluster_column = "microglia_clusters", 
    subplot_n = 9, 
    max_feature_plots = 11, 
    max_genes_plot_heatmap = 100,
    column_list = FALSE, 
    assay = "RNA") {
  
    # Messages 
    cat("Plotting markers from given df... \n")
    message(paste0("Parameters: markers_location: ", markers_location,
                   " - reduction_name: ", reduction_name,
                   " - name: ", name,
                   " - many_plot: ", many_plot,
                   " - feature_plot: ", feature_plot,
                   " - extension_plot: ", extension_plot,
                   " - heatmap_by_column: ", heatmap_by_column,
                   " - cluster_column: ", cluster_column,
                   " - subplot_n: ", subplot_n,
                   " - max_feature_plots: ", max_feature_plots,
                   " - max_genes_plot_heatmap: ", max_genes_plot_heatmap,
                   " - assay: ", assay,
                   " - column_list: ", column_list))

    # Set up output dir
    output_dir <- set_up_output(paste0(output_folder, "plot_markers_from_df_", name, "/"))
    
    # Load the markers
    markers_df <- readxl::read_excel(markers_location)
    
    # Saves the relative plots (either using many plots and a heatmap or the feature plots)
    if (feature_plot) {
        for (column in names(markers_df)) {
          
            # Skip if column list is defined and the column is in column_list
            if (is.character(column_list) && column %in% column_list) {
                message("skipping: ", column)
                next
            }
          
            # Extract the features that are present in the Seurat object for each column
            features <- dplyr::pull(markers_df, column) %>%
                .[. %in% rownames(seurat_object)] 
          
            # Message
            message("Source: ", column, " - Features: ", paste(features, collapse = ", "))
          
            # Iterate through the features in batches of subplot_n and plot (max 11 plots)
            if (length(features) == 0) {
                message("skipping, not enough features")
                next
            }
            lapply(seq(1, min(length(features), subplot_n * max_feature_plots), by = subplot_n), function(start_index) {
                end_index <- min(start_index + subplot_n - 1, length(features))
                save_plot(Seurat::FeaturePlot(seurat_object, features = features[start_index:end_index],
                                              cols = c("lightgrey", "blue"),
                                              reduction = reduction_name),
                          paste0(output_dir, "feature_plot_", column, "_", 
                                 start_index, "_", end_index, 
                                 extension_plot),
                          x = 10, y = 10) 
            })
        }
    }
    
    # Heatmap
    if (many_plot) {
        if (heatmap_by_column || is.character(column_list)) {
          
            # Iterate through the columns
            for (column in names(markers_df)) {
              
                # Skip if no data is present
                if (is.character(column_list) && !any(column %in% column_list)) {
                    message("skipping: ", column)
                    next
                }
              
                # Extract the features that are present in the Seurat object for each column
                features <- dplyr::pull(markers_df, column) %>%
                    .[. %in% rownames(seurat_object)] 
              
                # Create the heatmap, plots at most max_genes_plot_heatmap genes
                plot_heatmap(seurat_object, which = c("heatmap"), assay = assay,
                           cluster_column = cluster_column, name = paste0(output_dir, column), markers = features,
                           maxn_genes = max_genes_plot_heatmap)
            }
        } else {
          
            # Select all the features that are also present in the object
            features <- unlist(c(markers_df), use.names = FALSE) %>%
                .[. %in% rownames(seurat_object)]
  
            # Create the heatmap, plots at most max_genes_plot_heatmap genes
            plot_heatmap(seurat_object, which = c("heatmap"), assay = assay,
                       cluster_column = cluster_column, name = output_dir, markers = features,
                       maxn_genes = max_genes_plot_heatmap)
        }
    }
}

plot_markers_from_df_old <- function(seurat_object, 
    markers_location, 
    reduction_name = "umap_microglia_harmony_reduction", 
    name  = "",
    many_plot = TRUE, 
    feature_plot = FALSE, 
    extension_plot = ".png",
    heatmap_by_column = FALSE, 
    cluster_column = "microglia_clusters", 
    subplot_n = 9, 
    max_feature_plots = 11, 
    max_genes_plot_heatmap = 100,
    column_list = FALSE, 
    assay ="RNA") {
  
  
  # Messages 
  cat("Plotting markers from given df... \n")
  message(paste0("Parameters: markers_location: ", markers_location,
                 " - reduction_name: ", reduction_name,
                 " - name: ", name,
                 " - many_plot: ", many_plot,
                 " - feature_plot: ", feature_plot,
                 " - extension_plot: ", extension_plot,
                 " - heatmap_by_column: ", heatmap_by_column,
                 " - cluster_column: ", cluster_column,
                 " - subplot_n: ", subplot_n,
                 " - max_feature_plots: ", max_feature_plots,
                 " - max_genes_plot_heatmap: ", max_genes_plot_heatmap,
                 " - assay: ", assay,
                 " - column_list: ", column_list))

  # Set up output dir
  output_dir <- set_up_output(paste0(output_folder, "plot_markers_from_df_", name, "/"))
  
  # Load the markers
  markers_df <- readxl::read_excel(markers_location)
  
  # Saves the relative plots (either using many plots and a heatmap or the feature plots)
  if (feature_plot) {
    for (column in names(markers_df)) {
      
      # Skip if column list is defined and the column is in columnlist
      if (is.character(column_list) && column %in% column_list) {
        message("skipping: ", column)
        next
      }
      
      # Extract the features that are present in the seurat object for each column
      features <- markers_df %>%
        pull(column) %>%
        .[. %in% rownames(seurat_object)] 
      
      # Message
      message("Source: ", column, " - Features: ", paste(features, collapse = ", "))
      
      # Iterate thorugh the featres in batches of subplot_n and plot (max 11 plots)
      if (length(features) == 0) {
        message("skipping, not enough features")
        next
      }
      lapply(seq(1, min(length(features), subplot_n * max_feature_plots), by = subplot_n), function(start_index) {
        end_index <- min(start_index + subplot_n - 1, length(features))
        save_plot(FeaturePlot(seurat_object, features = features[start_index:end_index],
                              cols = c("lightgrey", "blue"),
                              reduction = reduction_name),
                  paste0(output_dir, "feature_plot_", column, "_", 
                         start_index, "_", end_index, 
                         extension_plot)
                  , x = 10, y = 10) 
      })
    }
  }
  
  
  # Heatmap
  if (many_plot) {
    if (heatmap_by_column || is.character(column_list)) {
      
      # Iterate through the columns
      for (column in names(markers_df)) {
        
        # Skip if no data is present
        if (is.character(column_list) && !any(column %in% column_list)) {
          message("skipping: ", column)
          next
        }
        
        # Extract the features that are present in the seurat object for each column
        features <- markers_df %>%
          pull(column) %>%
          .[. %in% rownames(seurat_object)] 
        
        # Create the heatmap, plots at most max_genes_plot_heatmap genes
        plot_heatmap(seurat_object, which = c("heatmap"), assay = assay,
                   cluster_column = cluster_column, name = paste0(output_dir, column), markers = features,
                   maxn_genes = max_genes_plot_heatmap)
      }
    } else {
      
      # Selects all the features that are also present in the object
      features <- unlist(c(markers_df), use.names = FALSE) %>%
        .[. %in% rownames(seurat_object)]

      # Create the heatmap, plots at most max_genes_plot_heatmap genes
      plot_heatmap(seurat_object, which = c("heatmap"), assay = assay,
                 cluster_column = cluster_column, name = output_dir, markers = features,
                 maxn_genes = max_genes_plot_heatmap)
      
    }
  }
}

#' Generate Plots for Publication from Seurat Object
#'
#' This function generates a variety of plots (bar plots, pie charts, feature plots, ridge plots, etc.)
#' from a Seurat object for use in publications. It can create cell type plots, subject-wise plots, and
#' other visualizations based on metadata from the Seurat object.
#'
#' @param seurat_object A Seurat object containing the single-cell RNA-seq data.
#' @param which A character vector specifying the types of plots to generate. Options include:
#'        "numberofcell_barplot", "numberofcell_pie_chart", "numberofcell_barplot_subject", 
#'        "numberofcell_pie_chart_subject", "numberofcell_pie_chart_cluster_subject", 
#'        "numberofcell_pie_chart_cluster_pathology", "feature_plots", "ridge_plots".
#' @param genes_to_plot A character vector specifying genes to plot in feature or ridge plots.
#' @param cluster_column A character string specifying the column name in the Seurat object's metadata 
#'        that contains the clustering information (default is "microglia_clusters").
#' @param extension_plot A character string specifying the file extension for the output plots (default is ".png").
#' @param name A character string that can be used to prefix the output plot filenames.
#' @param assay The assay to use for gene expression data (default is "RNA").
#' 
#' @return This function generates and saves the specified plots to a directory.
#' 
#' @details 
#' The function creates different types of visualizations based on the input Seurat object and metadata. 
#' It allows flexibility in terms of generating pie charts, bar plots, feature plots, and ridge plots. 
#' Additionally, the function saves these plots with user-defined names and extensions.
#'
#' Available plot types:
#' \itemize{
#'   \item \code{numberofcell_barplot}: Bar plot showing the number of cells per cluster.
#'   \item \code{numberofcell_pie_chart}: Pie chart showing the proportion of cells per cluster.
#'   \item \code{numberofcell_barplot_subject}: Bar plot showing the number of cells per cluster for each subject.
#'   \item \code{numberofcell_pie_chart_subject}: Pie chart showing the proportion of cells per cluster for each subject.
#'   \item \code{numberofcell_pie_chart_cluster_subject}: Pie chart showing the proportion of cells per cluster, with subjects in the legend.
#'   \item \code{numberofcell_pie_chart_cluster_pathology}: Pie chart showing the proportion of cells per pathology, with clusters in the legend.
#'   \item \code{feature_plots}: Feature plots for the specified genes, showing gene expression patterns across clusters.
#'   \item \code{ridge_plots}: Ridge plots for the specified genes, showing gene expression distribution across clusters.
#' }
#'
#' @examples
#' \dontrun{
#'   # Generate bar and pie charts for a Seurat object
#'   plots_for_paper(seurat_object, which = c("numberofcell_barplot", "numberofcell_pie_chart"))
#' 
#'   # Generate feature plots for specific genes
#'   plots_for_paper(seurat_object, which = "feature_plots", genes_to_plot = c("Gene1", "Gene2"))
#' }
#'
#' @export
plots_for_paper <- function(seurat_object, 
    which = c("numberofcell_barplot", 
      "numberofcell_pie_chart", 
      "numberofcell_barplot_subject", 
      "numberofcell_pie_chart_subject",
      "numberofcell_pie_chart_cluster_subject", 
      "numberofcell_pie_chart_cluster_pathology", 
      "feature_plots",
      "ridge_plots"), 
    genes_to_plot = FALSE, 
    cluster_column = "microglia_clusters", 
    extension_plot = ".png", 
    name = "", 
    assay = "RNA",
    reduction_name = "umap_harmony_microglia", 
    subplot_n = 9) {
  
  # Set up output dir
  output_dir <- set_up_output(paste0(output_folder, "plots_misc_", name, "/"), message)

  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  # Function definition to plot the pie plots
  pie_function <- function(name, subset_to_plot, categories_to_plot) {
    
    for (subset in unique(seurat_object@meta.data[subset_to_plot])) {
      
      # Calculate the frequency of each level
      level_counts <- table(seurat_object@meta.data[seurat_object@meta.data[subset_to_plot] == subset, categories_to_plot])
      
      # Create a data frame for ggplot
      pie_data <- data.frame(Level = names(level_counts),
                             Frequency = as.vector(level_counts))
      # Create the cake plot
      save_plot(pie_plot <- ggplot(pie_data, aes(x = "", y = Frequency, fill = Level)) +
                  geom_bar(width = 1, stat = "identity") +
                  coord_polar("y", start = 0) +
                  labs(title = paste0("Cell types ", subset), fill = "Level") +
                  theme_void() +
                  theme(legend.position = "right"),
                paste0(output_dir, name, "_pie_", subset_to_plot, extension_plot))
      
    }
  }
  
  # TODO: da fare function definiton for barplots
  bar_function <- function(name, bars, filling) {
    
    # Calculate the frequency of each level (for each cluster each subject)
    level_counts <- table(seurat_object@meta.data[[cluster_column]])
    subjects <- unique(seurat_object@meta.data$subject)
    names(subjects) <-  subjects
    
    # Etract cell type
    cell_types <- unique(seurat_object@meta.data[[cluster_column]])
    names(cell_types) <- cell_types
    level_counts <- lapply(cell_types, function(x) {
      
      counts <- lapply(subjects, function(y) {
        
        subject_counts <- nrow(seurat_object@meta.data[seurat_object@meta.data[[cluster_column]] == x & seurat_object@meta.data$subject == y, ])
        
        return(subject_counts)
      }) 
      
      return(counts)
    })
    frequency <- unlist(level_counts)
    
    # Create a data frame for ggplot
    pie_data <- data.frame(Bars = str_split_i(names(frequency), "\\.", 1),
                           Filling = str_split_i(names(frequency), "\\.", 2),
                           Frequency = as.vector(frequency))
    
    # Create the bar plot
    save_plot(ggplot(pie_data, aes(x = Bars, y = Frequency, fill = Filling)) +
                geom_bar(width = 1, stat = "identity", position = "fill") +
                labs(title = "Cell types", fill = "Subject") +
                theme_classic() +
                # scale_x_discrete(label = unique(pie_data$Cell)) +
                coord_flip() +
                theme(legend.position = "right"),
              paste0(output_dir, "barchart_subjects", extension_plot))
    
  }
  
  # Pie chart
  if ("numberofcell_pie_chart" %in% which) pie_function("pathology", "subject_pathology", cluster_column)
  # Cluster on legend, data from one subject
  if ("numberofcell_pie_chart_subject" %in% which) pie_function("subject", "subject", cluster_column)
  # Subject on legend, data from one cluster
  if ("numberofcell_pie_chart_cluster_subject" %in% which) pie_function("cluster_subject", cluster_column, "subject")
  # Pathology on legend, data from one cluster
  if ("numberofcell_pie_chart_cluster_pathology" %in% which) pie_function("cluster_subject_pathology", cluster_column, "subject_pathology")

  
  # Feature plots 
  if ("feature_plots" %in% which) {
    message("feature_plots")
    if (!reduction_name %in% names(seurat_object@reductions)) stop("plot_misc: the selected umap reduction does not exist")
    if (isFALSE(genes_to_plot)) stop("plot_misc: markers is a required_parameter for feature plots")
    genes <- genes_to_plot[genes_to_plot %in% Features(seurat_object[["RNA"]])]

    purrr::walk(seq(1, length(genes_to_plot), by = subplot_n), function(gene) {
       save_plot(
        FeaturePlot(seurat_object, 
          features = genes_to_plot[gene:ifelse(gene + subplot_n < length(genes_to_plot), gene + subplot_n, length(genes_to_plot))], 
          pt.size = 0.5, 
          reduction = reduction_name,
          cols = c("lightgrey", "#0026ff")), 
        paste0(output_dir, "feature_plots_genes_", 
          gene, "_", ifelse(gene + subplot_n < length(genes_to_plot), gene + subplot_n, length(genes_to_plot)), 
          extension_plot
          ),
        x = 10, y = 10
        )
    })
  }
  # Feature plots 
  if ("ridge_plots" %in% which) {
    message("ridge plots")
    if (isFALSE(genes_to_plot)) stop("plot_misc: markers is a required_parameter for ridge plots")

    genes <- genes_to_plot[genes_to_plot %in% Features(seurat_object[["RNA"]])]

    purrr::walk(seq(1, length(genes_to_plot), by = subplot_n), function(gene) {

      save_plot(
        RidgePlot(seurat_object, 
          features = genes_to_plot[gene:ifelse(gene + subplot_n < length(genes_to_plot), gene + subplot_n, length(genes_to_plot))],
          layer = "data",
          assay = assay
          ),
        paste0(output_dir, "ridge_plots_genes_", 
          gene, "_", ifelse(gene + subplot_n < length(genes_to_plot), gene + subplot_n, length(genes_to_plot)), 
          extension_plot
          ),
        x = 10, y = 10
        )
        # counts <- lapply(seq(10), function(gene) {
        #   if (genes_to_plot[gene] %in% rownames(seurat_object@assays[[assay]]$data))
        #     counts <- sum(seurat_object@assays[[assay]]$data[ 
        #       genes_to_plot[gene], # rownames -> genes
        #     ] > 0)
        # })
        # names(counts) <- genes_to_plot[seq(10)]

    })
  }
  
  # TODO: Pathology on legend, data from counts of one gene
  if ("numberofcell_pie_chart_gene_pathology" %in% which) {
    
    # One plot for each iteration (data coming from...)
    for (gene in genes_to_plot) {
      
      # if the gene total expression is 0 and it exists
      if (
        !(gene %in% rownames(seurat_object@assays[[assay]]$data)) 
        &&
        rowSums(t(seurat_object@assays[[assay]]$data[gene, ])) == 0
      ) 
        next
      
      # Iterate thrugh the pathologies to count how much that gene appears in each pathology
      # cant be an abs value because number o subjects for each path is different, it must be avg
      # problem: there are a lot fo nonzero elemtns, the mean does not work, the um does not work what can I use, if i can even use something)
      gene_count <- lapply(levels(factor(seurat_object@meta.data$subject_pathology)), function(pathology) {
        counts <- sum(seurat_object@assays[[assay]]$data[ 
          gene, # rownames -> genes
          rownames(seurat_object@meta.data[seurat_object@meta.data$subject_pathology == pathology, ]) # columnames -> cells
        ] 
        > 0) / # how many cells have a gene expression for this specific gene over 0? / tot cells of category
        length(rownames(seurat_object@meta.data[seurat_object@meta.data$subject_pathology == pathology, ]))
      })
      # for (pathology in seurat_object@meta.data$subject_pathology) {
      #   # Calculate the frequency of each level (... on legend)
      #   gene_counts <- mean(seurat_object@assays[[assay]]$data[ 
      #     gene, # rownames -> genes
      #     rownames(seurat_object@meta.data[seurat_object@meta.data$subject_pathology == pathology, ]) # columnames -> cells
      #   ])
      # }
      
      # Create a data frame for ggplot
      pie_data <- data.frame(Pathology = levels(factor(seurat_object@meta.data$subject_pathology)),
                             Frequency = unlist(gene_count))
      
      # Create the cake plot
      save_plot(ggplot(pie_data, aes(x = Pathology, y = Frequency, fill = Pathology)) +
                  geom_bar(stat = "identity") +
                  labs(title = paste0("Cell types ", gene), fill = "Pathology") +
                  theme(legend.position = "right"),
                paste0(output_dir, gene, "_bar_pathology", extension_plot))
      
    }
  }
  
  if ("numberofcell_pie_chart_cluster_pathology" %in% which) {
    
    for (cluster in unique(seurat_object@meta.data[[cluster_column]])) {
      
      # Calculate the frequency of each level
      level_counts <- table(seurat_object@meta.data[seurat_object@meta.data[[cluster_column]] == cluster, "subject_pathology"])
      
      # Create a data frame for ggplot
      pie_data <- data.frame(Level = names(level_counts),
                             Frequency = as.vector(level_counts))
      
      # Create the cake plot
      save_plot(pie_plot <- ggplot(pie_data, aes(x = "", y = Frequency, fill = Level)) +
                  geom_bar(width = 1, stat = "identity") +
                  coord_polar("y", start = 0) +
                  labs(title = paste0("Cell types ", cluster), fill = "Level") +
                  theme_void() +
                  theme(legend.position = "right"),
                paste0(output_dir, cluster, "_pie_pathology", extension_plot))
      
    }
    
  }
  # Bar plot
  if ("numberofcell_barplot" %in% which) {
    
    # Calculate the frequency of each level
    level_counts <- table(seurat_object@meta.data[[cluster_column]])
    pathologies <- unique(seurat_object@meta.data$subject_pathology)
    names(pathologies) <-  pathologies
    
    # Etract cell type
    cell_types <- unique(seurat_object@meta.data[[cluster_column]])
    names(cell_types) <- cell_types
    level_counts <- lapply(cell_types, function(x) {
      
      counts <- lapply(pathologies, function(y) {
        pathology_counts <- nrow(seurat_object@meta.data[
          seurat_object@meta.data[[cluster_column]] == x &
          seurat_object@meta.data$subject_pathology == y, 
        ])   
      }) 
      
      return(counts)
    })
    frequency <- unlist(level_counts)
    
    # Create a data frame for ggplot
    pie_data <- data.frame(Cell_type = str_split_i(names(frequency), "\\.", 1),
                           Pathology = str_split_i(names(frequency), "\\.", 2),
                           Frequency = as.vector(frequency))
    
    # Create the bar plot
    save_plot(ggplot(pie_data, aes(x = Cell_type, y = Frequency, fill = Pathology)) +
                geom_bar(width = 1, stat = "identity", position = "fill") +
                labs(title = "Cell types", fill = "Pathology") +
                theme_classic() +
                # scale_x_discrete(label = unique(pie_data$Cell)) +
                coord_flip() +
                theme(legend.position = "right"),
              paste0(output_dir, "barchart", extension_plot))

  
  }
  # Bar plot by subject
  if ("numberofcell_barplot_subject" %in% which) {
    
    
    # Calculate the frequency of each level (for each cluster each subject)
    level_counts <- table(seurat_object@meta.data[[cluster_column]])
    subjects <- unique(seurat_object@meta.data$subject)
    names(subjects) <-  subjects
    
    # Etract cell type
    cell_types <- unique(seurat_object@meta.data[[cluster_column]])
    names(cell_types) <- cell_types
    level_counts <- lapply(cell_types, function(x) {
      
      counts <- lapply(subjects, function(y) {
        
        subject_counts <- nrow(seurat_object@meta.data[seurat_object@meta.data[[cluster_column]] == x & seurat_object@meta.data$subject == y, ])
        
        return(subject_counts)
      }) 
      
      return(counts)
    })
    frequency <- unlist(level_counts)
    
    # Create a data frame for ggplot
    pie_data <- data.frame(Cell_type = str_split_i(names(frequency), "\\.", 1),
                           Subject = str_split_i(names(frequency), "\\.", 2),
                           Frequency = as.vector(frequency))
    
    # Create the bar plot
    save_plot(ggplot(pie_data, aes(x = Cell_type, y = Frequency, fill = Subject)) +
                geom_bar(width = 1, stat = "identity", position = "fill") +
                labs(title = "Cell types", fill = "Subject") +
                theme_classic() +
                # scale_x_discrete(label = unique(pie_data$Cell)) +
                coord_flip() +
                theme(legend.position = "right"),
              paste0(output_dir, "barchart_subjects", extension_plot))
    
    
  }
  # TODO: Dot plot
  if ("dotplot_plot" %in% which) {
    
  }

}
#### PCA and INTEGRATION ####

# Standard preprocessing pipeline (normalize, feature selection variable features, scale data)
#' Preprocess a Seurat Object: Normalization, Variable Feature Selection, and Scaling
#'
#' This function performs a standard preprocessing pipeline on a Seurat object or a list of Seurat objects. 
#' The pipeline includes data normalization, variable feature selection, and data scaling. Optionally, it can save 
#' the processed Seurat object and generate a plot of variable features.
#'
#' @param seurat_object A Seurat object or a list of Seurat objects to be preprocessed.
#' @param save A logical value indicating whether to save the processed Seurat object(s) to a file. 
#'        Default is \code{FALSE}.
#' @param vf_plot A logical value indicating whether to generate and save a plot of the top variable features. 
#'        Default is \code{TRUE}.
#' @param normalization A logical value indicating whether to perform data normalization. 
#'        Default is \code{TRUE}.
#' @param variable_features A logical value indicating whether to perform variable feature selection. 
#'        Default is \code{TRUE}.
#' @param scaling A logical value indicating whether to perform data scaling. 
#'        Default is \code{TRUE}.
#' @param check_for_norm A logical value indicating whether to check for normalization when scaling the data. 
#'        Default is \code{TRUE}.
#'
#' @return The processed Seurat object or list of Seurat objects.
#'
#' @details This function applies a series of preprocessing steps to a Seurat object or a list of Seurat objects. 
#' For a single object, the steps include normalization, selection of the top 2000 variable features, and scaling 
#' the data. For a list of objects, normalization and variable feature selection are applied to each object in the list.
#' If \code{vf_plot} is \code{TRUE}, a plot of the top variable features is saved. The processed object(s) can also 
#' be saved to a file if \code{save} is \code{TRUE}.
#'
#' @examples
#' # Preprocess a single Seurat object
#' processed_obj <- normalization_and_scaling(seurat_object = my_seurat, save = TRUE)
#'
#' # Preprocess a list of Seurat objects
#' processed_list <- normalization_and_scaling(seurat_object = list(seurat_obj1, seurat_obj2), save = FALSE)
#'
#' @export
normalization_and_scaling <- function(seurat_object, save = FALSE, vf_plot = TRUE, normalization = TRUE,
                                       variable_features = TRUE, scaling = TRUE, check_for_norm = TRUE) {
  
  # TODO: rendere questa funzione generic?
  # If the input is a list of objects and not a single object
  if (is.list(seurat_object)) {
    cat("Runing preprocessing on list, normalization and variable features... \n")

    for (i in seq_along(length(seurat_object))){
      
      # normalization, feature selection and scaling
      seurat_object[[i]] <- NormalizeData(seurat_object[[i]])
      seurat_object[[i]] <- FindVariableFeatures(seurat_object[[i]], selection.method = "vst", 
                                                 nfeatures = 2000, verbose = FALSE)
      
    }
    return(seurat_object)
  } 
  
  # If it is a single object
  else {
    cat("Runing normalization, scaling and variable features \n")
    
    # info message
    message(paste0("Parameters: normalization: ", normalization,
                 " - variable features: ", variable_features,
                 " - scaling: ", scaling,
                 " - variable features plot: ", vf_plot,
                 " - save: ", save,
                 " - check_for_norm: ", check_for_norm))
    
    # Normalization, feature selection and scaling
    if (normalization) seurat_object <- NormalizeData(seurat_object, verbose = FALSE)
    if (variable_features) seurat_object <- FindVariableFeatures(seurat_object, selection.method  = "vst", 
                                                                 nfeatures = 2000, verbose = FALSE)

    if (!check_for_norm) LayerData(seurat_object, assay = "counts", layer = "data") <-  
        GetAssayData(object = seurat_object, assay = "counts", layer = "counts")
    
    if (scaling) seurat_object <- ScaleData(seurat_object)#, check_for_norm=check_for_norm)
    
    if (vf_plot) save_plot(LabelPoints(VariableFeaturePlot(seurat_object),
                                    points = head(VariableFeatures(seurat_object)),
                                    repel = TRUE), 
                        paste0(output_folder, "quality_control/VariableFeaturePlot", extension_plot))
    
    # Save object
    if (save) saveRDS(seurat_object, file = paste0(output_folder, "after_normalization_and_varfeatures.rds"))
    
    return(seurat_object)
  }
  
}

# Runs pca and elbow plot if needed
#' Run PCA on a Seurat Object and Generate Elbow Plot
#'
#' This function performs Principal Component Analysis (PCA) on a Seurat object and optionally generates 
#' an elbow plot to help determine the number of significant principal components. The processed Seurat object 
#' can also be saved to a file.
#'
#' @param seurat_object A Seurat object containing the single-cell RNA-seq data.
#' @param plot A logical value indicating whether to generate and save an elbow plot. 
#'        Default is \code{FALSE}.
#' @param save A logical value indicating whether to save the Seurat object after running PCA. 
#'        Default is \code{FALSE}.
#' @param assay A character string specifying which assay to use for PCA. If \code{"default"}, the default assay is used. 
#'        Default is \code{"default"}.
#' @param reduction_name A character string specifying the name to assign to the PCA reduction. 
#'        Default is \code{"pca"}.
#' @param extension_plot A character string specifying the file extension for saving the elbow plot. 
#'        Default is \code{".png"}.
#'
#' @return The Seurat object with PCA results added.
#'
#' @details This function runs PCA on the specified Seurat object using either the default or specified assay. 
#' The PCA results are stored in the Seurat object under the name given by \code{reduction_name}. 
#' If \code{plot} is \code{TRUE}, an elbow plot is generated and saved. If \code{save} is \code{TRUE}, the updated 
#' Seurat object is saved to a file named \code{"after_pca.rds"} in the output folder.
#'
#' @examples
#' # Run PCA and save the Seurat object
#' seurat_obj <- run_pca(seurat_object = my_seurat, plot = TRUE, save = TRUE)
#'
#' @export
run_pca <- function(seurat_object, plot = FALSE, save = FALSE, assay = "default", reduction_name = "pca", extension_plot = ".png") {
  
  # Parameters
  cat("Running PCA... \n")
  message(paste0("reduction name: ", reduction_name,
               " - assay: ", assay,
               " - save plot: ", plot,
               " - save: ", save))
  
  # Run PCA
  if (assay == "default") seurat_object <- RunPCA(seurat_object, reduction.name = reduction_name, 
                                                  reduction.key = reduction_name, verbose = TRUE)
  else  seurat_object <- RunPCA(seurat_object, assay = assay, reduction.name = reduction_name, 
                            reduction.key = reduction_name)
  
  # Save plots and object
  if (plot) save_plot(ElbowPlot(seurat_object), paste0(output_folder, "elbow_plot_", assay, "_", reduction_name, extension_plot))
  if (save) saveRDS(seurat_object, file = paste0(output_folder, "after_pca.rds"))
  
  return(seurat_object)
}

# Layer integration using seurat
#' Integrate Layers of a Seurat Object Using Various Methods
#'
#' This function performs integration of layers within a Seurat object using various integration methods, 
#' such as Harmony, CCA, and PCA-based approaches. It can handle both single Seurat objects and lists of Seurat objects.
#'
#' @param seurat_object A Seurat object or a list of Seurat objects to be integrated.
#' @param save A logical value indicating whether to save the integrated Seurat object to a file. 
#'        Default is \code{FALSE}.
#' @param assay A character string specifying the assay to use for integration. Default is \code{"RNA"}.
#' @param new_reduction A character string specifying the name of the new reduction to add after integration. 
#'        Default is \code{"integrated.harmony"}.
#' @param reduction_method A character string specifying the method for integration. Options include \code{"harmony"},
#'        \code{"harmony_seurat"}, \code{"cca"}, \code{"jpca"}, and \code{"rpca"}. Default is \code{"harmony"}.
#' @param make_default A logical value indicating whether to make the new reduction the default reduction. 
#'        Default is \code{FALSE}.
#' @param orig_reduction A character string specifying the original reduction to use for integration. Default is \code{"pca"}.
#' @param new_assay A character string specifying a new assay to use. If \code{NA}, the original assay is used. Default is \code{NA}.
#'
#' @return The Seurat object with integrated layers.
#'
#' @details This function integrates layers of a Seurat object using the specified method. For a list of Seurat objects, 
#' integration is performed using \code{FindIntegrationAnchors} and \code{IntegrateData}. For a single Seurat object, 
#' integration is performed using one of the specified methods. The integrated layers are joined, and the updated 
#' Seurat object is returned. The function can also save the integrated Seurat object to a file if \code{save} is \code{TRUE}.
#'
#' @examples
#' # Integrate layers using Harmony method and save the object
#' integrated_obj <- layer_integration(seurat_object = my_seurat, reduction_method = "harmony", save = TRUE)
#'
#' # Integrate a list of Seurat objects
#' integrated_list <- layer_integration(seurat_object = list(seurat_obj1, seurat_obj2), reduction_method = "cca")
#'
#' @export
layer_integration <- function(seurat_object, save = FALSE, assay = "RNA",
                              new_reduction = "integrated.harmony", reduction_method = "harmony",
                              make_default = FALSE,  orig_reduction = "pca", new_assay = FALSE) {
  
  # Messages
  cat("Integrating layers...\n")
  message(paste0("Parameters: save: ", save,
                 " - assay: ", assay,
                 " - new_reduction: ", new_reduction,
                 " - reduction_method: ", reduction_method,
                 " - make_default: ", make_default,
                 " - orig_reduction: ", orig_reduction,
                 " - new_assay: ", new_assay))
  
  # Assign variables
  if (isFALSE(new_assay)) new_assay <- assay
  
  # If the input is a list, then the object are merged with integratedata
  if (is.list(seurat_object)) {
    cat("Runing integration on a list, normalization and variable features\n")
    
    anchors <- FindIntegrationAnchors(object.list = seurat_object)
    seurat_object <- IntegrateData(anchorset = anchors, normalization.method = "LogNormalize", )
    
    return(seurat_object)
  }
  
  # perform layer integration with selected method
  if (reduction_method == "harmony_seurat") {
    seurat_object  <- IntegrateLayers(
      object = seurat_object, method = HarmonyIntegration,
      orig.reduction = orig_reduction, new.reduction = new_reduction,
      verbose = TRUE
    )
  } else if (reduction_method == "cca") {
    seurat_object  <- IntegrateLayers(
      object = seurat_object, method = CCAIntegration,
      orig.reduction = orig_reduction, new.reduction = new_reduction,
      verbose = FALSE
    )
  } else if (reduction_method == "jpca") {
    seurat_object  <- IntegrateLayers(
      object = seurat_object, method = JointPCAIntegration,
      orig.reduction = orig_reduction, new.reduction = new_reduction,
      verbose = FALSE
    )
  } else if (reduction_method == "rpca") {
    seurat_object  <- IntegrateLayers(
      object = seurat_object, method = RPCAIntegration,
      orig.reduction = orig_reduction, new.reduction = new_reduction,
      verbose = FALSE
    )
  }  else if (reduction_method == "harmony") {
    check_packages("harmony")
    seurat_object <- RunHarmony(seurat_object, 
                                group.by.vars = c("subject"), 
                                reduction = orig_reduction, assay.use = assay, 
                                reduction.save = new_reduction)
  } else {
    stop("invalid parameter: reduction_method")
  }
  
  seurat_object <- JoinLayers(seurat_object)
  # Save object if needed
  if (save) saveRDS(seurat_object, file = paste0(output_folder, "after_layer_integration.rds"))
  
  return(seurat_object)
}

#### ANNOTATION ####

#' ScType reformatted (not tested) + docstring
#'
#' This function annotates a Seurat object using the scType method, which scores cells based on gene sets specific to cell types.
#'
#' @param seurat_object A Seurat object to be annotated.
#' @param assay The assay to be used from the Seurat object. Default is "RNA".
#' @param assignment_name The name for the new metadata column where annotations will be stored. Default is "annotation_scType".
#' @param clusters_column The metadata column in the Seurat object that contains clustering information. Default is "seurat_clusters".
#' @param name An optional name prefix for output files. Default is an empty string.
#'
#' @return A Seurat object with added annotations.
#' @export
#' 
scType_annotation <- function(seurat_object, 
                              assay = "RNA", 
                              assignment_name = "annotation_scType",
                              clusters_column = "seurat_clusters", 
                              name = "annotation_scType") {
    # Message
    cat("Running scType annotation...\n")
    message(paste0("Parameters: assay: ", assay,
                  " - assignment_name: ", assignment_name, 
                  " - clusters_column: ", clusters_column, 
                  " - name: ", name))
    
    # Dependencies
    library("HGNChelper")

    # Create directory
    output_dir <- paste0(output_folder, name, "/")
    if (!dir.exists(output_dir)) dir.create(output_dir)
    
    # Load gene set preparation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
    # Load cell type annotation function
    source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
    
    # DB file
    db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
    tissue <- "Brain" # e.g. Immune system, Pancreas, Liver, Eye, Kidney, Brain, etc.
    
    # Prepare gene sets
    gs_list <- gene_sets_prepare(db_, tissue)
    message("scType: data prepared")
    
    # Run scType
    es_max <- sctype_score(scRNAseqData = seurat_object[[assay]]$scale.data, 
                           scaled = TRUE, 
                           gs = gs_list$gs_positive, 
                           gs2 = gs_list$gs_negative)
    
    # Merge by cluster
    cl_results <- do.call("rbind", lapply(unique(seurat_object@meta.data[[clusters_column]]), function(cl) {
        es_max_cl <- sort(rowSums(es_max[, rownames(seurat_object@meta.data[seurat_object@meta.data[[clusters_column]] == cl, ])]), decreasing = TRUE)
        head(data.frame(cluster = cl, type = names(es_max_cl), scores = es_max_cl, ncells = sum(seurat_object@meta.data[[clusters_column]] == cl)), 10)
    }))

    sctype_scores <- cl_results[order(cl_results$cluster, -cl_results$scores), ]  # Sort by cluster and scores (descending)
    sctype_scores <- sctype_scores[!duplicated(sctype_scores$cluster), ]          # Keep only top scoring row per cluster
    
    # Set low-confident (low ScType score) clusters to "unknown"
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 6] <- "Unknown"
    
    # Save results in table
    openxlsx::write.xlsx(sctype_scores, paste0(output_dir, name, "results.xlsx"), sheetName = "results", append = TRUE)
    message(paste0("scType: annotation results saved in: ", output_dir))
    
    # Subset annotation results
    cluster_mapping <- as.data.frame(sctype_scores[c("cluster", "type")])
    
    # Copy metadata column
    metadata_assignments <- seurat_object@meta.data[[clusters_column]]
    
    # Update metadata column and add again to the Seurat object
    matching_rows <- match(metadata_assignments, cluster_mapping[["cluster"]])
    metadata_assignments <- cluster_mapping[["type"]][matching_rows]
    seurat_object <- AddMetaData(
        object = seurat_object,
        metadata = metadata_assignments,
        col.name = assignment_name
    )
    
    message("scType: metadata updated")
    
    return(seurat_object)
}

# Correct an annotation. changing specific values of a column
correct_annotation <- function(seurat_object, to_correct, annotation_column, 
                               new_annotation_column = NA, cluster_column = "seurat_clusters") {

  # if new name not given create it by adding corrected to the annotation column
  if (is.na(new_annotation_column)) new_annotation_column <- paste0(annotation_column, "_corrected")
  
  # copy te annotation column
  correction <- as.character(seurat_object@meta.data[[annotation_column]])
  
  # modify the name according to the information to correct
  for (item in names(to_correct)) {
    correction[seurat_object@meta.data[[cluster_column]] == item] <- to_correct[[item]]
  }
  
  # Add the new column to the seurat object
  seurat_object <- AddMetaData(
    object = seurat_object,
    metadata = correction,
    col.name = new_annotation_column
  )
  return(seurat_object)
}

#### OTHER ####

#' Perform Trajectory Analysis Using Monocle 3
#'
#' This function performs trajectory analysis on a Seurat object using the Monocle 3 package. It sets up the necessary environment,
#' converts the Seurat object to a Monocle CellDataSet object, performs clustering, and computes and visualizes the trajectory of cells.
#'
#' @param seurat_object A Seurat object containing the single-cell RNA-seq data. This object is converted to a Monocle CellDataSet for trajectory analysis.
#' @param extension_plot A character string specifying the file extension for saving the plots. Default is \code{".png"}.
#' @param name A character string specifying a name prefix for the output files. Default is \code{"test"}.
#'
#' @details
#' This function installs and loads necessary packages, including Monocle 3, and performs the following steps:
#' \itemize{
#'   \item Sets up the output directory for saving results.
#'   \item Converts the Seurat object to a Monocle CellDataSet object.
#'   \item Preprocesses the data by adding PCA and UMAP embeddings.
#'   \item Performs clustering on the cells.
#'   \item Computes the trajectory graph and orders the cells based on pseudotime.
#'   \item Saves plots of the clustering results and trajectory analysis.
#' }
#'
#' The function requires the installation of Monocle 3 and other dependencies, and uses functions from the SeuratWrappers package for conversion.
#'
#' @examples
#' # Perform trajectory analysis on a Seurat object and save the results
#' trajectory_analysis(seurat_object = my_seurat, extension_plot = ".png", name = "sample")
#'
#' @export
trajectory_analysis <- function(seurat_object, extension_plot = ".png", name = "test") {
  
  # Messages 
  cat("Performing trajectory analysis with monocle... \n")
  message(paste0("Parameters: name: ", name,
                 " - extension_plot: ", extension_plot))
  
  # Dependencies
  check_packages(c("BiocGenerics", "DelayedArray", "DelayedMatrixStats",
                       "limma", "lme4", "S4Vectors", "SingleCellExperiment",
                       "SummarizedExperiment", "batchelor", "HDF5Array",
                       "terra", "ggrastr", "devtools"))
  # Monocle
  devtools::install_github("cole-trapnell-lab/monocle3")
  library(monocle3)
  

  # Set up output dir
  output_dir <- set_up_output(paste0(output_folder, "monocle", name, "/"))
  
  # Seuart wraper 
  # devtools::install_github("satijalab/seurat-wrappers")
  remotes::install_github("satijalab/seurat-wrappers@community-vignette")
  
  # Automatic object creation
  # cds <- SeuratWrappers::as.cell_data_set(seurat_object)
  
  # Manual object creation
  cds <- new_cell_data_set(seurat_object@assays$RNA$data,
                              cell_metadata = seurat_object@meta.data)
  
  # Preprocessing
  reducedDim(cds, type = "PCA") <- seurat_object@reductions$pca@cell.embeddings 
  cds@reduce_dim_aux$prop_var_expl <- seurat_object@reductions$pca@stdev
  cds@int_colData@listData$reducedDims$UMAP <- seurat_object@reductions$umap_microglia_harmony_reduction@cell.embeddings
  
  # Clustering 
  cds <- cluster_cells(cds)
  
  # Save the plot of the identified clsters
  save_plot(plot_cells(cds, show_trajectory_graph = FALSE), 
            paste0(output_dir, "clustering", extension_plot))
  
  # Compute the connection graph
  cds <- learn_graph(cds, use_partition = FALSE)
  cds <- order_cells(cds)
  
  # Plot the results of the trajecotry analysis
  save_plot(plot_cells(cds, color_cells_by = "pseudotime"), 
            paste0(output_dir, "pseudotime_traj_graph", extension_plot))
  save_plot(plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = FALSE), 
            paste0(output_dir, "pseudotime", extension_plot))
            
  detach("package:openxlsx", unload=TRUE)
}

#' Run CellChat Analysis on Seurat Object
#'
#' This function performs a CellChat analysis on a Seurat object. It allows you to analyze specific clusters and subjects, save or load a precomputed CellChat object, and create plots related to the analysis.
#'
#' @param seurat_object A Seurat object containing single-cell RNA-seq data. This is the primary input for the CellChat analysis.
#' @param cluster_column Character. The name of the column in the Seurat object metadata that contains the cluster assignments. Default is "microglia_clusters".
#' @param clusters_to_analyze Character vector. A vector of cluster names or IDs to analyze. If left empty (`c()`), all clusters will be analyzed.
#' @param subject_column Character. The name of the column in the Seurat object metadata that contains the subject names. Default is "subject".
#' @param name Character. A name to identify the analysis, used for labeling and saving output files. Default is "test".
#' @param save_object Logical. If `TRUE`, the CellChat object will be saved to disk after analysis. Default is `TRUE`.
#' @param load_object Logical. If `TRUE`, a precomputed CellChat object will be loaded from a file instead of running the analysis. Default is `FALSE`.
#' @param obj_name Character. The filename of the saved or loaded CellChat object. Default is "cellchat_object.rds".
#' @param extension_plot Character. The file extension for saved plots (e.g., ".png", ".pdf"). Default is ".png".
#'
#' @return A CellChat object containing the results of the analysis. If `save_object` is `TRUE`, the object will also be saved to disk as `obj_name`.
#' 
#' @details 
#' This function performs a CellChat analysis to study cell-cell communication within specific clusters in a Seurat object. It allows customization of which clusters and subjects to analyze, whether to save or load a precomputed object, and the format of the output plots.
#' Documentation following: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
#'
#' @examples
#' \dontrun{
#'   # Example usage:
#'   cellchat_function(seurat_object, 
#'                     cluster_column = "seurat_clusters", 
#'                     clusters_to_analyze = c("Cluster_1", "Cluster_2"), 
#'                     subject_column = "subject", 
#'                     name = "CellChat_analysis", 
#'                     save_object = TRUE, 
#'                     load_object = FALSE, 
#'                     obj_name = "cellchat_results.rds", 
#'                     extension_plot = ".png")
#' }
#'
#' @export
cellchat_function <- function(seurat_object, cluster_column = "microglia_clusters", 
                              clusters_to_analyze = c(), subject_column = "subject",
                              name = "test", save_object = TRUE, load_object = FALSE,
                              obj_name = "cellchat_object.rds", extension_plot = ".png") {

  #devtools::install_github("jinworks/CellChat")
  #devtools::install_github("jokergoo/circlize")
  #devtools::install_github("jokergoo/ComplexHeatmap")
  #devtools::install_github("immunogenomics/presto")
  library("CellChat")
  # library("patchwork")
  # options(stringsAsFactors = FALSE)

  # Set op output dir
  output_dir <- set_up_output(file.path(output_folder, paste0("cellchat_", name, "/")), message)
  if (!dir.exists(file.path(output_dir, "circle_plots"))) dir.create(file.path(output_dir, "circle_plots"))
  

  
  ## Main ----
  main <- function() {
    

  }
  
  ## Functions ----
  # prepare the seurat object for the subsequent cell communication analysis
  prepare_data <- function() {
    
    # Preparing data
    assign("data_input", 
      seurat_object[["RNA"]]$data,
      envir = parent.frame()
    ) 
    seurat_object@meta.data[[cluster_column]]
    Idents(seurat_object) <- cluster_column
    
    # Add 1 to all idents if 0 is present 
    if ("0" %in% as.character(unique(seurat_object@meta.data[[cluster_column]]))) {
      
      new_idents <- as.character(
        as.numeric(
          levels(
            seurat_object@meta.data[[cluster_column]]))
        + 1)
      
      levels(seurat_object@meta.data[[cluster_column]]) <- new_idents
    }

    # update the cluster column in the 
    Idents(seurat_object) <- cluster_column
    labels <- Idents(seurat_object)
    
    # Create a dataframe of the cell labels 
    assign("meta", 
      data.frame(labels = factor(as.character(labels)),
        samples = factor(seurat_object@meta.data[[subject_column]]), 
        row.names = names(labels)),
      envir = parent.frame()
    )

  }
  
  # main loop for computing the cellchat object
  # source: https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
  compute_cellchat <- function() {
    # class: https://rdrr.io/github/sqjin/CellChat/man/CellChat-class.html
    cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")
    cell_chat_db <- CellChatDB.human 
    
    # For visual inspection
    # showDatabaseCategory(cell_chat_db)
    # dplyr::glimpse(cell_chat_db$interaction)
    
    # use a subset of cell_chat_db for cell-cell communication analysis (important, which one?)
    # set the used database in the object
    CellChatDB.use  <- subsetDB(cell_chat_db, search = "Secreted Signaling", key = "annotation") # use Secreted Signaling
    cellchat@DB <- CellChatDB.use

    # subset the expression data of signaling genes for saving computation cost
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database 
    
    # finds genes that are expressed in the dataset (subsetData)
    # future::plan("multisession", workers = 4) # do parallel
    cellchat <- identifyOverExpressedGenes(cellchat) # Identify over-expressed signaling genes associated with each cell group
    cellchat <- identifyOverExpressedInteractions(cellchat) # Identify over-expressed ligand-receptor interactions (pairs) within the used cell_chat_db

    # project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
    # cellchat <- projectData(cellchat, PPI.human)
    
    # Compute communication 
    cellchat <- computeCommunProb(cellchat, type = "triMean")
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    cellchat <- computeCommunProbPathway(cellchat) # Compute the communication probability/strength between any interacting cell groups

    # when everything is computed 
    # I can also use assign instead of the parent env
    # Calculate the aggregated network by counting the number of links or summarizing the communication probability
    cellchat <- aggregateNet(cellchat) 
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
    
    assign("cellchat", cellchat, envir = parent.frame())
    
    if (save_object) saveRDS(cellchat, file.path(output_dir, obj_name))
  }

  # function to plot all the results
  plot_results <- function() {
    
    plot_function <- function(p_function, name) {

      if (extension_plot == ".png") png(file = file.path(output_dir, paste0(name, ".png")))
      else if (extension_plot == ".svg") svg(file = file.path(output_dir, paste0(name, ".svg")))
      
      # Plot
      p_function
      dev.off()
      message("plot saved in:", file.path(output_dir,  paste0(name, extension_plot)))
    }
    
    
    plots_group <-  function() {
      
      group_size <- as.numeric(table(cellchat@idents))
      
      # Number of interactions
      # Strength of interactions
      png(file = file.path(output_dir,  paste0("number_of_interactions", extension_plot)))
      par(mfrow = c(1, 2), xpd = TRUE)
      # label edge and weight scale are for appearance 
      # https://rdrr.io/github/sqjin/CellChat/man/netVisual_circle.html
      netVisual_circle(cellchat@net$count, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
      netVisual_circle(cellchat@net$weight, vertex.weight = group_size, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")
      dev.off()
      
      message("plot saved in:", file.path(output_dir, paste0("number_of_interactions", extension_plot)))
      # png(file=paste0(output_dir, "interaction_strength", extension_plot))
      
      # dev.off()
      
      mat <- cellchat@net$weight
      
      
      png(file = file.path(output_dir, paste0("other_plot", extension_plot)))
      par(mfrow = c(3, 4), xpd = TRUE)
      for (i in seq_along(nrow(mat))) {
        mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
        mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, vertex.weight = group_size, weight.scale = TRUE, edge.weight.max = max(mat), title.name = rownames(mat)[i])
      }
      dev.off()
      
      message("plot saved in: ", file.path(output_dir, paste0("other_plot", extension_plot)))
      
      for (pathway in cellchat@netP$pathways) {
        pair_lr_cxcl <- extractEnrichedLR(cellchat, signaling = pathways_show, geneLR.return = FALSE)
        for (interaction_name in pair_lr_cxcl) {
          print(interaction_name)
          print(paste0("circle_", pathway, "_", interaction_name))
          plot_f(netVisual_individual(cellchat, signaling = pathways_show, pairLR.use = lr_show, layout = "circle"),
                 paste0("circle_", pathway, "_", interaction_name))
        }
      }
    }
    
    # plots_group()
    
    # signaling
    # pairLR
    # vertex_receiver
    # sources.use
    
    if (FALSE) {
      pathways_show <- cellchat@netP$pathways[[1]]
      pair_lr_cxcl <- extractEnrichedLR(cellchat, signaling = pathways_show, geneLR.return = FALSE)
      lr_show <- pair_lr_cxcl[1, ]
    }

    # pathway
    for (pathways_show in cellchat@netP$pathways) {

      pair_lr_cxcl <- extractEnrichedLR(cellchat, signaling = pathways_show, geneLR.return = FALSE)
      
      save_plot(netAnalysis_contribution(cellchat, signaling = pathways_show),
                    file.path(output_dir, paste0("contributions_", pathways_show, "", extension_plot)))
      
      for (lr_show in pair_lr_cxcl) {
        
        # LR show - ligand receptor
        # pathway show ->
        plot_function(netVisual_individual(cellchat, signaling = pathways_show, pairLR.use = lr_show, layout = "circle"),
                      paste0("circle_plots/", pathways_show, "_", lr_show))
        
        #
        #plot_function(netVisual_individual(cellchat, signaling = pathways_show, pairLR.use = lr_show, layout = "chord"),
        #              paste0("chord_plot_", pathways_show, "_", lr_show))

        #for (vertex_receiver in seq_along(levels(cellchat@idents))) {
        # isnide this loop chord plots or circle plots
        # plot_function(netVisual_aggregate(cellchat, signaling = pathways_show,  vertex.receiver = vertex_receiver), 
        #              paste0("circle_plot_pathway_vertex_", pathways_show, "_", lr_show, "_", vertex_receiver,"", extension_plot))
        #}
        # outside bubble plot o heatmap 
        
        
        #
        # (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
        # sources use -> a vector giving the index or the name of source cell groups
        # target use -> a vector giving the index or the name of target cell groups.
        
      } # questo va fatto per microglia 
      
      vertex_receiver <- seq(levels(cellchat@idents))
      plot_function(netVisual_aggregate(cellchat, signaling = pathways_show,  vertex.receiver = vertex_receiver), 
                    paste0("circle_plots/pathway_vertex_", pathways_show, "_", vertex_receiver))
      plot_function(netVisual_aggregate(cellchat, signaling = pathways_show), 
                    paste0("circle_plots/pathway_", pathways_show))

      # 

      name <- paste0("heatmap_plot_", pathways_show)
      if (extension_plot == ".png") png(file = file.path(output_dir, paste0(name, ".png")))
      else if (extension_plot == ".svg") svg(file = file.path(output_dir, paste0(name, ".svg")))
      
      # Plot
      netAnalysis_signalingRole_network(cellchat, signaling = pathways_show, width = 8, height = 2.5, font.size = 10)
      dev.off()
      message("plot saved in:", file.path(output_dir,  paste0(name, extension_plot)))
      # browser()
    }

    save_plot(netVisual_bubble(cellchat, sources.use = seq(levels(cellchat@idents)), targets.use = seq(levels(cellchat@idents)), remove.isolate = FALSE),
              paste0(output_dir, "/bubble_plot", extension_plot))
    
      # poi per altro invece va fatto la stessa cosa, devo solo costruire solo 3 clusters, due microglia infiammata + oligo. oppure 

    for (vertex_receiver in seq_along(levels(cellchat@idents))) {
      # non riesco a capire cosa ci andrebbe come input
      # save_plot(netVisual_hierarchy1(cellchat@data.signaling, vertex_receiver),
      #          paste0(output_dir, "hierarchy_ident_", vertex_receiver, "", extension_plot))
    }



    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")

    name <- "heatmap_signaling_role"
    if (extension_plot == ".png") png(file = file.path(output_dir, paste0(name, ".png")))
    else if (extension_plot == ".svg") svg(file = file.path(output_dir, paste0(name, ".svg")))
    
    # Plot
    ht1 + ht2
    dev.off()
    message("plot saved in:", file.path(output_dir,  paste0(name, extension_plot)))
    browser()
    # svg(file = file.path(output_dir, paste0("heatmap_signaling_role", ".svg")))
      
    # # Plot
    # p_function
    # dev.off()
    # message("plot saved in:", file.path(output_dir,  paste0(name, extension_plot)))
    # To see the connections between idents for all the pathways
    # target use -> microglia
    # sources use -> astrocytes  
    
    for (use in seq(levels(cellchat@idents))) {

      plot_function(netVisual_chord_gene(cellchat, sources.use = use, targets.use = seq(levels(cellchat@idents)),  small.gap = 1,
                                         big.gap = 10),
                    paste0("chord_gene_cluster_sources_", use, "_target_all"))
      plot_function(netVisual_chord_gene(cellchat, sources.use = seq(levels(cellchat@idents)), targets.use = use,  small.gap = 1,
                                         big.gap = 10),
                    paste0("chord_gene_cluster_sources_all_target_", use))
    }
    

    if (FALSE) {
      selectK(cellchat, pattern = "outgoing")
      cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 5)
      netAnalysis_river(cellchat, pattern = "outgoing")
      
      netAnalysis_river(cellchat, pattern = "outgoing")

    }

  }

  ## Script ----
  prepare_data()

  if (load_object) {
    tryCatch({
      cellchat <- readRDS(file.path(output_dir, obj_name))
      assign("cellchat", cellchat, envir = parent.frame())
      message(file.path(output_dir, obj_name))
    }, error = function(e) {
      if (grepl("gzfile", e$message)) {
        message("Error in opening gzfile: ", e$message, " Did you try to load a cellchat object that has not been yet created?")
        # Handle the error (e.g., skip, try another file, or terminate)
      } else {
        stop(e)  # Re-throw the error if it's not the expected one
      }
    })
  } else compute_cellchat()
  
  plot_results()
  openxlsx::write.xlsx(subsetCommunication(cellchat), file.path(output_dir, "interactions_dataframe.xlsx"))
  # prepare_data()
  
  # if (load_object) cellchat <- readRDS(paste0(output_dir, obj_name))
  # else compute_cellchat()
  
  # plot_results()



  
  # cellChat <- createCellChat(object = seurat_object, group.by = "ident", assay = "RNA")
  
}
