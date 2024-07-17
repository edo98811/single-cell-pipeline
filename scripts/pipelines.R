# To implement
preprocessing <- function(env_variables, clustering_settings) {    
}

# To implemennt
integration <- function(env_variables, clustering_settings) {
      
  ## Merge objects scaling, metadata, and pca ----
  # starting from a list of Seurat objects merges them 
  seurat_object <- merge_seurat_objects(seurat_objects, data_preparation$subjects_info$subject, save=FALSE)
  rm(seurat_objects)
  
  # Norm, Scaling, variable features
  seurat_object <- preprocessing_and_scaling(seurat_object, save=FALSE, normalization=FALSE,
                                             variable_features=TRUE, scaling=TRUE)
  
  # Prepare metadata
  seurat_object <- create_metadata(seurat_object, save=FALSE)
  
  # Run PCA before layer integration
  seurat_object <- run_pca(seurat_object, save=FALSE, plot=TRUE)
  
  ## new parameters ----
  
  # Before integration
  r <- "pca" # reduction to use for computations
  c <- "not_integrated_clusters" #  column in metadata
  # se no integration questi rimangono per dopo
  m <- "no_integration"# method for reduction
  a <- "RNA" # assay
  umap <- paste0("umap_", r)
  
  ## Clustering before integration ----
  
  if (FALSE){
    saveRDS(seurat_object, file = paste0(output_folder, "main_before_integration_and_clustering.rds"))
    seurat_object <- readRDS(paste0(output_folder, "main_before_integration_and_clustering.rds"))
  }
  
  # Select the number of dimension that expain more variance than a threshold
  PCs <- analyze_explained_variance(seurat_object, dstd, reduction_to_inspect = "pca") 
  if (!is.numeric(PCs)) {cat("Failed to identify correct number of dimensions, using default: 16"); PCs <- 16}
  
  # Perform Clustering
  seurat_object <- clustering(seurat_object, reduction=r, 
                              desired_resolution=desired_resolution, dimension=PCs,
                              save=FALSE, column_name = c)
  
  # Visualization
  message <- "visualisation of theclusters obtained before performing any type of integration"
  seurat_object <- visualization_UMAP(seurat_object, 
                                                 reduction_name=paste0("umap_", r), reduction=r, 
                                                 cluster_column=c, dimension=PCs,
                                                 save=FALSE, name="before_integration", message=message)
  ## new parameters ----
  # Integration and annotation parameters
  m <- "harmony"
  a <- "RNA"
  r <- paste0(m, "_reduction")
  c <- paste0(m, "_clusters")
  umap <- paste0("umap_", r)

  
  
  ## INTEGRATION ----
  if (FALSE){
    saveRDS(seurat_object, file = paste0(output_folder, "main_before_integration.rds"))
    seurat_object <- readRDS(paste0(output_folder, "main_before_integration.rds"))
    seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
  }
  
  
  seurat_object <- layer_integration(seurat_object, assay=a, make_default=TRUE,
                                     new_reduction = r, reduction_method = m)
  

}


clustering <- function(env_variables, clustering_settings){

    source("scripts/seurat_utils.r", local = TRUE)
    create_variables(env_variables)
      
    parameters <- .update_parameters(clustering_settings,  load_settings(settings_path)$clustering)
   
    if (parameters$create_subset) {
        if (length(parameters$subset) == 0) stop(" subsetting: invalit subset parameter")
        
        assign("seurat_object", create_object_from_cluster_id(seurat_object, 
                parameters$subset, clusters_column = parameters$cluster_column),
                envir = .GlobalEnv)
    }

    if (parameters$clustering) {
        PCs <- analyze_explained_variance(seurat_object, 
                                    dstd, 
                                    reduction_to_inspect = "pca")
                        # Perform Clustering
        seurat_object <- clustering(seurat_object, 
                                    reduction = r, 
                                    desired_resolution = d, 
                                    dimension = PCs,
                                    save = FALSE, 
                                    column_name = c)

        seurat_object <- visualization_UMAP(seurat_object, 
                                            reduction_name = umap, 
                                            cluster_column = c, 
                                            dimension = PCs,
                                            save = FALSE, 
                                            run_UMAP = FALSE, 
                                            name = parameters$name,
                                            extension_plot = extension_plot)
    }
    if (parameters$rename_clusters) {

        seurat_object <- correct_annotation(seurat_object, 
                                        parameters$to_correct, 
                                        c, 
                                        new_annotation_column = corrected_annotation, 
                                        cluster_column = c)
    
    
        seurat_object <- visualization_UMAP(seurat_object, 
                                        reduction_name = umap, 
                                        reduction = r, 
                                        cluster_column = corrected_annotation, 
                                        dimension = PCs,
                                        save = FALSE, 
                                        run_UMAP = FALSE, 
                                        name = paste0(parameters$name, "corrected"),
                                        extension_plot = extension_plot)
    }    

    if (parameters$save) saveRDS(seurat_object, parameters$save_path)
    return(seurat_object)

}


deg <- function(env_variables, deg_settings) {

    source("scripts/seurat_utils.r", local = TRUE)

    create_variables(env_variables)
    default_parameters <- load_settings(settings_path)$deg

    iterative_methods <- c(
        "condition_and_clusters",
        "condition_and_clusters_vf",
        "condition_and_clusters_DESeq2",
        "default",
        "default_vf",
        "default_LR")

    normal_methods <- c(
        "condition",
        "condition_vf",
        "condition_DESeq2")

    for (deg in seq_along(deg_settings)) {

        parameters <-  .update_parameters(deg_settings[[deg]], default_parameters)
        
        if(parameters$method %in% iterative_methods) {
            for (cluster in unique(seurat_object@meta.data[[parameters$cluster_column]])) {
                do.call(find_and_plot_markers, 
                        c(list(seurat_object, 
                            cluster_id=cluster,
                            reduction_name = umap, 
                            name = names(deg_settings)[deg]),
                            parameters[names(parameters) != "folder"]))
            }
        } else if (parameters$method %in% normal_methods) {
            do.call(find_and_plot_markers, 
                    c(list(seurat_object, 
                        reduction_name = umap, 
                        name = names(deg_settings)[deg]),
                        parameters[names(parameters) != "folder"]))

        } else if (parameters$method == "heatmap") {

            many_plots(seurat_object, 
                which = c("heatmap"), 
                assay = a,
                cluster_column = parameters$cluster_column, 
                name = names(deg_settings)[deg],
                extension_plot = parameters$extension_plot)
 
        } else if (parameters$method == "volcano") {
            if (isFALSE(parameters$folder)) stop("folder is a requuíred_parameter for volcano")

            volcano_plot(parameters$folder,
            extension_plot = parameters$extension_plot)
            
        } else stop("enrichment analysis: not valid method")
    }
}


wgcna <- function(env_variables, wgcna_settings){

    source("scripts/wgcna.r", local = TRUE)
    source("scripts/seurat_utils.r", local = TRUE)  

    create_variables(env_variables)
    default_parameters <- load_settings(settings_path)$wgcna

    ## WGCNA for clusters ---
    for (wgcna in seq_along(wgcna_settings)) {

        parameters <- .update_parameters(wgcna_settings[[wgcna]], default_parameters)
    
        if (!isFALSE(parameters$cluster))  {

            seurat_object_subset <- create_object_from_cluster_id(seurat_object, 
                                                                    parameters$cluster,
                                                                    assay = a,
                                                                    clusters_column = c)
        } else seurat_object_subset <-seurat_object
                
        do.call(wgcna_main, 
                c(list(seurat_object_subset, 
                name = names(wgcna_settings)[wgcna]), 
                parameters))
    }
}


enrichment <- function(env_variables, enrichment_settings){

    source("scripts/enrichment.r", local = TRUE)
    source("scripts/seurat_utils.r", local = TRUE)

    create_variables(env_variables)
    default_parameters <- load_settings(settings_path)$enrichment

    for (enrichment in seq_along(enrichment_settings)) {

        parameters <- .update_parameters(enrichment_settings[[enrichment]], default_parameters)
        
        do.call(enrichment_analysis, c(list(
                                    names(enrichment_settings)[enrichment], 
                                    parameters$markers_path), 
                                    parameters[names(parameters) != "markers_path"]))
    }

}

.update_parameters <- function(parameters, default) {

    message("setting parameters")
    for (parameter in names(default)) {
        
        # If the value is not present in the parameters given it is added from defaults
        if (!parameter %in% names(parameters)) {

            # If the default parameter is null stop the program (it is required)
            if (is.null(default[[parameter]])) stop(parameter, " required parameter")

            # Else add it 
            parameters[[parameter]] <- default[[parameter]]
        }
        

        message("   parameter ", parameter, ": ", parameters[[parameter]])
    }
    
    return(parameters)
}