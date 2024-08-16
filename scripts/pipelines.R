# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


preprocessing <- function(env_variables, preprocessing_settings) {

    source("scripts/seurat_utils.r", local = TRUE)
    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
      
    parameters <- .update_parameters(preprocessing_settings,  load_settings(settings_path)$preprocessing)

    # Preparation of subjects info
    data_preparation <- preparation_for_data_loading(data_folder, count_matrix_pattern, patient_info)

    # Seurat object creation and quality control
    seurat_object <- seurat_objects_and_quality_control(data_preparation$count_matrix_files,
                                                        data_preparation$subjects_info, save = FALSE,
                                                        normalization = TRUE, parts_to_remove = parameters$parts_to_remove)

    assign("seurat_object", seurat_object, envir = .GlobalEnv)

    if (parameters$save) saveRDS(seurat_object, file = paste0(output_folder, parameters$save_name))
}

integration <- function(env_variables, integration_settings) {

    source("scripts/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
      
    parameters <- .update_parameters(integration_settings,  load_settings(settings_path)$integration)
   
    data_preparation <- preparation_for_data_loading(data_folder, count_matrix_pattern, patient_info)
    seurat_object <- merge_seurat_objects(seurat_object, data_preparation$subjects_info$subject, save = FALSE)
    
    # Norm, Scaling, variable features
    seurat_object <- normalization_and_scaling(seurat_object, save = FALSE, normalization = FALSE,
                                                variable_features = TRUE, scaling = TRUE)
    
    # Prepare metadata
    seurat_object <- create_metadata(seurat_object, save = FALSE)
    
    # Run PCA before layer integration
    seurat_object <- run_pca(seurat_object, save = FALSE, plot = TRUE)
 
    # Select the number of dimension that expain more variance than a threshold
    PCs <- analyze_explained_variance(seurat_object, dstd, reduction_to_inspect = "pca") 
    
    # Perform Clustering
    seurat_object <- clustering(seurat_object, reduction = "pca", 
                                desired_resolution = d, dimension = PCs,
                                save = FALSE, column_name = "before_integration")
    
    # Visualization
    seurat_object <- visualization_UMAP(seurat_object, 
                                        reduction_name = "umap_before_integration", reduction = "pca", 
                                        cluster_column = "before_integration", dimension = PCs,
                                        save = FALSE, name = "before_integration",
                                        extension_plot = extension_plot)
    seurat_object <- layer_integration(seurat_object, assay = a, make_default = TRUE,
                                        new_reduction = r, reduction_method = parameters$method)


    # Select the number of dimension that expain more variance than a threshold
    PCs <- analyze_explained_variance(seurat_object, dstd, reduction_to_inspect = r) 

    # Perform Clustering
    seurat_object <- clustering(seurat_object, reduction = r, 
                                desired_resolution = d, dimension = PCs,
                                save = FALSE, column_name = c)

    # Visualization
    seurat_object <- visualization_UMAP(seurat_object, 
                                        reduction_name = umap, reduction = r, 
                                        cluster_column = c, dimension = PCs,
                                        save = FALSE, name = "after_integration",
                                        extension_plot = extension_plot)

    if (parameters$save) saveRDS(seurat_object, file = paste0(output_folder, parameters$save_name))
    assign("seurat_object", seurat_object, envir = .GlobalEnv)
}

annotation <- function(env_variables, annotation_settings) {

    source("scripts/seurat_utils.r", local = TRUE)
    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
      
    parameters <- .update_parameters(annotation_settings,  load_settings(settings_path)$annotation)
    
    # Annotation using scType
    if (parameters$annotate) seurat_object <- scType_annotation(seurat_object, assay = a,
                                        clusters_column = c,
                                        assignment_name = annotation)
    # Plot after annotation
        
    if (parameters$annotation_plots) {
        PCs <- analyze_explained_variance(seurat_object, dstd, reduction_to_inspect = r) 
        seurat_object <- visualization_UMAP(seurat_object, reduction_name = umap, 
                                    cluster_column = annotation, dimension = PCs,
                                    save = FALSE, run_UMAP = FALSE, name = "scType_annotation")
    }

    ## Correction of the annotation 
    if (parameters$correct) seurat_object <- correct_annotation(seurat_object, to_correct, 
                                        annotation, 
                                        new_annotation_column = corrected_annotation, 
                                        cluster_column = c)

    # Plot after correction
    if (parameters$corrected_annotation_plots) seurat_object <- visualization_UMAP(seurat_object, reduction_name=umap, 
                                                cluster_column=corrected_annotation, dimension=PCs,
                                                save=FALSE, run_UMAP=FALSE, name=paste0("scType_corrected_", m))

    # saveRDS(seurat_object, file = paste0(output_folder, "main_after_annotation.rds"))
    if (parameters$save) saveRDS(seurat_object, file = paste0(output_folder, parameters$save_name))
    assign("seurat_object", seurat_object, envir = .GlobalEnv)
    
}

clustering <- function(env_variables, clustering_settings) {

    source("scripts/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
      
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

    if (parameters$save) saveRDS(seurat_object, file = paste0(output_folder, parameters$save_name))
    assign("seurat_object", seurat_object, envir = .GlobalEnv)
    return(seurat_object)

}

deg <- function(env_variables, deg_settings) {

    source("scripts/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
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
                extension_plot = parameters$extension_plot,
                markers = parameters$markers)
 
        } else if (parameters$method == "volcano") {
            if (isFALSE(parameters$folder)) stop("folder is a required_parameter for volcano")

            volcano_plot(parameters$folder,
                extension_plot = parameters$extension_plot)
            
        } else if (parameters$method == "dimred") {
            if (isFALSE(parameters$markers)) stop("markers is a required_parameter for dimred")

            plots_for_paper(
                which = c("feature_plots"), 
                name = names(deg_settings)[deg],
                extension_plot = parameters$extension_plot,
                genes_to_plot = parameters$markers,
                assay = a,               
                cluster_column = parameters$cluster_column)
            
        } else stop("enrichment analysis: not valid method")
    }
}

wgcna <- function(env_variables, wgcna_settings){

    source("scripts/wgcna.r", local = TRUE)
    source("scripts/seurat_utils.r", local = TRUE)  

    create_variables(.update_parameters(env_variables,  load_settings(settings_path
    )$global_variables))
    default_parameters <- load_settings(settings_path)$wgcna

    ## WGCNA for clusters ---
    for (wgcna in seq_along(wgcna_settings)) {

        parameters <- .update_parameters(wgcna_settings[[wgcna]], default_parameters)

        if (parameters$method == "WGCNA") {
            if (!isFALSE(parameters$cluster))  {

                seurat_object_subset <- create_object_from_cluster_id(seurat_object, 
                                                                        parameters$cluster,
                                                                        assay = a,
                                                                        clusters_column = c)
            } else seurat_object_subset <- seurat_object
                    
            do.call(wgcna_main, 
                    c(list(seurat_object_subset, 
                    name = names(wgcna_settings)[wgcna]), 
                    parameters))
        }
        else if (parameters$method == "hdWGCNA") {
            source("scripts/hdwgcna.r", local = TRUE)

            do.call(hdwgcna, 
                    c(list(seurat_object, 
                    name = names(wgcna_settings)[wgcna]), 
                    parameters))
            
        }
    }
}


enrichment <- function(env_variables, enrichment_settings){

    source("scripts/enrichment.r", local = TRUE)
    source("scripts/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path
    )$global_variables))
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