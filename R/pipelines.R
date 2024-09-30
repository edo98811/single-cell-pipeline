# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


preprocessing <- function(env_variables, preprocessing_settings) {
    message("section: preprocessing")
    source("R/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
      
    parameters <- .update_parameters(preprocessing_settings,  load_settings(settings_path)$preprocessing)

    # Preparation of subjects info
    if (isFALSE(data_folder)) stop("Pipelines: Data folder not set and required to load data")
    data_preparation <- preparation_for_data_loading(data_folder, count_matrix_pattern, patient_info)

    # Seurat object creation and quality control
    seurat_object <- seurat_objects_and_quality_control(data_preparation$count_matrix_files,
                                                        data_preparation$subjects_info, save = FALSE,
                                                        normalization = TRUE, parts_to_remove = parameters$parts_to_remove)

    assign("seurat_object", seurat_object, envir = .GlobalEnv)

    if (parameters$save) {
        if (!isFALSE(parameters$save_name)) {
            saveRDS(seurat_object, file = parameters$save_name)
            message("object saved in: ", parameters$save_name)
        } else {
            warning(" main pipeline: If save is set to true and save_name to false the seurat object cannot be saved")
        }
    }
}

integration <- function(env_variables, integration_settings) {

    message("section: integration")
    source("R/seurat_utils.r", local = TRUE)

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
    pc_number <- analyze_explained_variance(seurat_object, dstd, reduction_to_inspect = "pca") 
    
    # Perform Clustering
    seurat_object <- clustering(seurat_object, reduction = "pca", 
                                desired_resolution = d, dimension = pc_number,
                                save = FALSE, column_name = "before_integration")
    
    # Visualization
    seurat_object <- visualization_UMAP(seurat_object, 
                                        reduction_name = "umap_before_integration", reduction = "pca", 
                                        cluster_column = "before_integration", dimension = pc_number,
                                        save = FALSE, name = "before_integration",
                                        extension_plot = extension_plot)
    seurat_object <- layer_integration(seurat_object, assay = a, make_default = TRUE,
                                        new_reduction = r, reduction_method = parameters$method)


    # Select the number of dimension that expain more variance than a threshold
    pc_number <- analyze_explained_variance(seurat_object, dstd, reduction_to_inspect = r) 

    # Perform Clustering
    seurat_object <- clustering(seurat_object, reduction = r, 
                                desired_resolution = d, dimension = pc_number,
                                save = FALSE, column_name = c)

    # Visualization
    seurat_object <- visualization_UMAP(seurat_object, 
                                        reduction_name = umap, reduction = r, 
                                        cluster_column = c, dimension = pc_number,
                                        save = FALSE, name = "after_integration",
                                        extension_plot = extension_plot)

    if (parameters$save) {
        if (!isFALSE(parameters$save_name)) {
            saveRDS(seurat_object, file = parameters$save_name)
            message("object saved in: ", parameters$save_name)
        } else {
            warning(" main pipeline: If save is set to true and save_name to false the seurat object cannot be saved")
        }
    }
    assign("seurat_object", seurat_object, envir = .GlobalEnv)
}

annotation <- function(env_variables, annotation_settings) {

    message("section: annotation")

    source("R/seurat_utils.r", local = TRUE)
    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
      
    parameters <- .update_parameters(annotation_settings,  load_settings(settings_path)$annotation)
    
    # Annotation using scType
    if (parameters$annotate) seurat_object <- scType_annotation(seurat_object, assay = a,
                                        clusters_column = c,
                                        assignment_name = annotation)
    # Plot after annotation
    if (parameters$annotation_plots) {
        pc_number <- analyze_explained_variance(seurat_object, dstd, reduction_to_inspect = r) 
        seurat_object <- visualization_UMAP(seurat_object, reduction_name = umap, 
                                    cluster_column = annotation, dimension = pc_number,
                                    save = FALSE, run_umap = FALSE, name = "scType_annotation")
    }

    ## Correction of the annotation 
    if (isFALSE(parameters$to_correct) && parameters$correct) stop("Pipelines: to_correct is a required_parameter for annotation correction")
    if (parameters$correct) seurat_object <- correct_annotation(seurat_object, parameters$to_correct, 
                                        annotation, 
                                        new_annotation_column = corrected_annotation, 
                                        cluster_column = c)

    # Plot after correction
    if (parameters$corrected_annotation_plots) seurat_object <- visualization_UMAP(seurat_object, reduction_name = umap, 
                                                cluster_column = corrected_annotation, dimension = pc_number,
                                                save = FALSE, run_umap = FALSE, name = file.path("scType_corrected_", m))

    # saveRDS(seurat_object, file = file.path(output_folder, "main_after_annotation.rds"))
    if (parameters$save) {
        if (!isFALSE(parameters$save_name)) {
            saveRDS(seurat_object, file = parameters$save_name)
            message("object saved in: ", parameters$save_name)
        } else {
            warning(" main pipeline: If save is set to true and save_name to false the seurat object cannot be saved")
        }
    }
    assign("seurat_object", seurat_object, envir = .GlobalEnv)
    
}

clustering <- function(env_variables, clustering_settings) {
    
    message("section: clustering")
    source("R/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
      
    parameters <- .update_parameters(clustering_settings,  load_settings(settings_path)$clustering)
   
    if (parameters$create_subset) {
        if (length(parameters$subset) == 0) stop("Pipelines:  subsetting: invalid subset parameter")

        assign("seurat_object", create_object_from_cluster_id(seurat_object, 
                parameters$subset, clusters_column = parameters$cluster_column),
                envir = .GlobalEnv)
    }

    if (parameters$before_clustering_plot) {
        pc_number <- analyze_explained_variance(seurat_object, 
                                    dstd, 
                                    reduction_to_inspect = r)


        seurat_object <- visualization_UMAP(seurat_object, 
                                            reduction_name = ifelse(isFALSE(parameters$umap_before_clustering), umap, parameters$umap_before_clustering),
                                            cluster_column = c, 
                                            reduction = r,
                                            save = FALSE, 
                                            run_umap = FALSE, 
                                            name = paste0(parameters$name, "_before_clustering"),
                                            extension_plot = extension_plot)
    }

    if (parameters$clustering || parameters$clustering_plot)
        pc_number <- analyze_explained_variance(seurat_object, 
                                    dstd, 
                                    reduction_to_inspect = r)

    if (parameters$clustering)  seurat_object <- clustering(seurat_object, 
                                    reduction = r, 
                                    desired_resolution = d, 
                                    dimension = pc_number,
                                    save = FALSE, 
                                    column_name = c)

    if (parameters$clustering_plot) seurat_object <- visualization_UMAP(seurat_object, 
                                            reduction_name = umap, 
                                            cluster_column = c, 
                                            reduction = r,
                                            dimension = pc_number,
                                            save = FALSE, 
                                            run_umap = TRUE, 
                                            name = parameters$name,
                                            extension_plot = extension_plot)

    if (parameters$rename_clusters) {
        if (isFALSE(parameters$to_correct)) stop("Pipelines: to_correct is a required_parameter for annotation correction")

        seurat_object <- correct_annotation(seurat_object, 
                                        parameters$to_correct, 
                                        c, 
                                        new_annotation_column = corrected_annotation, 
                                        cluster_column = c)
    
    
        seurat_object <- visualization_UMAP(seurat_object, 
                                        reduction_name = umap, 
                                        reduction = r, 
                                        cluster_column = corrected_annotation, 
                                        dimension = pc_number,
                                        save = FALSE, 
                                        run_umap = FALSE, 
                                        name = paste0(parameters$name, "corrected"),
                                        extension_plot = extension_plot)
    }    

    if (parameters$save) {
        if (!isFALSE(parameters$save_name)) {
            saveRDS(seurat_object, file = parameters$save_name)
            message("object saved in: ", parameters$save_name)
        } else {
            warning(" main pipeline: If save is set to true and save_name to false the seurat object cannot be saved")
        }
    }

    assign("seurat_object", seurat_object, envir = .GlobalEnv)

}

subsetting <- function(env_variables, subsetting_settings) {
        
    message("section: subsetting")
    source("R/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
    parameters <- .update_parameters(subsetting_settings,  load_settings(settings_path)$subsetting_settings)

    if (length(parameters$subset) == 0) stop("Pipelines:  subsetting: invalid subset parameter")

        assign("seurat_object", 
            create_object_from_cluster_id(seurat_object, 
                parameters$subset, 
                clusters_column = ifelse(isFALSE(parameters$cluster_column), c, parameters$cluster_column)),
            envir = .GlobalEnv)
    
    if (parameters$save) {
        if (!isFALSE(parameters$save_name)) {
            saveRDS(seurat_object, file = parameters$save_name)
            message("object saved in: ", parameters$save_name)
        } else {
            warning(" main pipeline: If save is set to true and save_name to false the seurat object cannot be saved")
        }
    }
}

main_pipeline <- function(env_variables, main_pipeline_settings) {

    message("section: main_pipeline")
    source("R/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path)$global_variables))
    default_parameters <- load_settings(settings_path)$main_pipeline
    
    for (name in seq_along(main_pipeline_settings)) {

        parameters <-  .update_parameters(main_pipeline_settings[[name]], default_parameters)

        if (parameters$method == "annotation") {

            if (parameters$automatic_annotation) seurat_object <- scType_annotation(seurat_object, 
                                                assay = a,
                                                clusters_column = ifelse(isFALSE(parameters$cluster_column), c, parameters$cluster_column),
                                                assignment_name = ifelse(isFALSE(parameters$annotation_column), annotation, parameters$annotation_column))


            ## Correction of the annotation 
            if (!isFALSE(parameters$to_correct)) 
                seurat_object <- correct_annotation(seurat_object, parameters$to_correct, 
                                                ifelse(isFALSE(parameters$annotation_column), annotation, parameters$annotation_column), 
                                                new_annotation_column = ifelse(isFALSE(parameters$corrected_annotation_column), 
                                                    ifelse(isFALSE(parameters$annotation_column), annotation, parameters$annotation_column), parameters$corrected_annotation_column), 
                                                cluster_column = ifelse(isFALSE(parameters$cluster_column), c, parameters$cluster_column))

            assign("seurat_object", seurat_object, envir = .GlobalEnv)

            if (parameters$save) {
                if (!isFALSE(parameters$save_name)) {
                    saveRDS(seurat_object, file = parameters$save_name)
                    message("object saved in: ", parameters$save_name)
                } else {
                    warning(" main pipeline: If save is set to true and save_name to false the seurat object cannot be saved")
                }
            }

        }
        if (parameters$method == "clustering") {
            seurat_object <- clustering(seurat_object, 
                                    reduction = ifelse(isFALSE(parameters$reduction), r, parameters$reduction), 
                                    desired_resolution = d, 
                                    dimension = pc_number,
                                    column_name = ifelse(isFALSE(parameters$cluster_column), c, parameters$cluster_column))

            assign("seurat_object", seurat_object, envir = .GlobalEnv)

            if (parameters$save) {
                if (!isFALSE(parameters$save_name)) {
                    saveRDS(seurat_object, file = parameters$save_name)
                    message("object saved in: ", parameters$save_name)
                } else {
                    warning(" main pipeline: If save is set to true and save_name to false the seurat object cannot be saved")
                }
            }
        }
        if (parameters$method == "plotting") {
            
            pc_number <- analyze_explained_variance(seurat_object, 
                            ifelse(isFALSE(parameters$explained_variance), dstd, parameters$explained_variance), 
                            reduction_to_inspect = ifelse(isFALSE(parameters$reduction), r, parameters$reduction))

            if (!isFALSE(parameters$compute_umap)) message("plotting: a new umap will be computed and saved in this slot: ", ifelse(isFALSE(parameters$umap_name), umap, parameters$umap_name))
            seurat_object <- visualization_UMAP(seurat_object, 
                                reduction_name = ifelse(isFALSE(parameters$umap_name), umap, parameters$umap_name), 
                                reduction = ifelse(isFALSE(parameters$reduction), r, parameters$reduction), 
                                cluster_column = ifelse(isFALSE(parameters$cluster_column), c, parameters$cluster_column), 
                                dimension = pc_number,
                                save = FALSE, 
                                run_umap = parameters$compute_umap, 
                                name = names(main_pipeline_settings)[name],
                                extension_plot = extension_plot)

            assign("seurat_object", seurat_object, envir = .GlobalEnv)
            if (parameters$save) {
                if (!isFALSE(parameters$save_name)) {
                    saveRDS(seurat_object, file = parameters$save_name)
                    message("object saved in: ", parameters$save_name)
                } else {
                    warning(" main pipeline: If save is set to true and save_name to false the seurat object cannot be saved")
                }
            }
        }
    }
}

deg <- function(env_variables, deg_settings) {

    message("section: differential gene expression analysis")
    source("R/seurat_utils.r", local = TRUE)

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

        if (parameters$method %in% iterative_methods) {
            for (cluster in unique(seurat_object@meta.data[[parameters$cluster_column]])) {
                try(do.call(find_and_plot_markers, 
                        c(list(seurat_object, 
                            cluster_id = cluster,
                            reduction_name = umap, 
                            name = names(deg_settings)[deg]),
                            parameters[names(parameters) != "folder"])))
            }
        } else if (parameters$method %in% normal_methods) {
            try(do.call(find_and_plot_markers, 
                    c(list(seurat_object, 
                        reduction_name = umap, 
                        name = names(deg_settings)[deg]),
                        parameters[names(parameters) != "folder"])))

        } else if (parameters$method == "heatmap") {

            try(plot_heatmap(seurat_object, 
                which = c("heatmap"), 
                assay = a,
                cluster_column = parameters$cluster_column, 
                name = names(deg_settings)[deg],
                extension_plot = extension_plot,
                markers = parameters$markers,
                sorting_method = parameters$sorting_method))
 
        } else if (parameters$method == "volcano") {
            if (isFALSE(parameters$folder)) stop("Pipelines: folder is a required_parameter for volcano")

            try(volcano_plot(parameters$folder,
                extension_plot = extension_plot))
            
        } else if (parameters$method == "plots_misc") {
            if (isFALSE(parameters$markers)) stop("Pipelines: markers is a required_parameter for plots_for_paper")

            try(plots_for_paper(seurat_object,
                which = parameters$which_other, 
                name = names(deg_settings)[deg],
                extension_plot = extension_plot,
                genes_to_plot = parameters$markers,
                assay = a,               
                cluster_column = parameters$cluster_column,
                reduction_name = ifelse(isFALSE(parameters$umap_name), umap, parameters$umap_name),
                subplot_n = parameters$subplot_n))
            
        } else if (parameters$method == "other_plots_from_df") {
            if (isFALSE(parameters$markers)) stop("Pipelines: markers is a required_parameter for plot_markers_from_df")
            if (!(is.character(parameters$markers) && grepl("\\.xlsx?$", parameters$markers, ignore.case = TRUE)))  stop("Pipelines: markers must point to an excel file for plot_markers_from_df")

            try(plot_markers_from_df(seurat_object,
                name = names(deg_settings)[deg],
                extension_plot = extension_plot,
                markers_location = parameters$markers,
                assay = a,               
                cluster_column = parameters$cluster_column,
                reduction_name = umap,
                subplot_n = parameters$subplot_n,
                many_plot = parameters$heatmap, 
                feature_plot = parameters$feature_plot,
                heatmap_by_column = parameters$heatmap_by_column, 
                max_feature_plots = parameters$max_feature_plots, 
                max_genes_plot_heatmap = parameters$max_genes_heatmap,
                column_list = parameters$column_list))
                
        } else stop("Pipelines: deg analysis: not valid method")
    }
}

wgcna <- function(env_variables, wgcna_settings) {

    message("section: wgcna")
    source("R/wgcna.r", local = TRUE)
    source("R/seurat_utils.r", local = TRUE)  

    create_variables(.update_parameters(env_variables,  load_settings(settings_path
    )$global_variables))
    default_parameters <- load_settings(settings_path)$wgcna

    ## WGCNA for clusters ---
    for (wgcna in seq_along(wgcna_settings)) {

        parameters <- .update_parameters(wgcna_settings[[wgcna]], default_parameters)
        
        if (!isFALSE(parameters$cluster))  {
            try(seurat_object_subset <- create_object_from_cluster_id(seurat_object, 
                                                                    parameters$cluster,
                                                                    assay = a,
                                                                    clusters_column = ifelse(isFALSE(parameters$cluster_column),
                                                                    c,
                                                                    parameters$cluster_column)))

            message("Seurat object subsetted according to required clusters")
        } else seurat_object_subset <- seurat_object
        
        if (!isFALSE(parameters$wgcna_subjects))  {
            tryCatch({seurat_object_subset <- create_object_from_cluster_id(seurat_object_subset, 
                                                                    names(parameters$wgcna_subjects),
                                                                    assay = a,
                                                                    clusters_column = parameters$subject_column)
                message("Seurat object subsetted according to required subjects")
            # {
            # "PD_001": "02_082",
            # "PD_002": "02_084",
            # "PD_005": "02_074",
            # "PD_007": "02_096",
            # "PD_008": "02_097",
            # "PD_009": "02_105",
            # "PD_016": "02_108",
            # "PD_017": "02_115"
            # }
            }, 
            error = function(e) {
                warning("Could not subset seurat object according to subjects given, cotinuing with complete object, error: ", e)
                seurat_object_subset <- seurat_object
            })
        }   
        if (parameters$method == "WGCNA") {

            try(do.call(wgcna_main, 
                    c(list(seurat_object_subset, 
                    name = names(wgcna_settings)[wgcna]), 
                    parameters)))
        }
        else if (parameters$method == "hdWGCNA") {
            source("R/hdwgcna.r", local = TRUE)

            try(do.call(hdwgcna, 
                    c(list(seurat_object_subset, 
                    name = names(wgcna_settings)[wgcna]), 
                    parameters)))
            
        }
    }
}

own_script <- function(env_variables, own_path) {
    create_variables(.update_parameters(env_variables,  load_settings(settings_path
    )$global_variables))

    message("opening personalized script in location: ", own_path)
    source(own_path)
}

enrichment <- function(env_variables, enrichment_settings) {

    message("section: enrichment")
    source("R/enrichment.r", local = TRUE)
    source("R/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path
    )$global_variables))
    default_parameters <- load_settings(settings_path)$enrichment 
 
    for (enrichment in seq_along(enrichment_settings)) { 

        parameters <- .update_parameters(enrichment_settings[[enrichment]], default_parameters)
        
        try(do.call(enrichment_analysis, c(list(
                                    names(enrichment_settings)[enrichment], 
                                    parameters$markers_path), 
                                    parameters[names(parameters) != "markers_path"])))
    } 
} 

cellchat_analysis <- function(env_variables, cellchat_settings) {

    message("section: cellchat")
    source("R/seurat_utils.r", local = TRUE)

    create_variables(.update_parameters(env_variables,  load_settings(settings_path
    )$global_variables))
    default_parameters <- load_settings(settings_path)$cellchat_analysis 

    #for (cellchat in seq_along(cellchat_settings)) { 

    # parameters <- .update_parameters(cellchat_settings[[cellchat]], default_parameters) 
    parameters <- .update_parameters(cellchat_settings, default_parameters) 

    try(do.call(cellchat_function, c(list(
                            seurat_object,
                            #names(cellchat_settings)[cellchat], 
                            extension_plot = extension_plot),
                            parameters)))
    #}
}

plot_info <- function() {

    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Pipelines: The Seurat package must be installed.")
    }

    # Print general information
    cat("\n\nLoaded Seurat Object Summary:\n")
    
    # Number of cells
    num_cells <- ncol(seurat_object)
    cat("Number of cells:", num_cells, "\n")
    
    # Number of features
    num_features <- nrow(seurat_object)
    cat("Number of features:", num_features, "\n")
    
    # Print metadata
    cat("Metadata columns:\n")
    print(colnames(seurat_object@meta.data))
    
    # Print identities
    cat("Identity classes (idents):\n")
    if (!is.null(seurat_object@active.ident)) {
        print(table(seurat_object@active.ident))
    } else {
        cat("No active identities found.\n")
    }
    
    # Print reductions
    cat("Reductions:\n")
    reductions <- names(seurat_object@reductions)
    if (length(reductions) > 0) {
        print(reductions)
    } else {
        cat("No reductions found.\n")
    }
    
    # Print reductions
    cat("subjects:\n")
    subjects <- unique(seurat_object@meta.data$subject)
    if (!is.null(subjects)) {
        print(subjects)
    } else {
        cat("No subjects found.\n")
    }

    # Print assays
    cat("Assays:\n")
    assays <- names(seurat_object@assays)
    if (length(assays) > 0) {
        print(assays)
    } else {
        cat("No assays found.\n")
    }
    
    # Print other fields
    cat("Other fields in the Seurat object:\n")
    other_fields <- setdiff(names(seurat_object@misc), c("active.ident", "reductions", "assays", "meta.data"))
    if (length(other_fields) > 0) {
        print(other_fields)
    } else {
        cat("No additional fields found.\n")
    }
    cat("\n\n")
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
        # message("   parameter ", parameter, ": ", parameters[[parameter]])
        message("   parameter ", parameter, ": ", 
            ifelse(is.vector(parameters[[parameter]]), 
            paste(parameters[[parameter]], collapse = ", "), 
            as.character(parameters[[parameter]])))
    }
    
    return(parameters)
} 