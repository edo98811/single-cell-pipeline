# can i make this iterative? (check parameters for deg)
check_parameters <- function(par_list, method) {
    if (method == "normal") {
        if (c %in% names(par_list)) {
            if (!is.character(par_list$c))  warning("non valid argument", par_list$c)
        } else {
            par_list$c <- default
            warning("parameter ", c, "not given, set to default", default) 
            }

    } else if (method == "normal") {
        if (is.par_list$name) {}

    } else if (method == "normal") {
        if (is.par_list$name) {}

    } else if (method == "normal") {
        if (is.par_list$name) {}

    }
}

# OLD deg pipeline
if (FALSE) {       
    message("old deg pipeline")
    # if (which_method$vf) {

    #     # For clusters considering only the variable features
    #     if(which_data$vf$clusters) {
    #         for (cluster in unique(seurat_object@meta.data[[c]])) {
    #             find_and_plot_markers(seurat_object, 
    #                                     cluster_id=cluster, 
    #                                     reduction_name=umap, 
    #                                     name="microglia_clusters_vf", 
    #                                     cluster_column=c, 
    #                                     method="default_vf")
    #         }
    #     }
    #     if (which_data$vf$within_clusters)
    #         # For pathology within cluster pd (only variable features)
    #         for (cluster in unique(seurat_object@meta.data[[c]])) {
    #         find_and_plot_markers(seurat_object, 
    #                                 cluster_id=cluster,
    #                                 reduction_name=umap, 
    #                                 name ="microglia_control_vs_pd_clusters_vf", 
    #                                 method = "condition_and_clusters_vf")
    #     }

    #     # For pathology dataset wide pd
    #     if (which_data$vf$all) {
    #         find_and_plot_markers(seurat_object, 
    #                                 reduction_name=umap,
    #                                 name ="microglia_control_vs_pd_nogenetic_vf", 
    #                                 method = "condition_vf")  
    #     }

    #     if (which_data$vf$all_genetic) {
    #         # For pathology dataset wide genetic pd
    #         # For pathology dataset wide genetic pd
    #         find_and_plot_markers(seurat_object, 
    #                             reduction_name=umap,
    #                             name ="microglia_control_vs_geneticpd", 
    #                             method = "condition", 
    #                             condition="genetic_PD")
        
    #     }
    # }

    # if (which_method$all) {
    
    #     if (which_data$all$clusters) {   
    #         # For pathology within cluster pd(all)
    #         for (cluster in unique(seurat_object@meta.data[[c]])) {
    #         find_and_plot_markers(seurat_object, 
    #                                 cluster_id=cluster,
    #                                 reduction_name=umap, 
    #                                 name ="microglia_control_vs_pd_clusters_all", 
    #                                 method = "condition_and_clusters")
    #         }
    #     }

    #     if (which_data$all$within_clusters) {
    #         # For pathology dataset wide pd(all)
    #         find_and_plot_markers(seurat_object, 
    #                             reduction_name=umap,
    #                             name ="microglia_control_vs_pd_nogenetic_all", 
    #                             method = "condition",
    #                             nothreshold=TRUE)
    #     }
    
    #     if (which_data$all$genetic_all) {
    #         # For pathology dataset wide genetic pd(all)
    #         find_and_plot_markers(seurat_object, 
    #                             reduction_name=umap,
    #                             name ="microglia_control_vs_geneticpd_all", 
    #                             method = "condition", 
    #                             condition="genetic_PD",
    #                             nothreshold=TRUE)     
    #     }

    #     if (which_data$all$all_genetic) {
    #         # For pathology dataset wide genetic pd
    #         find_and_plot_markers(seurat_object, 
    #                             reduction_name=umap,
    #                             name ="microglia_control_vs_geneticpd", 
    #                             method = "condition", 
    #                             condition="genetic_PD")
        
    #     }
    # }
    
    # if (which_method$normal) {
        
    #     if (which_data$normal$cluster) {
    #         # For clusters
    #         for (cluster in unique(seurat_object@meta.data[[c]])) {
    #         find_and_plot_markers(seurat_object, 
    #                                 cluster_id=cluster, 
    #                                 reduction_name=umap, 
    #                                 name="microglia_clusters", 
    #                                 cluster_column=c)
    #         }
    #     }

    #     if (which_data$normal$within_clusters) {
    #         # For pathology within cluster pd
    #         for (cluster in unique(seurat_object@meta.data[[c]])) {
    #             find_and_plot_markers(seurat_object, 
    #                                     cluster_id=cluster,
    #                                     reduction_name=umap, 
    #                                     name ="microglia_control_vs_pd_clusters", 
    #                                     method = "condition_and_clusters")
    #         }
    #     }
    #     if (which_data$normal$all) {
    #         # For pathology dataset wide pd
    #         find_and_plot_markers(seurat_object, 
    #                             reduction_name=umap,
    #                             name ="microglia_control_vs_pd_nogenetic", 
    #                             method = "condition")
            
    #         if (volcano) volcano_plot("microglia_control_vs_pd_nogenetic")
    #     }

    #     if (which_data$normal$all_genetic) {
    #         # For pathology dataset wide genetic pd
    #         find_and_plot_markers(seurat_object, 
    #                             reduction_name=umap,
    #                             name ="microglia_control_vs_geneticpd", 
    #                             method = "condition", 
    #                             condition="genetic_PD")
        
    #     }

        
    # }

    # if (which_method$many_plots) {
    
    #     if (which_data$many_plots$cluster) 
            
    #         many_plots(seurat_object, 
    #             which=c("heatmap"), 
    #             assay="RNA",
    #             cluster_column=c, 
    #             name ="microglia")
        
    #     if (which_data$many_plots$cluster_names) 
            
    #         many_plots(seurat_object, 
    #             which=c("heatmap"), 
    #             assay="RNA",
    #             cluster_column=corrected_annotation_microglia, 
    #             name ="microglia")
        
    #     if (which_data$many_plots$pathology) 
            
    #         many_plots(seurat_object, 
    #             which=c("heatmap"), 
    #             assay="RNA",
    #             cluster_column="subject_pathology", 
    #             name ="microglia")

    # }

}