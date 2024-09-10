# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de

# for debugging (error when using do.call and browser)
# https://stackoverflow.com/questions/44608323/crash-of-debugging-browser-in-r-rstudio-when-called-from-inside-do-call
options(deparse.max.lines = 10)
# install.packages("devtools")
devtools::install_github("immunogenomics/presto")
source("scripts/helper_functions.r")

#' Main Function to run the Pipeline
#'
#' It reads a configuration file in JSON format and runs the specified steps, including preprocessing, 
#' integration, clustering, annotation, differential expression analysis, WGCNA, and enrichment analysis.
#'
#' @param pipeline_file The name of the JSON file containing the pipeline configuration.
#'
#' @details
#' This function sources required scripts, loads the Seurat object if necessary, and runs the specified steps 
#' of the pipeline based on the settings provided in the JSON file. The pipeline can include the following steps:
#' * Preprocessing
#' * Integration
#' * Clustering
#' * Annotation
#' * Differential expression (DEG)
#' * WGCNA
#' * Enrichment analysis
#' 
#' Additionally, custom scripts can be executed as part of the pipeline.
#'
#' @return None. This function is used for its side effects, which include generating plots, saving Seurat objects, and writing output files.
#'
#' @examples
#' \dontrun{
#' main("pipeline_1.json")
#' }
#'
#' @export
main <- function(pipeline_file = "pipeline_wgcna.json") {
    source("scripts/pipelines.r", local = TRUE)

    # Load settings
    pipeline <- load_settings(paste0("pipelines/", pipeline_file))
    general_settings <- .update_parameters(pipeline$general,  load_settings(pipeline$general$settings_path)$general)

    setup_globals(general_settings$folder_destination, general_settings$data_folder, pipeline$general$settings_path)

    # Load microglia object if needed
    if (!isFALSE(general_settings$seurat_object) && isFALSE(general_settings$preprocessing))   load_seurat_object(general_settings$seurat_object)

    lapply(c("preprocessing", "integration", "clustering", "annotation", "deg", "wgcna", "enrichment", "own_script"), function(name) {
        if (name %in% names(pipeline$pipeline)) 
            pipeline$pipeline[[name]] <- FALSE
    })
    # Pipeline
    if (pipeline$pipeline$preprocessing)    preprocessing(pipeline$global_variables, pipeline$preprocessing)
    if (pipeline$pipeline$integration)      integration(pipeline$global_variables, pipeline$integration)
    if (pipeline$pipeline$clustering)       clustering(pipeline$global_variables, pipeline$clustering)
    if (pipeline$pipeline$annotation)       annotation(pipeline$global_variables, pipeline$annotation)
    if (pipeline$pipeline$deg)              deg(pipeline$global_variables, pipeline$deg)
    if (pipeline$pipeline$wgcna)            wgcna(pipeline$global_variables, pipeline$wgcna)
    if (pipeline$pipeline$enrichment)       enrichment(pipeline$global_variables, pipeline$enrichment)
    if (pipeline$pipeline$own_script)       own_script(pipeline$global_variables, pipeline$own_script)
}

check_packages <- function(list_of_packages) {

    new_packages <- setdiff(list_of_packages, rownames(installed.packages()))

    if (length(new_packages) > 0) {
        for (package in new_packages) {
            if (length(BiocManager::install(package, quietly = TRUE)) == 0) {
                print("Package installed via bioconductor")
            } else {
                install.packages(package)
                print("Package installed via package manager")
            }
        }
    }
    lapply(list_of_packages, library, character.only = TRUE)
}

load_settings <- function(file_path) {
    
    check_packages(c("jsonlite"))
    message("loading config file: ", file_path)
    # Check if the file exists
    if (!file.exists(file_path)) {
        stop(paste("Error: File does not exist at path:", file_path))
    }

    # Read the JSON file and convert it to a list
    json_data <- suppressWarnings(fromJSON(readLines(file_path)))
    
    return(json_data)
}

setup_globals <- function(folder, data_folder, settings_path) {
    library(Seurat)
    
    assign("project_folder", paste0(getwd(), "/"), envir = .GlobalEnv)
    assign("output_folder", folder, envir = .GlobalEnv)

    if (!dir.exists(folder)) {
        dir.create(folder, recursive = TRUE)
        message("Directory created: ", folder)
    }
    assign("data_folder", data_folder, envir = .GlobalEnv)
    if (isFALSE(data_folder)) message("data folder not set.")
    assign("settings_path", settings_path, envir = .GlobalEnv)
}


load_seurat_object <- function(file_path) {
    assign("seurat_object", readRDS(file_path), envir = .GlobalEnv)
}

create_variables <- function(var_list) {

    # Assign in the parent frame(caller env) all the variable passed as list
    for (i in seq(var_list)) {
        assign(names(var_list)[i], var_list[[i]], envir = parent.frame())
    }

}