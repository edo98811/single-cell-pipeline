# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de

# for debugging (error when using do.call and browser)
# https://stackoverflow.com/questions/44608323/crash-of-debugging-browser-in-r-rstudio-when-called-from-inside-do-call
options(deparse.max.lines = 10)
# install.packages("devtools")
devtools::install_github("immunogenomics/presto")
source("scripts/helper_functions.r")

main <- function(pipeline_file = "pipeline_wgcna.json") {
    source("scripts/pipelines.r", local = TRUE)

    # Load settings
    pipeline <- load_settings(paste0("pipelines/", pipeline_file))
    general_settings <- .update_parameters(pipeline$general,  load_settings(pipeline$general$settings_path)$general)

    setup_globals(general_settings$folder_destination, general_settings$data_folder, pipeline$general$settings_path)

    # Load microglia object if needed
    if (!pipeline$pipeline$preprocessing)   load_seurat_object(general_settings$seurat_object)

    # Pipeline
    if (pipeline$pipeline$preprocessing)    preprocessing(pipeline$global_variables, pipeline$preprocessing)
    if (pipeline$pipeline$integration)      integration(pipeline$global_variables, pipeline$integration)
    if (pipeline$pipeline$clustering)       clustering(pipeline$global_variables, pipeline$clustering)
    if (pipeline$pipeline$annotation)       annotation(pipeline$global_variables, pipeline$annotation)
    if (pipeline$pipeline$deg)              deg(pipeline$global_variables, pipeline$deg)
    if (pipeline$pipeline$wgcna)            wgcna(pipeline$global_variables, pipeline$wgcna)
    if (pipeline$pipeline$enrichment)       enrichment(pipeline$global_variables, pipeline$enrichment)

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