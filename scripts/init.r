# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


source("scripts/helper_functions.r")

main <- function(pipeline_file = "pipeline_wgcna.json") {
  source("scripts/pipelines.r", local = TRUE)

  pipeline <- load_settings(paste0("pipelines/", pipeline_file))

  setup_globals(pipeline$general$folder_destination)
  load_seurat_object(pipeline$general$microglia_object)

  if (pipeline$pipeline$clustering)       clustering(pipeline$global_variables, pipeline$clustering)
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

    # Read the JSON file and convert it to a list
    json_data <- fromJSON(readLines(file_path))
    
    return(json_data)
}

setup_globals <- function(folder = "microglia_correct") {

    assign("project_folder", paste0(getwd(), "/"), envir = .GlobalEnv)
    assign("output_folder", folder, envir = .GlobalEnv)

    if (!dir.exists(folder)) {
        dir.create(folder, recursive = TRUE)
        message("Directory created: ", folder)
    }
    
    assign("data_folder", paste0(get("project_folder", envir = .GlobalEnv), "input/"), envir = .GlobalEnv)
    assign("pattern", "filtered_feature_bc_matrix", envir = .GlobalEnv)
    assign("patient_info", "input/info_subj.txt", envir = .GlobalEnv)
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