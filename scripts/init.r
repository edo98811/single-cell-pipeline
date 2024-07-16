
# source("scripts/new/seurat_utils.r")
# source("scripts/new/helper_functions.r")
# source("scripts/new/enrichment.r")
# source("scripts/new/wgcna.r")

# here put a secret
Sys.setenv(GITHUB_PAT = '')

# comandi utili per installare pacchetti se mi servon
# rm(list=ls(all=TRUE)[sapply(mget(ls(all=TRUE)), class) == "data.frame"])

# takes as input a list containing packages names and if not present it installs them

check_packages <- function(list_of_packages){

  new_packages <- setdiff(list_of_packages, rownames(installed.packages()))

  if(length(new_packages) > 0) {
    for (package in new_packages){
      if(length(BiocManager::install(package, quietly = TRUE)) == 0) {
        print("Package installed via bioconductor")
      } else {
        install.packages(package)
        print("Package installed via package manager")
      }
    }
  }
  lapply(list_of_packages, library, character.only = TRUE)
}

check_packages("Seurat")

load_settings <- function(file_path) {
    
  check_packages(c("jsonlite"))

  # Read the JSON file and convert it to a list
  json_data <- fromJSON(readLines(file_path))
  
  return(json_data)
}

setup_globals <- function(folder = "microglia_correct"){

  # runID = "test_loading"
  assign("runID", folder, envir = .GlobalEnv)
  assign("project_folder", paste0(getwd(), "/"), envir = .GlobalEnv)
  assign("output_folder", define_output_folder(get("runID", envir = .GlobalEnv)), envir = .GlobalEnv)

  assign("data_folder", paste0(get("project_folder", envir = .GlobalEnv), "input/"), envir = .GlobalEnv)
  assign("pattern", "filtered_feature_bc_matrix", envir = .GlobalEnv)
  assign("patient_info", "input/info_subj.txt", envir = .GlobalEnv)
}

setup_globals_new <- function(folder = "microglia_correct"){

  # runID = "test_loading"
  
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