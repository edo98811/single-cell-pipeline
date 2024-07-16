source("scripts/new/init.r")
source("scripts/new/helper_functions.r")

main <- function(pipeline_file = "pipeline_wgcna.json") {
  source("scripts/new/pipelines.r", local = TRUE)

  pipeline <- load_settings(paste0("pipelines/", pipeline_file))

  setup_globals_new(pipeline$general$folder_destination)
  load_seurat_object(pipeline$general$microglia_object)

  if (pipeline$pipeline$clustering)       clustering(pipeline$global_variables, pipeline$clustering)
  if (pipeline$pipeline$deg)              deg(pipeline$global_variables, pipeline$deg)
  if (pipeline$pipeline$wgcna)            wgcna(pipeline$global_variables, pipeline$wgcna)
  if (pipeline$pipeline$enrichment)       enrichment(pipeline$global_variables, pipeline$enrichment)

}


pipeline <-  "useful/pipeline_wgcna_subset.json"
main(pipeline)
pipeline <-  "useful/pipeline_wgcna_subset_run.json"
main(pipeline)
pipeline <-  "useful/pipeline_enrich_wgcna_subset.json"
main(pipeline)
pipeline <-  "useful/pipeline_deg.json"
main(pipeline)
# main_with_switch <- function(option, pipeline, env_variables) {
#   switch(option,
#     "clustering" = {
#       clustering(env_variables, pipeline)
#     },
#     "deg" = {
#       deg(env_variables, pipeline)
#     },
#     "wgcna" = {
#       wgcna(env_variables, pipeline)
#     },
#     "enrichment" = {
#       enrichment(env_variables, pipeline)
#     },
#     {
#       print("Invalid option")
#     }
#   )
# }                                                            


# args = commandArgs(trailingOnly=TRUE)
# if (length(args)==0) {
#  #  main(pipeline)
# } else if (length(args)==1) {
#   # default output file
#   main(args[1])
# }


# "which": ["heatmap", "heatmapMRI", "histogram_plot", "heatmapMRI_zscore"],
# ["heatmap", "heatmapMRI", "heatmapMRI_zscore", "histogram_plot", "correlation_significance_scatter", "significance_membership_scatter"]



