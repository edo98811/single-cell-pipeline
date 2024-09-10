# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


source("scripts/init.r")

pipeline <-  "from_beginning.json"
main(pipeline)
pipeline <-  "from_beginning_s.json"
main(pipeline)
pipeline <-  "pipeline_complete_all.json"
main(pipeline)
pipeline <-  "pipeline_complete_microglia.json"
main(pipeline)
pipeline <-  "pipeline_complete_microglia_selection.json"
main(pipeline)
pipeline <-  "pipeline_hdwgcna_subset_run.json"
main(pipeline)
