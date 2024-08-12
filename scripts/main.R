# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


source("scripts/init.r")

pipeline <-  "from_beginning_n.json"
main(pipeline)
pipeline <-  "from_beginning_s.json"
main(pipeline)
pipeline <-  "pipeline_wgcna_subset_run.json"
main(pipeline)
pipeline <-  "useful/pipeline_deg.json"
main(pipeline)
