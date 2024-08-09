# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


source("scripts/init.r")

pipeline <-  "from_beginning_n.json"
main(pipeline)
pipeline <-  "from_beginning_s.json"
main(pipeline)
pipeline <-  "useful/pipeline_enrich_wgcna_subset.json"
main(pipeline)
pipeline <-  "useful/pipeline_deg.json"
main(pipeline)
