# Pipeline to simplify the single-cell processing with seurat
It works using pipeline files to describe the steps that are needed, the script is then simply launched from main like this: 

```
source("R/init.r")

main("path/pipeline_file.json")
```

# Subject info file description
To load the information about the subjects a file with the necessary info must be written, it is a text file with a very easy syntax. The path to the file must then be specified in the loading page. The spaces between the - must also be respected

example
```
subjectID - group
...
```

# Pipeline File

It specifies the parameters that are used in the pipeline, the settable parameters are described in settings file. 

## General information 
The pipeline sections are executed in the order in which they are described below. a "pipeline" section is necessary with at least one element set to true. The general settings section is also required, with at least the destination folder and the path to the settings file defined. The settings file contains all the default settings for the pipelines, if they are not overridden in the pipeline file, these are the settings that are used. The settings file needs to be created for each project, the pipeline files for each analysis that needs to be performed.

## Pipeline definition
The structure is provided in the examples. at the beginning the single sections that will be run need to be selected. by specifying it like this in a `pipeline` section:

- `preprocessing`: [true | false]
Preprocessing step loads the data from a given folder, then it performs the preprocessing sstep with doubletfiner, the seurat objects are then optionally saved with the save parameter.
- `integration`: [true | false]
The step after preprocessing, here the list of seurat objects created in the previous step can be merged and integrated using harmony. The seurat object can be saved at the end of this step with the save parameter.
- `clustering`: [true | false]
What the clustering section does is create a subset if the parameter is requested, and then work with this subset. Cluster the seurat iobject with the given parameters, the results are plotted. If rename_clusters is set to true and the to_correct parameter is given, then the annotation can also be corrected. The seurat object can be saved at the end of this step with the save parameter.
Te steps that can be performed by this function are: (by setting the correct marker to true or false in the settings section)
  - create subset: create a subset of the seurat object to be used for the whole pipeline, the subset parameter is mandatory
  - plot_before_clustering: wether to plot or not (after subsetting in this case) without computing umap
  - clustering: perform a new clustering 
  - clustering_plot: plot, performing a new clustering 
  - rename_clusters: the clusters are renamed 
All of these or only one of these can be ran. The order in which they are ran is the one defined above
- `annotation`: [true | false]
The annotation will be done with sctype, the parameters for the annotation are fixed on brain datasets but will be changed in the future. After annotation if correct is set to true and the to_correct parameter is given a manual correction can be performed. The annotation can also be competely manual if the parameter annoatate is set to false.  The seurat object can be saved at the end of this step with the save parameter.
- `deg`: [true | false]
The differentially expressed genes are found and the results saved, different commadns can be used, multiple analysis can be ran at the same time creating more than one field in the setting file
deg methods
    - `condition_and_clusters`: compute the markers betweeen the clusters
    - `condition_and_clusters_vf`: only considering variable features
    - `default`: compute the markers between the clusters 
    - `default_vf`: only considering variable features
    - `condition`: compute the markers beween the groups
    - `condition_vf`: only considering variable features
    other
    - `heatmap`: this creates an heatmap, either from a list of markers or from the differentially expressed markers (if the markers parameters is not defined) The needed parameters for this are: cluster_column, markers
    - `volcano`: this method is used to create a volcano plot in the folder of the deg analysis that is set as "folder" parameter, folder parameter is only used for this method
    - `paper`: plot different aspects of the seurat object
      parameters:  
        - `markers`: a list of markers
        - `cluster_column`: the column of metadata that is used for Idents
        - `subplot_n`: number of subplot per featurep plot
        - `which_other`: which plot to create, options
          - `numberofcell_barplot`:  
          - `numberofcell_pie_chart`: 
          - `numberofcell_barplot_subject`: 
          - `numberofcell_pie_chart_subject`: 
          - `numberofcell_pie_chart_cluster_subject`:
          - `numberofcell_pie_chart_cluster_pathology`: 
    - `plot_markers_from_df`: different plots from lists of markers that are in an excel file, possibilities:
      - feature plot (default)
      - heatmap, setting the plot_heatmap variable to true, this can be done with one heatmapp for each table column or for all in one heatmap
      parameters needed 
      - `markers`: in this case it must be an excel file
      - `cluster_column`: the column of metadata that is used for Idents
      - `subplot_n`: number of subplot per featurep plot
      - `heatmap`: wether to plot the heatmap
      - `feature_plot`: weather to plot the feature plot
      - `heatmap_by_column`: true if you want to create a different heatmap for each column of the dataframe. if false it printsall the genes in the excel file in the same heatmap
      - `subplot_n`: maximum number of subplots per feature plot (sugested 9 or 4)
      - `max_feature_plots`: maximum number of feature plots, if there are hundreds of genes in the dataframe the feature plots can be very long to generate, hence this parameter 
      - `max_genes_heatmap`: max genes per heatmap, if there are more genes another heatmap is generated
      - `column_list`: for which column of the dataframe are the plots created? if empty all the columns. if plotting the results of deg set this value to ["gene"]

- `wgcna`: [true | false]
To run the wgcna analysis, the same concept as the deg applies, the analysis can be ran multiple times with multiple different settings. Paraeters are described in the following section, 
possible plots or analysis are: 
  - `heatmap_pathology`: correlation heatmap for 
  - `TOM`: to develop
  - `dendro`: to develop
  - `heatmap_mri`: correlation heatmap for 
  - `heatmap_zscore`: correlation heatmap for 
  - `violin_plots`: violin plots for gene expression divided by condition
  - `histogram_plot`: histogram plot of gene expression by module
  - `histogram_plot_significance`: histogram plot of gene significance by module 
  - `significance_membership_scatter`: scatter plot of these two
  - `significance_log2fc_scatter`: scatter plot of these two
  - `correlation_avglog2fc_scatter`: scatter plot of these two
  - `corr_matrix`: to_develop
  - `significance_membership_model`: computes a linear model for each module between these two
if the cluster column is not defined the default (in general settings) is used 

- `enrichment`: [true | false]
to run enrichment analysis, here too the analysis can be ran multiple times. Enrichment analysis is ran in four different methods, depending on the parameters that are set
- on general deg results: for this te only parameter that needs to be set is the mraers path and eventually the cluster, cluster can be a single value or a list, or not set. If it is not set the enrichment is ran on all the files present in the folder given as deg results
- excluding a wgcna module: on all the genes excluding a single module the wgcna folder must be set, in the folder the module_genes.xlsx file (automatically created after running the wgcna) must be present
- on single wgcna modules: on a list of selected modules, the wgcna folder must be set, in the folder the module_genes.xlsx file (automatically created after running the wgcna) must be present 
- on selected wgcna modules according to the linear model computed in wgcna, for this selection the file needs to be specified, and the wgcna_folder as well. The modules selected are the ones that have a pvlaue less than 0.5 % in the model parameter estimation 

The settings
for `preprocessing`,`integration`,`clustering`,`annotation`: need to be specified in the apposite section 
for: `deg`,`wgcna`,`enrichment`: need to be specified in the subsection for the name of the analysis, for each pipeline can be run multiple times. The folder with the results will be called with the same name of the subsection in which the settings are in.

tip: if you create a section in the settings that does not have one of the names here descrbed it will not be run
example: 
```
   "enrichment": {
      "WGCNA_all": {
        "markers_path": "markers_microglia_control_vs_pd",
        "modules_significance_table":"linear_model_significance_membership_single_module_frontal_cortex_thickness.xlsx",
        "wgcna_folder": "WGCNA_all/"
      },
      "WGCNA_0": {
        "markers_path": "markers_microglia_control_vs_pd_clusters",
        "cluster": 2,
        "modules_significance_table":"linear_model_significance_membership_single_module_frontal_cortex_thickness.xlsx",
        "wgcna_folder": "WGCNA_cluster0/"
      },
      "WGCNA_2": {
        "markers_path": "markers_microglia_control_vs_pd_clusters",
        "cluster": 2,
        "modules_significance_table":"linear_model_significance_membership_single_module_frontal_cortex_thickness.xlsx",
        "wgcna_folder": "WGCNA_cluster2/"
      }
    }
```

# Settings File

The settings file stores the default parameters for each pipeline that is run. All the parameters specified in the pipeline files are overwritten when the pipeline is run. This file can serve as template when creating a pipeline file. 

## General

- **`seurat_object`** [boolean] - `false`
- **`folder_destination`** [string | null] - `null`
- **`data_folder`** [boolean] - `false`
- **`settings_path`** [string | null] - `null`

## Global Variables

- **`count_matrix_pattern`** [string] - `"filtered_feature_bc_matrix"`
- **`r_before_integration`** [string] - `"pca"`
- **`c_before_integration`** [string] - `"clusters"`
- **`patient_info`** [string] - `"input/info_subj.txt"`
- **`r`** [string] - `"harmony_reduction"`
- **`c`** [string] - `"harmony_clusters"`
- **`a`** [string] - `"RNA"`
- **`umap`** [string] - `"umap_microglia_harmony_reduction"`
- **`dstd`** [float] - `1.5`
- **`d`** [float] - `0.3`
- **`name`** [string] - `"only_omics_subjects"`
- **`annotation`** [string] - `"microglia_clusters"`
- **`corrected_annotation`** [string] - `"microglia_clusters_scType"`
- **`extension_plot`** [string] - `".png"`

## Preprocessing

- **`save`** [boolean] - `true`
- **`save_name`** [string] - `"after_preprocessing"`
- **`parts_to_remove`** [string] - `"_filtered_feature_bc_matrix"`

## Integration

- **`save`** [boolean] - `true`
- **`save_name`** [string] - `"after_integration.rds"`
- **`method`** [string] - `"harmony"`

## Annotation

- **`annotate`** [boolean] - `true`
- **`annotation_plots`** [boolean] - `true`
- **`annotation_column`** [string] - `"scType_annotation"`
- **`corrected_annotation`** [object] - `{ "13": "Border macrophages" }`
- **`corrected_annotation_column`** [string] - `"scType_annotation_corrected"`
- **`correct`** [boolean] - `false`
- **`corrected_annotation_plots`** [boolean] - `true`
- **`save`** [boolean] - `true`
- **`save_name`** [string] - `"after_annotation.rds"`

## Clustering

- **`create_subset`** [boolean] - `false`
- **`cluster_column`** [string] - `"clusters"`
- **`subset`** [array] - `[]`
- **`clustering`** [boolean] - `false`
- **`clustering_plot`** [boolean] - `true`
- **`rename_clusters`** [boolean] - `false`
- **`to_correct`** [array] - `[]`
- **`name`** [string] - `"open_data"`
- **`save`** [boolean] - `true`
- **`save_name`** [string] - `"after_clustering"`

## Enrichment

- **`markers_path`** [string | null] - `null`
- **`count_threshold`** [integer] - `0`
- **`wgcna_folder`** [boolean] - `false`
- **`wgcna_module`** [boolean] - `false`
- **`modules_significance_table`** [boolean] - `false`
- **`module_threshold`** [integer] - `500`
- **`wgcna_exclude`** [boolean] - `false`
- **`cluster`** [boolean] - `false`
- **`minGSsize`** [integer] - `15`
- **`maxGSsize`** [integer] - `500`
- **`organism`** [string] - `"9606"`
- **`num_tries`** [integer] - `3`
- **`raw`** [boolean] - `false`
- **`extension_plot`** [string] - `".png"`
- **`scoring`** [string] - `"log2FC"`
- **`n_gene_enrich`** [integer] - `300`

## WGCNA

- **`method`** [string] - `"WGCNA"`
- **`cluster`** [list] - `false`
- **`wgcna_file`** [string] - `"bwnet.rds"`
- **`save_net`** [boolean] - `false`
- **`load_net`** [boolean] - `false`
- **`soft_power`** [boolean] - `false`
- **`subject_pathology_column`** [string] - `"subject_pathology"`
- **`subject_column`** [string] - `"subject"`
- **`extension_plot`** [string] - `".png"`
- **`markers_analysis_pd`** [string] - `"microglia_control_vs_pd_clusters"`
- **`markers_analysis_gpd`** [string] - `false`
- **`hub_gene_threshold`** [array] - `[0.3, 1]`
- **`which`** [array] - `[]`
- **`type`** [string] - `"unsigned"`
- **`TOMType`** [string] - `"unsigned"`
- **`mergeCutHeight`** [float] - `0.25`
- **`minModuleSize`** [integer] - `50`
- **`wgcna_subjects`** [object] - `{ "PD_001": "02_082",}`
- **`regions`** [array] - `["superiorfrontal", "caudalmiddlefrontal", "rostralmiddlefrontal"]`
- **`regions_plot`** [array] - `[]`
- **`data_source_mri`** [string] - `"aparc"`

### Available WGCNA Plots

- `heatmap_pathology`
- `TOM`
- `dendro`
- `heatmap_mri`
- `heatmap_zscore`
- `violin_plots`
- `histogram_plot`
- `histogram_plot_significance`
- `significance_membership_scatter`
- `significance_log2fc_scatter`
- `correlation_avglog2fc_scatter`
- `corr_matrix`
- `significance_membership_model`

## DEG (Differential Expression Gene) Analysis

- **`save_data`** [boolean] - `true`
- **`cluster_column`** [string] - `"harmony_clusters"`
- **`assay`** [string] - `"RNA"`
- **`method`** [string] - `"default"`
- **`subset_id`** [string] - `"control"`
- **`o2`** [string] - `"NA"`
- **`condition`** [string] - `"PD"`
- **`nothreshold`** [boolean] - `false`
- **`control`** [string] - `"non_PD"`
- **`condition_column`** [string] - `"subject_pathology"`
- **`extension_plot`** [string] - `".png"`
- **`folder`** [boolean] - `false`

## Own Script

- **`path`** [string | null] - `null`



## To write your own script

to write your own script, the own script parametr needs to be set to true, the global settings can be set and also the specific settings (but in this case there no practical difference between the two, apart from the way they need to be accessed)
The global settings are accessed by simply writing the variable name how it is written in the pipeline file, the sepcific settings are accessed with parameters$variable_name.
If you want the seurat object to be preloaded for you you can just set it in the gloab variables, otherwie just dont set anything or leave it to false. 