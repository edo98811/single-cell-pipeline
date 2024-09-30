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



## Pipeline Definition

This document outlines the structure and functionality of a pipeline for data preprocessing, integration, and analysis. The pipeline consists of multiple steps, each of which can be executed independently by setting corresponding parameters in the configuration file.

---

### 1. Preprocessing
- **Enabled**: `[true | false]`

The preprocessing step loads data from a specified folder, processes it with DoubletFinder, and creates Seurat objects. If the `save` parameter is enabled, these objects can be saved in the `save_name`location.

### 2. Integration
- **Enabled**: `[true | false]`

Following preprocessing, Seurat objects can be merged and integrated using Harmony. As with preprocessing, the final Seurat object can be saved in the `save_name`location if the `save` parameter is set.

---

## Main Pipeline

This section covers basic operations such as clustering, dimensionality reduction plotting, and cluster annotation. Each operation is controlled by specific parameters and can be run multiple times.

### 1. Plotting: 
method: `plotting`
Visualize the results of clustering or annotation using UMAP dimensionality reduction.
- **umap**: Specifies the name of the reduction slot for UMAP. Defaults to `umap` if not provided.
- **compute_UMAP**: `[true | false]` — Whether to compute a new UMAP. If `false`, a pre-existing UMAP is required.
- **reduction**: Specifies the reduction method for UMAP. Defaults to the global `r` parameter.
- **cluster_column**: Specifies the column to use for cluster names. Defaults to the global `c` parameter.
- **save**: Whether to save the resulting object.
- **save_name**: The name under which the object will be saved.

### 2. Annotation
method: `annotation`
Clusters can be annotated either automatically using `scType` or manually using a defined mapping.
- **auto**: `[true | false]` — When `true`, clusters are automatically annotated with `scType`. If no manual annotation is provided via the `to_annotate` parameter, this must be set to `true`.
- **to_annotate**: Mapping of existing cluster names to new names, saved in the `annotation_column`.  Defined as r list, (like python dictionary), in JSON as a standard object.
- **annotation_column**: The column where annotations are saved (applies to both automatic and manual annotations).
- **cluster_column**: The source column for clustering results, used in both automatic and manual annotation. Defaults to the global `c` parameter.
- **save**: Whether to save the resulting object.
- **save_name**: The name under which the object will be saved.

### 3. Clustering
method: `clustering`
 Executes Louvain clustering with Seurat’s algorithm, storing results in the metadata in the given column
- **desired_resolution**: The resolution for Louvain clustering. Defaults to the global `dsrd` parameter.
- **reduction**: The reduction method used for clustering. Defaults to the global `r` parameter.
- **cluster_column**: Column where clustering results will be stored. Defaults to the global `c` parameter.
- **save**: Whether to save the resulting object.
- **save_name**: The name under which the object will be saved.
---

## Subsetting

The Seurat object can be subset based on specific clusters, with the option for manual correction. The parameters for this step are:

- **save**: Whether to save the resulting object.
- **save_name**: The name under which the object will be saved.
- **cluster_column**: The cluster column used for subsetting.
- **subset**: Specifies which clusters to keep in the subset, defined as array.

---

## Differential Expression (DEG)

Different methods are available for identifying differentially expressed genes (DEGs). Each method can be run multiple times with different settings.

### DEG Parameters:
- **feature_plots_top9_deg**: make a plot of the top 9 genes according to adjusted pvalue, the umap used is the one defined in the global variables
- **method**: the method to use, an option for a normal markers analysis can be selected or one of the additional methods (volcano, paper, other_plots_from_df, heatmap)

## Other possibilities (visualisation)
- **heatmap**: Generates a heatmap from a list of markers or from DEGs if no markers are provided.
  - **cluster_column**: The column containing the annotation for the plot.
  - **markers**: A list of markers for the heatmap.
  - **clusters**: Specifies clusters to visualize.
  - **maxn_genes**: Maximum number of genes to plot (default 100).
  - **n_genes**: Number of genes per cluster (default 25).
  - **maxn_genes_per_plot**: Maximum genes per plot (default 100).
  - **sorting_method**: method to sort the found markers for filtering, only the top n are plotted according to this metric.

- **volcano**: Creates a volcano plot from the table in the `folder` specified for DEG analysis.
- **plots_misc**: Generates plots for various Seurat object properties.
  - **markers**: A list of markers for feature and ridge plots.
  - **cluster_column**: Metadata column for the plot.
  - **subplot_n**: Number of subplots for the feature plot.
  - **which_other**: Specifies additional plots such as bar or pie charts. (names to be given in an r array)
  - **umap_name**: umap reduction to use for visualisation, if not given uses default in global parameters
- **other_plots_from_df**: Plots markers from a dataframe. Possibilities are:
      - `feature plot` plot, setting the **feature_plot** variable to true, also the number of genes per feature plots can be selected
      - `heatmap`, setting the **plot_heatmap** variable to true, this can be done with one heatmapp for each table column or for all in one heatmap
      parameters needed 
  - **plot_heatmap**: To be set to true or false
  - **markers**: Excel file containing markers.
  - **cluster_column**: Metadata column for plotting.
  - **subplot_n**: Number of subplots for feature plots.
  - **heatmap**: Option to plot a heatmap.
  - **feature_plot**: Option to plot feature plots.
  - **heatmap_by_column**: If `true`, generates one heatmap per dataframe column.
  - **subplot_n**: maximum number of subplots per feature plot (suggested 9 or 4)
  - **max_feature_plots**: maximum number of feature plots, if there are hundreds of genes in the dataframe the feature plots can be very long to generate, hence this parameter 
  - **max_genes_heatmap**: max genes per heatmap, if there are more genes another heatmap is generated
  - **column_list**: for which column of the dataframe are the plots created? if empty all the columns. if plotting the results of deg set this value to ["gene"]
---

## WGCNA
- **Enabled**: `[true | false]`

WGCNA analysis can be run multiple times with different settings. If the cluster column is not defined the default c (in global settings) is used. Parameters for WGCNA are:

- **heatmap_pathology**: Correlation heatmap.
- **TOM**: Topological overlap matrix (to be developed).
- **dendro**: Dendrogram (to be developed).
- **heatmap_mri**: MRI correlation heatmap.
- **heatmap_zscore**: Z-score correlation heatmap.
- **violin_plots**: Violin plots for gene expression by condition.
- **histogram_plot**: Histogram of gene expression by module.
- **significance_membership_scatter**: Scatter plot of significance and membership.
- **significance_log2fc_scatter**: Scatter plot of significance and log2 fold change.
- **significance_membership_model**: Linear model between significance and membership.

---

## Enrichment

Enrichment analysis can be performed using one of four methods, depending on the parameters set:

1. **General DEG Results**: Enrichment based on all DEGs. If only the folder of the DEG analysis is provided.
2. **Excluding a WGCNA Module**: Runs enrichment analysis excluding a specific WGCNA module. It is set by proiding the parameter **wgcna_exclude**
3. **Single WGCNA Module**: Runs enrichment on selected WGCNA modules. If the selected wgcna modules are provided it is run on those, if they are over the module threshold. If the wgcna modules are not provided  it is run only on those which respect the condition of having a sufficient number of genes in the module (the parameter `module_threshold` decides this, the default value is 500).
4. **Linear Model in WGCNA**: Runs enrichment on WGCNA modules selected via a linear model. If a module significance table computed with the wgcna section is provided.

The enrichment can also be done for a single cluster by setting the parameter **cluester** and if necessary **cluster_column**

The possible methods are: 
1. **GSEA**: using clusterProfiler
2. **ORA**: using clusterProfiler
3. **enrichr**: using rbioapi
4. **panther**: using rbioapi

ClusterProfiler: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html
rbioapi: https://cran.r-project.org/web/packages/rbioapi/vignettes/rbioapi.html

## The settings
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

### Available metods for deg

  - `condition_and_clusters`: Compares markers between clusters.
  - `condition_and_clusters_vf`: Similar, but only considers variable features.
  - `default`: Finds markers between clusters.
  - `default_vf`: Finds markers between clusters but only for variable features.
  - `condition`: Compares markers between different groups.
  - `condition_vf`: Compares markers between groups but only for variable features.

### Available other methods

  - `paper`: 
  - `heatmap`: 
  - `other_plots_from_df`: 
  - `volcano`: 

### Available plotting options in paper

  - `numberofcell_barplot`:  
  - `numberofcell_pie_chart`: 
  - `numberofcell_barplot_subject`: 
  - `numberofcell_pie_chart_subject`: 
  - `numberofcell_pie_chart_cluster_subject`:
  - `numberofcell_pie_chart_cluster_pathology`: 
  - `feature_plot`: plots the feature plots with the given markers

## Own Script

- **`path`** [string | null] - `null`



## To write your own script

to write your own script, the own script parametr needs to be set to true, the global settings can be set and also the specific settings (but in this case there no practical difference between the two, apart from the way they need to be accessed)
The global settings are accessed by simply writing the variable name how it is written in the pipeline file, the sepcific settings are accessed with parameters$variable_name.
If you want the seurat object to be preloaded for you you can just set it in the gloab variables, otherwie just dont set anything or leave it to false. 


---

### Additional Features:
- **Cluster Renaming:** If `rename_clusters` is set to `true` and the `to_correct` parameter is provided, cluster names can be corrected.
- **Saving the Seurat Object:** The Seurat object can be saved at the end of this step with the `save` parameter.
  
Each of these steps can be run independently or together based on the required settings.

