# settings 

This document provides an overview of the configuration parameters used in the pipeline. Each section and parameter is briefly described below. The setting file must contain all the default parameters (so that in the pipeline files you need to specify only the ones that change). 
A required parameter is assined to null in the settings file, if the null parameters in settings are not assigned in the pipeline files an error is thrown.


## Clustering

- **create_subset**: `boolean` - Indicates whether to create a subset of the data.
- **cluster_column**: `string` - The column name used for clustering.
- **subset**: `list` - List of subsets to be created.
- **clustering**: `boolean` - Flag to enable or disable clustering.
- **rename_clusters**: `boolean` - Indicates whether to rename clusters.
- **to_correct**: `list` - List of clusters to correct.
- **name**: `string` - Name of the clustering.
- **save**: `boolean` - Flag to save the clustering results.
- **save_path**: `string` - Path to save the clustering results.

## Enrichment

- **markers_path**: `string (required)` - Path to the markers file.
- **count_threshold**: `integer` - Threshold for the count.
- **wgcna_folder**: `boolean` - Indicates if WGCNA folder is used.
- **wgcna_module**: `boolean` - Indicates if WGCNA module is used.
- **modules_significance_table**: `boolean` - Flag to use modules significance table.
- **module_threshold**: `integer` - Threshold for WGCNA modules.
- **wgcna_exclude**: `boolean` - Exclude WGCNA modules.
- **cluster**: `boolean` - Flag to enable clustering for enrichment.
- **minGSsize**: `integer` - Minimum gene set size for enrichment analysis.
- **maxGSsize**: `integer` - Maximum gene set size for enrichment analysis.
- **organism**: `string` - NCBI taxonomy ID for the organism (9606 for Homo sapiens).
- **num_tries**: `integer` - Number of tries for enrichment analysis.
- **raw**: `boolean` - Flag to use raw data.
- **extension_plot**: `string` - File extension for plots.
- **scoring**: `string` - Scoring method for enrichment analysis.
- **n_gene_enrich**: `integer` - Number of genes for enrichment.

## WGCNA (Weighted Gene Co-expression Network Analysis)

- **cluster**: `boolean` - Flag to enable clustering for WGCNA.
- **wgcna_file**: `string` - Path to the WGCNA file.
- **save_net**: `boolean` - Flag to save the WGCNA network.
- **load_net**: `boolean` - Flag to load the WGCNA network.
- **soft_power**: `boolean` - Soft power threshold for WGCNA.
- **subject_pathology_column**: `string` - Column for subject pathology.
- **subject_column**: `string` - Column for subject information.
- **extension_plot**: `string` - File extension for plots.
- **markers_analysis_pd**: `string` - Analysis markers for PD.
- **markers_analysis_gpd**: `string` - Analysis markers for GPD.
- **hub_gene_threshold**: `list` - Threshold for hub genes.
- **which**: `list` - Specifies which analyses to run.
- **type**: `string` - Network type for WGCNA.
- **TOMType**: `string` - Topological Overlap Matrix type.
- **mergeCutHeight**: `float` - Cut height for merging modules.
- **minModuleSize**: `integer` - Minimum module size.
- **wgcna_subjects**: `dictionary` - Mapping of subject IDs to study IDs:
- **regions**: `list[string]` - List of brain regions:
- **regions_plot**: `list` - Regions to plot.
- **data_source_mri**: `string` - Source of MRI data.

## DEG (Differential Expression Genes)

- **save_data**: `boolean` - Flag to save the DEG data.
- **cluster_column**: `string` - Column for clustering.
- **assay**: `string` - Assay type for DEG.
- **method**: `string` - Method for DEG analysis.
- **subset_id**: `string` - Subset identifier.
- **o2**: `string` - Additional parameter for DEG.
- **condition**: `string` - Condition for DEG analysis.
- **nothreshold**: `string` - Flag to disable threshold.
- **control**: `string` - Control group for DEG analysis.
- **condition_column**: `string` - Column for condition.
- **extension_plot**: `string` - File extension for plots.
- **folder**: `boolean` - Flag to use folder for DEG.