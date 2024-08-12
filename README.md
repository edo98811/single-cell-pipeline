# Pipeline to simplify the single-cell processing with seurat
It works using pipeline files to describe the steps that are needed, the script is then simply launched from main.

# Subject info file description
subjectID - group
...

# Pipeline File

It specifies the parameters that are used in the pipeline, the settable parameters are described in settings file. 

The structure is provided in the examples. at the beginning the single sections that will be run need to be selected. by specifying it like this in a `pipeline` section:

- `preprocessing`: [true | false]
- `integration`: [true | false]
- `clustering`: [true | false]
- `annotation`: [true | false]
- `deg`: [true | false]
- `wgcna`: [true | false]
- `enrichment`: [true | false]

The settings
for `preprocessing`,`integration`,`clustering`,`annotation` -> need to be specified in the apposite section 
for: `deg`,`wgcna`,`enrichment`: need to be specified in the subsection for the name of the analysis, for each pipeline these can be run multiple times, the folder with the results will be called with the same name of the subsection in which the settings are in.

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

- `seurat_object` [string] - Path to the Seurat object file. Required.
- `folder_destination` [string | null] - Destination folder for outputs. Required.
- `data_folder` [string] - Path to the folder containing input data.
- `settings_path` [string | null] - Path to the settings file. Required.

## Global Variables

- `count_matrix_pattern` [string] - Pattern used to identify count matrix files.
- `r_before_integration` [string] - Reduction method before integration.
- `c_before_integration` [string] - Clustering method before integration.
- `patient_info` [string] - Path to the patient information file.
- `r` [string] - Reduction method used for integration.
- `c` [string] - Clustering method used after integration.
- `a` [string] - Assay used for analysis.
- `umap` [string] - UMAP embedding used for visualization.
- `dstd` [float] - Standard deviation threshold for filtering.
- `d` [float] - Distance threshold for clustering.
- `name` [string] - Name of the dataset or analysis.
- `annotation` [string] - Column name for initial annotation.
- `corrected_annotation` [string] - Column name for corrected annotation.
- `extension_plot` [string] - File extension for saved plots.

## Preprocessing

- `save` [boolean] - Whether to save the preprocessing output.
- `save_name` [string] - Name of the file to save after preprocessing.
- `parts_to_remove` [string] - Parts of file names to remove during preprocessing.

## Integration

- `save` [boolean] - Whether to save the integration output.
- `save_name` [string] - Name of the file to save after integration.
- `method` [string] - Method used for integration (e.g., "harmony").

## Annotation

- `annotate` [boolean] - Whether to perform annotation.
- `annotation_plots` [boolean] - Whether to generate plots for annotation.
- `annotation_column` [string] - Column name for annotation.
- `corrected_annotation` [object] - Dictionary of corrected annotations.
- `correctet_annotation_column` [string] - Column name for corrected annotations.
- `correct` [boolean] - Whether to apply corrections to annotations.
- `corrected_annotation_plots` [boolean] - Whether to generate plots for corrected annotations.
- `save` [boolean] - Whether to save the annotation output.
- `save_name` [string] - Name of the file to save after annotation.

## Clustering

- `create_subset` [boolean] - Whether to create a subset of the data.
- `cluster_column` [string] - Column name for clustering.
- `subset` [array] - List of clusters to subset.
- `clustering` [boolean] - Whether to perform clustering.
- `rename_clusters` [boolean] - Whether to rename clusters.
- `to_correct` [array] - List of clusters to correct.
- `name` [string] - Name of the clustering analysis.
- `save` [boolean] - Whether to save the clustering output.
- `save_name` [string] - Name of the file to save after clustering.

## Enrichment

- `markers_path` [string | null] - Path to markers file. Required.
- `count_threshold` [integer] - Threshold for counting markers.
- `wgcna_folder` [boolean] - Whether to use WGCNA folder.
- `wgcna_module` [boolean] - Whether to use WGCNA modules.
- `modules_significance_table` [boolean] - Whether to generate a significance table for modules.
- `module_threshold` [integer] - Threshold for modules.
- `wgcna_exclude` [boolean] - Whether to exclude certain WGCNA results.
- `cluster` [boolean] - Whether to perform clustering in enrichment analysis.
- `minGSsize` [integer] - Minimum gene set size for enrichment.
- `maxGSsize` [integer] - Maximum gene set size for enrichment.
- `organism` [string] - Organism for enrichment analysis (e.g., "9606" for human).
- `num_tries` [integer] - Number of tries for enrichment analysis.
- `raw` [boolean] - Whether to use raw data for enrichment analysis.
- `extension_plot` [string] - File extension for saved enrichment plots.
- `scoring` [string] - Scoring method for enrichment analysis (e.g., "log2FC").
- `n_gene_enrich` [integer] - Number of genes to use for enrichment.

## WGCNA

- `method` [string] - Method used for WGCNA analysis.
- `cluster` [boolean] - Whether to perform clustering in WGCNA.
- `wgcna_file` [string] - Path to the WGCNA file.
- `save_net` [boolean] - Whether to save the network generated by WGCNA.
- `load_net` [boolean] - Whether to load an existing WGCNA network.
- `soft_power` [boolean] - Whether to calculate soft power for WGCNA.
- `subject_pathology_column` [string] - Column name for subject pathology.
- `subject_column` [string] - Column name for subject identification.
- `extension_plot` [string] - File extension for saved WGCNA plots.
- `markers_analysis_pd` [string] - Name of the markers analysis for PD.
- `markers_analysis_gpd` [string] - Name of the markers analysis for GPD.
- `hub_gene_threshold` [array] - Thresholds for hub gene identification.
- `which` [array] - List of modules to include in WGCNA analysis.
- `type` [string] - Type of network (e.g., "unsigned").
- `TOMType` [string] - Type of Topological Overlap Matrix (TOM).
- `mergeCutHeight` [float] - Height cutoff for merging modules.
- `minModuleSize` [integer] - Minimum module size for WGCNA.
- `wgcna_subjects` [object] - Dictionary mapping subject IDs to WGCNA IDs.
- `regions` [array] - List of regions included in the analysis.
- `regions_plot` [array] - List of regions to include in plots.
- `data_source_mri` [string] - Data source for MRI-based analysis.

## DEG (Differential Expression Gene) Analysis

- `save_data` [boolean] - Whether to save the DEG analysis output.
- `cluster_column` [string] - Column name for clustering in DEG analysis.
- `assay` [string] - Assay used for DEG analysis.
- `method` [string] - Method used for DEG analysis (e.g., "default").
- `subset_id` [string] - ID of the subset for DEG analysis.
- `o2` [string] - Placeholder parameter, currently set to "NA".
- `condition` [string] - Condition to compare in DEG analysis (e.g., "PD").
- `nothreshold` [string] - Whether to apply thresholds in DEG analysis (set to "false").
- `control` [string] - Control group in DEG analysis (e.g., "non_PD").
- `condition_column` [string] - Column name for conditions in DEG analysis.
- `extension_plot` [string] - File extension for saved DEG plots.
- `folder` [boolean] - Whether to save DEG results in a separate folder.


# Pipeline file Documentation 

The pipeline file is the one that defined the 