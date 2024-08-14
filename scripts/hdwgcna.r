
hdwgcna <- function(seurat_obj, name = "test", ...) {



    # devtools::install_github("NightingaleHealth/ggforestplot")
    # devtools::install_github("smorabit/hdWGCNA", ref="dev")

    library(hdWGCNA)
    library(cowplot)
    library(patchwork)
    # enableWGCNAThreads(nThreads = 8)
    # tutorial: https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html
    # Set folders
    output_dir <- set_up_output(paste0(output_folder, "wgcna_", name, "/"), message)
    
    browser()

    # Select the genes that are expressed in at least 5% of cells in the test field in misc seurat obkect
    seurat_obj <- SetupForWGCNA(
        seurat_obj,
        gene_select = "fraction", # the gene selection approach
        fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
        wgcna_name = "test" # the name of the hdWGCNA experiment
    )
    # construct metacells  in each group ( group by to keep the cells coming from these groups in the same metacells)
    seurat_obj <- MetacellsByGroups(
        seurat_obj = seurat_obj,
        group.by = c("subject", "microglia_clusters"), # specify the columns in seurat_obj@meta.data to group by
        reduction = "harmony_reduction", # select the dimensionality reduction to perform KNN on
        k = 25, # nearest-neighbors parameter
        max_shared = 10, # maximum number of shared cells between two metacells
        ident.group = "microglia_clusters" # set the Idents of the metacell seurat object
    )

    # The group.by parameter determines which groups metacells will be constructed in. We 
    # only want to construct metacells from cells that came from the same biological sample of origin,
    # so it is critical to pass that information to hdWGCNA via the group.by parameter. 
    # Additionally, we usually construct metacells for each cell type separately. 
    # Thus, in this example, we are grouping by Sample and cell_type to achieve the 
    # desired result.
    
    # normalize metacell expression matrix:
    seurat_obj <- NormalizeMetacells(seurat_obj)

    # Questo posso anche non metterlo, serve se voglio fare il wgcna only per un subset delle cellule
    # Pero questo é molto utile prché riduce il lavoro che ho fatto per la wgcna normale 
    # quando ad esempio voglio lavorare solo con un subset
    seurat_obj <- SetDatExpr(
        seurat_obj,
        group_name = c(
            "PD_001",
            "PD_002",
            "PD_005",
            "PD_007",
            "PD_008",
            "PD_009",
            "PD_012",
            "PD_016",
            "PD_017"
        ),
        group.by = "subject",
        assay = "RNA", # using RNA assay
        slot = "data" # using normalized data
    )

    # Test different soft powers:
    seurat_obj <- TestSoftPowers(
        seurat_obj,
        networkType = "unsigned" # you can also use "unsigned" or "signed hybrid"
    )
    
    # plot the results:
    plot_list <- PlotSoftPowers(seurat_obj)

    # assemble with patchwork (to save with save plot)
    save_plot(wrap_plots(plot_list, ncol = 2), paste0(output_dir, "/module_feature_plot1.png"))

    # construct co-expression network:
    seurat_obj <- ConstructNetwork(
        seurat_obj,
        soft_power = 10,
        tom_name = "test", # name of the topoligical overlap matrix written to disk
        networkType = "unsigned",
        overwrite_tom = "true"
    )

    save_plot(wrap_plots(PlotDendrogram(seurat_obj, main = "TEST_1"), paste0(output_dir, "/dendro_test.png")))
    # need to run ScaleData first or else harmony throws an error:
    #seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

    # compute all MEs in the full single-cell dataset
    seurat_obj <- ModuleEigengenes(
        seurat_obj,
        group.by.vars = "Sample",
    )
    # harmonized module eigengenes:
    hMEs <- GetMEs(seurat_obj)

    # module eigengenes:
    MEs <- GetMEs(seurat_obj, harmonized=FALSE)

    # make a featureplot of hMEs for each module
    plot_list <- ModuleFeaturePlot(
        seurat_obj,
        features="hMEs", # plot the hMEs
        order=TRUE # order so the points with highest hMEs are on top
    )

    # stitch together with patchwork
    save_plot(wrap_plots(plot_list, ncol=6), paste0(output_dir, "/module_feature_plot1.png"))
    
    # make a featureplot of hub scores for each module
    plot_list <- ModuleFeaturePlot(
        seurat_obj,
        features="scores", # plot the hub gene scores
        order="shuffle", # order so cells are shuffled
        ucell = TRUE # depending on Seurat vs UCell for gene scoring
    )

    # stitch together with patchwork
    save_plot(wrap_plots(plot_list, ncol=6), paste0(output_dir, "/module_feature_plot2.png"))

    correlation_heatmaps()
}

correlation_heatmaps <- function() {
    wgcna_subjects <- list(
        "PD_001" = "02_082",
        "PD_002" = "02_084",
        "PD_005" = "02_074",
        "PD_007" = "02_096",
        "PD_008" = "02_097",
        "PD_009" = "02_105",
        "PD_012" = "02_104",
        "PD_016" = "02_108",
        "PD_017" = "02_115"
        )

    regions <- c(
        "superiorfrontal",
        "caudalmiddlefrontal", 
        "rostralmiddlefrontal"
        )

    source("scripts/wgcna_plots_helpers.r", local = TRUE)

    traits <- load_mri_traits(wgcna_subjects = wgcna_subjects, regions = regions)

    detailed_regions <- colnames(traits)
    traits$subject <- rownames(traits)

    seurat_obj@meta.data <- merge(seurat_obj@meta.data, traits, by = "subject", all.x = TRUE)
    seurat_obj@meta.data <- seurat_obj@meta.data[, 1:14]
    # colnames(seurat_obj@meta.data)[1:14]

    seurat_obj <- ModuleTraitCorrelation(
        seurat_obj,
        traits = detailed_regions,
        group.by = "subject"
    )

    mt_cor <- GetModuleTraitCorrelation(seurat_obj)
    excel_filename <- paste0(output_dir, "test.xlsx")

    wb <- write_on_excel("correlation", as.data.frame(mt_cor$cor$all_cells), mode = "colorscale")
    wb <- write_on_excel("p.value", as.data.frame(mt_cor$p$all_cells), wb = wb)
    wb <- write_on_excel("t.statistic", as.data.frame(mt_cor$t$all_cells), wb = wb)
    openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)

    make_heatmap(cor_data)
}