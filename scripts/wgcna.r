# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de


# WGCNA analysis Arguments:
#' @param  name
#' @param  wgcna_file
#' @param  save_net
#' @param  load_net
#' @param  soft_power
#' @param  subject_pathology_column
#' @param  subject_column
#' @param  extension_plot
#' @param  markers_analysisGPD
#' @param  markers_analysisPD
#' @param  hub_genes_list
#' @param  hub_genes_threshold
#' @param  module_of_interest
#' @param  which
#' @return seurat_object


# Analysis using the WGCNA package
WGCNA_main <- function(seurat_object, name = "test", wgcna_file = "bwnet.rds", save_net = TRUE, 
                        load_net = FALSE, soft_power = NA, ...) {
    {

        # Check types of arguments
        if (!inherits(seurat_object, "Seurat")) {
            stop("seurat_object must be a Seurat object")
        }

        if (!is.character(name)) {
            stop("name argument must be a character")
        }

        if (!is.character(wgcna_file)) {
            stop("wgcna_file argument must be a character")
        }

        if (!is.logical(save_net)) {
            stop("save_net argument must be a logical value")
        }

        if (!is.logical(load_net)) {
            stop("load_net argument must be a logical value")
        }

        if (!is.na(soft_power) && !is.numeric(soft_power)) {
            stop("soft_power argument must be a numeric value or NA")
        }


    }

    cat("Running WGCNA...\n")
    message(paste0("Parameters: name: ", name, " - wgcna_file: ", wgcna_file, " - save_net: ",
        save_net, " - load_net: ", load_net, " - soft_power: ", soft_power))

    # Dependencies
    # check_packages(c("Seurat", "WGCNA", "GEOquery", "tidyverse", "gridExtra", "dplyr",
    #     "readxl", "openxlsx", "DESeq2", "ggpmisc"))
    library(magrittr)

    # Set folders
    output_dir <- set_up_output(paste0(output_folder, "WGCNA_", name, "/"), message)

    # Main loop Data preparation
    column_data <- prepare_column_data(seurat_object, ...)
    norm_counts <- prepare_data(seurat_object, column_data, ...)

    # Soft power
    if (is.na(soft_power) & !load_net)
        soft_power <- soft_power_intuition(norm_counts, ...)

    # Matrix computation
    if (load_net) bwnet <- readRDS(paste0(output_dir, wgcna_file)) 
    else bwnet <- matrix_computation(norm_counts, soft_power, save_net, output_dir, wgcna_file, ...)

    # Plots and other functions
    source("scripts/new/wgcna_plots.r", local = TRUE)
    plot_functions(bwnet, norm_counts, column_data, output_dir, ...)
    save_module_genes(bwnet, norm_counts, column_data, output_dir)

    return(seurat_object)
}

# Prepare info on subjects
prepare_column_data <- function(seurat_object, subject_pathology_column = "subject_pathology", subject_column = "subject", ...) {

    # Extract metadata from seurat object
    meta.data <- seurat_object@meta.data

    # Create subject meta data
    column_data <- data.frame(subject = unique(meta.data[[subject_column]]), row.names = unique(meta.data[[subject_column]]))
    column_data$pathology <- lapply(unique(meta.data[[subject_column]]), function(x) unique(meta.data[meta.data[[subject_column]] ==
        x, ][[subject_pathology_column]]))

    return(column_data)
}

# Prepare the count matrix
prepare_data <- function(seurat_object, column_data, subject_column = "subject", ...) {

    library(DESeq2)

    # Get data from seurat object (transpose to have genes on columns and cells
    # on rows)
    counts_data <- t(GetAssayData(object = seurat_object, assay = "RNA", layer = "counts"))
    meta.data <- seurat_object@meta.data

    # Group counts by subject
    subject_data <- data.frame(matrix(ncol = 0, nrow = length(colnames(counts_data))),
        row.names = colnames(counts_data))
    for (subject in unique(meta.data[[subject_column]])) {
        subject_data <- cbind(subject_data, colSums(counts_data[row.names(meta.data[meta.data[[subject_column]] ==
            subject, ]), ]))
        colnames(subject_data)[length(colnames(subject_data))] <- subject
    }

    # Create dds and preprocessing + filtering (come giustifico questo input)
    dds <- DESeqDataSetFromMatrix(countData = subject_data, colData = column_data,
        design = ~1)  # not spcifying model 

    dds75 <- dds[rowSums(counts(dds) >= 15) >= nrow(column_data) * 0.75, ]
    message("number of genes: ", nrow(dds75))  # 13284 genes

    # Perform variance stabilization
    dds_norm <- vst(dds75, fitType = "local")

    # Get normalized counts
    norm_counts <- assay(dds_norm) %>%
        t()

    return(norm_counts)
}
    
# Compute WGCNA network matrix
matrix_computation <- function(norm_counts, soft_power, save_net, output_dir, wgcna_file, type = "unsigned", TOMType = "signed", minModuleSize = 0, mergeCutHeight = 0.25, ...) {
    
    library(WGCNA)

    allowWGCNAThreads()

    # convert matrix to numeric
    norm_counts[] <- sapply(norm_counts, as.numeric)

    # Create temp variable for cor (WGCNA uses a different one)
    temp_cor <- cor
    cor <- WGCNA::cor

    # memory estimate w.r.t blocksize
    bwnet <- blockwiseModules(norm_counts, 
        maxBlockSize = 12000, 
        TOMType = TOMType,
        networkType = type,
        saveTOMs = FALSE, 
        power = soft_power, 
        mergeCutHeight = mergeCutHeight, 
        numericLabels = FALSE,
        randomSeed = 1234,
        minModuleSize = minModuleSize,
        verbose = 3)

    cor <- temp_cor
    
    # Save in file
    if (save_net)
        saveRDS(bwnet, paste0(output_dir, wgcna_file))

    # Return data
    return(bwnet)
}
  
# Make plots and intuition for WGCNA
soft_power_intuition <- function(norm_counts, extension_plot = ".png", ...) {

    library(WGCNA)
    power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

    # Call the network topology analysis function input -> rows - subj, cols -
    # geneid
    sft <- pickSoftThreshold(as.data.frame(norm_counts), powerVector = power, networkType = "signed",
        verbose = 5)
    sft.data <- sft$fitIndices

    # Visualization to pick power
    a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) + geom_point() +
        geom_text(nudge_y = 0.1) + geom_hline(yintercept = 0.8, color = "red") +
        labs(x = "Power", y = "Scale free topology model fit, signed R^2") + theme_classic()
    a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) + geom_point() + geom_text(nudge_y = 0.1) +
        labs(x = "Power", y = "Mean Connectivity") + theme_classic()
    save_plot(grid.arrange(a1, a2, nrow = 2), paste0(output_dir, "networkinfo", extension_plot))

    # Ask for user input
    cat("#### INPUT REQUESTED #### \n")
    soft_power <- readline("Which power do you want to select? ")
    soft_power <- as.numeric(unlist(strsplit(soft_power, ",")))

    # Return value
    return(soft_power)
}

# Per caricare la lista di markers dalla cartella, dovrebbe essere ok caricare
# anche senza fare in modo che funzioni per i clusters, magari pensare ad un
# mod di aggiungerlo, ma non priorita
load_log2fc_old <- function(
    markers_analysisPD = "microglia_control_vs_pd_nogenetic",
    markers_analysisGPD = "microglia_control_vs_geneticpd", 
    cluster = FALSE, ...) {

    library("openxlsx")
    library("tibble")
    
    if (isFALSE(cluster)) {
        message(paste0("loading... ", output_folder, "markers_", markers_analysisPD,
            "/expressed_markers_all_", markers_analysisPD, ".xlxs"))
        markers_tablePD <- read.xlsx(paste0(output_folder, "markers_", markers_analysisPD,
            "/expressed_markers_all_", markers_analysisPD, ".xlsx"))

        message(paste0("loading... ", output_folder, "markers_", markers_analysisGPD,
            "/expressed_markers_all_", markers_analysisGPD, ".xlxs"))
        markers_tableGPD <- read.xlsx(paste0(output_folder, "markers_", markers_analysisGPD,
            "/expressed_markers_all_", markers_analysisGPD, ".xlsx"))
    } else {
        message(paste0("loading... ", output_folder, "markers_", markers_analysisPD,
            "/expressed_markers_", cluster, "_", markers_analysisPD, ".xlxs"))
        markers_tablePD <- read.xlsx(paste0(output_folder, "markers_", markers_analysisPD,
            "/expressed_markers_", cluster, "_", markers_analysisPD, ".xlsx"))

        message(paste0("loading... ", output_folder, "markers_", markers_analysisGPD,
            "/expressed_markers_", cluster, "_", markers_analysisGPD, ".xlxs"))
        markers_tableGPD <- read.xlsx(paste0(output_folder, "markers_", markers_analysisGPD,
            "/expressed_markers_", cluster, "_", markers_analysisGPD, ".xlsx"))
    }
    
    return(list(GPD = column_to_rownames(as.data.frame(markers_tableGPD), var = "gene"),
        PD = column_to_rownames(as.data.frame(markers_tablePD), var = "gene")))
}


plot_functions_old <- function(bwnet, norm_counts, column_data, output_dir, extension_plot = ".png",
    reduction_name = "umap_microglia_harmony_reduction", module_of_interest = c(""),
    hub_genes_threshold = c(0.3, 1), which = c(""), cluster = cluster, ...) {

    
    # Library
    library(WGCNA)
    library(magrittr)
    library(ggplot2)
    library(tibble)
    library(openxlsx)
    library(ggpmisc)
    

    write_on_excel <- function(sheet_name, data, wb = NULL, mode = "pvalue") {

        library(openxlsx)
        
        # Create a new workbook and add a worksheet
        if (is.null(wb)) wb <- openxlsx::createWorkbook()
        # detach("package:openxlsx", unload=TRUE)
        openxlsx::addWorksheet(wb, sheet_name)

        # Write the data to the worksheet
        openxlsx::writeDataTable(wb, sheet_name, data, rowNames = TRUE)


        # Define the style for conditional formatting
        if (mode == "colorscale")
            openxlsx::conditionalFormatting(
                wb, 
                sheet = sheet_name, 
                cols = 1:ncol(data)+1, 
                rows = 1:nrow(data)+1, 
                type = "colorScale", 
                rule = c(1, 0, -1),
                style = c("#6386be", "#FFFFFF", "#ec4d50"),
                gradient = TRUE
            )
        else if (mode == "pvalue") {
            openxlsx::conditionalFormatting(
                wb, 
                sheet = sheet_name, 
                cols = 1:ncol(data)+1, 
                rows = 1:nrow(data)+1, 
                type = "expression", 
                rule = "<=0.05",
                style = createStyle(bgFill = "yellow"),
                gradient = TRUE
            )
        }
            
        # Save the workbook
        return(wb)
    }
    make_heatmap <- function(excel_filename, traits) {
        
        # remotes::install_github('kevinblighe/CorLevelPlot')
        # library(CorLevelPlot)
        library(tibble)
        library(openxlsx)

        module_eigengenes <- bwnet$MEs

        # traits_binary <- column_data %>% mutate(across(pathology, ~
        # as.integer(. == levels(pathology))))

        colnames(traits) <- unique(column_data$pathology)
        row.names(traits) <- row.names(column_data)

        heatmap_data <- merge(module_eigengenes, traits, by = "row.names")
        heatmap_data <- heatmap_data %>%
            column_to_rownames(var = "Row.names")

        # save correlation results
        res <- corAndPvalue(x = heatmap_data[, (length(heatmap_data) - ncol(traits) + 1):length(heatmap_data)],
            y = heatmap_data[, 1:(length(heatmap_data) - ncol(traits))])

        excel_filename <- paste0(output_dir, "trait_correlation_heatmap.xlsx")

        wb <- write_on_excel("correllation", as.data.frame(t(res$cor)), mode = "colorscale")
        wb <- write_on_excel("p.value", as.data.frame(t(res$p)), wb = wb)
        wb <- write_on_excel("t.statistic", as.data.frame(t(res$t)), wb = wb)
        openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)
        # head(heatmap_data)


        # # create plot fror wgcna
        # try(dev.off())
        # png(file = paste0(output_dir, "traits_correlation_heatmap", extension_plot))

        # CorLevelPlot(heatmap_data, x = names(heatmap_data)[(length(heatmap_data) -
        #     2):length(heatmap_data)], y = names(heatmap_data)[1:(length(heatmap_data) -
        #     3)], col = c("blue1", "skyblue", "white", "pink", "red"))
        # dev.off()
        message("heatmap saved in: ", output_dir, excel_filename)
    }
    # Define the function to handle different options
    handle_analysis <- function(option, ...) {
    result <- switch(option,
                    "heatmap" = {
                        print("Executing heatmap analysis")
                        # Add your heatmap analysis code here
                        heatmap_result <- "Heatmap result"
                        heatmap_result
                    },
                    "dendro" = {
                        print("Executing dendrogram analysis")
                        # Add your dendrogram analysis code here
                        dendro_result <- "Dendrogram result"
                        dendro_result
                    },
                    "heatmapMRI" = {
                        print("Executing heatmapMRI analysis")
                        # Add your heatmapMRI analysis code here
                        heatmapMRI_result <- "heatmapMRI result"
                        heatmapMRI_result
                    },
                    "violin_plots" = {
                        print("Executing violin plots analysis")
                        # Add your violin plots analysis code here
                        violin_plots_result <- "Violin plots result"
                        violin_plots_result
                    },
                    "normalized_expr" = {
                        print("Executing normalized expression analysis")
                        # Add your normalized expression analysis code here
                        normalized_expr_result <- "Normalized expression result"
                        normalized_expr_result
                    },
                    "module_trait" = {
                        print("Executing module trait analysis")
                        # Add your module trait analysis code here
                        module_trait_result <- "Module trait result"
                        module_trait_result
                    },
                    "histogram_plot" = {
                        print("Executing histogram plot analysis")
                        # Add your histogram plot analysis code here
                        histogram_plot_result <- "Histogram plot result"
                        histogram_plot_result
                    },
                    "significance_membership_scatter" = {
                        print("Executing significance membership scatter plot analysis")
                        # Add your significance membership scatter plot analysis code here
                        significance_membership_scatter_result <- "Significance membership scatter plot result"
                        significance_membership_scatter_result
                    },
                    "corr_matrix" = {
                        print("Executing correlation matrix analysis")
                        # Add your correlation matrix analysis code here
                        corr_matrix_result <- "Correlation matrix result"
                        corr_matrix_result
                    },
                    {
                        print("Invalid option")
                        NA
                    }
    )
    return(result)
    }

    if (!is.character(which) || !is.vector(which)) 
        # stop("which argument must be a character vector containingone or more of these values: 
        # heatmap, TOM, dendro, heatmapMRI, violin_plots, normalized_expr, module_trait, histogram_plot, significance_membership_scatter, corr_matrix")

    # Plotting function
    if ("heatmap" %in% which) {

        # remotes::install_github('kevinblighe/CorLevelPlot')
        library(CorLevelPlot)
        library(tibble)
        library(openxlsx)

        module_eigengenes <- bwnet$MEs

        # traits_binary <- column_data %>% mutate(across(pathology, ~
        # as.integer(. == levels(pathology))))

        traits_binary <- data.frame(lapply(unique(column_data$pathology), function(x) {
            as.numeric(column_data$pathology == x)
        }))

        colnames(traits_binary) <- unique(column_data$pathology)
        row.names(traits_binary) <- row.names(column_data)

        heatmap_data <- merge(module_eigengenes, traits_binary, by = "row.names")
        heatmap_data <- heatmap_data %>%
            column_to_rownames(var = "Row.names")

        # save correlation results
        res <- corAndPvalue(x = heatmap_data[, (length(heatmap_data) - ncol(traits_binary) + 1):length(heatmap_data)],
            y = heatmap_data[, 1:(length(heatmap_data) - ncol(traits_binary))])

        excel_filename <- paste0(output_dir, "trait_correlation_heatmap.xlsx")

        wb <- write_on_excel("correllation", as.data.frame(t(res$cor)), mode = "colorscale")
        wb <- write_on_excel("p.value", as.data.frame(t(res$p)), wb = wb)
        wb <- write_on_excel("t.statistic", as.data.frame(t(res$t)), wb = wb)
        openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)
        # head(heatmap_data)

        # create plot fror wgcna
        # try(dev.off())
        # png(file = paste0(output_dir, "traits_correlation_heatmap", extension_plot))

        # CorLevelPlot(heatmap_data, x = names(heatmap_data)[(length(heatmap_data) -
        #     2):length(heatmap_data)], y = names(heatmap_data)[1:(length(heatmap_data) -
        #     3)], col = c("blue1", "skyblue", "white", "pink", "red"))
        # dev.off()
        

        message("plot saved in: ", output_dir, excel_filename)
    }
    if ("TOM" %in% which) {
        library(WGCNA)

        adjacency_matrix <- adjacency(norm_counts, power = soft_power, type = "signed")
        TOM <- 1 - TOMsimilarity(adjacency_matrix)
        saveRDS(TOM, paste0(output_dir, "TOM.rds"))
        save_plot(TOMplot(TOM, bwnet$dendrograms[[1]]))
        # 6B. Intramodular analysis: Identifying driver genes ---------------

        # Calculating the module dissimilarity eigengenes
        MEDiss = 1 - cor(MEs)

        # Clustering the eigengenes modules
        METree = hclust(as.dist(MEDiss), method = "average")
        MEDissThres = 0.25
        # Plotting the result sizeGrWindow(7, 6)
        # save_plot(plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "") + 
        #     abline(h = MEDissThres, col = "red"))

    }
    if ("dendro" %in% which) {
        
        library(WGCNA)

        png(file = paste0(output_dir, "modules_dendrogram", extension_plot))
        plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
            c("unmerged", "merged"), dendroLabels = FALSE, addGuide = TRUE, hang = 0.03,
            guideHang = 0.05)
        dev.off()
        message("plot saved in: ", output_dir, "modules_dendrogram", extension_plot)

    }
    if ("heatmapMRI" %in% which) {

        # remotes::install_github('kevinblighe/CorLevelPlot')
        library(CorLevelPlot)
        library(tibble)
        library(openxlsx)
        # traits <- load_mri_data_for_wgcna()
        source("scripts/new/mri_analysis.r", local = TRUE)
        traits <- mri_for_wgcna()

        # compute the correlation of the Meigengene with the
        heatmap_data <- bwnet$MEs %>%
            merge(traits, by = "row.names") %>%
            column_to_rownames(var = "Row.names")

        ## Heatmap ---- Check its existence
        excel_filename <- paste0(output_dir, "mri_correlation_heatmap_new.xlsx")
        # if (file.exists(excel_filename))
        #     file.remove(excel_filename)

        res <- corAndPvalue(y = heatmap_data[, (length(heatmap_data) - ncol(traits) + 1):length(heatmap_data)],
            x = heatmap_data[, 1:(length(heatmap_data) - ncol(traits))])


        wb <- write_on_excel("correlation", as.data.frame(res$cor), mode = "colorscale")
        wb <- write_on_excel("p.value", as.data.frame(res$p), wb = wb)
        wb <- write_on_excel("t.statistic", as.data.frame(res$t), wb = wb)
        openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)

        # saveWorkbook(write_on_excel("correllation", as.data.frame(t(res$cor))), excel_filename, overwrite)
        # # write.xlsx(, file = excel_filename, sheetName = "corr",
        #     append = TRUE)
        # write.xlsx(as.data.frame(t(res$p)), file = excel_filename, sheetName = "p.values",
        #     append = TRUE)
        # write.xlsx(as.data.frame(t(res$t)), file = excel_filename, sheetName = "t.statistic",
        #     append = TRUE)
        # library(openxlsx)


        # png(file = paste0(output_dir, "mri_correlation_heatmap", extension_plot))
        # CorLevelPlot(heatmap_data, x = names(heatmap_data)[(length(heatmap_data) -
        #     2):length(heatmap_data)], y = names(heatmap_data)[1:(length(heatmap_data) -
        #     9)], col = c("blue1", "skyblue", "white", "pink", "red"))
        # dev.off()
        
        message("plot saved in: ", output_dir, excel_filename)

    }
    if ("heatmapMRI_zscore" %in% which) {

        # remotes::install_github('kevinblighe/CorLevelPlot')
        library(CorLevelPlot)
        library(tibble)
        library(openxlsx)
        # traits <- load_mri_data_for_wgcna()
        source("scripts/new/mri_analysis.r", local = TRUE)
        traits <- zscore_for_wgcna()

        # compute the correlation of the Meigengene with the
        heatmap_data <- bwnet$MEs %>%
            merge(traits, by = "row.names") %>%
            column_to_rownames(var = "Row.names")
        
        ## Heatmap ---- Check its existence
        excel_filename <- paste0(output_dir, "mri_correlation_heatmap_zscore.xlsx")
        # if (file.exists(excel_filename))
        #     file.remove(excel_filename)

        res <- corAndPvalue(
            y = heatmap_data[, (length(heatmap_data) - ncol(traits) + 1):length(heatmap_data)],
            x = heatmap_data[, 1:(length(heatmap_data) - ncol(traits))])

        wb <- write_on_excel("correlation", as.data.frame(res$cor), mode = "colorscale")
        wb <- write_on_excel("p.value", as.data.frame(res$p), wb = wb)
        wb <- write_on_excel("t.statistic", as.data.frame(res$t), wb = wb)
        openxlsx::saveWorkbook(wb, excel_filename, overwrite = TRUE)

        # saveWorkbook(write_on_excel("correllation", as.data.frame(t(res$cor))), excel_filename, overwrite)
        # # write.xlsx(, file = excel_filename, sheetName = "corr",
        #     append = TRUE)
        # write.xlsx(as.data.frame(t(res$p)), file = excel_filename, sheetName = "p.values",
        #     append = TRUE)
        # write.xlsx(as.data.frame(t(res$t)), file = excel_filename, sheetName = "t.statistic",
        #     append = TRUE)
        # library(openxlsx)


        # png(file = paste0(output_dir, "mri_correlation_heatmap", extension_plot))
        # CorLevelPlot(heatmap_data, x = names(heatmap_data)[(length(heatmap_data) -
        #     2):length(heatmap_data)], y = names(heatmap_data)[1:(length(heatmap_data) -
        #     9)], col = c("blue1", "skyblue", "white", "pink", "red"))
        # dev.off()
        
        message("heatmap saved in: ", output_dir, excel_filename)

    }
    # Create violin plots of the modules
    if ("violin_plots" %in% which) {

        module_df <- data.frame(gene_id = names(bwnet$colors), colors =  WGCNA::labels2colors(bwnet$colors))

        traits_binary <- data.frame(lapply(unique(column_data$pathology), function(x) {
            as.numeric(column_data$pathology == x)
        }))
        colnames(traits_binary) <- unique(column_data$pathology)
        row.names(traits_binary) <- row.names(column_data)

        module_gene_mapping <- as.data.frame(bwnet$colors)
        for (color in unique(module_gene_mapping[[1]])) {

            module_genes <- module_gene_mapping %>%
                filter(bwnet$colors == color) %>%
                rownames()

            norm_counts_subset <- norm_counts[, colnames(norm_counts) %in% module_genes]

            # Calculate the gene significance and associated p-values
            # controllare qui come e stato fatto il testo eh (significance)
            gene_signf_corr <- cor(norm_counts_subset, traits_binary, use = "p")
            gene_signf_corr_pvals <- corPvalueStudent(gene_signf_corr, n_samples)

            # For each module creates te violin plots for the 10 genes with the highest correlation to the traits?
            # TODO: ceck + source
            for (column in colnames(gene_signf_corr_pvals)) {

                genes <- gene_signf_corr_pvals %>%
                  as.data.frame() %>%
                  arrange(column) %>%
                  head(10) %>%
                  row.names()
                
                # Calls the violin plot function from seurat_utiles
                violin_plot(seurat_object, genes, name = paste0("WGCNA/test_2", color,
                  column), markers_analysisPD = "microglia_control_vs_pd_nogenetic_all",
                  markers_analysisGPD = "microglia_control_vs_geneticpd_all")

            }

            # Using the gene significance you can identify genes that have a
            # high significance for trait of interest Using the module
            # membership measures you can identify genes with high module
            # membership in interesting modules.
        }
    }

    if ("normalized_expr" %in% which) {
        message("normalized expression NOT IMPLEMENTED")

        next

        # pick out a few modules of interest here
        modules_of_interest = c("green", "turquoise", "tan")

        # df of module colors
        module_df <- data.frame(gene_id = names(bwnet$colors), colors =  WGCNA::labels2colors(bwnet$colors))

        # Pull out list of genes in that module
        submod = module_df %>%
            subset(colors %in% modules_of_interest)

        row.names(module_df) <- module_df$gene_id

        subexpr = t(norm_counts[, submod$gene_id])

        submod_df = data.frame(subexpr) %>%
            mutate(gene_id = row.names(.)) %>%
            pivot_longer(-gene_id) %>%
            mutate(module = module_df[gene_id, ]$color)

        submod_df %>%
            ggplot(., aes(x = name, y = value, group = gene_id)) + geom_line(aes(color = module),
            alpha = 0.2) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) +
            facet_grid(rows = vars(module)) + labs(x = "subject", y = "normalized expression") %>%
            save_plot(., paste0(output_dir, "normalized_expression", extension_plot))

    }
    if ("module_trait" %in% which) {

        message("module trait NOT implemented")

        next

        # Convert labels to colors for plotting
        mergedColors =  WGCNA::labels2colors(bwnet$colors)

        # Get Module Eigengenes per cluster
        MEs0 <- moduleEigengenes(norm_counts, mergedColors)$eigengenes

        # Reorder modules so similar modules are next to each other
        MEs0 <- orderMEs(MEs0)
        module_order = names(MEs0) %>%
            gsub("ME", "", .)

        # Add treatment names
        MEs0$patient = row.names(MEs0)

        # tidy & plot data
        mME = MEs0 %>%
            pivot_longer(-patient) %>%
            mutate(name = gsub("ME", "", name), name = factor(name, levels = module_order))

        mME %>%
            ggplot(., aes(x = patient, y = name, fill = value)) + geom_tile() + theme_bw() +
            scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0,
                limit = c(-1, 1)) + theme(axis.text.x = element_text(angle = 90)) +
            labs(title = "Module-trait Relationships", y = "Modules", fill = "corr") %>%
                save_plot(., paste0(output_dir, "Module-trait_rel", extension_plot))
    }
    if ("histogram_plot_significance" %in% which) {


            # markers <- markers_list[[markers_name]]
        source("scripts/new/mri_analysis.r", local = TRUE)

        traits <- zscore_for_wgcna()
        gene_significance <- cor(norm_counts[rownames(traits),], traits, use = "p")
        # Creating merged dataset for plotting
        plot_df <- merge(data.frame(bwnet$colors), gene_significance,
            by = "row.names", all = FALSE)


        plot_df$bwnet.colors <- factor(plot_df$bwnet.colors)
        for (column in c("superiorfrontal_thickness", "frontal_cortex_thickness")) {
            # Plot (with correct colors)
            color_mapping <- setNames(c(levels(plot_df$bwnet.colors)), levels(plot_df$bwnet.colors))
            save_plot(
                ggplot(
                    plot_df, aes(x = .data[["bwnet.colors"]], y = abs(.data[[column]]), color = .data[["bwnet.colors"]]))
                     + geom_boxplot()
                     + scale_color_manual(values = color_mapping) 
                     + theme(
                        legend.position = "none",            # Remove the legend
                        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
                        )
                , paste0(output_dir,column, "_boxplot", extension_plot))

            save_plot(
                ggplot(
                    plot_df, aes(x = .data[["bwnet.colors"]], y = .data[[column]], color = .data[["bwnet.colors"]]))
                     + geom_violin(alpha = 0.2) 
                     + geom_jitter(alpha = 0.8, size = 0.2) 
                     + scale_color_manual(values = color_mapping) 
                     + theme(
                        legend.position = "none",            # Remove the legend
                        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
                        ),
                paste0(output_dir, column, "_violin", extension_plot))
        }
        #}  # scale_fill_manual(values = c('red', 'blue', 'green')) +
    }
    if ("histogram_plot" %in% which) {

        markers_list <- load_log2fc(cluster = cluster , ...)
        #for (markers_name in names(markers_list)) {
            markers_name <- "PD"

            # markers <- markers_list[[markers_name]]

            # Creating merged dataset for plotting
            plot_df <- merge(data.frame(bwnet$colors), markers_list[[markers_name]],
                by = "row.names", all = FALSE)

            plot_df$abs_avg_log2FC <- abs(plot_df$avg_log2FC)
            plot_df$bwnet.colors <- factor(plot_df$bwnet.colors)

            # Plot (with correct colors)
            color_mapping <- setNames(c(levels(plot_df$bwnet.colors)), levels(plot_df$bwnet.colors))
            save_plot(
                ggplot(
                    plot_df, aes(x = bwnet.colors, y = abs_avg_log2FC, color = bwnet.colors))
                     + geom_boxplot()
                     + scale_color_manual(values = color_mapping) 
                     + theme(
                        legend.position = "none",            # Remove the legend
                        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
                        )
                , paste0(output_dir,markers_name, "_boxplot", extension_plot))

            save_plot(
                ggplot(
                    plot_df, aes(x = bwnet.colors, y = avg_log2FC, color = bwnet.colors))
                     + geom_violin(alpha = 0.2) 
                     + geom_jitter(alpha = 0.8, size = 0.2) 
                     + scale_color_manual(values = color_mapping) 
                     + theme(
                        legend.position = "none",            # Remove the legend
                        axis.text.x = element_text(angle = 45, hjust = 1)  # Rotate x-axis labels
                        ),
                paste0(output_dir, markers_name, "_violin", extension_plot))
        #}  # scale_fill_manual(values = c('red', 'blue', 'green')) +
    }
    # il significance membership precedente (one plot per module)
    if ("membership_log2fc_scatter" %in% which) {
        library(tibble)

        # Compute membershio measure
        n_samples <- nrow(norm_counts)
        module_membership_measure <- cor(norm_counts, bwnet$MEs, use = "p")  # ahh ho invertito e quindo ho una direzione divers
        # module_membership_measure.pvals <-
        # corPvalueStudent(module_membership_measure, n_samples)

        markers_list <- load_log2fc(cluster = cluster, ...)
        
        # Iterate though gpd and pd
    # for (markers_name in names(markers_list)) {
        markers_name <- "PD"
        # Creating merged dataset for plotting (merge colors with
        # expression values with membership measures)
        plot_df <- column_to_rownames(merge(data.frame(bwnet$colors), markers_list[[markers_name]],
            by = "row.names", all = FALSE), var = "Row.names")
        plot_df <- column_to_rownames(merge(plot_df, module_membership_measure,
            by = "row.names", all = FALSE), var = "Row.names")
        
        plot_df$abs_avg_log2FC <- abs(plot_df$avg_log2FC)

        for (module in unique(bwnet$colors)) {
                # markers <- markers_list[[markers_name]] as
                library("ggrepel")

                # Pull out list of genes in that module
                submod <- plot_df %>%
                    subset(bwnet.colors == module)
                message(module, " module size: ", nrow(submod))

                # membership score
                module_membership_column_abs <- paste0("ME", module, "_abs")
                module_membership_column <- paste0("ME", module)
                submod[[module_membership_column_abs]] <- abs(submod[[module_membership_column]])
                plot_df[[module_membership_column_abs]] <- abs(plot_df[[module_membership_column]])

                # Skip modules that are too small
                if (nrow(submod) < 9)
                  next

                # Function to visualize and plot genes that are higher than
                # cutoff parameters
                plot_hub_genes <- function(df) {

                  FCcutoff <- hub_genes_threshold[[2]]
                  MMcutoff <- hub_genes_threshold[[1]]

                  save_plot(ggplot(df, aes(x = .data[[module_membership_column]],
                    y = .data[["avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                    geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
                      x) + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                    after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5,
                    color = "black", label.x = "right", label.y = "top") + theme_minimal() +
                    geom_hline(yintercept = c(-FCcutoff, FCcutoff), colour = "black") +
                    geom_vline(xintercept = c(-MMcutoff, MMcutoff), colour = "black") +
                    geom_text_repel(data = subset(df, abs(df[["avg_log2FC"]]) > FCcutoff &
                      abs(df[[module_membership_column]]) > MMcutoff), aes(label = row.names(subset(df,
                      abs(df[["avg_log2FC"]]) > FCcutoff & abs(df[[module_membership_column]]) >
                        MMcutoff))), xlim = c(NA, NA), ylim = c(NA, NA)), paste0(output_dir,
                    markers_name, "_scatter_hub/", module, "_scatter_plot", extension_plot))

                  if (module %in% module_of_interest) {
                    hub_genes_list <- row.names(subset(submod, abs(df[["avg_log2FC"]]) >
                      FCcutoff & abs(df[[module_membership_column]]) > MMcutoff))

                    # feature_plots(seurat_object_microglia, hub_genes_list,
                    # name = paste0('WGCNA_', name, '/feature_plot_hubgene_',
                    # markers_name, '_', module, '/'), reduction_name =
                    # reduction_name)

                    write.xlsx(data.frame(hub_genes_list), paste0(output_dir, markers_name,
                      module, "_hubgene.xlsx"))

                    many_plots(seurat_object_microglia, which = c("dotplot", "heatmap"),
                      cluster_column = "microglia_clusters", name = paste0("WGCNA_modules/",
                        module, markers_name), markers = hub_genes_list)
                  }
                }

                # scatter plot of the module genes
                save_plot(ggplot(submod, aes(x = .data[[module_membership_column]],
                  y = .data[["avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
                    x) + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                  after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                  label.x = "right", label.y = "top") + theme_minimal(), paste0(output_dir,
                  markers_name, "_scatter/", module, "_scatter_plot", extension_plot))

                # scatter plot of all the genes
                save_plot(ggplot(plot_df, aes(x = .data[[module_membership_column]],
                  y = .data[["avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                    fill = "white") + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                  after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                  label.x = "right", label.y = "top") + theme_minimal(), paste0(output_dir,
                  markers_name, "_scatter_all/", module, "_scatter_plot", extension_plot))

                # scatter plot of the module genes
                save_plot(ggplot(submod, aes(x = .data[[module_membership_column_abs]],
                  y = .data[["abs_avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                  geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
                    x) + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                  after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                  label.x = "right", label.y = "top") + theme_minimal(), paste0(output_dir,
                  markers_name, "_scatter_abs/", module, "_scatter_plot", extension_plot))
                
                # scatter plot of all the genes
                save_plot(ggplot(plot_df, aes(x = .data[[module_membership_column_abs]],
                  y = .data[["abs_avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                  geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                    fill = "white") + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                  after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                  label.x = "right", label.y = "top") + theme_minimal(), paste0(output_dir,
                  markers_name, "_scatter_all_abs/", module, "_scatter_plot", extension_plot))

                if ("hub_genes" %in% which)
                  plot_hub_genes(submod)

            }
        # }

        # plots[[module]] <- ggplot()
    }
    # one plot per region (significance)
    if ("significance_log2fc_scatter" %in% which) {
        library(tibble)

        # Compute membershio measure
        source("scripts/new/mri_analysis.r", local = TRUE)
        traits <- zscore_for_wgcna()

        gene_significance <- cor(norm_counts, traits, use = "p")
        # module_membership_measure.pvals <-
        # corPvalueStudent(module_membership_measure, n_samples)

        markers_list <- load_log2fc(cluster = cluster, ...)
        
        # Iterate though gpd and pd
        # for (markers_name in names(markers_list)) {
            markers_name <- "PD"
            # Creating merged dataset for plotting (merge colors with
            # expression values with membership measures)
            plot_df <- column_to_rownames(merge(data.frame(bwnet$colors), markers_list[[markers_name]],
                by = "row.names", all = FALSE), var = "Row.names")
            plot_df <- column_to_rownames(merge(plot_df, gene_significance,
                by = "row.names", all = FALSE), var = "Row.names")

            plot_df$bwnet.colors <- factor(plot_df$bwnet.colors)
            plot_df$abs_avg_log2FC <- abs(plot_df$avg_log2FC)


            # for (region_significance_column in unique(colnames(gene_significance))) {
            for (region_significance_column in  c("superiorfrontal_thickness", "frontal_cortex_thickness")) {

                region_significance_column_abs <- paste0(region_significance_column, "_abs")
                plot_df[[region_significance_column_abs]] <- abs(plot_df[[region_significance_column]])

                for (module in levels(plot_df$bwnet.colors)) {

                    submod <- plot_df %>%
                        subset(bwnet.colors == module)
                    # scatter plot of all the genes
                    save_plot(ggplot(submod, aes(x = .data[[region_significance_column]],
                        y = .data[["avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                        geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                        fill = "white") + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                        after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                        label.x = "right", label.y = "top") + theme_minimal() + theme(legend.position = "none"), paste0(output_dir,
                        markers_name, "_scatter_significance_log2fc_all/", region_significance_column, "_", module, "_scatter_plot", extension_plot))
                    
                    # scatter plot of all the genes
                    save_plot(ggplot(submod, aes(x = .data[[region_significance_column_abs]],
                        y = .data[["abs_avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                        geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                        fill = "white") + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                        after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                        label.x = "right", label.y = "top") + theme_minimal(), paste0(output_dir,
                        markers_name, "_scatter_significance_log2fc_abs/", region_significance_column, "_", module, "_scatter_plot", extension_plot))

                }

                # membership score
    
            }
        # }

        # plots[[module]] <- ggplot()
    }
    # significance membersip (one plot per module)
    if ("significance_membership_scatter" %in% which) {
        library(tibble)

        # Compute membershio measure
        n_samples <- nrow(norm_counts)
        module_membership_measure <- cor(norm_counts, bwnet$MEs, use = "p")  # ahh ho invertito e quindo ho una direzione divers
        # module_membership_measure.pvals <-
        # corPvalueStudent(module_membership_measure, n_samples)
        source("scripts/new/mri_analysis.r", local = TRUE)
        traits <- zscore_for_wgcna()

        gene_significance <- cor(norm_counts, traits, use = "p")
        # Iterate though gpd and pd
        # for (markers_name in names(markers_list)) {
            markers_name <- "PD"
            # Creating merged dataset for plotting (merge colors with
            # expression values with membership measures)
            plot_df <- column_to_rownames(merge(data.frame(bwnet$colors), gene_significance,
                by = "row.names", all = FALSE), var = "Row.names")
            plot_df <- column_to_rownames(merge(plot_df, module_membership_measure,
                by = "row.names", all = FALSE), var = "Row.names")
            
            # plot_df$abs_avg_log2FC <- abs(plot_df$avg_log2FC)

            for (module in unique(bwnet$colors)) {
                # markers <- markers_list[[markers_name]] as
                library("ggrepel")

                # Pull out list of genes in that module
                submod <- plot_df %>%
                    subset(bwnet.colors == module)
                # message(module, " module size: ", nrow(submod))

                # membership score
                module_membership_column_abs <- paste0("ME", module, "_abs")
                module_membership_column <- paste0("ME", module)
                submod[[module_membership_column_abs]] <- abs(submod[[module_membership_column]])
                plot_df[[module_membership_column_abs]] <- abs(plot_df[[module_membership_column]])

                # Skip modules that are too small
                if (nrow(submod) < 9)
                  next

                # Function to visualize and plot genes that are higher than
                # cutoff parameters
                plot_hub_genes <- function(df) {

                  FCcutoff <- hub_genes_threshold[[2]]
                  MMcutoff <- hub_genes_threshold[[1]]

                  save_plot(ggplot(df, aes(x = .data[[module_membership_column]],
                    y = .data[["avg_log2FC"]])) + geom_point(color = module, alpha = 0.5) +
                    geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
                      x) + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                    after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5,
                    color = "black", label.x = "right", label.y = "top") + theme_minimal() +
                    geom_hline(yintercept = c(-FCcutoff, FCcutoff), colour = "black") +
                    geom_vline(xintercept = c(-MMcutoff, MMcutoff), colour = "black") +
                    geom_text_repel(data = subset(df, abs(df[["avg_log2FC"]]) > FCcutoff &
                      abs(df[[module_membership_column]]) > MMcutoff), aes(label = row.names(subset(df,
                      abs(df[["avg_log2FC"]]) > FCcutoff & abs(df[[module_membership_column]]) >
                        MMcutoff))), xlim = c(NA, NA), ylim = c(NA, NA)), paste0(output_dir,
                    markers_name, "_scatter_hub/", module, "_scatter_plot", extension_plot))

                  if (module %in% module_of_interest) {
                    hub_genes_list <- row.names(subset(submod, abs(df[["avg_log2FC"]]) >
                      FCcutoff & abs(df[[module_membership_column]]) > MMcutoff))

                    # feature_plots(seurat_object_microglia, hub_genes_list,
                    # name = paste0('WGCNA_', name, '/feature_plot_hubgene_',
                    # markers_name, '_', module, '/'), reduction_name =
                    # reduction_name)

                    write.xlsx(data.frame(hub_genes_list), paste0(output_dir, markers_name,
                      module, "_hubgene.xlsx"))

                    many_plots(seurat_object_microglia, which = c("dotplot", "heatmap"),
                      cluster_column = "microglia_clusters", name = paste0("WGCNA_modules/",
                        module, markers_name), markers = hub_genes_list)
                  }
                }
                for (column in c("superiorfrontal_thickness", "frontal_cortex_thickness")) {

                    save_plot(ggplot(submod, aes(x = .data[[module_membership_column]],
                    y = .data[[column]])) + geom_point(color = module, alpha = 0.5) +
                    geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
                        x) + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                    after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                    label.x = "right", label.y = "top") + theme_minimal(), paste0(output_dir,
                    markers_name, "_scatter/", module, "_", column, "_scatter_plot", extension_plot))

                    # scatter plot of all the genes
                    save_plot(ggplot(plot_df, aes(x = .data[[module_membership_column]],
                    y = .data[[column]])) + geom_point(color = module, alpha = 0.5) +
                    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                        fill = "white") + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                    after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                    label.x = "right", label.y = "top") + theme_minimal(), paste0(output_dir,
                    markers_name, "_scatter_all/", module, "_", column, "_scatter_plot", extension_plot))

                    # scatter plot of the module genes
                    save_plot(ggplot(submod, aes(x = .data[[module_membership_column_abs]],
                    y = abs(.data[[column]]))) + geom_point(color = module, alpha = 0.5) +
                    geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~
                        x) + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                    after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                    label.x = "right", label.y = "top") + theme_minimal(), paste0(output_dir,
                    markers_name, "_scatter_abs/", module, "_", column, "_scatter_plot", extension_plot))
                    
                    # scatter plot of all the genes
                    save_plot(ggplot(plot_df, aes(x = .data[[module_membership_column_abs]],
                    y = abs(.data[[column]]))) + geom_point(color = module, alpha = 0.5) +
                    geom_smooth(method = "lm", se = FALSE, formula = y ~ x, color = "black",
                        fill = "white") + stat_poly_eq(formula = y ~ x, aes(label = paste(after_stat(eq.label),
                    after_stat(rr.label), sep = "~~~~")), parse = TRUE, size = 5, color = "black",
                    label.x = "right", label.y = "top") + theme_minimal(), paste0(output_dir,
                    markers_name, "_scatter_all_abs/", module, "_", column, "_scatter_plot", extension_plot))
                }
                # scatter plot of the module genes 
 

                if ("hub_genes" %in% which)
                  plot_hub_genes(submod)

            }
        # }

        # plots[[module]] <- ggplot()
    }
    # scatter of module correlation / avg significance
    if ("correlation_avglog2fc_scatter" %in% which){
        library(ggrepel)
        # library(tibble)

        source("scripts/new/mri_analysis.r", local = TRUE)
        # Compute membershio measure - source?
        traits <- zscore_for_wgcna()
        # compute the correlation of the Meigengene with the
        heatmap_data <- bwnet$MEs %>%
            merge(traits, by = "row.names") %>%
            column_to_rownames(var = "Row.names")

        correlation <- t(corAndPvalue(x = heatmap_data[, (length(heatmap_data) - ncol(traits) + 1):length(heatmap_data)],
            y = heatmap_data[, 1:(length(heatmap_data) - ncol(traits))])$cor)

        markers_list <- load_log2fc(cluster = cluster, ...)

        #for (markers_name in names(markers_list)) {
        markers_name <- "PD"

            # Creating merged dataset for plotting (merge colors with
        # expression values with membership measures)
        plot_df <- column_to_rownames(merge(data.frame(bwnet$colors), markers_list[[markers_name]],
            by = "row.names", all = FALSE), var = "Row.names")
        plot_df[["gene"]] <- rownames(plot_df)

        significance <- sapply(unique(bwnet$colors), 
            function(module) {

                submod <- plot_df %>%
                        subset(bwnet.colors == module)
                        
                submod <- submod[, c("avg_log2FC", "gene")]

                # Skip modules that are too small 
                # if (nrow(submod) < 9) next
                # else 
                return(mean(abs(submod$avg_log2FC)))
                }
            )
        for (region in colnames(correlation)) {
            submod <- data.frame(abs(correlation[, region]), significance)
            colnames(submod)[1] <- "gene_significance"
            colnames(submod)[2] <- "abs_log2fc"
            submod$color <- substring(rownames(submod), 3)
            save_plot(
                ggplot(
                    submod, aes(x = gene_significance, y = abs_log2fc)) +
                    geom_point(alpha = 0.5) +
                    geom_smooth(method = "lm", se = FALSE, color = "black", formula = y ~ x) +
                    stat_poly_eq(
                        formula = y ~ x,
                        aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~~")),
                        parse = TRUE,
                        size = 5,
                        color = "black",
                        label.x = "right",
                        label.y = "top"
                    ) +
                    geom_text_repel(aes(label = color), size = 3) +
                    theme_minimal() + 
                    theme(legend.position = "none"),
                paste0(output_dir, markers_name, "correlation_significance_scatter/", region, "_scatter_plot", extension_plot))
        }
    }

    if ("corr_matrix" %in% which) {

        library(openxlsx)
        # TODO: check if implemented
        cor_matrix <- cor(bwnet$MEs)
        # from gpt, to test: Install and load the necessary packages
        # install.packages('gplots') library(gplots)


        wb <- write_on_excel("correllation", as.data.frame(cor_matrix))
        openxlsx::saveWorkbook(wb, paste0(output_dir, "module_correlation.xlsx"), overwrite = TRUE)

        # Define color palette
        my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

        # Create the heatmap


        png(file = paste0(output_dir, "modules_correlation_heatmap", extension_plot))
        heatmap(cor_matrix)
        dev.off()
        message("plot saved in: ", output_dir, "modules_correlation_heatmap", extension_plot)
    }
}


# To save list of he module genes in excel files
save_module_genes <- function(bwnet, norm_counts, column_data, output_dir) {
    library("tibble")

    # Compute module membership
    module_membership_measure <- cor(bwnet$MEs, norm_counts, use = "p")

    # Create a df with modules
    module_df <- data.frame(gene_id = names(bwnet$colors), colors =  WGCNA::labels2colors(bwnet$colors),
        row.names = names(bwnet$colors))

    # Create a df that contains the membership measure
    module_df_membership <- column_to_rownames(merge(data.frame(bwnet$colors), t(module_membership_measure),
        by = "row.names", all = FALSE), var = "Row.names")
    module_df_membership$gene_id <- row.names(module_df_membership)  # because otherwise the selection afterwards doesnt include rownames

    list_of_dfs <- list()

    for (module in unique(bwnet$colors)) {

        list_of_dfs[[module]] <- module_df_membership[, c("gene_id", paste0("ME",
            module))] %>%
            subset(module_df_membership$bwnet.colors %in% module)
        message(module, " module size: ", nrow(list_of_dfs[[module]]))

    }
    # tranform this in for cycle list_of_dfs <-
    # lapply(sort(unique(bwnet$colors)), function(x) {
    # data.frame(module_df_membership[bwnet$colors %in% x, c('gene_id',
    # paste0('ME', x))]) })

    # names(list_of_dfs) <- sort(unique(bwnet$colors))

    wb <- openxlsx::createWorkbook()

    # Loop through the list of dataframes
    for (name in names(list_of_dfs)) {

        # Add dataframe to a new sheet in the Excel workbook
        openxlsx::addWorksheet(wb, sheetName = name) # openxlsx
        # sheet <- createSheet(wb, sheetName = name)
        openxlsx::writeData(wb, sheet = name, x = list_of_dfs[[name]], colNames = FALSE) # openxlsx
        # addDataFrame(list_of_dfs[[name]], sheet, row.names = FALSE)
    }

    # Save the Excel workbook
    openxlsx::saveWorkbook(wb, file = paste0(output_dir, "module_genes.xlsx"), overwrite = TRUE)
    message("module genes saved in: ", paste0(output_dir, "module_genes.xlsx"))

    write.xlsx(table(bwnet$colors), paste0(output_dir, "module_recap.xlsx"))
    message("module recap saved in: ", paste0(output_dir, "module_recap.xlsx"))

}

