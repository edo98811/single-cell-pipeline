# Author: Edoardo Filippi 
# mail: efilippi@uni-mainz.de

#' @details To run the WGCNA on complex dta coming from a seurat oject
# 
# Arguments:
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


# Analysis using the wgcna package
wgcna_main <- function(seurat_object, name = "test", wgcna_file = "bwnet.rds", save_net = TRUE, 
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

        if (!isFALSE(soft_power) && !is.numeric(soft_power)) {
            stop("soft_power argument must be a numeric value or FALSE")
        }


    }

    cat("Running wgcna...\n")
    message(paste0("Parameters: name: ", name, " - wgcna_file: ", wgcna_file, " - save_net: ",
        save_net, " - load_net: ", load_net, " - soft_power: ", soft_power))

    # Dependencies
    # check_packages(c("Seurat", "wgcna", "GEOquery", "tidyverse", "gridExtra", "dplyr",
    #     "readxl", "openxlsx", "DESeq2", "ggpmisc"))
    library(magrittr)

    # Set folders
    output_dir <- set_up_output(paste0(output_folder, "wgcna_", name, "/"), message)

    # Main loop Data preparation
    column_data <- prepare_column_data(seurat_object, ...)
    norm_counts <- prepare_data(seurat_object, column_data, ...)

    # Soft power
    if (isFALSE(soft_power) & !load_net)
        soft_power <- soft_power_intuition(norm_counts, output_dir, ...)

    # Matrix computation
    if (load_net) bwnet <- readRDS(paste0(output_dir, wgcna_file)) 
    else bwnet <- matrix_computation(norm_counts, soft_power, save_net, output_dir, wgcna_file, ...)

    # Plots and other functions
    source("scripts/wgcna_plots.r", local = TRUE)
    plot_functions(bwnet, norm_counts, column_data, output_dir, ...)
    save_module_genes(bwnet, norm_counts, column_data, output_dir)

    return(seurat_object)
}

# Prepare info on subjects
prepare_column_data <- function(seurat_object, subject_pathology_column = "subject_pathology", subject_column = "subject", ...) {

    # Extract metadata from seurat object
    meta_data <- seurat_object@meta.data

    # Create subject meta data
    column_data <- data.frame(subject = unique(meta_data[[subject_column]]), row.names = unique(meta_data[[subject_column]]))
    column_data$pathology <- lapply(unique(meta_data[[subject_column]]), function(x) unique(meta_data[meta_data[[subject_column]] ==
        x, ][[subject_pathology_column]]))

    return(column_data)
}

# Prepare the count matrix
prepare_data <- function(seurat_object, column_data, subject_column = "subject", ...) {

    library(DESeq2)

    # Get data from seurat object (transpose to have genes on columns and cells
    # on rows)
    counts_data <- t(GetAssayData(object = seurat_object, assay = "RNA", layer = "counts"))
    meta_data <- seurat_object@meta.data

    # Group counts by subject
    subject_data <- data.frame(matrix(ncol = 0, nrow = length(colnames(counts_data))),
        row.names = colnames(counts_data))
    for (subject in unique(meta_data[[subject_column]])) {
        subject_data <- cbind(subject_data, colSums(counts_data[row.names(meta_data[meta_data[[subject_column]] ==
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
    
# Compute wgcna network matrix
matrix_computation <- function(norm_counts, soft_power, save_net, output_dir, wgcna_file, type = "unsigned", TOMType = "signed", minModuleSize = 0, mergeCutHeight = 0.25, ...) {
    
    library(WGCNA)

    allowWGCNAThreads()

    # convert matrix to numeric
    norm_counts[] <- sapply(norm_counts, as.numeric)

    # Create temp variable for cor (wgcna uses a different one)
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
  
# Make plots and intuition for wgcna
soft_power_intuition <- function(norm_counts, output_dir, extension_plot = ".png", ...) {

    library(WGCNA)
    library(ggplot2)
    library(gridExtra)

    power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

    # Call the network topology analysis function input -> rows - subj, cols -
    # geneid
    sft <- pickSoftThreshold(as.data.frame(norm_counts), powerVector = power, networkType = "signed",
        verbose = 5)
    sft_data <- sft$fitIndices

    # Visualization to pick power
    a1 <- ggplot(sft_data, aes(Power, SFT.R.sq, label = Power)) + geom_point() +
        geom_text(nudge_y = 0.1) + geom_hline(yintercept = 0.8, color = "red") +
        labs(x = "Power", y = "Scale free topology model fit, signed R^2") + theme_classic()
    a2 <- ggplot(sft_data, aes(Power, mean.k., label = Power)) + geom_point() + geom_text(nudge_y = 0.1) +
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

